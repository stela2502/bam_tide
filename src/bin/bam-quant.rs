// src/bin/bam-quant.rs
//
// Quantify a 10x-style BAM against a splice index (GTF-derived) into scdata,
// with optional genome-refined read cleanup and optional SNP ref/alt side-channel.
//
// Pipeline:
// 1) Load splice index
// 2) Optionally load genome FASTA
// 3) Optionally load SNP VCF
// 4) Stream BAM:
//      - collect CB/UB
//      - build AlignedRead if genome and/or VCF were supplied
//      - optionally refine AlignedRead against Genome
//      - build SplicedRead from the BAM record for gene/transcript quantification
//        (or from refined AlignedRead once that conversion exists)
//      - keep optional AlignedRead for SNP matching
// 5) Rayon over chunks:
//      - normal gene/transcript quantification
//      - optional SNP ref/alt quantification
// 6) Export:
//      - normal gene/transcript 10x-style matrix
//      - optional intronic matrix
//      - optional SNP ref/alt matrices
//
// Important:
// Genome refinement is an addition to the normal quantification path, not a separate mode.
// SNPs are collected as a side-channel when --vcf is supplied.
//
// Expected optional usage:
//   bam-quant \
//     --quant-mode transcript \
//     --bam possorted_genome_bam.bam \
//     --index splice.idx \
//     --genome genome.fa \
//     --vcf variants.vcf.gz \
//     --outpath quant_out
//
// This writes:
//   quant_out/
//   quant_out_ref/
//   quant_out_alt/
// and optionally:
//   quant_out_intronic/

use anyhow::{anyhow, Context, Result};
use clap::{Parser, ValueEnum};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader, Record};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::collections::HashSet;

use bam_tide::core::ref_block::record_to_blocks;
use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};

use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex, SplicedRead, Strand};

use int_to_str::int_to_str::IntToStr;
use mapping_info::MappingInfo;

use scdata::cell_data::GeneUmiHash;
use scdata::{MatrixValueType, Scdata};

use snp_index::{AlignedRead, Genome, RefineOptions, SnpIndex, VcfReadOptions };

const CHUNK: usize = 2_000_000;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum QuantMode {
    Gene,
    Transcript,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "bam-quant",
    about = "Quantify 10x BAM against splice index into scdata, optionally collecting SNP ref/alt matrices"
)]
pub struct QuantCli {
    /// Input BAM (10x / CellRanger BAM)
    #[arg(long, short)]
    pub bam: PathBuf,

    /// Splice index path (built from GTF beforehand)
    #[arg(long, short)]
    pub index: PathBuf,

    /// Outpath for the 10x mtx-formatted outfiles
    #[arg(long, short)]
    pub outpath: PathBuf,

    /// Split Intronic from rest.
    ///
    /// This is currently not recommended as exon/intron detection seems to be too strict
    /// for normal sequencing data.
    #[arg(long, short, default_value_t = false)]
    pub split_intronic: bool,

    /// Minimum MAPQ
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Use only read1 (recommended for 10x; reduces duplicate mate-counting noise)
    #[arg(long, default_value_t = false)]
    pub read1_only: bool,

    /// Rayon thread count (0 = default)
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Collect Gene or Transcript names
    #[arg(long, value_enum, default_value_t = QuantMode::Gene)]
    pub quant_mode: QuantMode,

    /// Max reads to process (debug/dev)
    #[arg(long)]
    pub max_reads: Option<usize>,

    /// Min read counts per reported cell (debug/dev)
    #[arg(long, default_value_t = 400)]
    pub min_cell_counts: usize,

    /// Optional reference genome FASTA.
    ///
    /// If supplied, BAM-derived AlignedRead objects are refined against the genome.
    /// This is useful for SNP matching and can also improve downstream splice/transcript
    /// detection once the refined alignment is used to build the SplicedRead.
    #[arg(long)]
    pub genome: Option<PathBuf>,

    /// Optional SNP VCF.
    ///
    /// If supplied, SNP ref/alt matrices are written in addition to the normal
    /// gene/transcript matrix. Requires --genome.
    #[arg(long)]
    pub vcf: Option<PathBuf>,

    /// Minimum SNP anchor/support passed to snp_index.match_read().
    #[arg(long, default_value_t = 20)]
    pub snp_min_anchor: u32,

    /// Disable genome-based AlignedRead refinement even if --genome is supplied.
    #[arg(long, default_value_t = false)]
    pub no_genome_refine: bool,

    // ------------------------------
    // MatchOptions exposed to user
    // ------------------------------
    /// If true, require read blocks to be on a compatible strand.
    #[arg(long, default_value_t = false)]
    pub require_strand: bool,

    /// If true, require the read to have the exact same splice junction chain as the transcript.
    #[arg(long, default_value_t = false)]
    pub require_exact_junction_chain: bool,

    /// Maximum allowed 5′ overhang (bp). If exceeded -> OverhangTooLarge.
    #[arg(long, default_value_t = 100)]
    pub max_5p_overhang_bp: u32,

    /// Maximum allowed 3′ overhang (bp). If exceeded -> OverhangTooLarge.
    #[arg(long, default_value_t = 100)]
    pub max_3p_overhang_bp: u32,

    /// Allowed sequencing error gap. If exceeded -> JunctionMismatch.
    #[arg(long, default_value_t = 5)]
    pub allowed_intronic_gap_size: u32,
}

#[derive(Clone)]
struct Job {
    cell: u64,
    umi: u64,
    spliced: SplicedRead,

    /// Refined/cleaned read used for SNP matching.
    ///
    /// This is currently created when --genome or --vcf is supplied.
    /// If --genome is supplied and --no-genome-refine is not set, it is refined
    /// against the reference genome before being stored here.
    aligned: Option<AlignedRead>,
}

struct OptionalSnp {
    index: SnpIndex,
}

fn aux_tag_str<'a>(rec: &'a Record, tag: [u8; 2]) -> Option<&'a str> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(&tag).ok()? {
        Aux::String(s) => Some(s),
        _ => None,
    }
}

/// Many 10x barcodes look like "AAAC...-1". Strip the suffix before encoding.
fn normalize_10x_barcode(cb: &str) -> &str {
    match cb.split_once('-') {
        Some((core, _)) => core,
        None => cb,
    }
}

/// Convert a DNA string (ACGT only) into a packed u64 via 2-bit encoding.
fn dna_to_u64(seq: &str) -> Option<u64> {
    if !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    let tool = IntToStr::new(seq.as_bytes());
    Some(tool.into_u64())
}

/// Build SplicedRead + union span [start0, end0) from a BAM record.
///
/// This is still the canonical splice path.
/// The optional genome-refined AlignedRead is currently used for SNP matching.
/// Once snp-index exposes a stable "AlignedRead -> RefBlock/SplicedRead" conversion,
/// this function can be mirrored with a refined-read based version to make the
/// transcript matching directly benefit from the cleaned alignment representation.
fn record_to_spliced_read(rec: &Record, chr_id: usize) -> Option<(SplicedRead, u32, u32)> {
    let blocks: Vec<RefBlock> = record_to_blocks(rec);
    if blocks.is_empty() {
        return None;
    }

    let strand = if rec.is_reverse() {
        Strand::Minus
    } else {
        Strand::Plus
    };

    let mut spliced = SplicedRead::new(chr_id, strand, blocks);
    let (start0, end0) = spliced.finalize();
    Some((spliced, start0, end0))
}

fn main() -> Result<()> {
    let args = QuantCli::parse();

    if args.threads > 0 {
        // If build_global fails because Rayon was already initialized, ignore.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading index {}", args.index.display()))?;
    println!("{idx}");

    let chr_map = build_chr_map_fuzzy(&idx);

    let genome = match &args.genome {
        Some(path) => Some(
            Genome::from_fasta(path)
                .with_context(|| format!("reading genome FASTA {}", path.display()))?,
        ),
        None => None,
    };

    if args.vcf.is_some() && genome.is_none() {
        return Err(anyhow!("--vcf requires --genome"));
    }

    let mut reader = Reader::from_path(&args.bam)
        .with_context(|| format!("bam file could not be read: {}", args.bam.display()))?;
    let header = reader.header().clone();

    let chr_names: Vec<String> = (0..header.target_count())
        .map(|i| {
            std::str::from_utf8(header.tid2name(i))
                .unwrap()
                .to_string()
        })
        .collect();

    let chr_lengths: Vec<u32> = (0..header.target_count())
        .map(|i| header.target_len(i).unwrap_or(0) as u32)
        .collect();

    let vcf_opts = VcfReadOptions::default();
    let snp = match &args.vcf {
        Some(path) => {
            println!("Loading SNP index");
            let snp_index = SnpIndex::from_vcf_path(
                path,
                chr_names.clone(),
                chr_lengths.clone(),
                10_000,
                &vcf_opts,
            )
            .with_context(|| format!("reading SNP VCF {}", path.display()))?;

            eprintln!("{snp_index}");

            Some(OptionalSnp { index: snp_index })
        }
        None => None,
    };

    if genome.is_some() {
        eprintln!(
            "[INFO] genome refinement enabled: reads will be represented as AlignedRead and refined before SNP matching"
        );
        eprintln!(
            "[INFO] this cleanup layer can also be used to improve splice/transcript matching once SplicedRead construction from refined AlignedRead is enabled"
        );
    }

    if snp.is_some() {
        eprintln!("[INFO] SNP side-channel enabled: writing additional ref/alt matrices");
    }

    let match_opts = MatchOptions {
        require_strand: args.require_strand,
        require_exact_junction_chain: args.require_exact_junction_chain,
        max_5p_overhang_bp: args.max_5p_overhang_bp,
        max_3p_overhang_bp: args.max_3p_overhang_bp,
        allowed_intronic_gap_size: args.allowed_intronic_gap_size,
    };



    // ----------------------------
    // Stage 1: stream BAM -> jobs
    // ----------------------------
    let mut jobs: Vec<Job> = Vec::new();
    let mut n_seen: usize = 0;

    // ----------------------------
    // Merge partial results
    // ----------------------------
    let mut merged = Scdata::new(1, MatrixValueType::Real);
    let mut merged_intron = Scdata::new(1, MatrixValueType::Real);

    let mut merged_snp_ref = Scdata::new(1, MatrixValueType::Real);
    let mut merged_snp_alt = Scdata::new(1, MatrixValueType::Real);

    let mut merged_report = MappingInfo::new(
        None,
        args.min_mapq as f32,
        args.max_reads.unwrap_or(usize::MAX),
    );
    merged_report.start_counter();

    #[cfg(debug_assertions)]
    let mut i = 0usize;

    for r in reader.records() {
        #[cfg(debug_assertions)]
        {
            i += 1;
            if i - n_seen > 1_000_0000 {
                panic!(
                    "We have found more than 1 mio reads that did not pass initial filters.\n read {i}, processed {n_seen}, Wrong gene model?\n{merged_report}"
                );
            }
        }

        let rec = r.context("BAM read error")?;

        if rec.is_unmapped() {
            merged_report.report("unmapped");
            continue;
        }
        if rec.mapq() < args.min_mapq {
            merged_report.report("mapq failed");
            continue;
        }
        if rec.is_secondary() || rec.is_supplementary() {
            merged_report.report("secondary or supplementary");
            continue;
        }
        if args.read1_only && !rec.is_first_in_template() {
            merged_report.report("read!=1");
            continue;
        }

        let cb_raw = match aux_tag_str(&rec, *b"CB") {
            Some(v) => v,
            None => {
                merged_report.report("no CB tag");
                continue;
            }
        };
        let ub = match aux_tag_str(&rec, *b"UB") {
            Some(v) => v,
            None => {
                merged_report.report("no UB tag");
                continue;
            }
        };

        let cb = normalize_10x_barcode(cb_raw);

        let cell = match dna_to_u64(cb) {
            Some(v) => v,
            None => {
                merged_report.report("invalid CB");
                continue;
            }
        };
        let umi = match dna_to_u64(ub) {
            Some(v) => v,
            None => {
                merged_report.report("invalid UB");
                continue;
            }
        };

        let tid = rec.tid();
        if tid < 0 {
            merged_report.report("tid below 0");
            continue;
        }

        let chr_name = std::str::from_utf8(header.tid2name(tid as u32))
            .context("Invalid chromosome name in BAM header")?;

        let chr_id = match chr_map.get(chr_name) {
            Some(&id) => id,
            None => {
                merged_report.report("contig not in index");
                continue;
            }
        };

        // Optional genome/SNP read representation.
        // If --genome was supplied, this becomes the cleaned/refined read representation.
        // If only --vcf were supplied, we reject above because SNP matching needs genome.
        let aligned = if genome.is_some() || snp.is_some() {
            let mut read = AlignedRead::from_record(&rec, chr_id);

            if let Some(genome) = &genome {
                if !args.no_genome_refine {
                    read.refine_against_genome(genome, RefineOptions::default());
                }
            }

            Some(read)
        } else {
            None
        };

        let (spliced, _start0, _end0) = match record_to_spliced_read(&rec, chr_id) {
            Some(v) => v,
            None => {
                merged_report.report("no ref blocks");
                continue;
            }
        };

        jobs.push(Job {
            cell,
            umi,
            spliced,
            aligned,
        });

        n_seen += 1;
        if let Some(maxr) = args.max_reads {
            if n_seen >= maxr {
                break;
            }
        }

        if jobs.len() >= CHUNK {
            println!("Processing chunk of size {}", CHUNK);
            merged_report.stop_file_io_time();

            match args.quant_mode {
                QuantMode::Gene => process_chunk_gene(
                    &jobs,
                    &idx,
                    snp.as_ref(),
                    match_opts,
                    &mut merged,
                    &mut merged_intron,
                    &mut merged_snp_ref,
                    &mut merged_snp_alt,
                    &mut merged_report,
                    args.min_mapq,
                    args.snp_min_anchor,
                )?,
                QuantMode::Transcript => process_chunk_transcript(
                    &jobs,
                    &idx,
                    snp.as_ref(),
                    match_opts,
                    &mut merged,
                    &mut merged_intron,
                    &mut merged_snp_ref,
                    &mut merged_snp_alt,
                    &mut merged_report,
                    args.min_mapq,
                    args.snp_min_anchor,
                )?,
            }

            jobs.clear();
        }
    }

    // ----------------------------------------
    // Stage 2: final chunk
    // ----------------------------------------
    if !jobs.is_empty() {
        println!("Processing final chunk of size {}", jobs.len());
        merged_report.stop_file_io_time();

        match args.quant_mode {
            QuantMode::Gene => process_chunk_gene(
                &jobs,
                &idx,
                snp.as_ref(),
                match_opts,
                &mut merged,
                &mut merged_intron,
                &mut merged_snp_ref,
                &mut merged_snp_alt,
                &mut merged_report,
                args.min_mapq,
                args.snp_min_anchor,
            )?,
            QuantMode::Transcript => process_chunk_transcript(
                &jobs,
                &idx,
                snp.as_ref(),
                match_opts,
                &mut merged,
                &mut merged_intron,
                &mut merged_snp_ref,
                &mut merged_snp_alt,
                &mut merged_report,
                args.min_mapq,
                args.snp_min_anchor,
            )?,
        }

        jobs.clear();
    }

    println!("Writing outfiles");

    match args.quant_mode {
        QuantMode::Gene => {
            let features = GeneFeatureIndex::new(&idx);

            if args.split_intronic {
                merged.finalize_for_export(args.min_cell_counts, &features);

                let pass: std::collections::HashSet<u64> =
                    merged.export_cell_ids().iter().copied().collect();
                merged_intron.finalize_for_cells(&pass, &features);

                merged_report.stop_multi_processor_time();

                println!("Writing matrix files");
                merged
                    .write_sparse(&args.outpath, &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing gene matrix")?;

                println!("Writing intronic matrix files");
                merged_intron
                    .write_sparse(&add_suffix(&args.outpath, "_intronic"), &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing intronic gene matrix")?;
            } else {
                merged.merge(&merged_intron);
                merged.finalize_for_export(args.min_cell_counts, &features);

                merged_report.stop_multi_processor_time();

                println!("Writing matrix files");
                merged
                    .write_sparse(&args.outpath, &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing gene matrix")?;
            }
        }

        QuantMode::Transcript => {
            let features = TranscriptFeatureIndex::new(&idx);

            if args.split_intronic {
                merged.finalize_for_export(args.min_cell_counts, &features);

                let pass: std::collections::HashSet<u64> =
                    merged.export_cell_ids().iter().copied().collect();
                merged_intron.finalize_for_cells(&pass, &features);

                merged_report.stop_multi_processor_time();

                println!("Writing matrix files");
                merged
                    .write_sparse(&args.outpath, &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing transcript matrix")?;

                println!("Writing intronic matrix files");
                merged_intron
                    .write_sparse(&add_suffix(&args.outpath, "_intronic"), &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing intronic transcript matrix")?;
            } else {
                merged.merge(&merged_intron);
                merged.finalize_for_export(args.min_cell_counts, &features);

                merged_report.stop_multi_processor_time();

                println!("Writing matrix files");
                merged
                    .write_sparse(&args.outpath, &features)
                    .map_err(|e| anyhow!(e))
                    .context("writing transcript matrix")?;
            }
        }
    }

    if let Some(snp) = &snp {
        println!("Writing SNP ref/alt matrix files");

        // Keep SNP matrices aligned to the same cells as the main matrix.
        // This uses the main exported cells as whitelist.
        let pass: std::collections::HashSet<u64> =
            merged.export_cell_ids().iter().copied().collect();




        if merged_snp_ref.is_empty(){
            eprintln!("No SNPs collected");
        }else {


            merged_snp_alt.finalize_for_cells(&pass, &snp.index);

            merged_snp_ref.retain_features( &merged_snp_alt.observed_feature_ids() );

            merged_snp_ref.finalize_for_cells(&pass, &snp.index); 

            merged_snp_alt
                .write_sparse(&add_suffix(&args.outpath, "_alt"), &snp.index)
                .map_err(|e| anyhow!(e))
                .context("writing SNP alt matrix")?;    

            merged_snp_ref
                .write_sparse(&add_suffix(&args.outpath, "_ref"), &snp.index)
                .map_err(|e| anyhow!(e))
                .context("writing SNP ref matrix")?;
        }
    }

    merged_report.stop_file_io_time();

    println!("{merged_report}");

    let log_str = format!("{merged_report}");
    let log_path = args.outpath.with_extension("log");
    let mut file = File::create(&log_path)
        .with_context(|| format!("failed to create log file {}", log_path.display()))?;
    file.write_all(log_str.as_bytes())
        .with_context(|| format!("failed to write log file {}", log_path.display()))?;

    Ok(())
}

fn add_suffix(path: &PathBuf, suffix: &str) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| std::path::Path::new(""));
    let stem = path.file_stem().unwrap().to_string_lossy();
    let ext = path.extension().map(|e| e.to_string_lossy());

    let new_name = match ext {
        Some(e) => format!("{stem}{suffix}.{e}"),
        None => format!("{stem}{suffix}"),
    };

    parent.join(new_name)
}

#[allow(clippy::too_many_arguments)]
fn process_chunk_gene(
    jobs: &[Job],
    idx: &SpliceIndex,
    snp: Option<&OptionalSnp>,
    match_opts: MatchOptions,
    merged: &mut Scdata,
    merged_intron: &mut Scdata,
    merged_snp_ref: &mut Scdata,
    merged_snp_alt: &mut Scdata,
    merged_report: &mut MappingInfo,
    min_mapq: u8,
    snp_min_anchor: u32,
) -> Result<()> {
    let threads = rayon::current_num_threads().max(1);
    let chunk_size = (jobs.len() / threads).max(10_000);

    let partials: Vec<(Scdata, Scdata, Scdata, Scdata, MappingInfo)> = jobs
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut sc = Scdata::new(1, MatrixValueType::Real);
            let mut sc_intron = Scdata::new(1, MatrixValueType::Real);
            let mut sc_snp_ref = Scdata::new(1, MatrixValueType::Real);
            let mut sc_snp_alt = Scdata::new(1, MatrixValueType::Real);

            let mut rep = MappingInfo::new(None, min_mapq as f32, usize::MAX);

            for job in chunk {
                let gene_hits = idx.match_genes(&job.spliced, match_opts);
                if gene_hits.is_empty() {
                    rep.report("no hit");
                } else {
                    let g = &gene_hits[0];
                    rep.report(g.best_hit.class.to_string());

                    let gid = g.gene_id;

                    if g.best_hit.class == MatchClass::Intronic {
                        sc_intron.try_insert(
                            &job.cell,
                            GeneUmiHash(gid as u64, job.umi),
                            1.0,
                            &mut rep,
                        );
                    } else {
                        sc.try_insert(
                            &job.cell,
                            GeneUmiHash(gid as u64, job.umi),
                            1.0,
                            &mut rep,
                        );
                    }
                }

                add_snp_hits(
                    snp,
                    &job.aligned,
                    job.cell,
                    job.umi,
                    snp_min_anchor,
                    &mut sc_snp_ref,
                    &mut sc_snp_alt,
                    &mut rep,
                );
            }

            (sc, sc_intron, sc_snp_ref, sc_snp_alt, rep)
        })
        .collect();

    merged_report.stop_multi_processor_time();

    for (sc, sc_intron, sc_snp_ref, sc_snp_alt, rep) in partials {
        merged.merge(&sc);
        merged_intron.merge(&sc_intron);
        merged_snp_ref.merge(&sc_snp_ref);
        merged_snp_alt.merge(&sc_snp_alt);
        merged_report.merge(&rep);
    }

    merged_report.stop_single_processor_time();

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_chunk_transcript(
    jobs: &[Job],
    idx: &SpliceIndex,
    snp: Option<&OptionalSnp>,
    match_opts: MatchOptions,
    merged: &mut Scdata,
    merged_intron: &mut Scdata,
    merged_snp_ref: &mut Scdata,
    merged_snp_alt: &mut Scdata,
    merged_report: &mut MappingInfo,
    min_mapq: u8,
    snp_min_anchor: u32,
) -> Result<()> {
    let threads = rayon::current_num_threads().max(1);
    let chunk_size = (jobs.len() / threads).max(10_000);

    let partials: Vec<(Scdata, Scdata, Scdata, Scdata, MappingInfo)> = jobs
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut sc = Scdata::new(1, MatrixValueType::Real);
            let mut sc_intron = Scdata::new(1, MatrixValueType::Real);
            let mut sc_snp_ref = Scdata::new(1, MatrixValueType::Real);
            let mut sc_snp_alt = Scdata::new(1, MatrixValueType::Real);

            let mut rep = MappingInfo::new(None, min_mapq as f32, usize::MAX);

            for job in chunk {
                let tx_hits = idx.match_transcripts(&job.spliced, match_opts);
                if tx_hits.is_empty() {
                    rep.report("no hit");
                } else {
                    let h = &tx_hits[0];
                    rep.report(h.hit.class.to_string());

                    let tid = h.transcript_id;

                    if h.hit.class == MatchClass::Intronic {
                        sc_intron.try_insert(
                            &job.cell,
                            GeneUmiHash(tid as u64, job.umi),
                            1.0,
                            &mut rep,
                        );
                    } else {
                        sc.try_insert(
                            &job.cell,
                            GeneUmiHash(tid as u64, job.umi),
                            1.0,
                            &mut rep,
                        );
                    }
                }

                add_snp_hits(
                    snp,
                    &job.aligned,
                    job.cell,
                    job.umi,
                    snp_min_anchor,
                    &mut sc_snp_ref,
                    &mut sc_snp_alt,
                    &mut rep,
                );
            }

            (sc, sc_intron, sc_snp_ref, sc_snp_alt, rep)
        })
        .collect();

    merged_report.stop_multi_processor_time();

    for (sc, sc_intron, sc_snp_ref, sc_snp_alt, rep) in partials {
        merged.merge(&sc);
        merged_intron.merge(&sc_intron);
        merged_snp_ref.merge(&sc_snp_ref);
        merged_snp_alt.merge(&sc_snp_alt);
        merged_report.merge(&rep);
    }

    merged_report.stop_single_processor_time();

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn add_snp_hits(
    snp: Option<&OptionalSnp>,
    aligned: &Option<AlignedRead>,
    cell: u64,
    umi: u64,
    snp_min_anchor: u32,
    sc_snp_ref: &mut Scdata,
    sc_snp_alt: &mut Scdata,
    rep: &mut MappingInfo,
) {
    let Some(snp) = snp else {
        return;
    };

    let Some(read) = aligned else {
        rep.report("snp requested but no aligned read");
        return;
    };

    let hits = snp.index.match_read(read, snp_min_anchor.try_into().unwrap());

    if hits.ref_ids.is_empty() && hits.alt_ids.is_empty() {
        rep.report("no snp hit");
        return;
    }

    for snp_id in hits.ref_ids {
        sc_snp_ref.try_insert(
            &cell,
            GeneUmiHash(snp_id as u64, umi),
            1.0,
            rep,
        );
    }

    for snp_id in hits.alt_ids {
        sc_snp_alt.try_insert(
            &cell,
            GeneUmiHash(snp_id as u64, umi),
            1.0,
            rep,
        );
    }

    rep.report("snp hit");
}

fn build_chr_map_fuzzy(idx: &SpliceIndex) -> std::collections::HashMap<String, usize> {
    let mut map = std::collections::HashMap::new();

    for (i, n) in idx.chr_names.iter().enumerate() {
        // exact
        map.entry(n.clone()).or_insert(i);

        // without chr prefix
        let no_chr = n.strip_prefix("chr").unwrap_or(n).to_string();
        map.entry(no_chr).or_insert(i);

        // with chr prefix
        let with_chr = if n.starts_with("chr") {
            n.clone()
        } else {
            format!("chr{n}")
        };
        map.entry(with_chr).or_insert(i);

        // mito aliases
        if n == "MT" {
            map.entry("chrM".to_string()).or_insert(i);
            map.entry("M".to_string()).or_insert(i);
        }
        if n == "chrM" {
            map.entry("MT".to_string()).or_insert(i);
            map.entry("M".to_string()).or_insert(i);
        }
    }

    map
}