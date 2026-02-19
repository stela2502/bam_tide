// src/bin/bam-quant.rs
//
// Quantify a 10x-style BAM against a splice index (GTF-derived) into scdata.
//
// Pipeline:
// 1) Stream BAM -> collect jobs (cell_u64, umi_u64, spliced_read, span)
// 2) Rayon over chunks -> local Scdata + MappingInfo
// 3) Merge locals into one Scdata + MappingInfo
//
// Assumptions (adjust if your crate differs):
// - Tags: CB (cell barcode), UB (UMI). CB may end with "-1" etc.
// - record_to_blocks(&Record) exists and returns Vec<RefBlock> 0-based half-open.
// - SpliceIndex has:
//     - transcripts: Vec<Transcript>
//     - chrom_names or equivalent to build chr_name -> chr_id mapping
//     - candidates_for_span_union(chr_id, start0, end0) -> Vec<TranscriptId>
// - Transcript has: gene_id: usize and match_spliced_read(&SplicedRead, MatchOptions) -> MatchHit
// - Strand enum has Plus/Minus (or similar)
// - MappingInfo has fields ok_reads, pcr_duplicates, local_dup, and ideally a merge method.
//   If no merge method exists, see `merge_mapping_info_fallback()` below.

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader, Record};
use std::path::PathBuf;

use scdata::cell_data::GeneUmiHash;
use scdata::{IndexedGenes, MatrixValueType, Scdata};

use mapping_info::MappingInfo;

use int_to_str::int_to_str::IntToStr;

const CHUNK: usize = 2_000_000;

// ---- Your splice index / transcript-matching crate ----
// Adjust these paths to your actual crate/module names.
use gtf_splice_index::model::MatchHit;
use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex, SplicedRead, Strand, TranscriptId};

// ---- Your ref-block conversion ----
// Adjust these paths to where RefBlock + record_to_blocks live in bam_tide.
use bam_tide::core::ref_block::record_to_blocks;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "bam-quant",
    about = "Quantify 10x BAM against splice index into scdata"
)]
pub struct QuantCli {
    /// Input BAM (10x / CellRanger BAM)
    #[arg(long, short)]
    pub bam: std::path::PathBuf,

    /// Splice index path (built from GTF beforehand)
    #[arg(long, short)]
    pub index: std::path::PathBuf,

    /// Outpath for the 10x mmx formated outfiles
    #[arg(long, short)]
    pub outpath: std::path::PathBuf,

    /// Minimum MAPQ
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Use only read1 (recommended for 10x; reduces duplicate mate counting noise)
    #[arg(long, default_value_t = false)]
    pub read1_only: bool,

    /// Rayon thread count (0 = default)
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Max reads to process (debug/dev)
    #[arg(long)]
    pub max_reads: Option<usize>,

    /// Min read counts per reported cell (debug/dev)
    #[arg(long, default_value_t = 400)]
    pub min_cell_counts: usize,

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
    chr_id: usize,
    start0: u32,
    end0: u32,
    spliced: SplicedRead,
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
fn record_to_spliced_read(rec: &Record, chr_id: usize) -> Option<(SplicedRead, u32, u32)> {
    let mut blocks: Vec<RefBlock> = record_to_blocks(rec);
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

/// Choose the "best" transcript among candidates, using:
/// 1) higher MatchClass rank
/// 2) smaller total overhang (5p + 3p)
///
/// Only accept MatchClass ∈ {Compatible, ExactJunctionChain, Intronic}.
fn pick_best_transcript(
    idx: &SpliceIndex,
    candidates: &[TranscriptId],
    read: &SplicedRead,
    opts: MatchOptions,
) -> Option<(TranscriptId, MatchHit)> {
    let mut best: Option<(TranscriptId, MatchHit)> = None;

    for &tid in candidates {
        let tr = &idx.transcripts[tid];
        let hit = tr.match_spliced_read(read, opts);

        if hit.class.rank() < 2 {
            //
            continue;
        }

        match best {
            None => best = Some((tid, hit)),
            Some((_btid, bhit)) => {
                if hit.class > bhit.class {
                    best = Some((tid, hit));
                } else if hit.class == bhit.class {
                    let bo = bhit.overhang_5p_bp + bhit.overhang_3p_bp;
                    let ho = hit.overhang_5p_bp + hit.overhang_3p_bp;
                    if ho < bo {
                        best = Some((tid, hit));
                    }
                }
            }
        }
    }

    best
}

/// Build a chr_name -> chr_id lookup.
///
/// Adjust this if your SpliceIndex already provides a method, e.g.
/// `idx.chr_name_to_id(&str) -> Option<usize>`.
fn build_chr_map(idx: &SpliceIndex) -> std::collections::HashMap<String, usize> {
    // Common patterns:
    // - idx.chr_names: Vec<String>
    // - idx.chrom_names: Vec<String>
    //
    // Change the field name here to match your SpliceIndex.
    let names: &Vec<String> = &idx.chr_names;

    let mut map = std::collections::HashMap::with_capacity(names.len());
    for (i, n) in names.iter().enumerate() {
        map.insert(n.clone(), i);
    }
    map
}

fn main() -> Result<()> {
    let args = QuantCli::parse();

    if args.threads > 0 {
        // If build_global fails (already initialized), ignore.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading index {}", args.index.display()))?;
    println!("{idx}");

    //let chr_map = build_chr_map(&idx);
    let chr_map = build_chr_map_fuzzy(&idx);

    let match_opts = MatchOptions {
        require_strand: args.require_strand,
        require_exact_junction_chain: args.require_exact_junction_chain,
        max_5p_overhang_bp: args.max_5p_overhang_bp,
        max_3p_overhang_bp: args.max_3p_overhang_bp,
        allowed_intronic_gap_size: args.allowed_intronic_gap_size,
    };

    let mut reader = Reader::from_path(&args.bam)
        .with_context(|| format!("bam file could not be read: {}", args.bam.display()))?;
    let header = reader.header().clone();

    // ----------------------------
    // Stage 1: stream BAM -> jobs
    // ----------------------------
    let mut jobs: Vec<Job> = Vec::new();
    let mut n_seen: usize = 0;

    // ----------------------------
    //    merge partial results
    // ----------------------------
    let mut merged = Scdata::new(1, MatrixValueType::Real);
    let mut merged_intron = Scdata::new(1, MatrixValueType::Real);
    let mut merged_report = MappingInfo::new(
        None,
        args.min_mapq as f32,
        args.max_reads.unwrap_or(usize::MAX),
    );
    #[cfg(debug_assertions)]
    let mut i = 0;
    for r in reader.records() {
        #[cfg(debug_assertions)]
        {
            i += 1;
            if i - n_seen > 1_000_0000 {
                panic!(
                    "We have found more than 1 mio reads that did not pass initial filters.\n read {i}, processed {n_seen},  Wrong gene model?\n{merged_report}"
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
            merged_report.report("secondary or supplemantary");
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
            None => continue,
        };
        let umi = match dna_to_u64(ub) {
            Some(v) => v,
            None => continue,
        };

        let tid = rec.tid();
        if tid < 0 {
            merged_report.report("tid below 1");
            continue;
        }
        let chr_name = std::str::from_utf8(header.tid2name(tid as u32))
            .context("Invalid chromosome name in BAM header")?;

        let chr_id = match chr_map.get(chr_name) {
            Some(&id) => id,
            None => {
                merged_report.report("contig not in index");
                continue; // contig not in splice index
            }
        };

        let (spliced, start0, end0) = match record_to_spliced_read(&rec, chr_id) {
            Some(v) => v,
            None => continue,
        };

        jobs.push(Job {
            cell,
            umi,
            chr_id,
            start0,
            end0,
            spliced,
        });

        n_seen += 1;
        if let Some(maxr) = args.max_reads {
            if n_seen >= maxr {
                break;
            }
        }

        if jobs.len() >= CHUNK {
            println!("Processing chunk of size {}", CHUNK);
            process_chunk(
                &jobs,
                &idx,
                match_opts,
                &mut merged,
                &mut merged_intron,
                &mut merged_report,
                args.min_mapq,
            )?;
            jobs.clear();
        }
    }

    // ----------------------------------------
    // Stage 2: parallel chunks -> partial scdata
    // ----------------------------------------

    if jobs.len() > 0 {
        println!("Processing final chunk of size {}", jobs.len());
        process_chunk(
            &jobs,
            &idx,
            match_opts,
            &mut merged,
            &mut merged_intron,
            &mut merged_report,
            args.min_mapq,
        )?;
        jobs.clear();
    }

    println!("Writing outfiles");
    // TODO: persist merged scdata to disk in your preferred format.
    let genes = IndexedGenes::from_names(&idx.chr_names);
    let _ = merged.write_sparse(&args.outpath, &genes, args.min_cell_counts);
    let _ = merged_intron.write_sparse(
        &add_suffix(&args.outpath, "_intronic"),
        &genes,
        args.min_cell_counts,
    );
    println!("{merged_report}");

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

fn process_chunk(
    jobs: &[Job],
    idx: &SpliceIndex,
    match_opts: MatchOptions,
    merged: &mut Scdata,
    merged_intron: &mut Scdata,
    merged_report: &mut MappingInfo,
    min_mapq: u8,
) -> Result<()> {
    let threads = rayon::current_num_threads().max(1);
    let chunk_size = (jobs.len() / threads).max(10_000);

    let partials: Vec<(Scdata, Scdata, MappingInfo)> = jobs
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut sc = Scdata::new(1, MatrixValueType::Real);
            let mut sc_intron = Scdata::new(1, MatrixValueType::Real);
            let mut rep = MappingInfo::new(None, min_mapq as f32, usize::MAX);

            for job in chunk {
                let cands = idx.candidates_for_span_union(job.chr_id, job.start0, job.end0);
                if cands.is_empty() {
                    rep.report("no match");
                    continue;
                }

                let best = pick_best_transcript(idx, &cands, &job.spliced, match_opts);
                let (tid, hit) = match best {
                    Some(v) => v,
                    None => {
                        rep.report("no hit");
                        continue;
                    }
                };

                rep.report(hit.class.to_string());
                let gene_id = idx.transcripts[tid].gene_id;
                rep.report(gene_id.to_string());
                println!("Found gene_id {gene_id}");

                if hit.class == MatchClass::Intronic {
                    sc_intron.try_insert(&job.cell, GeneUmiHash(gene_id, job.umi), 1.0, &mut rep);
                } else {
                    sc.try_insert(&job.cell, GeneUmiHash(gene_id, job.umi), 1.0, &mut rep);
                }
            }

            (sc, sc_intron, rep)
        })
        .collect();

    for (sc, sc_intron, rep) in partials {
        merged.merge(&sc);
        merged_intron.merge(&sc_intron);
        merged_report.merge(&rep);
    }

    Ok(())
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
