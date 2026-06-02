use crate::ont_normalizer::cli::Cli;
crate::fastq::{record::FastqRecord, writer::FastqWriter};

use crate::read_tag_table::{
    ReadTagTableWriter,
    ReadTagWriteRecord,
};
use sc_primer::{Orientation, PrimerDetector, PrimerMatch};

use use crate::fastq::record::::FastqRecord;

use mapping_info::MappingInfo;

use anyhow::{Context, Result};
use rust_htslib::bam::{Read, Reader};
use std::fs::File;
use std::io::{BufWriter};
use std::path::PathBuf;
use rayon::prelude::*;
use std::io::Write;
use scdata::Scdata;
use scdata::cell_data::GeneUmiHash;


const CHUNK_SIZE :usize = 1_000_000;

#[derive(Debug, Clone)]
struct NormalizedPrimerRecord {
    insert: FastqRecord,
    original_read_id: String,
    orientation: &'static str,
    raw_cb: String,
    quality_cb: String,
    raw_umi: String,
    quality_umi: String,
}

struct ChunkPartial {
    fastq_records: Vec<FastqRecord>,
    read_tags: ReadTagTable,
    tag_counts: Scdata,
    stats: MappingInfo,
}

#[derive(Debug, Clone)]
pub struct OntNormalizerConfig {
    pub bam: PathBuf,
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub min_transcript_len: usize,
    pub primer: sc_primer::PrimerDetector,
    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: MappingInfo,
    feature_tag_table: Scdata,
}

impl OntNormalizer {
    pub fn new(config: OntNormalizerConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
            feature_tag_table: Scdata::new(1, MatrixValueType::Real),
        }
    }

    pub fn write_feature_tag_table_if_present(&mut self) -> Result<()> {
        if self.feature_tag_table.is_empty() {
            return Ok(());
        }

        let mapper = match self.config.feature_tag_mapper.as_ref() {
            Some(mapper) => mapper,
            None => return Ok(()),
        };

        let out_dir = self
            .config
            .out
            .with_extension("")
            .join("feature_tag_table_unfiltered");

        std::fs::create_dir_all(&out_dir)?;

        self.feature_tag_table
            .write_sparse(&out_dir, mapper)
            .map_err(|e| anyhow::anyhow!("writing feature tag table failed: {e}"))?;

        Ok(())
    }

    pub fn from_cli(cli: Cli) -> Self {
        Ok(
            Self::new(OntNormalizerConfig {
                bam: cli.bam,
                out: cli.out,
                read_tags: cli.read_tags,
                min_transcript_len: cli.min_transcript_len,
                primer: cli.primer.detector()?,
                tags: cli.tags.mapper()?,
                threads: cli.threads,
                gzip_level: cli.gzip_level,
                gzip: !cli.no_gzip,
            })
        )
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &OntNormalizerConfig {
        &self.config
    }

    pub fn run(&mut self) -> Result<()> {
        let mut bam = Reader::from_path(&self.config.bam)
            .with_context(|| format!("failed to open BAM: {}", self.config.bam.display()))?;

        if self.config.threads > 1 {
            bam.set_threads(self.config.threads.saturating_sub(1).max(1))
                .context("failed to set BAM reader threads")?;
        }

        let mut fastq =
            FastqWriter::new(&self.config.out, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!("failed to create FASTQ: {}", self.config.out.display())
                })?;

        let mut tags = self.create_tag_writer()?;

        let detector = PrimerDetector::tenx_v3();

        let mut chunk: Vec<FastqRecord> = Vec::with_capacity(CHUNK_SIZE);

        for rec_result in bam.records() {
            let rec = rec_result?;
            chunk.push(FastqRecord::from_bam_record(&rec));

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &detector, &mut fastq, &mut tags).expect("Chunk was not prcessed correctly");
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &detector, &mut fastq, &mut tags).expect("Chunk was not prcessed correctly");
        }

        fastq.finish()?;
        tags.flush()?;

        Ok(())
    }

    fn process_chunk<W: Write>(
        &mut self,
        input: &[FastqRecord],
        detector: &PrimerDetector,
        tag_mapper: Option<&FastTagMapper>,
        fastq: &mut FastqWriter,
        tags: &mut ReadTagTableWriter<W>,
    ) -> Result<()> {
        let partials: Vec<ChunkPartial> = input
            .par_iter()
            .map(|read| {
                let mut out = ChunkPartial {
                    molecules: Vec::new(),
                    tag_counts: Scdata::new(1, MatrixValueType::Real),
                    stats: MappingInfo::new(None, 0.0, 0),
                };

                let matches = match detector.detect_all(&read.seq, &read.qual) {
                    Ok(matches) => matches,
                    Err(_) => {
                        out.stats.report("no_primer_match");
                        return out;
                    }
                };

                match matches.len() {
                    0 => out.stats.report("no_primer_match"),
                    1 => out.stats.report("one_primer_match"),
                    _ => out.stats.report("multi_primer_match"),
                }

                for (match_index, primer_match) in matches.iter().enumerate() {
                    let insert = match primer_match.get_insert(&read.seq, &read.qual) {
                        Ok(x) => x,
                        Err(_) => {
                            out.stats.report("bad_insert_slice");
                            continue;
                        }
                    };

                    if insert.seq.len() < self.config.min_insert_len {
                        out.stats.report("short_insert");
                        continue;
                    }

                    let cell = match primer_match.get_cell(&read.seq, &read.qual) {
                        Ok(x) => x,
                        Err(_) => {
                            out.stats.report("bad_cell_slice");
                            continue;
                        }
                    };

                    let umi = match primer_match.get_umi(&read.seq, &read.qual) {
                        Ok(x) => x,
                        Err(_) => {
                            out.stats.report("bad_umi_slice");
                            continue;
                        }
                    };

                    let cell_id = encode_cell(cell.seq);
                    let umi_id = encode_umi(umi.seq);

                    let orientation = match primer_match.orientation {
                        Orientation::Forward => "forward",
                        Orientation::ReverseComplement => "reverse_complement",
                    };

                    let insert_id = format!("{}_{}", read.id, match_index);

                    let mut insert_record = FastqRecord {
                        id: insert_id,
                        desc: None,
                        seq: insert.seq.to_vec(),
                        qual: insert.qual.to_vec(),
                    };

                    if primer_match.orientation == Orientation::ReverseComplement {
                        insert_record = insert_record.revcomp();
                    }

                    if let Some(mapper) = tag_mapper {
                        let tag_call = mapper.call(&insert_record.seq);

                        if let Some(tag_id) = tag_call.best_tag_id() {
                            let feature_umi = GeneUmiHash(tag_id, umi_id);

                            out.tag_counts.try_insert(
                                GeneUmiHash {&cell_id,feature_umi},
                                1.0,
                                &mut out.stats,
                            );

                            out.stats.report("feature_tag_match");
                            continue;
                        } else {
                            out.stats.report(tag_call.status());
                        }
                    }
                    let read_tag_record = ReadTagRecord {
                        read_id: insert_record.id.clone(),
                        original_read_id: Some(read.id.clone()),
                        cell: String::from_utf8_lossy(cell.seq).to_string(),
                        cell_qual: Some(String::from_utf8_lossy(cell.qual).to_string()),
                        umi: String::from_utf8_lossy(umi.seq).to_string(),
                        umi_qual: Some(String::from_utf8_lossy(umi.qual).to_string()),
                    };

                    out.read_tags.insert(read_tag_record);
                    out.fastq_records.push(insert_record);
                }

                out
            })
            .collect();
        
        for partial in partials {
            self.stats.merge(&partial.stats);
            self.read_tags.merge(partial.read_tags);
            self.feature_tag_counts.merge(&partial.tag_counts);

            for record in partial.fastq_records {
                fastq.write(&record)?;
            }
        }
        Ok(())
    }

    pub fn stats_report(&self) -> String {
        let info = &self.stats;

        let total = info.get_issue_count("total_records");
        let zero = info.get_issue_count("zero_cassette");
        let one = info.get_issue_count("one_cassette");
        let multi = info.get_issue_count("multi_cassette");

        let emitted = info.get_issue_count("emitted_molecules");
        let written = info.get_issue_count("fastq_reads_written");

        let fwd = info.get_issue_count("forward_molecules");
        let rev = info.get_issue_count("reverse_molecules");

        let failed_poly_t = info.get_issue_count("failed_poly_t");
        let too_short = info.get_issue_count("too_short_after_adapter");

        let with_cassette = one + multi;

        let pct = |n: usize, d: usize| {
            if d == 0 {
                0.0
            } else {
                (n as f64 / d as f64) * 100.0
            }
        };

        let mean = if total == 0 {
            0.0
        } else {
            emitted as f64 / total as f64
        };

        format!(
            r#"bam-ont-normalizer summary
============================

Input
-----
reads processed        : {total}
reads w/o cassette     : {zero} ({zero_pct:.2}%)
reads w/ cassette      : {with_cassette} ({cassette_pct:.2}%)
  ├─ one cassette      : {one}
  └─ multi cassette    : {multi} ({multi_pct:.2}%)

Output
------
FASTQ reads written    : {written}
molecules emitted      : {emitted}
mean molecules/read    : {mean:.3}

Orientation
-----------
forward                : {fwd} ({fwd_pct:.2}%)
reverse                : {rev} ({rev_pct:.2}%)

Rejected
--------
failed polyT           : {failed_poly_t}
too short after adapter: {too_short}
"#,
            total = total,
            zero = zero,
            zero_pct = pct(zero, total),
            with_cassette = with_cassette,
            cassette_pct = pct(with_cassette, total),
            one = one,
            multi = multi,
            multi_pct = pct(multi, total),
            written = written,
            emitted = emitted,
            mean = mean,
            fwd = fwd,
            rev = rev,
            fwd_pct = pct(fwd, emitted),
            rev_pct = pct(rev, emitted),
            failed_poly_t = failed_poly_t,
            too_short = too_short,
        )
    }

    fn create_tag_writer(&self) -> Result<ReadTagTableWriter<BufWriter<File>>> {
        let file = File::create(&self.config.tags).with_context(|| {
            format!(
                "failed to create tags TSV: {}",
                self.config.tags.display()
            )
        })?;

        let writer = BufWriter::new(file);

        ReadTagTableWriter::new(writer)
    }


}
