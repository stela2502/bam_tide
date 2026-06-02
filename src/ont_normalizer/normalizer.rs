use crate::ont_normalizer::cli::Cli;
use use crate::fastq::writer::::FastqWriter;
use crate::read_tag_table::{
    ReadTagTableWriter,
    ReadTagWriteRecord,
};
use crate::primer::{PrimerDetector, PrimerSplit};
use use crate::fastq::record::::FastqRecord;

use mapping_info::MappingInfo;

use anyhow::{Context, Result};
use rust_htslib::bam::{Read, Reader};
use std::fs::File;
use std::io::{BufWriter};
use std::path::PathBuf;
use rayon::prelude::*;
use std::io::Write;

const CHUNK_SIZE :usize = 1_000_000;

#[derive(Debug, Clone)]
pub struct OntNormalizerConfig {
    pub bam: PathBuf,
    pub out: PathBuf,
    pub tags: PathBuf,
    pub adapter: Vec<u8>,
    pub min_adapter_match: usize,
    pub min_transcript_len: usize,
    pub max_adapter_mismatches: usize,
    pub cb_len: usize,
    pub umi_len: usize,
    pub poly_t_min: usize,
    pub poly_t_window: usize,
    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: MappingInfo,
}

impl OntNormalizer {
    pub fn new(config: OntNormalizerConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn from_cli(cli: Cli) -> Self {
        Self::new(OntNormalizerConfig {
            bam: cli.bam,
            out: cli.out,
            tags: cli.tags,
            adapter: cli.adapter.into_bytes(),
            min_adapter_match: cli.min_adapter_match,
            min_transcript_len: cli.min_transcript_len,
            max_adapter_mismatches: cli.max_adapter_mismatches,
            cb_len: cli.cb_len,
            umi_len: cli.umi_len,
            poly_t_min: cli.poly_t_min,
            poly_t_window: cli.poly_t_window,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        })
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
        fastq: &mut FastqWriter,
        tags: &mut ReadTagTableWriter<W>,
    ) -> Result<()>  {

        let partials: Vec<(Vec<PrimerSplit>, MappingInfo)> = input
            .par_iter()
            .map(|read| {
                let mut local_stats = MappingInfo::new(None, 0.0, 0);

                let splits = match detector.split_read(read, &mut local_stats) {
                    Ok(splits) => splits,
                    Err(_) => {
                        local_stats.report("no_primer_split");
                        Vec::new()
                    }
                };

                (splits, local_stats)
            })
            .collect();

        for (splits, local_stats) in partials {
            self.stats.merge(&local_stats);

            for split in splits {
                fastq.write(&split.insert)?;
                let rec = ReadTagWriteRecord {
                    read_id: &split.insert.id,
                    original_read_id: None,
                    orientation: None,
                    status:"written",
                    raw_cb: split
                        .cell_id
                        .as_ref()
                        .map(|r| String::from_utf8_lossy(&r.seq).to_string())
                        .unwrap_or_default(),
                    quality_cb: split
                        .cell_id
                        .as_ref()
                        .map(|r| r.qual_string() )
                        .unwrap_or_default(),
                    raw_umi: split
                        .umi
                        .as_ref()
                        .map(|r| String::from_utf8_lossy(&r.seq).to_string())
                        .unwrap_or_default(),
                    quality_umi: split
                        .umi
                        .as_ref()
                        .map(|r| r.qual_string() )
                        .unwrap_or_default(),
                };

                tags.write_record(&rec)?;
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
