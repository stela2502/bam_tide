use crate::illumina_normalizer::cli::{Chemistry, Cli};
use crate::illumina_normalizer::fastq_pair_reader::FastqPairReader;
use crate::read_tag_table::{ReadTagTableWriter, ReadTagWriteRecord};
use use crate::fastq::record::::FastqRecord;
use use crate::fastq::writer::::FastqWriter;
use crate::primer::rhapsody_cellids::RhapsodyCellIds;

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

const CHUNK_SIZE: usize = 1_000_000;

#[derive(Debug, Clone)]
pub struct IlluminaNormalizerConfig {
    pub r1: PathBuf,
    pub r2: PathBuf,
    pub out_r1: PathBuf,
    pub out_r2: PathBuf,
    pub tags: PathBuf,
    pub chemistry: Chemistry,
    pub cb_len: usize,
    pub umi_len: usize,
    pub poly_t_min: usize,
    pub sample_species: String,
    pub sample_kmer_size: usize,
    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

pub struct IlluminaNormalizer {
    config: IlluminaNormalizerConfig,
    stats: MappingInfo,
}

#[derive(Debug)]
struct NormalizedPair {
    artificial_r1: FastqRecord,
    r2: FastqRecord,
    raw_cb: String,
    quality_cb: String,
    raw_umi: String,
    quality_umi: String,
    original_r1_id: String,
    original_r2_id: String,
}

impl IlluminaNormalizer {
    pub fn new(config: IlluminaNormalizerConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn from_cli(cli: Cli) -> Self {
        let cb_len = cli.cb_len.unwrap_or_else(|| cli.chemistry.default_cb_len());
        let umi_len = cli.umi_len.unwrap_or_else(|| cli.chemistry.default_umi_len());

        Self::new(IlluminaNormalizerConfig {
            r1: cli.r1,
            r2: cli.r2,
            out_r1: cli.out_r1,
            out_r2: cli.out_r2,
            tags: cli.tags,
            chemistry: cli.chemistry,
            cb_len,
            umi_len,
            poly_t_min: cli.poly_t_min,
            sample_species: cli.sample_species,
            sample_kmer_size: cli.sample_kmer_size,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        })
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &IlluminaNormalizerConfig {
        &self.config
    }

    pub fn run(&mut self) -> Result<()> {
        if self.config.threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.config.threads)
                .build_global()
                .ok();
        }

        let mut reader = FastqPairReader::from_paths(&self.config.r1, &self.config.r2)?;

        let mut out_r1 = FastqWriter::new(
            &self.config.out_r1,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create FASTQ: {}", self.config.out_r1.display()))?;

        let mut out_r2 = FastqWriter::new(
            &self.config.out_r2,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create FASTQ: {}", self.config.out_r2.display()))?;

        let mut tags = self.create_tag_writer()?;
        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK_SIZE);

        while let Some(pair) = reader.next_pair()? {
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &mut out_r1, &mut out_r2, &mut tags)?;
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut out_r1, &mut out_r2, &mut tags)?;
        }

        out_r1.finish()?;
        out_r2.finish()?;
        tags.flush()?;

        Ok(())
    }

    fn process_chunk<W: Write>(
        &mut self,
        input: &[(FastqRecord, FastqRecord)],
        out_r1: &mut FastqWriter,
        out_r2: &mut FastqWriter,
        tags: &mut ReadTagTableWriter<W>,
    ) -> Result<()> {
        let config = self.config.clone();

        let partials: Vec<(Option<NormalizedPair>, MappingInfo)> = input
            .par_iter()
            .map(|(r1, r2)| {
                let mut local_stats = MappingInfo::new(None, 0.0, 0);
                local_stats.report("total_pairs");

                let result = match config.chemistry {
                    Chemistry::TenxV2 | Chemistry::TenxV3 | Chemistry::TenxV4 => {
                        normalize_10x_pair(r1, r2, &config, &mut local_stats)
                    }
                    Chemistry::BdV1 | Chemistry::BdV2_96 | Chemistry::BdV2_384 => {
                        normalize_bd_pair(r1, r2, &config, &mut local_stats)
                    }
                };

                match result {
                    Ok(pair) => (Some(pair), local_stats),
                    Err(_) => {
                        local_stats.report("failed_pairs");
                        (None, local_stats)
                    }
                }
            })
            .collect();

        for (pair, local_stats) in partials {
            self.stats.merge(&local_stats);

            let Some(pair) = pair else {
                continue;
            };

            out_r1.write(&pair.artificial_r1)?;
            out_r2.write(&pair.r2)?;

            let rec = ReadTagWriteRecord {
                read_id: &pair.artificial_r1.id,
                original_read_id: Some(&pair.original_r1_id),
                orientation: Some(self.config.chemistry.as_str()),
                raw_cb: pair.raw_cb,
                quality_cb: pair.quality_cb,
                raw_umi: pair.raw_umi,
                quality_umi: pair.quality_umi,
                status: "written",
            };

            tags.write_record(&rec)?;
            self.stats.report("fastq_pairs_written");
        }

        Ok(())
    }

    pub fn stats_report(&self) -> String {
        let info = &self.stats;

        let total = info.get_issue_count("total_pairs");
        let written = info.get_issue_count("fastq_pairs_written");
        let failed = info.get_issue_count("failed_pairs");
        let failed_bd = info.get_issue_count("failed_bd_cell_id");
        let failed_tenx = info.get_issue_count("failed_10x_barcode");
        let failed_polyt = info.get_issue_count("failed_poly_t");

        let pct = |n: usize, d: usize| {
            if d == 0 { 0.0 } else { (n as f64 / d as f64) * 100.0 }
        };

        format!(
            r#"illumina-normalizer summary
===========================

Input
-----
FASTQ pairs processed   : {total}
chemistry               : {chemistry}

Output
------
FASTQ pairs written     : {written} ({written_pct:.2}%)
failed pairs            : {failed} ({failed_pct:.2}%)

Rejected
--------
failed BD cell id       : {failed_bd}
failed 10x barcode      : {failed_tenx}
failed polyT            : {failed_polyt}
"#,
            total = total,
            chemistry = self.config.chemistry.as_str(),
            written = written,
            written_pct = pct(written, total),
            failed = failed,
            failed_pct = pct(failed, total),
            failed_bd = failed_bd,
            failed_tenx = failed_tenx,
            failed_polyt = failed_polyt,
        )
    }

    fn create_tag_writer(&self) -> Result<ReadTagTableWriter<BufWriter<File>>> {
        let file = File::create(&self.config.tags).with_context(|| {
            format!("failed to create tags TSV: {}", self.config.tags.display())
        })?;

        ReadTagTableWriter::new(BufWriter::new(file))
    }
}

fn normalize_10x_pair(
    r1: &FastqRecord,
    r2: &FastqRecord,
    config: &IlluminaNormalizerConfig,
    stats: &mut MappingInfo,
) -> Result<NormalizedPair> {
    let needed = config.cb_len + config.umi_len;

    if r1.len() < needed {
        stats.report("failed_10x_barcode");
        bail!("R1 too short for 10x barcode + UMI");
    }

    let cb_seq = &r1.seq[..config.cb_len];
    let cb_qual = &r1.qual[..config.cb_len];
    let umi_seq = &r1.seq[config.cb_len..needed];
    let umi_qual = &r1.qual[config.cb_len..needed];

    if has_n(cb_seq) || has_n(umi_seq) {
        stats.report("failed_10x_barcode");
        bail!("10x barcode or UMI contains N");
    }

    let artificial_id = normalized_id(&r1.id, 0);
    let artificial_seq = [cb_seq, umi_seq].concat();
    let artificial_qual = [cb_qual, umi_qual].concat();

    stats.report("accepted_10x_pairs");

    Ok(NormalizedPair {
        artificial_r1: FastqRecord::new(artificial_id.clone(), &artificial_seq, &artificial_qual),
        r2: rename_fastq_record(r2, artificial_id),
        raw_cb: String::from_utf8_lossy(cb_seq).to_string(),
        quality_cb: ascii_qual(cb_qual),
        raw_umi: String::from_utf8_lossy(umi_seq).to_string(),
        quality_umi: ascii_qual(umi_qual),
        original_r1_id: r1.id.clone(),
        original_r2_id: r2.id.clone(),
    })
}

fn normalize_bd_pair(
    r1: &FastqRecord,
    r2: &FastqRecord,
    config: &IlluminaNormalizerConfig,
    stats: &mut MappingInfo,
) -> Result<NormalizedPair> {
    let version = config
        .chemistry
        .rhapsody_version()
        .expect("BD chemistry must provide a Rhapsody version");

    let cellids = RhapsodyCellIds::new(version);
    let hit = match cellids.to_cellid_from_slice(&r1.seq) {
        Ok(hit) => hit,
        Err(err) => {
            stats.report("failed_bd_cell_id");
            bail!("failed BD cell id detection: {err:?}");
        }
    };

    if hit.umi_seq.len() != config.umi_len {
        stats.report("failed_bd_cell_id");
        bail!("unexpected BD UMI length: {}", hit.umi_seq.len());
    }

    if !has_polyt_after(&r1.seq, hit.umi_end, config.poly_t_min) {
        stats.report("failed_poly_t");
        bail!("failed polyT after BD UMI");
    }

    let cb_qual = r1.qual[hit.cell_start..hit.cell_end].to_vec();
    let umi_qual = r1.qual[hit.umi_start..hit.umi_end].to_vec();

    let artificial_id = normalized_id(&r1.id, 0);
    let artificial_seq = [hit.cell_seq.as_slice(), hit.umi_seq.as_slice()].concat();
    let artificial_qual = [cb_qual.as_slice(), umi_qual.as_slice()].concat();

    stats.report("accepted_bd_pairs");

    Ok(NormalizedPair {
        artificial_r1: FastqRecord::new(artificial_id.clone(), &artificial_seq, &artificial_qual),
        r2: rename_fastq_record(r2, artificial_id),
        raw_cb: String::from_utf8_lossy(&hit.cell_seq).to_string(),
        quality_cb: ascii_qual(&cb_qual),
        raw_umi: String::from_utf8_lossy(&hit.umi_seq).to_string(),
        quality_umi: ascii_qual(&umi_qual),
        original_r1_id: r1.id.clone(),
        original_r2_id: r2.id.clone(),
    })
}

fn normalized_id(original: &str, index: usize) -> String {
    format!("{original}/mol{index}")
}

fn rename_fastq_record(rec: &FastqRecord, id: String) -> FastqRecord {
    FastqRecord::new(id, &rec.seq, &rec.qual)
}

fn ascii_qual(qual: &[u8]) -> String {
    qual.iter().map(|q| (q + 33) as char).collect()
}

fn has_n(seq: &[u8]) -> bool {
    seq.iter().any(|b| b.to_ascii_uppercase() == b'N')
}

fn has_polyt_after(seq: &[u8], start: usize, min_t: usize) -> bool {
    if start >= seq.len() {
        return false;
    }

    let mut t_count = 0usize;
    for b in &seq[start..] {
        let b = b.to_ascii_uppercase();
        if b == b'T' || b == b'N' {
            t_count += 1;
        } else {
            break;
        }
    }

    t_count >= min_t
}

