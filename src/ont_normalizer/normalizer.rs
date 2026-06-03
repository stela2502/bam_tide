use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::ont_normalizer::cli::Cli;
use crate::read_tag_table::{ReadTagTable, ReadTagRecord, ReadTagTableWriter};
use crate::tags::FastTagMapper;
use crate::index::FastTagFeatureIndex;

use anyhow::{Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader};
use sc_primer::{Orientation, PrimerDetector};
use scdata::{MatrixValueType, Scdata};
use scdata::cell_data::GeneUmiHash;

use std::path::PathBuf;

const CHUNK_SIZE: usize = 1_000_000;

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

    pub primer: PrimerDetector,
    pub feature_tag_mapper: Option<FastTagMapper>,

    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

impl OntNormalizerConfig {
    fn process_read(&self, read: &FastqRecord) -> ChunkPartial {
        let detector = &self.primer;
        let tag_mapper = self.feature_tag_mapper.as_ref();
        let min_transcript_len = self.min_transcript_len;
        let mut out = ChunkPartial {
            fastq_records: Vec::new(),
            read_tags: ReadTagTable::new(),
            tag_counts: Scdata::new(1, MatrixValueType::Real),
            stats: MappingInfo::new(None, 0.0, 0),
        };

        out.stats.report("total_records");

        let matches = match detector.detect_all(&read.seq, &read.qual) {
            Ok(matches) => matches,
            Err(_) => {
                out.stats.report("zero_cassette");
                out.stats.report("no_primer_match");
                return out;
            }
        };

        match matches.len() {
            0 => {
                out.stats.report("zero_cassette");
                return out;
            }
            1 => out.stats.report("one_cassette"),
            _ => out.stats.report("multi_cassette"),
        }

        for (match_index, primer_match) in matches.iter().enumerate() {
            let insert = match primer_match.get_insert(&read.seq, &read.qual) {
                Ok(x) => x,
                Err(_) => {
                    out.stats.report("bad_insert_slice");
                    continue;
                }
            };

            if insert.seq.len() < min_transcript_len {
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

            let cell_id = encode_sequence_id(&cell.seq);
            let umi_id = encode_sequence_id(&umi.seq);

            let orientation = match primer_match.orientation {
                Orientation::Forward => {
                    out.stats.report("forward_molecules");
                    "forward"
                }
                Orientation::ReverseComplement => {
                    out.stats.report("reverse_molecules");
                    "reverse_complement"
                }
            };

            let insert_id = format!("{}/mol{}", read.id, match_index);

            let mut insert_record = FastqRecord {
                id: insert_id,
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

                    out.tag_counts
                        .try_insert(&cell_id, feature_umi, 1.0, &mut out.stats);

                    out.stats.report("feature_tag_match");
                    out.stats.report("emitted_molecules");
                    continue;
                }

                out.stats.report(tag_call.status());
            }

            out.read_tags.insert(ReadTagRecord::new(
                insert_record.id.clone(),
                Some(read.id.clone()),
                &cell.seq,
                &cell.qual,
                &umi.seq,
                &umi.qual,
            ));
            out.fastq_records.push(insert_record);
            out.stats.report("emitted_molecules");
        }

        out
    }
}
 

pub struct OntNormalizer {
    config: OntNormalizerConfig,
    stats: MappingInfo,
    read_tags: ReadTagTable,
    feature_tag_table: Scdata,
}

impl OntNormalizer {
    pub fn new(config: OntNormalizerConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
            read_tags: ReadTagTable::new(),
            feature_tag_table: Scdata::new(1, MatrixValueType::Real),
        }
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        Ok(Self::new(OntNormalizerConfig {
            bam: cli.bam,
            out: cli.out,
            read_tags: cli.read_tags,
            min_transcript_len: cli.min_transcript_len,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            feature_tag_mapper: cli
                .feature_tags
                .mapper()
                .map_err(anyhow::Error::msg)?,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        }))
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &OntNormalizerConfig {
        &self.config
    }

    pub fn write_feature_tag_table_if_present(&mut self) -> Result<()> {
        if self.feature_tag_table.is_empty() {
            return Ok(());
        }

        let Some(mapper) = self.config.feature_tag_mapper.as_ref() else {
            return Ok(());
        };

        let out_dir = self
            .config
            .out
            .with_extension("")
            .join("feature_tag_table_unfiltered");

        std::fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;

        let ftfi = FastTagFeatureIndex::new( mapper );

        self.feature_tag_table
            .write_sparse(&out_dir, &ftfi)
            .map_err(|e| anyhow::anyhow!("writing feature tag table failed: {e}"))?;

        Ok(())
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

        let mut chunk: Vec<FastqRecord> = Vec::with_capacity(CHUNK_SIZE);

        for rec_result in bam.records() {
            let rec = rec_result.context("failed to read BAM record")?;
            chunk.push(FastqRecord::from_bam_record(&rec));

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &mut fastq)?;
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut fastq)?;
        }

        fastq.finish()?;

        self.read_tags
            .save(&self.config.read_tags)
            .with_context(|| {
                format!(
                    "failed to write read-tag table {}",
                    self.config.read_tags.display()
                )
            })?;

        self.write_feature_tag_table_if_present()?;

        Ok(())
    }

    fn process_chunk(
        &mut self,
        input: &[FastqRecord],
        fastq: &mut FastqWriter,
    ) -> Result<()> {

        let partials: Vec<ChunkPartial> = input
            .par_iter()
            .map(|read| self.config.process_read(read ))
            .collect();

        for partial in partials {
            self.stats.merge(&partial.stats);
            self.read_tags.merge(partial.read_tags);
            self.feature_tag_table.merge(&partial.tag_counts);

            for record in partial.fastq_records {
                fastq.write(&record)?;
                self.stats.report("fastq_reads_written");
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

        let feature_tags = info.get_issue_count("feature_tag_match");
        let bad_insert = info.get_issue_count("bad_insert_slice");
        let bad_cell = info.get_issue_count("bad_cell_slice");
        let bad_umi = info.get_issue_count("bad_umi_slice");
        let too_short = info.get_issue_count("short_insert");

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
feature-tag molecules  : {feature_tags}

Orientation
-----------
forward                : {fwd} ({fwd_pct:.2}%)
reverse                : {rev} ({rev_pct:.2}%)

Rejected
--------
bad insert slice       : {bad_insert}
bad cell slice         : {bad_cell}
bad UMI slice          : {bad_umi}
too short insert       : {too_short}
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
            feature_tags = feature_tags,
            fwd = fwd,
            rev = rev,
            fwd_pct = pct(fwd, emitted),
            rev_pct = pct(rev, emitted),
            bad_insert = bad_insert,
            bad_cell = bad_cell,
            bad_umi = bad_umi,
            too_short = too_short,
        )
    }
}

    
fn bytes_to_string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).to_string()
}

fn qual_to_string(qual: &[u8]) -> String {
    qual.iter()
        .map(|q| q.saturating_add(33) as char)
        .collect()
}

fn encode_sequence_id(seq: &[u8]) -> u64 {
    let mut out = 0_u64;

    for &base in seq.iter().take(32) {
        out <<= 2;
        out |= match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
    }

    out
}
