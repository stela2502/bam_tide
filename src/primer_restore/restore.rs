use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;
use sc_primer::{Orientation, PrimerDetector};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use crate::ont_normalizer::fastq_record::FastqRecord;
use crate::ont_normalizer::fastq_writer::FastqWriter;
use crate::ont_normalizer::read_tag_table::{
    ReadTagTable,
    ReadTagTableConfig,
    ReadTagRecord,
};

#[derive(Debug, Clone)]
pub struct PrimerRestoreConfig {
    pub fastq: PathBuf,
    pub out: PathBuf,
    pub read_tags: ReadTagTableConfig,
    pub primer: PrimerDetector,
    pub gzip: bool,
    pub gzip_level: u32,
    pub die_on_error: bool,
}

pub struct PrimerRestore {
    config: PrimerRestoreConfig,
    stats: MappingInfo,
}

impl PrimerRestore {
    pub fn new(config: PrimerRestoreConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
        }
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn run(&mut self) -> Result<()> {
        let read_tags = ReadTagTable::from_config(&self.config.read_tags)
            .with_context(|| {
                format!(
                    "failed to load read-tag table: {}",
                    self.config.read_tags.path.display()
                )
            })?;

        let mut reader = FastqReader::from_path(&self.config.fastq)
            .with_context(|| {
                format!("failed to open FASTQ: {}", self.config.fastq.display())
            })?;

        let mut writer = FastqWriter::new(
            &self.config.out,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| {
            format!("failed to create restored FASTQ: {}", self.config.out.display())
        })?;

        while let Some(record) = reader.next_record()? {
            self.stats.report("total_reads");

            let Some(tag) = read_tags.get(record.id.as_str()) else {
                self.stats.report("missing_read_tag");
                if self.config.die_on_error {
                    bail!(
                        "primer-restore encountered {errors} restore errors \
                         and --die-on-error is active"
                    );
                }
                continue;
            };

            let restored = match self.restore_one(&record, tag) {
                Ok(restored) => restored,
                Err(err) => {
                    self.stats.report("primer_restore_failed");
                    eprintln!(
                        "primer-restore: failed to restore read '{}': {err}",
                        record.id
                    );
                    if self.config.die_on_error {
                        bail!(
                            "primer-restore encountered {errors} restore errors \
                             and --die-on-error is active"
                        );
                    }
                    continue;
                }
            };

            writer.write(&restored)?;
            self.stats.report("restored_reads");
        }

        writer.finish()?;
        self.finish_or_fail()
    }

    fn restore_one(
        &self,
        record: &FastqRecord,
        tag: &ReadTagRecord,
    ) -> Result<FastqRecord> {
        let (_target_cell, primer) = self
            .config
            .primer
            .generate(tag)
            .with_context(|| format!("could not generate primer for {}", tag.read_id))?;

        let mut seq = Vec::with_capacity(primer.seq.len() + record.seq.len());
        seq.extend_from_slice(&primer.seq);
        seq.extend_from_slice(&record.seq);

        let mut qual = Vec::with_capacity(primer.qual.len() + record.qual.len());
        qual.extend_from_slice(&primer.qual);
        qual.extend_from_slice(&record.qual);

        Ok(FastqRecord::new(record.id.clone(), &seq, &qual))
    }

    fn finish_or_fail(&self) -> Result<()> {
        let total = self.stats.get_issue_count("total_reads");
        let restored = self.stats.get_issue_count("restored_reads");
        let missing = self.stats.get_issue_count("missing_read_tag");
        let failed = self.stats.get_issue_count("primer_restore_failed");
        let errors = missing + failed;

        eprintln!(
            "\
primer-restore summary
======================
FASTQ reads processed      : {total}
FASTQ reads restored       : {restored}
reads without tag record   : {missing}
primer generation failures : {failed}
"
        );

        if errors > 0 {
            eprintln!(
                "WARNING: primer-restore did not restore {errors} reads. \
                 Output FASTQ is incomplete."
            );
        }

        Ok(())
    }
}

struct FastqReader {
    reader: Box<dyn BufRead>,
    line: String,
}

impl FastqReader {
    fn from_path(path: &Path) -> Result<Self> {
        let reader = open_maybe_gz(path)?;
        Ok(Self {
            reader: Box::new(BufReader::new(reader)),
            line: String::new(),
        })
    }

    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        let id = match self.read_line()? {
            None => return Ok(None),
            Some(line) => line,
        };

        if !id.starts_with('@') {
            bail!("invalid FASTQ record: id line does not start with '@': {id}");
        }

        let seq = self
            .read_line()?
            .context("truncated FASTQ record: missing sequence line")?;

        let plus = self
            .read_line()?
            .context("truncated FASTQ record: missing plus line")?;

        if !plus.starts_with('+') {
            bail!("invalid FASTQ record for {id}: plus line does not start with '+'");
        }

        let qual = self
            .read_line()?
            .context("truncated FASTQ record: missing quality line")?;

        if seq.len() != qual.len() {
            bail!(
                "invalid FASTQ record for {id}: sequence length {} != quality length {}",
                seq.len(),
                qual.len()
            );
        }

        let clean_id = id[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        let qual_phred: Vec<u8> = qual
            .as_bytes()
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect();

        Ok(Some(FastqRecord::new(
            clean_id,
            seq.as_bytes(),
            &qual_phred,
        )))
    }

    fn read_line(&mut self) -> Result<Option<String>> {
        self.line.clear();
        let n = self.reader.read_line(&mut self.line)?;

        if n == 0 {
            return Ok(None);
        }

        while self.line.ends_with('\n') || self.line.ends_with('\r') {
            self.line.pop();
        }

        Ok(Some(self.line.clone()))
    }
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}