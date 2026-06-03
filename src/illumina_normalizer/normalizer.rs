use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::illumina_normalizer::cli::{Cli, InsertRead, PrimerRead};
use crate::ngs_normalizer::{
    FastqPairReader, NgsNormalizerSupport, NormalizedMolecule, NormalizerPartial, CHUNK_SIZE,
};
use crate::read_tag_table::ReadTagTable;

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use sc_primer::{Orientation, PrimerDetector};
use scdata::Scdata;

use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct IlluminaNormalizerConfig {
    pub r1: PathBuf,
    pub r2: PathBuf,
    pub out: PathBuf,
    pub read_tags: PathBuf,
    pub primer_read: PrimerRead,
    pub insert_read: InsertRead,
    pub primer: PrimerDetector,
    pub feature_tags: Option<crate::tags::FastTagMapper>,
    pub min_insert_len: usize,
    pub threads: usize,
    pub gzip_level: u32,
    pub gzip: bool,
}

impl IlluminaNormalizerConfig {
    fn normalize_pair(
        &self,
        r1: &FastqRecord,
        r2: &FastqRecord,
        stats: &mut MappingInfo,
        tag_counts: &mut Scdata,
    ) -> Result<Option<NormalizedMolecule>> {
        stats.report("total_pairs");

        let primer_source = match self.primer_read {
            PrimerRead::R1 => r1,
            PrimerRead::R2 => r2,
        };

        let primer_match = match self.primer.detect_first(&primer_source.seq, &primer_source.qual) {
            Ok(Some(x)) => x,
            Ok(None) => {
                stats.report("no_primer_match");
                bail!("no primer match");
            }
            Err(err) => {
                stats.report("no_primer_match");
                bail!("primer detection failed: {err}");
            }
        };

        let cell = primer_match
            .get_cell(&primer_source.seq, &primer_source.qual)
            .map_err(|err| {
                stats.report("bad_cell_slice");
                anyhow::anyhow!(err)
            })?;

        let umi = primer_match
            .get_umi(&primer_source.seq, &primer_source.qual)
            .map_err(|err| {
                stats.report("bad_umi_slice");
                anyhow::anyhow!(err)
            })?;

        let mut emitted = match self.insert_read {
            InsertRead::Detected => primer_match
                .get_insert(&primer_source.seq, &primer_source.qual)
                .map(|insert| FastqRecord::new("tmp", &insert.seq, &insert.qual))
                .map_err(|err| {
                    stats.report("bad_insert_slice");
                    anyhow::anyhow!(err)
                })?,
            InsertRead::R1 => r1.clone(),
            InsertRead::R2 => r2.clone(),
        };

        if emitted.seq.len() < self.min_insert_len {
            stats.report("short_insert");
            bail!("emitted read is shorter than --min-insert-len");
        }

        emitted.id = NgsNormalizerSupport::normalized_molecule_id(&r1.id, 0);

        if self.insert_read == InsertRead::Detected
            && primer_match.orientation == Orientation::ReverseComplement
        {
            emitted = emitted.revcomp();
        }

        NgsNormalizerSupport::report_orientation(stats, primer_match.orientation);

        if NgsNormalizerSupport::maybe_collect_feature_tag(
            self.feature_tags.as_ref(),
            &emitted.seq,
            &cell.seq,
            &umi.seq,
            tag_counts,
            stats,
        ) {
            stats.report("accepted_pairs");
            return Ok(None);
        }

        stats.report("accepted_pairs");

        Ok(Some(NormalizedMolecule {
            fastq: emitted,
            original_read_id: Some(Self::original_read_id(r1, r2, self.primer_read)),
            orientation: primer_match.orientation,
            cell_seq: cell.seq.to_vec(),
            cell_qual: cell.qual.to_vec(),
            umi_seq: umi.seq.to_vec(),
            umi_qual: umi.qual.to_vec(),
        }))
    }

    fn original_read_id(r1: &FastqRecord, r2: &FastqRecord, primer_read: PrimerRead) -> String {
        match primer_read {
            PrimerRead::R1 => r1.id.clone(),
            PrimerRead::R2 => r2.id.clone(),
        }
    }
}

pub struct IlluminaNormalizer {
    config: IlluminaNormalizerConfig,
    stats: MappingInfo,
    read_tags: ReadTagTable,
    feature_tag_table: Scdata,
}

impl IlluminaNormalizer {
    pub fn new(config: IlluminaNormalizerConfig) -> Self {
        Self {
            config,
            stats: NgsNormalizerSupport::new_stats(),
            read_tags: ReadTagTable::new(),
            feature_tag_table: NgsNormalizerSupport::new_feature_tag_table(),
        }
    }

    pub fn from_cli(cli: Cli) -> Result<Self> {
        Ok(Self::new(IlluminaNormalizerConfig {
            r1: cli.r1,
            r2: cli.r2,
            out: cli.out,
            read_tags: cli.read_tags,
            primer_read: cli.primer_read,
            insert_read: cli.insert_read,
            primer: cli.primer.detector().map_err(anyhow::Error::msg)?,
            feature_tags: cli.feature_tags.mapper().map_err(anyhow::Error::msg)?,
            min_insert_len: cli.min_insert_len,
            threads: cli.threads,
            gzip_level: cli.gzip_level,
            gzip: !cli.no_gzip,
        }))
    }

    pub fn stats(&self) -> &MappingInfo {
        &self.stats
    }

    pub fn config(&self) -> &IlluminaNormalizerConfig {
        &self.config
    }

    pub fn run(&mut self) -> Result<()> {
        NgsNormalizerSupport::configure_rayon_threads(self.config.threads);

        let mut reader = FastqPairReader::from_paths(&self.config.r1, &self.config.r2)
            .with_context(|| {
                format!(
                    "failed to open FASTQ pair: {} and {}",
                    self.config.r1.display(),
                    self.config.r2.display()
                )
            })?;

        let mut fastq = FastqWriter::new(
            &self.config.out,
            self.config.gzip,
            self.config.gzip_level,
        )
        .with_context(|| format!("failed to create FASTQ: {}", self.config.out.display()))?;

        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK_SIZE);
        let mut processed_pairs = 100;
        
        while let Some(pair) = reader.next_pair()? {
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                processed_pairs += CHUNK_SIZE;
                self.process_chunk(&chunk, &mut fastq)?;
                chunk.clear();

                eprintln!(
                    "processed {} read pairs; written {} FASTQ reads; failed {} pairs; sample-tag reads {}",
                    self.stats.get_issue_count("total_pairs"),
                    self.stats.get_issue_count("fastq_reads_written"),
                    self.stats.get_issue_count("failed_pairs"),
                    self.stats.get_issue_count("feature_tag_match"),
                );

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

        NgsNormalizerSupport::write_feature_tag_table_if_present(
            &mut self.feature_tag_table,
            self.config.feature_tags.as_ref(),
            &self.config.out,
        )?;

        Ok(())
    }

    fn process_chunk(
        &mut self,
        input: &[(FastqRecord, FastqRecord)],
        fastq: &mut FastqWriter,
    ) -> Result<()> {
        let config = self.config.clone();

        let partials: Vec<NormalizerPartial> = input
            .par_iter()
            .map(|(r1, r2)| {
                let mut out = NormalizerPartial::new();

                match config.normalize_pair(r1, r2, &mut out.stats, &mut out.feature_tag_table) {
                    Ok(Some(molecule)) => out.push_molecule(molecule),
                    Ok(None) => {}
                    Err(_) => out.stats.report("failed_pairs"),
                }

                out
            })
            .collect();

        for partial in partials {
            partial.merge_into(
                &mut self.stats,
                &mut self.read_tags,
                &mut self.feature_tag_table,
                fastq,
            )?;
        }

        Ok(())
    }

    pub fn stats_report(&self) -> String {
        let info = &self.stats;

        let total = info.get_issue_count("total_pairs");
        let written = info.get_issue_count("fastq_reads_written");
        let failed = info.get_issue_count("failed_pairs");
        let accepted = info.get_issue_count("accepted_pairs");
        let no_primer = info.get_issue_count("no_primer_match");
        let short_insert = info.get_issue_count("short_insert");
        let bad_insert = info.get_issue_count("bad_insert_slice");
        let bad_cell = info.get_issue_count("bad_cell_slice");
        let bad_umi = info.get_issue_count("bad_umi_slice");
        let feature_tags = info.get_issue_count("feature_tag_match");
        let fwd = info.get_issue_count("forward_molecules");
        let rev = info.get_issue_count("reverse_molecules");

        let pct = |n: usize, d: usize| {
            if d == 0 {
                0.0
            } else {
                (n as f64 / d as f64) * 100.0
            }
        };

        format!(
            r#"bam-illumina-normalizer summary
=================================

Input
-----
FASTQ pairs processed : {total}
primer read           : {primer_read}
insert read           : {insert_read}

Output
------
accepted pairs        : {accepted} ({accepted_pct:.2}%)
FASTQ reads written   : {written} ({written_pct:.2}%)
feature-tag molecules : {feature_tags}
failed pairs          : {failed} ({failed_pct:.2}%)

Orientation
-----------
forward               : {fwd}
reverse               : {rev}

Rejected
--------
no primer match       : {no_primer}
short insert/read     : {short_insert}
bad insert slice      : {bad_insert}
bad cell slice        : {bad_cell}
bad UMI slice         : {bad_umi}
"#,
            total = total,
            primer_read = self.config.primer_read.as_str(),
            insert_read = self.config.insert_read.as_str(),
            accepted = accepted,
            accepted_pct = pct(accepted, total),
            written = written,
            written_pct = pct(written, total),
            feature_tags = feature_tags,
            failed = failed,
            failed_pct = pct(failed, total),
            fwd = fwd,
            rev = rev,
            no_primer = no_primer,
            short_insert = short_insert,
            bad_insert = bad_insert,
            bad_cell = bad_cell,
            bad_umi = bad_umi,
        )
    }
}
