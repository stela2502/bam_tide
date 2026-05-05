use crate::ont_normalizer::cassette::{Cassette, CassetteExtractor, Orientation};
use crate::ont_normalizer::cli::Cli;
use crate::ont_normalizer::fastq_writer::FastqWriter;
use crate::ont_normalizer::stats::NormalizeStats;

use mapping_info::MappingInfo;

use anyhow::{Context, Result};
use rust_htslib::bam::{Read, Reader};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

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
            stats: MappingInfo::new(None,0.0,0),
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
            bam.set_threads(self.config.threads.max(3))
                .context("failed to set BAM reader threads")?;
        }

        let mut fastq =
            FastqWriter::new(&self.config.out, self.config.gzip, self.config.gzip_level)
                .with_context(|| {
                    format!("failed to create FASTQ: {}", self.config.out.display())
                })?;

        let mut tags = self.create_tag_writer()?;

        let extractor = CassetteExtractor::new(
            self.config.adapter.clone(),
            self.config.min_adapter_match,
            self.config.cb_len,
            self.config.umi_len,
            self.config.poly_t_min,
            self.config.poly_t_window,
            self.config.min_transcript_len,
            self.config.max_adapter_mismatches,
        );

        for rec_result in bam.records() {
            let rec = rec_result?;
            self.process_record(&rec, &extractor, &mut fastq, &mut tags)?;
        }

        fastq.finish()?;
        tags.flush()?;

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
            if d == 0 { 0.0 } else { (n as f64 / d as f64) * 100.0 }
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

    fn create_tag_writer(&self) -> Result<BufWriter<File>> {
        let file = File::create(&self.config.tags).with_context(|| {
            format!("failed to create tags TSV: {}", self.config.tags.display())
        })?;

        let mut writer = BufWriter::new(file);

        writeln!(
            writer,
            "output_read_id\toriginal_read_id\tmolecule_index\torientation\tadapter_start\tadapter_end\tsegment_start\tsegment_end\traw_cb\tquality_cb\traw_umi\tquality_umi\tpoly_t_start\tpoly_t_len\tstatus"
        )?;

        Ok(writer)
    }

    fn process_record(
        &mut self,
        rec: &rust_htslib::bam::Record,
        extractor: &CassetteExtractor,
        fastq: &mut FastqWriter,
        tags: &mut BufWriter<File>,
    ) -> Result<()> {
        self.stats.report("total_records");

        let read_name = String::from_utf8_lossy(rec.qname()).to_string();
        let seq = rec.seq().as_bytes();
        let qual = rec.qual().to_vec();

        let cassettes = extractor.extract_both_orientations(
            &seq,
            &qual,
            &mut self.stats,
        );

        match cassettes.len() {
            0 => self.stats.report("zero_cassette"),
            1 => self.stats.report("one_cassette"),
            _ => self.stats.report("multi_cassette"),
        }

        let (rc_seq, rc_qual) = extractor.revcomp_with_qual(&seq, &qual);

        for (idx, cassette) in cassettes.iter().enumerate() {
            let mol_index = idx + 1;
            let out_id = format!("{read_name}/mol{mol_index}");

            let (oriented_seq, oriented_qual) = match cassette.orientation {
                Orientation::Forward => (&seq[..], &qual[..]),
                Orientation::ReverseComplement => (&rc_seq[..], &rc_qual[..]),
            };

            self.write_molecule(
                fastq,
                tags,
                &out_id,
                &read_name,
                mol_index,
                cassette,
                oriented_seq,
                oriented_qual,
            )?;
        }

        Ok(())
    }

    fn write_molecule(
        &mut self,
        fastq: &mut FastqWriter,
        tags: &mut BufWriter<File>,
        out_id: &str,
        original_id: &str,
        mol_index: usize,
        cassette: &Cassette,
        seq: &[u8],
        qual: &[u8],
    ) -> Result<()> {
        let mol_seq = &seq[cassette.segment_start..cassette.segment_end];
        let mol_qual = &qual[cassette.segment_start..cassette.segment_end];

        fastq.write_record(out_id, mol_seq, mol_qual)?;

        writeln!(
            tags,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tPASS",
            out_id,
            original_id,
            mol_index,
            cassette.orientation.as_str(),
            cassette.adapter_start,
            cassette.adapter_end,
            cassette.segment_start,
            cassette.segment_end,
            Self::bytes_to_string(&cassette.cb),
            Self::qual_to_string(&cassette.cb_qual),
            Self::bytes_to_string(&cassette.umi),
            Self::qual_to_string(&cassette.umi_qual),
            cassette.poly_t_start,
            cassette.poly_t_len,
        )?;

        self.stats.report("emitted_molecules");
        self.stats.report("fastq_reads_written");

        match cassette.orientation {
            Orientation::Forward => self.stats.report("forward_molecules"),
            Orientation::ReverseComplement => self.stats.report("reverse_molecules"),
        }

        Ok(())
    }

    fn bytes_to_string(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).to_string()
    }

    fn qual_to_string(qual: &[u8]) -> String {
        qual.iter()
            .map(|q| char::from(q.saturating_add(33)))
            .collect()
    }
}
