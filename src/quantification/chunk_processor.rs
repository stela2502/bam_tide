//chunk_processor.rs
use anyhow::Result;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex};
use rayon::prelude::*;
use scdata::cell_data::GeneUmiHash;

use crate::quantification::cli::QuantMode;
use crate::quantification::job::Job;
use crate::quantification::snp::SnpSideChannel;
use crate::results::QuantData;

pub struct ChunkProcessor<'a> {
    idx: &'a SpliceIndex,
    snp: Option<&'a SnpSideChannel>,
    match_opts: MatchOptions,
    #[allow(dead_code)]
    min_mapq: u8,
}

impl<'a> ChunkProcessor<'a> {
    pub fn new(
        idx: &'a SpliceIndex,
        snp: Option<&'a SnpSideChannel>,
        match_opts: MatchOptions,
        min_mapq: u8,
    ) -> Self {
        Self {
            idx,
            snp,
            match_opts,
            min_mapq,
        }
    }

    pub fn process_into(
        &self,
        quant_mode: QuantMode,
        jobs: &[Job],
        merged: &mut QuantData,
    ) -> Result<()> {
        if jobs.is_empty() {
            return Ok(());
        }

        let threads = rayon::current_num_threads().max(1);
        let chunk_size = (jobs.len() / threads).max(10_000);

        let partials: Vec<QuantData> = jobs
            .par_chunks(chunk_size)
            .map(|chunk| self.process_partial(quant_mode, chunk))
            .collect();

        merged.report.stop_multi_processor_time();

        for partial in partials {
            merged.merge(&partial);
        }

        merged.report.stop_single_processor_time();

        Ok(())
    }

    fn process_partial(&self, quant_mode: QuantMode, jobs: &[Job]) -> QuantData {
        let mut out = QuantData::new();

        for job in jobs {
            match quant_mode {
                QuantMode::Gene => self.add_gene_hit(job, &mut out),
                QuantMode::Transcript => self.add_transcript_hit(job, &mut out),
            }

            self.add_snp_hits(job, &mut out);
        }

        out
    }

    fn add_gene_hit(&self, job: &Job, out: &mut QuantData) {
        let hits = self.idx.match_genes(&job.spliced, self.match_opts);

        if hits.is_empty() {
            out.report.report("no hit");
            return;
        }

        let hit = &hits[0];
        out.report.report(hit.best_hit.class.to_string());

        let feature_id = hit.gene_id as u64;
        let feature_umi = GeneUmiHash(feature_id, job.umi);

        if hit.best_hit.class == MatchClass::Intronic {
            out.intron
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        } else {
            out.gene
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        }
    }

    fn add_transcript_hit(&self, job: &Job, out: &mut QuantData) {
        let hits = self.idx.match_transcripts(&job.spliced, self.match_opts);

        if hits.is_empty() {
            out.report.report("no hit");
            return;
        }

        let hit = &hits[0];
        out.report.report(hit.hit.class.to_string());

        let feature_id = hit.transcript_id as u64;
        let feature_umi = GeneUmiHash(feature_id, job.umi);

        if hit.hit.class == MatchClass::Intronic {
            out.intron
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        } else {
            out.gene
                .try_insert(&job.cell, feature_umi, 1.0, &mut out.report);
        }
    }

    fn add_snp_hits(&self, job: &Job, out: &mut QuantData) {
        let Some(snp) = self.snp else {
            return;
        };

        snp.add_hits(
            job.aligned.as_ref(),
            job.cell,
            job.umi,
            &mut out.snp_ref,
            &mut out.snp_alt,
            &mut out.report,
        );
    }
}
