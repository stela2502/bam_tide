//job.rs
use anyhow::{Context, Result};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{HeaderView, Record};

use crate::core::ref_block::record_to_blocks;
use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{SplicedRead, Strand};

use int_to_str::int_to_str::IntToStr;

use mapping_info::MappingInfo;
use snp_index::{AlignedRead, Genome, RefineOptions, SnpIndex};

#[derive(Clone)]
pub struct Job {
    pub cell: u64,
    pub umi: u64,
    pub spliced: SplicedRead,

    /// Refined/cleaned read used for SNP matching.
    ///
    /// Currently created when genome/SNP support is active.
    pub aligned: Option<AlignedRead>,
}

pub struct JobBuilder<'a> {
    header: &'a HeaderView,
    chr_map: &'a std::collections::HashMap<String, usize>,
    genome: Option<&'a Genome>,
    snp: Option<&'a SnpIndex>,
    min_mapq: u8,
    read1_only: bool,
    refine_against_genome: bool,
}

impl<'a> JobBuilder<'a> {
    pub fn new(
        header: &'a HeaderView,
        chr_map: &'a std::collections::HashMap<String, usize>,
    ) -> Self {
        Self {
            header,
            chr_map,
            genome: None,
            snp: None,
            min_mapq: 0,
            read1_only: false,
            refine_against_genome: false,
        }
    }

    pub fn with_genome(mut self, genome: Option<&'a Genome>, refine_against_genome: bool) -> Self {
        self.genome = genome;
        self.refine_against_genome = refine_against_genome;
        self
    }

    pub fn with_snp_index(mut self, snp: Option<&'a SnpIndex>) -> Self {
        self.snp = snp;
        self
    }

    pub fn with_min_mapq(mut self, min_mapq: u8) -> Self {
        self.min_mapq = min_mapq;
        self
    }

    pub fn read1_only(mut self, read1_only: bool) -> Self {
        self.read1_only = read1_only;
        self
    }

    pub fn build(&self, rec: &Record, report: &mut MappingInfo) -> Result<Option<Job>> {
        if rec.is_unmapped() {
            report.report("unmapped");
            return Ok(None);
        }

        if rec.mapq() < self.min_mapq {
            report.report("mapq failed");
            return Ok(None);
        }

        if rec.is_secondary() || rec.is_supplementary() {
            report.report("secondary or supplementary");
            return Ok(None);
        }

        if self.read1_only && !rec.is_first_in_template() {
            report.report("read!=1");
            return Ok(None);
        }

        let cb_raw = match Self::aux_tag_str(rec, *b"CB") {
            Some(v) => v,
            None => {
                report.report("no CB tag");
                return Ok(None);
            }
        };

        let ub = match Self::aux_tag_str(rec, *b"UB") {
            Some(v) => v,
            None => {
                report.report("no UB tag");
                return Ok(None);
            }
        };

        let cb = Self::normalize_10x_barcode(cb_raw);

        let cell = match Self::dna_to_u64(cb) {
            Some(v) => v,
            None => {
                report.report("invalid CB");
                return Ok(None);
            }
        };

        let umi = match Self::dna_to_u64(ub) {
            Some(v) => v,
            None => {
                report.report("invalid UB");
                return Ok(None);
            }
        };

        let tid = rec.tid();

        if tid < 0 {
            report.report("tid below 0");
            return Ok(None);
        }

        let chr_name = std::str::from_utf8(self.header.tid2name(tid as u32))
            .context("Invalid chromosome name in BAM header")?;

        let chr_id = match self.chr_map.get(chr_name) {
            Some(&id) => id,
            None => {
                report.report("contig not in index");
                return Ok(None);
            }
        };

        let spliced = match Self::record_to_spliced_read(rec, chr_id) {
            Some(v) => v,
            None => {
                report.report("no ref blocks");
                return Ok(None);
            }
        };

        let aligned = self.build_aligned_read(rec, chr_name);
        /*
        if rec.qname() == b"795ade5c-7cba-4e99-a3a3-f105872d4b18_0" {
            eprintln!("FOUND DEBUG READ: 795ade5c-7cba-4e99-a3a3-f105872d4b18_0");
            eprintln!("  chr_id={chr_id}");
            eprintln!("  bam_pos0={}", rec.pos());
            eprintln!("  cigar={}", rec.cigar());

            match &aligned {
                Some(read) => {
                    eprintln!("  aligned.ref_span={:?}", read.ref_span());
                    eprintln!("  aligned.strand={:?}", read.strand);
                    eprintln!("  aligned.seq_len={}", read.seq.len());
                    eprintln!("  aligned.qual_len={:?}", read.qual.as_ref().map(|q| q.len()));
                    eprintln!("  aligned.ops_len={}", read.ops.len());

                    for pos0 in [7_674_894 - 1, 7_674_953 - 1, 7_675_994 - 1] {
                        eprintln!(
                            "  probe pos1={} base_at_ref_pos={:?}",
                            pos0 + 1,
                            read.base_at_ref_pos(pos0)
                        );
                    }

                    match self.snp{
                        Some(snp) => {
                            let (ref_ids, alt_ids, other_ids) = snp.get_ref_alt_other_ids_for_read(&read, 0);
                            eprintln!(
                                "Detected SNP ids:\n\tref: {:?}\n\talt: {:?}\n\tother: {:?}",
                                ref_ids,alt_ids, other_ids
                            );
                        },
                        None => {
                            eprintln!("The sno index was not stored correctly!");
                        }
                    }
                }
                None => {
                    eprintln!("  aligned=None");
                }
            }
        }
        */

        Ok(Some(Job {
            cell,
            umi,
            spliced,
            aligned,
        }))
    }

    fn build_aligned_read(&self, rec: &Record, chr_name: &str) -> Option<AlignedRead> {
        if self.genome.is_none() && self.snp.is_none() {
            return None;
        }

        let chr_id = if let Some(snp) = self.snp {
            snp.chr_id(chr_name)?
        } else {
            *self.chr_map.get(chr_name)?
        };

        let mut read = AlignedRead::from_record(rec, chr_id);

        if let Some(genome) = self.genome
            && self.refine_against_genome
        {
            read.refine_against_genome(genome, RefineOptions::default());
        }

        Some(read)
    }

    fn aux_tag_str(rec: &Record, tag: [u8; 2]) -> Option<&str> {
        match rec.aux(&tag).ok()? {
            Aux::String(s) => Some(s),
            _ => None,
        }
    }

    fn normalize_10x_barcode(cb: &str) -> &str {
        match cb.split_once('-') {
            Some((core, _)) => core,
            None => cb,
        }
    }

    fn dna_to_u64(seq: &str) -> Option<u64> {
        if !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            return None;
        }

        let tool = IntToStr::new(seq.as_bytes());
        Some(tool.into_u64())
    }

    fn record_to_spliced_read(rec: &Record, chr_id: usize) -> Option<SplicedRead> {
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
        spliced.finalize();

        Some(spliced)
    }
}
