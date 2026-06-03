use crate::fastq::{record::FastqRecord, writer::FastqWriter};
use crate::illumina_normalizer::cli::{Cli, InsertRead, PrimerRead};
use crate::read_tag_table::{ReadTagTable, ReadTagRecord, ReadTagTableWriter};
use crate::index::FastTagFeatureIndex;

use scdata::{Scdata, MatrixValueType};

use anyhow::{bail, Context, Result};
use mapping_info::MappingInfo;
use rayon::prelude::*;
use sc_primer::{Orientation, PrimerDetector};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

const CHUNK_SIZE: usize = 1_000_000;

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

impl IlluminaNormalizerConfig{
    fn normalize_pair(
        &self,
        r1: &FastqRecord,
        r2: &FastqRecord,
        stats: &mut MappingInfo,
    ) -> Result<NormalizedRead> {
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

        let cell = primer_match.get_cell(&primer_source.seq, &primer_source.qual).map_err(|err| {
            stats.report("bad_cell_slice");
            anyhow::anyhow!(err)
        })?;

        let umi = primer_match.get_umi(&primer_source.seq, &primer_source.qual).map_err(|err| {
            stats.report("bad_umi_slice");
            anyhow::anyhow!(err)
        })?;

        let orientation = match primer_match.orientation {
            Orientation::Forward => "forward",
            Orientation::ReverseComplement => "reverse_complement",
        };

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

        let read_id = normalized_id(&r1.id);
        emitted.id = read_id;

        if self.insert_read == InsertRead::Detected && primer_match.orientation == Orientation::ReverseComplement {
            emitted = emitted.revcomp();
        }

        stats.report("accepted_pairs");

        Ok(NormalizedRead {
            fastq: emitted,
            original_read_id: original_read_id(r1, r2, self.primer_read),
            orientation,
            cell_seq: cell.seq.to_vec(),
            cell_qual: cell.qual.to_vec(),
            umi_seq: umi.seq.to_vec(),
            umi_qual: umi.qual.to_vec(),
        })
    }
}

#[derive(Debug, Clone)]
struct NormalizedRead {
    fastq: FastqRecord,
    original_read_id: String,
    orientation: &'static str,
    cell_seq: Vec<u8>,
    cell_qual: Vec<u8>,
    umi_seq: Vec<u8>,
    umi_qual: Vec<u8>,
}

struct ChunkPartial {
    reads: Vec<NormalizedRead>,
    stats: MappingInfo,
}

pub struct IlluminaNormalizer {
    config: IlluminaNormalizerConfig,
    stats: MappingInfo,
    tag_counts: Scdata,
    read_tags: ReadTagTable,
}

impl IlluminaNormalizer {
    pub fn new(config: IlluminaNormalizerConfig) -> Self {
        Self {
            config,
            stats: MappingInfo::new(None, 0.0, 0),
            read_tags: ReadTagTable::new(),
            tag_counts: Scdata::new(1, MatrixValueType::Real),
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
            primer: cli
                    .primer
                    .detector()
                    .map_err(anyhow::Error::msg)?,
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
        if self.config.threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.config.threads)
                .build_global()
                .ok();
        }

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

        while let Some(pair) = reader.next_pair()? {
            chunk.push(pair);

            if chunk.len() >= CHUNK_SIZE {
                self.process_chunk(&chunk, &mut fastq)?;
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            self.process_chunk(&chunk, &mut fastq)?;
        }

        fastq.finish()?;

        self.read_tags.save(&self.config.read_tags)?;

        Ok(())
    }

    fn process_chunk(
        &mut self,
        input: &[(FastqRecord, FastqRecord)],
        fastq: &mut FastqWriter,
    ) -> Result<()> {
        let config = self.config.clone();

        let partials: Vec<ChunkPartial> = input
            .par_iter()
            .map(|(r1, r2)| {
                let mut out = ChunkPartial {
                    reads: Vec::new(),
                    stats: MappingInfo::new(None, 0.0, 0),
                };

                out.stats.report("total_pairs");

                match self.config.normalize_pair(r1, r2, &mut out.stats) {
                    Ok(read) => out.reads.push(read),
                    Err(_) => out.stats.report("failed_pairs"),
                }

                out
            })
            .collect();

        for partial in partials {
            self.stats.merge(&partial.stats);

            for read in partial.reads {
                fastq.write(&read.fastq)?;

                self.read_tags.insert(ReadTagRecord::new(
                    read.fastq.id.clone(),
                    Some(read.original_read_id.clone()),
                    &read.cell_seq,
                    &read.cell_qual,
                    &read.umi_seq,
                    &read.umi_qual,
                ));

                self.stats.report("fastq_reads_written");
            }
        }

        Ok(())
    }

    pub fn stats_report(&self) -> String {
        let info = &self.stats;

        let total = info.get_issue_count("total_pairs");
        let written = info.get_issue_count("fastq_reads_written");
        let failed = info.get_issue_count("failed_pairs");
        let no_primer = info.get_issue_count("no_primer_match");
        let short_insert = info.get_issue_count("short_insert");
        let bad_cell = info.get_issue_count("bad_cell_slice");
        let bad_umi = info.get_issue_count("bad_umi_slice");

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
FASTQ reads written   : {written} ({written_pct:.2}%)
failed pairs          : {failed} ({failed_pct:.2}%)

Rejected
--------
no primer match       : {no_primer}
short insert/read     : {short_insert}
bad cell slice        : {bad_cell}
bad UMI slice         : {bad_umi}
"#,
            total = total,
            primer_read = self.config.primer_read.as_str(),
            insert_read = self.config.insert_read.as_str(),
            written = written,
            written_pct = pct(written, total),
            failed = failed,
            failed_pct = pct(failed, total),
            no_primer = no_primer,
            short_insert = short_insert,
            bad_cell = bad_cell,
            bad_umi = bad_umi,
        )
    }

    fn create_read_tag_writer(&self) -> Result<ReadTagTableWriter<BufWriter<File>>> {
        let file = File::create(&self.config.read_tags).with_context(|| {
            format!(
                "failed to create read-tags TSV: {}",
                self.config.read_tags.display()
            )
        })?;

        ReadTagTableWriter::new(BufWriter::new(file))
    }
}



fn normalized_id(original_r1_id: &str) -> String {
    format!("{original_r1_id}/mol0")
}

fn original_read_id(r1: &FastqRecord, r2: &FastqRecord, primer_read: PrimerRead) -> String {
    match primer_read {
        PrimerRead::R1 => r1.id.clone(),
        PrimerRead::R2 => r2.id.clone(),
    }
}

fn qual_string(qual: &[u8]) -> String {
    qual.iter().map(|q| q.saturating_add(33) as char).collect()
}

/// Minimal paired FASTQ reader.
///
/// If you already have a crate-level FASTQ reader, replace this with that one.
/// This local implementation is intentionally simple and keeps this normalizer
/// independent from the old BD-only Illumina code.
struct FastqPairReader {
    r1: Box<dyn FastqRead>,
    r2: Box<dyn FastqRead>,
}

impl FastqPairReader {
    fn from_paths(r1: &PathBuf, r2: &PathBuf) -> Result<Self> {
        Ok(Self {
            r1: open_fastq_reader(r1)?,
            r2: open_fastq_reader(r2)?,
        })
    }

    fn next_pair(&mut self) -> Result<Option<(FastqRecord, FastqRecord)>> {
        let r1 = self.r1.next_record()?;
        let r2 = self.r2.next_record()?;

        match (r1, r2) {
            (Some(r1), Some(r2)) => Ok(Some((r1, r2))),
            (None, None) => Ok(None),
            (Some(_), None) => bail!("R1 has more records than R2"),
            (None, Some(_)) => bail!("R2 has more records than R1"),
        }
    }
}

trait FastqRead: Send {
    fn next_record(&mut self) -> Result<Option<FastqRecord>>;
}

fn open_fastq_reader(path: &PathBuf) -> Result<Box<dyn FastqRead>> {
    use flate2::read::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead + Send> = if path.extension().is_some_and(|x| x == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(Box::new(SimpleFastqReader { reader, line: String::new() }))
}

struct SimpleFastqReader {
    reader: Box<dyn std::io::BufRead + Send>,
    line: String,
}

impl FastqRead for SimpleFastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        self.line.clear();

        if self.reader.read_line(&mut self.line)? == 0 {
            return Ok(None);
        }

        let id_line = self.line.trim_end().to_string();
        if !id_line.starts_with('@') {
            bail!("invalid FASTQ record: expected @ header, got {id_line}");
        }

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let seq = self.line.trim_end().as_bytes().to_vec();

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let plus = self.line.trim_end().to_string();
        if !plus.starts_with('+') {
            bail!("invalid FASTQ record: expected + line, got {plus}");
        }

        self.line.clear();
        self.reader.read_line(&mut self.line)?;
        let qual_ascii = self.line.trim_end().as_bytes().to_vec();

        if seq.len() != qual_ascii.len() {
            bail!(
                "invalid FASTQ record {}: sequence length {} != quality length {}",
                id_line,
                seq.len(),
                qual_ascii.len()
            );
        }

        let qual: Vec<u8> = qual_ascii
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect();

        Ok(Some(FastqRecord::new(
            id_line.trim_start_matches('@').to_string(),
            &seq,
            &qual,
        )))
    }
}
