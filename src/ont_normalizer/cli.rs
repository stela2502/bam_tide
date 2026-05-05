use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Clone, Parser)]
#[command(
    author,
    version,
    about = "Normalize messy ONT/Dorado BAM reads into one 10x-style molecule per FASTQ read",
    long_about = "\
Normalize messy ONT/Dorado BAM reads into clean FASTQ records.

The tool searches each BAM read in both orientations for 10x-style 3' barcode
cassettes:

    adapter + cell barcode + UMI + polyT + transcript

Each detected molecule is written as one FASTQ record:

    original_read_name/mol<N>

The output sequence is normalized to the expected molecule orientation. Reads
without a valid cassette are not emitted, but they are counted in the final
summary. A TSV sidecar file records the extracted CB/UMI, qualities, coordinates,
orientation, and status for every emitted molecule.

The normalizer does NOT perform whitelist correction. CB and UMI are raw slices
from the read. Barcode correction can be done later.

Typical use:

    bam-ont-normalizer \\
      --bam dorado.bam \\
      --out normalized.fastq.gz \\
      --tags molecule_tags.tsv \\
      --threads 8 \\
      --gzip-level 1
"
)]
pub struct Cli {
    #[arg(
        long, short,
        value_name = "BAM",
        help = "Input Dorado/ONT BAM file. The BAM may be unmapped; query sequence and qualities are used."
    )]
    pub bam: PathBuf,

    #[arg(
        long, short,
        value_name = "FASTQ[.GZ]",
        help = "Output normalized FASTQ file. By default this is gzip-compressed unless --no-gzip is set."
    )]
    pub out: PathBuf,

    #[arg(
        long, short,
        value_name = "TSV",
        help = "Output molecule metadata TSV. Contains one row per emitted molecule with CB/UMI, qualities, coordinates, orientation, and status."
    )]
    pub tags: PathBuf,

    #[arg(
        long,
        default_value = "CTACACGACGCTCTTCCGATCT",
        value_name = "SEQ",
        help = "Expected 10x adapter/primer sequence before CB+UMI+polyT."
    )]
    pub adapter: String,

    #[arg(
        long,
        default_value_t = 13,
        value_name = "BP",
        help = "Minimum adapter suffix length allowed for cassette detection. Useful when ONT reads contain only the adapter suffix before CB/UMI."
    )]
    pub min_adapter_match: usize,

    #[arg(
        long,
        default_value_t = 2,
        value_name = "N",
        help = "Maximum mismatches allowed in the adapter/suffix match. Bases N/n are treated as unknown and do not count as mismatches."
    )]
    pub max_adapter_mismatches: usize,

    #[arg(
        long,
        default_value_t = 16,
        value_name = "BP",
        help = "Cell barcode length immediately after the adapter."
    )]
    pub cb_len: usize,

    #[arg(
        long,
        default_value_t = 12,
        value_name = "BP",
        help = "UMI length immediately after the cell barcode."
    )]
    pub umi_len: usize,

    #[arg(
        long,
        default_value_t = 10,
        value_name = "N",
        help = "Minimum number of T bases required inside the polyT detection window."
    )]
    pub poly_t_min: usize,

    #[arg(
        long,
        default_value_t = 14,
        value_name = "BP",
        help = "Window size used to validate the expected polyT start after adapter + CB + UMI."
    )]
    pub poly_t_window: usize,

    #[arg(
        long,
        default_value_t = 20,
        value_name = "BP",
        help = "Minimum transcript sequence length required after the polyT stretch. Molecules shorter than this are counted but not emitted."
    )]
    pub min_transcript_len: usize,

    #[arg(
        long,
        default_value_t = 4,
        value_name = "N",
        help = "Number of HTSlib threads used for BAM reading/decompression."
    )]
    pub threads: usize,

    #[arg(
        long,
        default_value_t = 1,
        value_name = "0-9",
        help = "Gzip compression level for FASTQ output. Level 1 is fast and recommended for large ONT files."
    )]
    pub gzip_level: u32,

    #[arg(
        long,
        default_value_t = false,
        help = "Write plain FASTQ instead of gzip-compressed FASTQ. Useful for speed benchmarking."
    )]
    pub no_gzip: bool,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
