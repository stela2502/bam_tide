use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Chemistry {
    #[value(name = "tenx-v2")]
    TenxV2,

    #[value(name = "tenx-v3")]
    TenxV3,

    #[value(name = "tenx-v4")]
    TenxV4,

    #[value(name = "bd-v1")]
    BdV1,

    #[value(name = "bd-v2-96")]
    BdV2_96,

    #[value(name = "bd-v2-384")]
    BdV2_384,
}

impl Chemistry {
    pub fn is_bd(self) -> bool {
        matches!(
            self,
            Self::BdV1 | Self::BdV2_96 | Self::BdV2_384
        )
    }

    pub fn is_10x(self) -> bool {
        matches!(
            self,
            Self::TenxV2 | Self::TenxV3 | Self::TenxV4
        )
    }

    pub fn rhapsody_version(self) -> Option<&'static str> {
        match self {
            Self::BdV1 => Some("v1"),
            Self::BdV2_96 => Some("v2.96"),
            Self::BdV2_384 => Some("v2.384"),
            _ => None,
        }
    }

    pub fn default_cb_len(self) -> usize {
        match self {
            Self::TenxV2 | Self::TenxV3 | Self::TenxV4 => 16,
            Self::BdV1 | Self::BdV2_96 | Self::BdV2_384 => 0,
        }
    }

    pub fn default_umi_len(self) -> usize {
        match self {
            Self::TenxV2 => 10,
            Self::TenxV3 | Self::TenxV4 => 12,
            Self::BdV1 => 8,
            Self::BdV2_96 | Self::BdV2_384 => 6,
        }
    }
}

#[derive(Debug, Clone, Parser)]
#[command(
    author,
    version,
    about = "Normalize BD Rhapsody Illumina FASTQ pairs into STAR-friendly reads plus barcode metadata",
    long_about = "\
Normalize BD Rhapsody Illumina paired FASTQ reads.

The tool reads R1/R2 FASTQ pairs. R1 is parsed as a BD Rhapsody
cell-label + UMI read:

    C1 + linker + C2 + linker + C3 + UMI + polyT/trailing sequence

The normalizer decodes the BD cell id, extracts the UMI and qualities, and
writes an artificial barcode read plus the original R2 read.

The output is designed for workflows where STAR performs mapping only and
bam-quant performs cell/UMI accounting from the sidecar table.

Typical use:

    bd-illumina-normalizer \\
      --r1 sample_R1.fastq.gz \\
      --r2 sample_R2.fastq.gz \\
      --out-r1 normalized_R1.fastq.gz \\
      --out-r2 normalized_R2.fastq.gz \\
      --tags molecule_tags.tsv \\
      --rhapsody-version v2.384 \\
      --threads 8 \\
      --gzip-level 1
"
)]
pub struct Cli {
    #[arg(long, value_name = "FASTQ[.GZ]", help = "Input BD Rhapsody R1 FASTQ containing cell labels and UMI.")]
    pub r1: PathBuf,

    #[arg(long, value_name = "FASTQ[.GZ]", help = "Input BD Rhapsody R2 FASTQ containing biological sequence.")]
    pub r2: PathBuf,

    #[arg(long, value_name = "FASTQ[.GZ]", help = "Output artificial normalized R1 FASTQ.")]
    pub out_r1: PathBuf,

    #[arg(long, value_name = "FASTQ[.GZ]", help = "Output normalized R2 FASTQ. Usually the original R2 sequence with synchronized read names.")]
    pub out_r2: PathBuf,

    #[arg(long, short, value_name = "TSV", help = "Output molecule metadata TSV with BD cell id, UMI, qualities, shift, sample id, and status.")]
    pub tags: PathBuf,

    #[arg(
        long,
        value_enum,
        default_value_t = Chemistry::BdV2_384,
        value_name = "CHEMISTRY",
        help = "Input single-cell chemistry. Supported: tenx-v2, tenx-v3, tenx-v4, bd-v1, bd-v2-96, bd-v2-384."
    )]
    pub chemistry: Chemistry,

    #[arg(
        long,
        default_value = "none",
        value_name = "none|human|mouse",
        help = "Optional BD sample-tag primer set for sample id detection."
    )]
    pub sample_species: String,

    #[arg(
        long,
        default_value_t = 9,
        value_name = "BP",
        help = "K-mer size used for BD sample id matching."
    )]
    pub sample_kmer_size: usize,

    #[arg(
        long,
        default_value_t = 6,
        value_name = "BP",
        help = "Expected BD UMI length."
    )]
    pub umi_len: usize,

    #[arg(
        long,
        default_value_t = 10,
        value_name = "N",
        help = "Minimum T bases required after UMI to accept trailing polyT."
    )]
    pub poly_t_min: usize,

    #[arg(
        long,
        default_value_t = 4,
        value_name = "N",
        help = "Worker threads for paired FASTQ processing."
    )]
    pub threads: usize,

    #[arg(
        long,
        default_value_t = 1,
        value_name = "0-9",
        help = "Gzip compression level for FASTQ output."
    )]
    pub gzip_level: u32,

    #[arg(
        long,
        default_value_t = false,
        help = "Write plain FASTQ instead of gzip-compressed FASTQ."
    )]
    pub no_gzip: bool,

}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}