use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Clone, Parser)]
#[command(author, version, about = "Normalize messy ONT/Dorado BAM reads into one 10x-style molecule per FASTQ read")]
pub struct Cli {
    #[arg(long)]
    pub bam: PathBuf,

    #[arg(long)]
    pub out: PathBuf,

    #[arg(long)]
    pub tags: PathBuf,

    #[arg(long, default_value = "CTACACGACGCTCTTCCGATCT")]
    pub adapter: String,

    #[arg(long, default_value_t = 13)]
    pub min_adapter_match: usize,

    #[arg(long, default_value_t = 16)]
    pub cb_len: usize,

    #[arg(long, default_value_t = 12)]
    pub umi_len: usize,

    #[arg(long, default_value_t = 10)]
    pub poly_t_min: usize,

    #[arg(long, default_value_t = 14)]
    pub poly_t_window: usize,

    #[arg(long, default_value_t = 4)]
    pub threads: usize,

    #[arg(long, default_value_t = 1)]
    pub gzip_level: u32,

    #[arg(long, default_value_t = false)]
    pub no_gzip: bool,

    #[arg(long, default_value_t = 20)]
    pub min_transcript_len: usize,

    #[arg(long, default_value_t = 2)]
    pub max_adapter_mismatches: usize,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
