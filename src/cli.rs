//cli.rs
// src/cli.rs
use clap::{Parser, ValueEnum};
use crate::bed_data::Normalize;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum AnalysisType {
    Bulk,
    #[value(name = "single-cell")]
    SingleCell,
}


/// Shared CLI options for coverage exporters (bedGraph / bigWig)
#[derive(Parser, Debug, Clone)]
pub struct CoverageCli {
    /// Input BAM file (sorted by chromosome position)
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Output file (bedGraph or BigWig depending on binary)
    #[arg(short = 'o', long)]
    pub outfile: String,

    /// tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
    #[arg(short = 'c', long, default_value = "CB")]
    pub cell_tag: String,

    /// tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
    #[arg(short = 'u', long, default_value = "UB")]
    pub umi_tag: String,

    /// Collect single cell info or bulk
    #[arg(short = 'a', long, value_enum, default_value_t = AnalysisType::Bulk)]
    pub analysis_type: AnalysisType,

    /// Normalize the data somehow
    #[arg(short = 'n', long, value_enum, default_value_t = Normalize::Not)]
    pub normalize: Normalize,

    /// Bin width for coverage calculation
    #[arg(short = 'w', long, default_value_t = 50)]
    pub width: u32,

    /// Collect only R1 areas
    #[arg(long)]
    pub only_r1: bool,

    /// Minimum mapping quality to include a read
    #[arg(long, default_value_t = 0)]
    pub min_mapping_quality: u8,

    /// Include secondary alignments
    #[arg(long, default_value_t = false)]
    pub include_secondary: bool,

    /// Include supplementary alignments
    #[arg(long, default_value_t = false)]
    pub include_supplementary: bool,

    /// Include duplicate-marked reads
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,

    /// Exclude reads with ANY of these SAM flag bits set (deeptools --samFlagExclude).
    /// Example: 2816 = secondary(256) + QC-fail(512) + supplementary(2048)
    #[arg(long)]
    pub sam_flag_exclude: Option<u16>,


}
