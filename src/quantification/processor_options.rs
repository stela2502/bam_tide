use crate::quantification::cli::{QuantCli};
use crate::quantification::QuantMode;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct ProcessorOptions {
    pub min_mapq: u8,
    pub read1_only: bool,
    pub require_strand: bool,
    pub quant_mode: QuantMode,
    pub read_tag_table: Option<PathBuf>,
    pub rt_cell_column: String,
    pub rt_cell_qual_column: String,
    pub rt_umi_column: String,
    pub rt_umi_qual_column: String,
}

impl From<&QuantCli> for ProcessorOptions {
    fn from(args: &QuantCli) -> Self {
        Self {
            min_mapq: args.min_mapq,
            read1_only: args.read1_only,
            require_strand: args.require_strand,
            quant_mode: args.quant_mode,
            read_tag_table: args.read_tags.read_tag_table.clone(),
            rt_cell_column: args.read_tags.rt_cell_column.clone(),
            rt_cell_qual_column: args.read_tags.rt_cell_qual_column.clone(),
            rt_umi_column: args.read_tags.rt_umi_column.clone(),
            rt_umi_qual_column: args.read_tags.rt_umi_qual_column.clone(),
        }
    }
}