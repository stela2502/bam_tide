#![allow(clippy::too_many_arguments)]

pub mod bed_data;
pub mod cli;
pub mod compare_report;
pub mod core;
pub mod data_iter;
pub mod index;
pub mod ont_normalizer;
pub mod results;
pub mod multi_subset_bam;

pub mod encoder;
pub mod quantification;

pub use results::QuantData;
