//bam-quant-testdata.rs
use std::fs;
use std::process::ExitCode;

use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{MatchClass, MatchOptions, SpliceIndex, SplicedRead, Strand};
use snp_index::{Genome, SnpIndex};

use bam_tide::encoder::{SamEncoder, TestDataCli};
use clap::Parser;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let cli = TestDataCli::parse();

    let outfile = cli.outfile.clone();

    let mut encoder = SamEncoder::new(cli)?;
    let lines = encoder.generate()?;

    std::fs::write(&outfile, lines.join("\n") + "\n")
        .map_err(|e| format!("failed to write SAM {:?}: {e}", outfile))?;

    Ok(())
}
