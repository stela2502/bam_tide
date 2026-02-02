//bam2bigwig.rs
use clap::Parser;
use bam_tide::cli::CoverageCli;

fn main() {
    let opts = CoverageCli::parse();
    // next: bam_tide::coverage::run_bigwig(opts)
    println!("{opts:?}");
}
