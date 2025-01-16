use bam_tide::bam_to_bedgraph;
use clap::Parser;
use clap::Command;

/// A command-line tool for converting BAM coverage data to BigWig format.
#[derive(Parser, Debug)]
#[clap(version = "0.0.1", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Args {
    /// Input BAM file (sorted by chromosome position).
    #[arg(short, long)]
    bam: String,

    /// Output BigWig file.
    #[arg(short, long)]
    bigwig: String,

    /// Bin width for coverage calculation (default: 50bp).
    #[arg(short, long, default_value_t = 50)]
    bin_width: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    println!("BamTide CLI");
    println!("Processing BAM: {}", args.bam);
    println!("Output BigWig: {}", args.bigwig);
    println!("Bin Width: {} bp", args.bin_width);

    bam_to_bedgraph(&args.bam, &args.bigwig, args.bin_width)?;

    println!("Conversion completed successfully.");
    Ok(())
}
