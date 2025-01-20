use bam_tide::bed_data::BedData;
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
    outfile: String,

    /// Bin width for coverage calculation (default: 50bp).
    #[arg(short, long, default_value_t = 50)]
    width: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    println!("BamTide CLI");
    println!("Processing BAM: {}", args.bam);
    println!("Output BigWig: {}", args.outfile);
    println!("Bin Width: {} bp", args.width);

    let bed_data = BedData::new( &args.bam, args.width, 1 );

    match bed_data.write_bedgraph(&args.outfile, ){
        Ok(_) => {
            println!("Conversion completed successfully.");
        },
        Err(e) => {
            println!("Some error occured - check if that is a problem: {e:?}")
        }
    }

    
    Ok(())
}
