use clap::Parser;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write, BufRead};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::path::PathBuf;

use bam_tide::ATACtoRNAMapper; // adjust based on your project structure

use rustody::int_to_str::IntToStr;

#[derive(Parser)]
#[command(name = "translate_barcodes")]
#[command(version = "1.0")]
#[command(about = "Translates 10x barcodes.tsv.gz using HashMaps created by generate_atac_to_rna_hashmap.", long_about = None)]
struct Cli {
    /// Optional, a path to the ATAC-to-RNA translation binary file
    #[arg(short, long)]
    map: Option<PathBuf>,

    /// Path to input barcodes.tsv.gz
    #[arg(short, long)]
    input: PathBuf,

    /// Path to output translated barcodes.tsv.gz
    #[arg(short, long)]
    output: PathBuf,
}



fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    let map = match &args.map{
        Some(path) => path,
        None => &PathBuf::new( ),
    };
    let mapper = ATACtoRNAMapper::from_binary_file(&map)
        .map_err(|e| format!("Failed to load mapping: {}\nThe file is available in the github repo!", e))?;

    let input = File::open(&args.input)?;
    let decoder = GzDecoder::new(input);
    let reader = BufReader::new(decoder);

    let output = File::create(&args.output)?;
    let encoder = GzEncoder::new(output, Compression::default());
    let mut writer = BufWriter::new(encoder);

    let mut i2s = IntToStr::new(b"AAA".to_vec(), 32).unwrap();

    for line in reader.lines() {
        let barcode = line?;
        let atac_id = barcode_to_u32(&barcode, &mut i2s); 

        if let Some(rna_id) = mapper.translate(atac_id) {
            let rna_barcode = u32_to_barcode(rna_id, &mut i2s);
            writeln!(writer, "{}", rna_barcode)?;
        } else {
            eprintln!("Warning: no mapping found for barcode '{}'", barcode);
            writeln!(writer, "{}_unmapped", barcode)?; // fallback: write original + note
        }
    }

    Ok(())
}

// These are placeholder signatures; you will implement the actual logic
fn barcode_to_u32(barcode: &str, i2s: &mut IntToStr) -> u32 {
    // Your encoding logic
    i2s.str_to_u32( barcode )
}

fn u32_to_barcode(encoded: u32, i2s: &mut IntToStr) -> String {
    // Your decoding logic
    i2s.u64_to_string(16, &(encoded.into()) )
}
