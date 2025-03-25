use rust_htslib::bam::{Reader,Read};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs::File;

use std::time::SystemTime;


use clap::{Parser};
#[derive(Parser)]
#[clap(version = "0.4.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,

    /// the gtf file fitting to the Bam file (text or gzipped)
    #[clap( long)]
    fq: String,

    /// the outpath
    #[clap( long)]
    fq2: String,
}


fn main(){
	let now = SystemTime::now();
    let opts: Opts = Opts::parse();

    match Path::new(&opts.fq)?.exists() {
    	Some(_) => {},
    	None => panic!("Missing bam file option!")
    }

    // Create output directory if needed
    if let Some(parent) = Path::new(&opts.fq).parent() {
        if ! parent.exists() {
        	if let Err(err) = fs::create_dir_all(&opts.outpath) {
	            panic!("Error creating directory {}: {}", &opts.outpath, err);
	        } else {
	            println!("New output directory created successfully!");
	        }
        }
    }

    let mut reader = match Reader::from_path(bam_file) {
    	Ok(r) => r,
    	Err(e) => panic!("Error opening BAM file: {}", e),
    };

    // Gzipped FASTQ output files
    let fq1_writer = GzEncoder::new(
    	File::create(fq1_path).expect(&format!("Failed to create ofile {} output file", &opts.fq )),
    	Compression::default()
    );

    let fq2_writer = GzEncoder::new(
    	File::create(fq2_path).expect(&format!("Failed to create ofile {} output file2", &opts.fq2 )), 
    	Compression::default()
    );
    let mut fq1_writer = BufWriter::new(fq1_writer);
    let mut fq2_writer = BufWriter::new(fq2_writer);

    let temp_storage = HashMap::<String, rust_htslib::bam::Record>::new();

    for r in reader.records() {
    	let mut record = read?;
        let qname = String::from_utf8_lossy(record.qname()).to_string();

        if let Some(mate) = temp_storage.remove(&qname) {
            // Write both mates to FASTQ
            write_fastq(&mut fq1_writer, &mate)?;
            write_fastq(&mut fq2_writer, &record)?;
        } else {
            // Store single-end read until its mate is found
            temp_storage.insert(qname, record.clone());
        }
    }

    println!("Paired-end FASTQ files written and gzipped.");

    if temp_storage.len() > 0 {
    	eprintln!("We have unpaired reads in the bam file (n={}})- is it truncated?", teml_storage.len() )
    }
    Ok(())
}

// Function to write a single BAM record in FASTQ format
fn write_fastq<W: Write>(writer: &mut W, record: &Record) -> std::io::Result<()> {
    let qname = String::from_utf8_lossy(record.qname());
    let seq = record.seq().as_bytes();
    let qual: String = record.qual()
	    .iter()
	    .map(|&q| (q + 33 ) as char)
	    .collect();

    writeln!(writer, "@{}\n{}\n+\n{}", qname, String::from_utf8_lossy(&seq), qual)
}