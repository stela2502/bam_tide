use rust_htslib::bam::{Reader, Read, Record};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs::{self, File, read};
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;

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

    if ! Path::new(&opts.bam).exists() {
        panic!("I need a bam file to work on!")
    }
    

    // Create output directory if needed
    if let Some(parent) = Path::new(&opts.fq).parent() {
        if ! parent.exists() {
        	if let Err(err) = fs::create_dir_all(&parent) {
	            panic!("Error creating directory {}: {}", &opts.fq, err);
	        } else {
	            println!("New output directory created successfully!");
	        }
        }
    }

    let mut reader = match Reader::from_path(&opts.bam) {
    	Ok(r) => r,
    	Err(e) => panic!("Error opening BAM file: {}", e),
    };

    // Gzipped FASTQ output files
    let f1 = File::create(&opts.fq).expect(&format!("Failed to create ofile {} output file", &opts.fq ));
    let buff_f1 = BufWriter::new( f1 );
    let mut fq1_writer = GzEncoder::new( buff_f1, Compression::default() );
    
    
    let f2 = File::create(&opts.fq2).expect(&format!("Failed to create ofile {} output file2", &opts.fq2 ));
    let buff_f2 = BufWriter::new(f2);
    let mut fq2_writer = GzEncoder::new( buff_f2, Compression::default() );

    let mut temp_storage = HashMap::<String, rust_htslib::bam::Record>::new();

    for r in reader.records() {
    	let mut record = r.expect("Bam read could not be collected");
        let qname = String::from_utf8_lossy(record.qname()).to_string();

        if let Some(mate) = temp_storage.remove(&qname) {
            // Write both mates to FASTQ
            writeln!(fq1_writer, "{}", to_string(&record) ).expect("Faled to write data to read1");
            writeln!(fq2_writer, "{}",to_string(&mate) ).expect("Faled to write data to read2");
        } else {
            // Store single-end read until its mate is found
            temp_storage.insert(qname, record.clone());
        }
    }

    println!("Paired-end FASTQ files written and gzipped.");

    if temp_storage.len() > 0 {
    	eprintln!("We have unpaired reads in the bam file (n={})- is it truncated?", temp_storage.len() )
    }
}

// Function to write a single BAM record in FASTQ format
fn to_string( record: &Record) -> String {
    let qname = String::from_utf8_lossy(record.qname());
    let seq = record.seq().as_bytes();
    let qual: String = record.qual()
	    .iter()
	    .map(|&q| (q + 33 ) as char)
	    .collect();

    format!( "@{}\n{}\n+\n{}", qname, String::from_utf8_lossy(&seq), qual)
}
