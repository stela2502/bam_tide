use rust_htslib::bam::{Reader, Read, Record};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs::{self, File, read};
use std::path::Path;
use flate2::write::GzEncoder;
use flate2::Compression;
use bam_tide::read_data::ReadData;

use std::time::SystemTime;

/// Used to regenerate the R1 and R2 fastq files from a CellRanger bam file.
use clap::{Parser};
#[derive(Parser)]
#[clap(version = "0.4.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,

    /// the prefix for the R1 R2 and I1 fastq.gz files
    #[clap(short, long)]
    prefix: String,

}


fn main(){
	let now = SystemTime::now();
    let opts: Opts = Opts::parse();

    panic!("Do not use this program! use https://github.com/10XGenomics/bamtofastq\nIf you are on COSMOS-SENS bamtofastq is installed in the Rustody/1.4 module.");

    if ! Path::new(&opts.bam).exists() {
        panic!("I need a bam file to work on!")
    }
    

    // Create output directory if needed
    if let Some(parent) = Path::new(&opts.prefix).parent() {
        if ! parent.exists() {
        	if let Err(err) = fs::create_dir_all(&parent) {
	            panic!("Error creating directory {}: {}", &opts.prefix, err);
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
    // CellRanger whnts this file name structure: SampleName_S1_L001_R[12]_001.fastq.gz
    let f1_str = format!("{}_{}", &opts.prefix, "S1_L001_R1_001.fastq.gz");
    let f1 = File::create(&f1_str).expect(&format!("Failed to create ofile {} output file", &f1_str ));
    let buff_f1 = BufWriter::new( f1 );
    let mut fq1_writer = GzEncoder::new( buff_f1, Compression::default() );
    
    let f2_str = format!("{}_{}", &opts.prefix, "S1_L001_R2_001.fastq.gz");
    let f2 = File::create(&f2_str).expect(&format!("Failed to create ofile {} output file", &f2_str ));
    let buff_f2 = BufWriter::new(f2);
    let mut fq2_writer = GzEncoder::new( buff_f2, Compression::default() );

    /*
    let i1_str = format!("{}_{}", &opts.prefix, "S1_L001_I1_001.fastq.gz");
    let i1 = File::create(&i1_str).expect(&format!("Failed to create ofile {} output file", &i1_str ));
    let buff_i1 = BufWriter::new(i1);
    let mut fi1_writer = GzEncoder::new( buff_i1, Compression::default() );
    */

    for r in reader.records() {
    	let mut record = r.expect("Bam read could not be collected");
        let qname = String::from_utf8_lossy(record.qname()).to_string();
        let qual = &record.qual()
	    .iter()
	    .map(|&q| (q + 33 ))
	    .collect::<Vec<u8>>(); 
        writeln!(fq1_writer, "{}", as_fastq(&qname, &record.seq().as_bytes(), qual ) ).expect("Faled to write data to read1");

        writeln!(fq2_writer, "{}", as_fastq(
            &qname, 
            get_tag_value( &record, b"CR" ).unwrap_or_else(|| panic!("CR tag not found?! {:?}", &qname ) ).as_bytes(),
            get_tag_value( &record, b"CY" ).unwrap_or_else(|| panic!("CY tag not found!? {:?}", &qname ) ).as_bytes()
            /*
            ReadData::get_tag_value( &record, b"CR" ).unwrap_or( continue ).as_bytes(),
            ReadData::get_tag_value( &record, b"CY" ).unwrap_or( continue ).as_bytes() 
            */
            ) 
        ).expect("Faled to write data to read2");

        /*
        // The I1 fastq file does only contain the sample sequencing tag and is not necessary any more.
        // The data for that is also not sotored in every bam file any more.
        // I can not regenerate that here.
        writeln!(fi1_writer, "{}", as_fastq(
            &qname, 
            get_tag_value( &record, b"UR" ).unwrap_or_else(|| panic!("UR tag not found?! {:?}", &qname )  ).as_bytes(),
            get_tag_value( &record, b"UY" ).unwrap_or_else(|| panic!("UY tag not found? {:?}", &qname ) ).as_bytes()
            /*
            ReadData::get_tag_value( &record, b"BC" ).unwrap_or( continue ).as_bytes(),
            ReadData::get_tag_value( &record, b"QT" ).unwrap_or( continue ).as_bytes() 
            */
            ) 
        ).expect("Faled to write data to read2");
        */
    }

    println!("R1/R2/I1 triplet FASTQ files written and gzipped.");

}

// Function to write a single BAM record in FASTQ format
fn as_fastq( qname:&str, seq:&[u8], qual:&[u8] ) -> String {
    format!( "@{}\n{}\n+\n{}", qname, String::from_utf8_lossy(&seq), String::from_utf8_lossy(qual) )
}


fn get_tag_value(record: &Record, tag: &[u8; 2]) -> Option<String> {
        if let Some(value) = record.aux(tag).ok() {
            match value {
                rust_htslib::bam::record::Aux::String(s) => Some(s.to_string()),
                _ => None,
            }
        } else {
            None
        }
    }
