use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use bam::{BamReader, Record, Header, RecordReader};  // assuming bam crate for reading BAM files

use bigtools::{BigWigWrite, ChromInfo};


pub struct BedData {
    genome_info: Vec<(String, usize, usize)>, // (chromosome name, length, bin offset)
    coverage_data: Vec<u32>, // coverage data array for bins
    bin_width: usize, // bin width for coverage
}

impl BedData {

    // Constructor to initialize a BedData instance
    pub fn new(bam_file: &str , bin_width: usize) -> Self {

    	let mut bam_reader = BamReader::from_path(&bam_file, 1).unwrap();
	    let header = bam_reader.header().to_owned();
	    let genome_info:Vec<(String, usize, usize)> = Self::create_ref_id_to_name_map( &header,bin_width );

	    let num_bins = genome_info
            .iter()
            .map(|(_, length, _)| (length + bin_width - 1) / bin_width)
            .sum::<usize>();

	    let mut coverage_data = vec![0u32; num_bins];

	    let mut record = bam::Record::new();
	    let mut last_start = 0; 

	    // Process BAM records
	    loop {
	        match bam_reader.read_into(&mut record) {
	            Ok(true) => {},
	            Ok(false) => break,
	            Err(e) => panic!("{}", e),
	        }
	        let chrom_idx:usize = record.ref_id().try_into().unwrap();
	        let (chrom_name, chrom_length, chrom_offset) = &genome_info[chrom_idx];

	        let start:usize = record.start().try_into().unwrap();
	        let end = start + record.sequence().raw().len();
	        println!("Trying start to end {start} to {end} length {}", {end - start});

	        // Update bins
	        for pos in (start..end).step_by(bin_width) {
	            if chrom_length < &pos {
	                let index= chrom_offset+ pos / bin_width;
	                coverage_data[index] += 1;
	            }
	            
	        }
	    }

        BedData {
            genome_info,
            coverage_data,
            bin_width,
        }
    }

    /// Create a mapping of reference ID to reference name, chromosome length, and bin offset.
	/// `bin_width` is the size of each bin.
	fn create_ref_id_to_name_map( header_view: &Header, bin_width: usize) -> Vec<(String, usize, usize)> {
	    let reference_names = header_view.reference_names();
	    let reference_lengths = header_view.reference_lengths();

	    let mut total_bins = 0; // This will store the cumulative bin count (offset)

	    // Zip the reference names and lengths together and compute the bin offsets
	    reference_names
	        .iter()
	        .zip(reference_lengths.iter())
	        .map(|(name, &length)| {
	            // Calculate number of bins for this chromosome
	            let bins = (length as usize + bin_width - 1) / bin_width;

	            // Create the tuple: (chromosome name, length, bin offset)
	            let result = (name.to_string(), length as usize, total_bins);

	            // Update the total_bins to include this chromosome's bins
	            total_bins += bins;

	            result
	        })
	        .collect()
	}


	/// Write the coverage data to a BEDGraph file.
	fn write_bedgraph(
		&self,
	    file_path: &str, // Output file path
	) -> std::io::Result<()> {
	    // Open the output file
	    let mut file = File::create(file_path)?;

	    // Iterate over the genome info and coverage data
	    let mut bin_start: usize = 0; // Start position of each bin
	    for (chrom_index, (chrom_name, chrom_length, bin_offset)) in self.genome_info.iter().enumerate() {
	        let mut bin_end = 0; // End position of each bin
	        
	        // Iterate through the coverage data for the current chromosome
	        for bin_index in 0..(chrom_length / self.bin_width) {
	            // Calculate start and end positions for each bin
	            bin_start = bin_index * self.bin_width;
	            bin_end = (bin_index + 1) * self.bin_width;

	            if bin_end > *chrom_length {
	                bin_end = *chrom_length; // Adjust for last bin that may not fill up to full width
	            }

	            // Get the coverage value for this bin (you'll need to map the correct coverage)
	            let coverage_value = self.coverage_data.get(bin_offset + bin_index).unwrap_or(&0);

	            // Write to the BedGraph file
	            if *coverage_value > 0 {
	                writeln!(
	                    file,
	                    "{}\t{}\t{}\t{}",
	                    chrom_name,   // Chromosome name
	                    bin_start,    // Start position (0-based)
	                    bin_end,      // End position (1-based)
	                    coverage_value // Coverage value
	                ).unwrap();
	            }
	            
	        }
	    }

	    Ok(())
	}


	fn write_bigwig_data(
		&self,
	    file: &mut File,
	) -> Result<(), io::Error> {
	    // Iterate over the genome info and coverage data for each chromosome
	    let mut bin_start: usize = 0;
	    let mut bin_end: usize = 0;
	    let chr_sizes : HashMap<String, u32> =  self.genome_info
	    	.iter()
	    	.map(|(name, length, _)| (name.clone(), *length as u32)) // Clone the name and cast the length to u32
    		.collect(); // Collect the result into a HashMap
		let bw_writer = BigWigWrite::new(file, chr_sizes )?;

	    for (chrom_index, (chrom_name, chrom_length, bin_offset)) in self.genome_info.iter().enumerate() {
	        // Create a vector to store data for this chromosome
	        let mut chrom_data: Vec<u32> = Vec::new();
	        
	        for bin_index in 0..(chrom_length / self.bin_width) {
	            bin_start = bin_index * self.bin_width;
	            bin_end = (bin_index + 1) * self.bin_width;

	            if bin_end > *chrom_length {
	                bin_end = *chrom_length;
	            }

	            // Get the coverage value for this bin
	            let coverage_value = self.coverage_data.get(bin_offset + bin_index).unwrap_or(&0);

	            // Store the coverage value (you may want to use a different data type for BigWig)
	            chrom_data.push(*coverage_value);
	        }

	        // Write the chromosome's data into BigWig format

	        

	        bw_writer.write_data(chrom_name, &chrom_data)?;

	        // Reset the vector for the next chromosome
	        chrom_data.clear();
	    }

	    Ok(())
	}


	fn write_bigwig_header(
		&self,
	    file: &mut File,
	) -> Result<(), io::Error > {
	    // Create a list of ChromInfo to hold chromosome names and their sizes
	    let chrom_info: Vec<ChromInfo> = self.genome_info.iter().map(|(name, length, _)| {
	        ChromInfo {
	            name: name.clone(),
	            length: *length as u32,  // Length in u32 (as required for BigWig format)
	        }
	    }).collect();

	    // Create the BigWig writer with the header
	    let bw_writer = BigWigWrite::new(file, &chrom_info, self.bin_width as u32)?;

	    // Write the header (metadata about chromosomes and bins)
	    bw_writer.write_header()?;

	    Ok(())
	}

	pub fn write_bigwig(&self, bigwig_file:&str ) -> Result<(),  io::Error>{
		let mut file = File::create( bigwig_file )?;
	    // Write the header
	    self.write_bigwig_header( &mut file )?;

	    // Write the coverage data
	    self.write_bigwig_data(&mut file )?;
	}

	pub fn bam_to_bigwig(
	    bam_file: &str,
	    bigwig_file: &str,
	    bin_width: usize,
	) -> Result<(),  io::Error> {
	    // Open BAM file and create a reader
	    let me = Self::new( bam_file, bin_width );
	    me.write_bigwig( bigwig_file );


	    Ok(())
	}
}