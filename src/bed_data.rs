use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use bam::{BamReader, Record, Header, RecordReader};  // assuming bam crate for reading BAM files
use tokio::runtime::Runtime;
use std::thread::Builder;

use bigtools::BigWigWrite;
use bigtools::beddata::BedParserStreamingIterator;
use std::path::Path;


use crate::data_iter::DataIter;


/// Represents a single value in a bigWig file
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "write", derive(Serialize, Deserialize))]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl Value{
	pub fn flat(&self) -> ( u32, u32, f32) {
		( self.start, self.end, self.value )
	}
}


pub struct BedData {
    pub genome_info: Vec<(String, usize, usize)>, // (chromosome name, length, bin offset)
    pub search: HashMap<String, usize>, // get id for chr
    pub coverage_data: Vec<u32>, // coverage data array for bins
    pub bin_width: usize, // bin width for coverage
    pub threads: usize, // how many worker threads should we use here?
}

impl BedData {

    // Constructor to initialize a BedData instance
    pub fn new(bam_file: &str , bin_width: usize, threads:usize) -> Self {

    	let mut bam_reader = match BamReader::from_path(&bam_file, 1){
    		Ok(f) => f,
    		Err(e) => {
    			panic!("We hit an erro in Bam file reading: {e:?}");
    		}
    	};
	    let header = bam_reader.header().clone();
	    let genome_info:Vec<(String, usize, usize)> = Self::create_ref_id_to_name_map( &header,bin_width );
	    let search = Self::genome_info_to_search( &genome_info );

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
	            Ok(false) => {
	            	println!("bam file file read completely");
	               	break
	              },
	            Err(e) => panic!("{}", e),
	        }
	        let chrom_idx:usize = record.ref_id().try_into().unwrap();
	        let (chrom_name, chrom_length, chrom_offset) = &genome_info[chrom_idx];

	        let start:usize = record.start().try_into().unwrap();
	        let end = start + record.sequence().raw().len();

	        // Update bins
	        for pos in (start..end).step_by(bin_width) {
	        	//println!("Trying to add {chrom_name}:{pos} - 1");
	            if chrom_length >= &pos {
	                let index= chrom_offset+ pos / bin_width;
	                coverage_data[index] += 1;
	                //println!("Adding a one to the data in index {index}");
	            }else {
	            	//println!( "pos {pos} is outside of the chromsome?! {chrom_length}")
	            }
	            
	        }
	    }

        BedData {
            genome_info,
            search,
            coverage_data,
            bin_width,
            threads,
        }
    }

    pub fn genome_info_to_search(genome_info: &Vec< (String, usize, usize) > ) -> HashMap<String, usize> {
        genome_info
            .iter()
            .enumerate()
            .map(|(index, (name, _, _))| (name.clone(), index))
            .collect()
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

	pub fn id_for_chr_start( &self, chr:&str, start:usize ) -> Option<usize> {
		match self.search.get( chr ){
			Some(id) => {
				// Some ( offset + id)
				Some( self.genome_info[*id].2 + start / self.bin_width )
			},
			None => None
		}
	}

	/// get the chromosme name for and value id.
	/// This checks the offsets of all chromosomes and returns None if the id exeeds the range.
	pub fn current_chr_for_id(&self, id:usize) -> Option<( String, usize, usize)> {
		for ( chr, length, offset ) in &self.genome_info{
			if id >= *offset &&  (id - offset) * self.bin_width < *length {
				return Some( (chr.clone(), *length, *offset) );
			}
		}
		// the id did not map to any entry here!
		return None
	}




	/// Write the coverage data to a BEDGraph file.
	pub fn write_bedgraph(
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

	pub fn write_bigwig( &self, file: &str) -> Result<(),String>{	

		let runtime = tokio::runtime::Builder::new_multi_thread()
		     .worker_threads( self.threads )
		    .build()
		     .expect("Unable to create runtime.");

		let iter = DataIter::new( self );
		let data = BedParserStreamingIterator::wrap_infallible_iter(iter, true);

		let chrom_map: HashMap<String, u32> = self.genome_info.iter().map(|(chrom, len, _)| (chrom.clone(), *len as u32)).collect();

		let outfile = Path::new(file);

		// Create the BigWig writer
        let outb = BigWigWrite::create_file(outfile, chrom_map)
            .map_err(|e| format!("Failed to create BigWig file: {}", e))?;

        // Write data using the runtime
        outb.write(data, runtime)
            .map_err(|e| format!("Failed to write BigWig file: {}", e))?;

        Ok(())
	}


}