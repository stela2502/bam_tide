use std::collections::HashMap;
use std::fs::File;
use std::io::{Write};
//use tokio::runtime::Runtime;
//use std::thread::Builder;

use bigtools::BigWigWrite;
use bigtools::beddata::BedParserStreamingIterator;
use std::path::Path;

use rust_htslib::bam::{Reader,Read};
use rust_htslib::bam::Header;

use crate::gtf_logics::{ create_ref_id_to_name_hashmap, AnalysisType };
use crate::read_data::ReadData;
use crate::data_iter::DataIter;

use rustody::genes_mapper::cigar::Cigar;

//const BUFFER_SIZE: usize = 1_000_000;
use clap::ValueEnum;


/// Represents a single value in a bigWig file
#[derive(Copy, Clone, Debug, PartialEq)]
//#[cfg_attr(feature = "write", derive(Serialize, Deserialize))]
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


/// The 'supported' normalization options
#[derive(ValueEnum, Clone, Debug)]
pub enum Normalize {
    Not,
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
    pub fn new(bam_file: &str , bin_width: usize, threads:usize, analysis_type: &AnalysisType, cell_tag: &[u8;2], umi_tag: &[u8;2] ) -> Self {


	    let mut reader = match Reader::from_path(bam_file) {
	    	Ok(r) => r,
	    	Err(e) => panic!("Error opening BAM file: {}", e),
	    };

	    let cigar = Cigar::new( "" );

	    let header = Header::from_template(reader.header());
    	// (chr name, length, offset)
	    let genome_info:Vec<(String, usize, usize)> = Self::create_ref_id_to_name_vec( &header,bin_width );

	    // i32 id to chr name
	    let ref_id_to_name = create_ref_id_to_name_hashmap(  &header );

	    let search = Self::genome_info_to_search( &genome_info );

	    let mut singlets = HashMap::<String, ReadData>::new();

	    let num_bins = genome_info
            .iter()
            .map(|(_, length, _)| (length + bin_width - 1) / bin_width)
            .sum::<usize>();

	    let mut coverage_data = vec![0u32; num_bins];

	    let mut lines = 0;

	    // Process BAM records
	    for r in reader.records() {
	        // Read a record from BAM file
	        let record = match r {
	            Ok(r) => r,
	            Err(e) => panic!("I could not collect a read: {e:?}"),
	        };
			// Choose the correct function to extract the data based on the AnalysisType
	        let data_tuple = match *analysis_type {
	            AnalysisType::SingleCell => ReadData::from_single_cell(&record, &ref_id_to_name, &cell_tag, &umi_tag),
	            AnalysisType::Bulk => ReadData::from_bulk(&record, &ref_id_to_name, &umi_tag, lines, "1"),
	        };

	        // Here we really only need start and end at the moment.
	        let region = match data_tuple {
	            Ok(ref res) => {
	                let qname = &res.0; // Cell ID or read name as key
	                let read_data = &res.1;

	                if read_data.is("paired") {
	                    match singlets.remove(qname) {
	                        Some(first_read) => {
	                            // Mate found! Process the pair
	                            //println!("Found paired reads: {:?} <-> {:?}", first_read, read_data);
	                            if first_read.chromosome != read_data.chromosome {
	                            	eprintln!( "chr missmatch!\n{}\n{}",first_read, read_data );
	                            	None
	                            }
	                            // add a +1 to the whole region.
	                            else if first_read.start < read_data.start {
	                            	Some((first_read.start, first_read.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ), )
	                            }else {
	                            	Some((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &first_read.cigar ).1 as i32 ))
	                            }
	                        }
	                        None => {
	                            // Handle orphaned reads (mate unmapped)
	                            if read_data.is("mate_unmapped") {
	                                //println!("Orphaned read (mate unmapped): {:?}", read_data);
	                                // Process it as a single read
	                                Some((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ))
	                            } else {
	                                // No mate found, store this read for future pairing
	                                singlets.insert(qname.to_string(), read_data.clone());
	                                //println!("Storing read for future pairing: {:?}", read_data);
	                                None
	                            }
	                        }
	                    }
	                } else {
	                    // Unpaired read (not part of a pair at all)
	                    //println!("Unpaired read: {}", read_data);
	                    Some ((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ))
	                }
	            }
	            Err("missing_Chromosome") => {
	                //eprintln!("Missing chromosome for BAM entry - assuming end of usable data.\n{:?}", record);
	                None
	            }
	            Err(_err) => {
	                //mapping_info.report(err);
	                None
	            }
	        };

		    /*if lines % BUFFER_SIZE == 0 {
		        pb.set_message(format!("{} million reads processed", lines / BUFFER_SIZE));
		        pb.inc(1);
		    }*/

		    if let Some( (start, end) ) = region {

		    	let (chrom_name, chrom_length, chrom_offset) = &genome_info[ record.tid()as usize ];
				//println!("I am filling {}:{}-{} with +1",chrom_name, start, end );
			    // Update bins
		        for pos in (start ..end).step_by(bin_width) {
		        	//println!("Trying to add {chrom_name}:{pos} - 1");
		            if chrom_length >= &(pos as usize) {
		                let index= chrom_offset+ pos as usize / bin_width;
		                coverage_data[index] += 1;
		                //println!("Adding a one to the data in index {index}");
		            }else {
		            	panic!( "pos {pos} is outside of the chromsome?! {chrom_length}")
		            }
		            
		        }
		    }
		    
	 		lines += 1;

	    }

        BedData {
            genome_info,
            search,
            coverage_data,
            bin_width,
            threads,
        }
    }

    /// Here we have implemented all the normalization features.
    pub fn normalize(&mut self, by: &Normalize ) {
    	match by {
    		Normalize::Not => {
    			// do just notthing ;-)
    		},
    	}
    }

    /// Get the hashMap to locate the position of the chromosome in the annotation table
    pub fn genome_info_to_search(genome_info: &Vec< (String, usize, usize) > ) -> HashMap<String, usize> {
        genome_info
            .iter()
            .enumerate()
            .map(|(index, (name, _, _))| (name.clone(), index))
            .collect()
    }

    /// Create a mapping of reference ID to reference name, chromosome length, and bin offset.
	/// `bin_width` is the size of each bin.
	/// returns a vector with (chromsome name, chr length, chr offset - for this chromosome )
	/// later called "annotation table"
	pub fn create_ref_id_to_name_vec( header: &Header, bin_width: usize) -> Vec<(String, usize, usize)> {

		let header_map = header.to_hashmap();
		let mut ref_id_to_name = Vec::with_capacity(30);
		let mut total_bins = 0; 

		if let Some(reference_info) = header_map.get("SQ") {
	        for record in reference_info {
	            if let Some(length) = record.get("LN") {
	                // The "LN" tag holds the length of the chromosome
	                if let Ok(length_int) = length.parse::<usize>() {
	                    if let Some(name) = record.get("SN") {
	                    	// Calculate number of bins for this chromosome
	                    	let bins = (length_int as usize + bin_width - 1) / bin_width;

	                    	// Create the tuple: (chromosome name, length, bin offset)
	            			let result = (name.to_string(), length_int as usize, total_bins);

	            			total_bins += bins;

	                        ref_id_to_name.push(result);
	                    }
	                }
	            }
	        }
	    }
	    ref_id_to_name
	}

	/// an example of how to calculate the id from the chromosome + position
	pub fn id_for_chr_start( &self, chr:&str, start:usize ) -> Option<usize> {
		match self.search.get( chr ){
			Some(id) => {
				// Some ( offset + id)
				Some( self.genome_info[*id].2 + start / self.bin_width )
			},
			None => None
		}
	}

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
	    let mut bin_start: usize; // Start position of each bin
	    for (_chrom_index, (chrom_name, chrom_length, bin_offset)) in self.genome_info.iter().enumerate() {
	    	#[allow(unused_assignments)]
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