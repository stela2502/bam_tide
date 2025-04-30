use std::collections::HashMap;
use std::fs::File;
use std::io::{Write};
//use tokio::runtime::Runtime;
//use std::thread::Builder;

use bigtools::BigWigWrite;
use bigtools::beddata::BedParserStreamingIterator;

use std::path::Path;
use std::collections::HashSet;

use rust_htslib::bam::{Reader,Read};
use rust_htslib::bam::Header;

use crate::data_iter::DataIter;
use crate::feature_matcher::*;
//use crate::gtf::{gene::Gene, gene::RegionStatus};
use crate::gtf::{ExonIterator,SplicedRead, exon_iterator::ReadResult};
use crate::gtf_logics::{ create_ref_id_to_name_hashmap, AnalysisType, MatchType };
use crate::mutation_processor::MutationProcessor;
use crate::read_data::ReadData;
use crate::bed_data::ChrArea;

use rustody::genes_mapper::cigar::Cigar;
use rustody::mapping_info::MappingInfo;
use rustody::singlecelldata::{SingleCellData, IndexedGenes, cell_data::GeneUmiHash};


//const BUFFER_SIZE: usize = 1_000_000;
use clap::ValueEnum;

use rayon::prelude::*;

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
	// just collect - no normalization whatsoever
	Not,
    // Reads Per Kilobase per Million mapped reads;
    // number of reads per bin / (number of mapped reads (in millions) * bin length (kb))
    Rpkm,
    // Counts Per Million mapped reads
    // number of reads per bin / number of mapped reads (in millions)
    Cpm,
    // Bins Per Million mapped reads
    // number of reads per bin / sum of all reads per bin (in millions)
    Bpm,
    // reads per genomic content (1x normalization)
    // number of reads per bin / scaling factor for 1x average coverage.
    // This scaling factor, in turn, is determined from the sequencing depth: 
    // (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of
    // the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize.
    Rpgc,
}


pub struct BedData {
    pub genome_info: Vec<(String, usize, usize)>, // (chromosome name, length, bin offset)
    pub search: HashMap<String, usize>, // get id for chr
    pub coverage_data: Vec<f32>, // coverage data array for bins
    pub bin_width: usize, // bin width for coverage
    pub threads: usize, // how many worker threads should we use here?
    pub nreads: usize,
}


impl FeatureMatcher for BedData{

	fn extract_gene_ids(
		&self,
		read_result: &Option<Vec<ReadResult>>,
		_data: &ReadData,
		_mapping_info: &mut MappingInfo,
		) -> Vec<String> {
		if let Some(results) = read_result {
			results
			.iter()
			.map(|result| result.gene.clone() )
			.collect()
		}else {
			vec![]
		}
	}

	fn init_search(
		&self,
		chr: &str,
		_start: usize,
		_iterator: &mut ExonIterator,
		) -> Result<(), QueryErrors>{
		match self.search.get(chr) {
			Some(_) => Ok(()),
			None => Err( QueryErrors::ChrNotFound )
		}
	}

	fn process_feature(
		&self,
		data: &(ReadData, Option<ReadData>),
		mutations: &Option<MutationProcessor>,
		_iterator: &mut ExonIterator,
		exp_gex: &mut SingleCellData,
		exp_idx: &mut IndexedGenes,
		mut_gex: &mut SingleCellData,
		mut_idx: &mut IndexedGenes,
		mapping_info: &mut MappingInfo,
		_match_type: &MatchType,
		){
    	let primary_read = &data.0; // Always present
        let mate_read = data.1.as_ref(); // Optional paired read

        // Parse cell ID from the primary read
        let cell_id = match Self::parse_cell_id( &primary_read.cell_id ) {
        	Ok(id) => id,
        	Err(_) => return,
        };

        let mut guhs = HashSet::<GeneUmiHash>::new();
        
        if let Some(mate) = mate_read {
        	if primary_read.chromosome != mate.chromosome {
        		return // TODO treat translocations?!
        	}

        	self.collect_ghus( &mate, &mut guhs, exp_idx );
        	if let Some(processor) = mutations {
	             processor.handle_mutations( primary_read, "unknown", mut_idx, mut_gex, mapping_info, &cell_id, self.collect_gene_areas(primary_read) );
	        }
        }
        self.collect_ghus( &primary_read, &mut guhs, exp_idx);
        for guh in guhs.iter(){
        	//println!("Inserting GHU {}",guh);
        	exp_gex.try_insert( &cell_id, *guh, mapping_info );
        }
        if let Some(processor) = mutations {
             processor.handle_mutations( primary_read, "unknown", mut_idx, mut_gex, mapping_info, &cell_id, self.collect_gene_areas(primary_read) );
        }
    }

}

impl BedData {

	fn collect_gene_areas( &self, read:&ReadData ) -> Option<HashSet<ChrArea>>{

		let mut gene_areas = HashSet::<ChrArea>::new();
		let ( _chr, _length, _offset) = &self.genome_info[ *self.search.get(&read.chromosome).unwrap() ];
		let (spliced_read, _final_position) = SplicedRead::new_mut( &read.cigar, read.start.try_into().unwrap() );
		for exon in spliced_read.exons {
			for pos in (exon.start..exon.end).step_by(self.bin_width) {
				let id = pos / self.bin_width;
				let gene_area = ChrArea{
					start:  id*self.bin_width,
					end: (id+1) * self.bin_width,
					chr: read.chromosome.to_string(),
				};
				gene_areas.insert( gene_area );
			}
		}
		Some(gene_areas)
	}

	fn collect_ghus( &self, read: &ReadData, ghus: &mut HashSet<GeneUmiHash>, idx: &mut IndexedGenes, ) {

		let (_chr, _length, _offset) = &self.genome_info[ *self.search.get(&read.chromosome).unwrap() ];

		let (spliced_read, _final_position) = SplicedRead::new_mut( &read.cigar, read.start.try_into().unwrap() );
		//println!("found a spliced read: {:?}", spliced_read);
		for exon in spliced_read.exons {
			//println!("found an exon {}..{}", exon.start, exon.end );
			for pos in (exon.start..exon.end).step_by(self.bin_width) {

				let id = pos / self.bin_width;
                let gene= format!( "{}:{}-{}", read.chromosome, id*self.bin_width, (id+1) * self.bin_width );
                //let intern_id = id + offset;
                //println!("This should be the matching bed area: {}", &gene);
				let guh = GeneUmiHash( idx.get_gene_id( &gene ) , read.umi );
				
				if ghus.insert( guh ) {
					//self.coverage_data[intern_id] +=1.0;
				}
				
			}
		}
	}

	pub fn init( bam_file: &str, bin_width:usize, threads:usize )-> Self{

		let reader = match Reader::from_path(bam_file) {
    		Ok(r) => r,
    		Err(e) => panic!("Error opening BAM file: {}", e),
    	};

    	let header = Header::from_template(reader.header());
    	// (chr name, length, offset)

		let genome_info:Vec<(String, usize, usize)> = BedData::create_ref_id_to_name_vec( &header,bin_width );

		let search = Self::genome_info_to_search( &genome_info );	  
		let num_bins = genome_info
		.iter()
		.map(|(_, length, _)| (length + bin_width - 1) / bin_width)
		.sum::<usize>();

		let coverage_data = vec![0.0_f32; num_bins];

		Self {
			genome_info,
			search,
			coverage_data,
			bin_width,
			threads,
			nreads: 0,
		}
	}
// Constructor to initialize a BedData instance
    pub fn new(bam_file: &str , bin_width: usize, threads:usize, analysis_type: &AnalysisType, cell_tag: &[u8;2], umi_tag: &[u8;2], add_introns:bool ) -> Self {
    	
    	let mut ret = Self::init( bam_file, bin_width, threads );

    	ret.process_bam(bam_file, analysis_type, cell_tag, umi_tag, add_introns );
    	ret
    }

	pub fn process_bam( &mut self, bam_file: &str, analysis_type: &AnalysisType, cell_tag: &[u8;2], umi_tag: &[u8;2], add_introns:bool ) {


		let mut reader = match Reader::from_path(bam_file) {
    		Ok(r) => r,
    		Err(e) => panic!("Error opening BAM file: {}", e),
    	};

    	let cigar = Cigar::new( "" );
		// i32 id to chr name
    	let header = Header::from_template(reader.header());

		let ref_id_to_name = create_ref_id_to_name_hashmap(  &header );

    	let mut singlets = HashMap::<String, ReadData>::new();
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
	        let region: Option<Vec<(usize, usize)>> = match data_tuple {
	        	Ok(ref res) => {
	                let qname = &res.0; // Cell ID or read name as key
	                let read_data = &res.1;

	                if read_data.is("paired") {
	                	match singlets.remove(qname) {
	                		Some(first_read) => {
	                            // Mate found! Process the pair
	                            #[cfg(debug_assertions)]
	                            println!("Found paired reads: {:?} <-> {:?}", first_read, read_data);
	                            if first_read.chromosome != read_data.chromosome {
	                            	eprintln!( "chr missmatch!\n{}\n{}",first_read, read_data );
	                            	break
	                            }

	                            let mut covered_regions = Vec::new();

	                            let first_read_regions = cigar.read_on_database_matching_positions( &first_read.cigar, first_read.start, add_introns );
	                            let second_read_regions = cigar.read_on_database_matching_positions( &read_data.cigar, read_data.start, add_introns);

	                            covered_regions.extend(first_read_regions);
	                            covered_regions.extend(second_read_regions);

	                            Some(covered_regions)
                        		/*
	                            // add a +1 to the whole region.
	                            else if first_read.start < read_data.start {
	                            	Some((first_read.start, first_read.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ), )
	                            }else {
	                            	Some((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &first_read.cigar ).1 as i32 ))
	                            }*/

	                        }
	                        None => {
	                            // Handle orphaned reads (mate unmapped)
	                            if read_data.is("mate_unmapped") {
	                            	#[cfg(debug_assertions)]
	                            	println!("Orphaned read (mate unmapped): {:?}", read_data);
	                                // Process it as a single read

	                                Some (cigar.read_on_database_matching_positions(&read_data.cigar, read_data.start, add_introns))

	                                // Some((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ))
	                            } else {
	                                // No mate found, store this read for future pairing
	                                singlets.insert(qname.to_string(), read_data.clone());
	                                #[cfg(debug_assertions)]
	                                println!("Storing read for future pairing: {:?}", read_data);
	                                None
	                            }
	                        }
	                    }
	                } else {
	                    // Unpaired read (not part of a pair at all)
	                    #[cfg(debug_assertions)]
	                    println!("Unpaired read: {}", read_data);

	                    Some (cigar.read_on_database_matching_positions( &read_data.cigar, read_data.start, add_introns ))
                        //Some ((read_data.start, read_data.start + cigar.calculate_covered_nucleotides( &read_data.cigar ).1 as i32 ))
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

		    if let Some(vec) = region {
		    	//println!("I got {} regions to process!", vec.len());
		    	for (start, end) in vec {
		    		let start_window = (start) / self.bin_width;
		    		let end_window = (end )/ self.bin_width;
		    		self.nreads +=1;

		    		#[allow(unused_variables)]
		    		let (chrom_name, chrom_length, chrom_offset) = &self.genome_info[ record.tid()as usize ];

		    		for id in start_window..=end_window {
		    			#[cfg(debug_assertions)]
		    			println!("add to {chrom_name}:{} - +50", id* self.bin_width );
		    			if chrom_length / self.bin_width >= id {
		    				let index= chrom_offset+ id;
		    				self.coverage_data[index] += 1.0;
			                //println!("Adding a one to the data in index {index}");
			            }else {
			            	panic!( "pos {} is outside of the chromsome?! {chrom_length}", id as usize * self.bin_width )
			            }
			            
			        }


			    }
			}

			lines += 1;
		}
	    // clean up singlets
	    #[cfg(debug_assertions)]
	    println!( "clean up {} still still unpaired reads", singlets.len());
	    for region in &singlets {

	    	for (start, end) in cigar.read_on_database_matching_positions( &region.1.cigar, region.1.start, add_introns ){
	    		let start_window = (start +1) / self.bin_width;
	    		let end_window = (end +1 )/ self.bin_width;
	    		self.nreads +=1;

	    		#[allow(unused_variables)]
	    		let (chrom_name, chrom_length, chrom_offset) = match self.search.get ( &region.1.chromosome){
	    			Some(ret) => &self.genome_info[*ret],
	    			None => continue,
	    		};

	    		self.nreads +=1;

	    		for id in start_window..=end_window {
	    			#[cfg(debug_assertions)]
	    			println!("add to {chrom_name}:{} - +50", id* self.bin_width );
	    			if chrom_length / self.bin_width >= id {
	    				let index= chrom_offset+ id;
	    				self.coverage_data[index] += 1.0;
		                //println!("Adding a one to the data in index {index}");
		            }else {
		            	panic!( "pos {} is outside of the chromsome?! {chrom_length}", id as usize * self.bin_width )
		            }
		            
		        }
		    }
		    lines +=1;
		}
		
		//panic!("Finished with parsing the file");

	}

    /// Here we have implemented all the normalization features.
   	/// Not: just collect - no normalization whatsoever
   	///
    /// Rpkm: Reads Per Kilobase per Million mapped reads;
    /// number of reads per bin / (number of mapped reads (in millions) * bin length (kb))
    ///
    /// Cpm: Counts Per Million mapped reads
    /// number of reads per bin / number of mapped reads (in millions)
    ///
    /// Bmp: Bins Per Million mapped reads
    /// number of reads per bin / sum of all reads per bin (in millions)
    ///
    /// Rpgc: reads per genomic content (1x normalization)
    /// number of reads per bin / scaling factor for 1x average coverage.
    /// This scaling factor, in turn, is determined from the sequencing depth: 
    /// (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of
    /// the sequencing depth computed for the sample to match the 1x coverage. This option requires --effectiveGenomeSize.
    pub fn normalize(&mut self, by: &Normalize ) {
    	match by {
    		Normalize::Not => {
    			// do just notthing ;-)
    		},
    		Normalize::Rpkm => {
    			// number of mapped reads * bin length / (in millions) / (kb)
    			let div: f32 = (self.nreads * self.bin_width )  as f32 / 1_000_000.0 / 1_000.0 ; 
    			// number of reads per bin / div
    			self.coverage_data.par_iter_mut().for_each(|x| *x /= div);
    		},
    		Normalize::Cpm => {
    			//  number of mapped reads (in millions)
    			let div: f32 = self.nreads as f32 / 1_000_000.0;
    			// number of reads per bin / div
    			self.coverage_data.par_iter_mut().for_each(|x| *x /= div);
    		},
    		Normalize::Bpm => {
    			// sum of all reads per bin (in millions)
    			let div: f32 = self.coverage_data.par_iter().sum();
    			// number of reads per bin / div
    			self.coverage_data.par_iter_mut().for_each(|x| *x /= div);
    		},
    		Normalize::Rpgc => {
    			panic!("Sorry - not implemented at the moment");
    		}
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
	    let mut iter = DataIter::new( self );

	    while let Some(values) = &iter.next(){
	    	writeln!(
	    		file,
	    		"{}\t{}\t{}\t{}",
	            values.0,   // Chromosome name
	            values.1.start, values.1.end, values.1.value
	            ).unwrap();
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