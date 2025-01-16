
use bam::{BamReader, Record, Header, RecordReader};
use bam::record::tags::{StringType, TagValue};

use std::collections::HashMap;
use std::fs::File;
use std::io::{Write};
use std::process::exit;

/// Processes a BAM file to calculate coverage and writes the result to a BigWig file.
pub fn bam_to_bedgraph(
    bam_file: &str,
    bedgraph_file: &str,
    bin_width: usize,
) -> Result<(), String > {
    // Open BAM file
    let mut bam_reader = BamReader::from_path(&bam_file, 1).unwrap();
    let header = bam_reader.header().to_owned();
    let genome_info:Vec<(String, usize, usize)> = create_ref_id_to_name_map( &header,bin_width );
    let mut geneome_sizes: HashMap<String, usize> = genome_info.iter()
        .map(|(name, length, _offset)| (name.to_string(), *length))
        .collect();

    let num_bins = calculate_total_bins(&genome_info, bin_width);
    println!("The total number of bins we need: {num_bins}for these genome sizes: {:?}",geneome_sizes);
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

    println!("Collected the data - writing file");

    // here is where you can implement the normalization steps.

    // Write data to text
    write_bedgraph(bedgraph_file, &genome_info, &coverage_data, bin_width).unwrap();

    Ok(())
}

/// Create a mapping of reference ID to reference name, chromosome length, and bin offset.
/// `bin_width` is the size of each bin.
fn create_ref_id_to_name_map(header_view: &Header, bin_width: usize) -> Vec<(String, usize, usize)> {
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

/// Calculate total number of bins needed.
pub fn calculate_total_bins(genome_info: &Vec<(String, usize, usize)>, bin_width: usize) -> usize {
    genome_info
        .iter() // We want to iterate over the vector
        .map(|(_, size, _)| (size + bin_width - 1) / bin_width) // We only care about the chromosome length
        .sum() // Sum all the individual bin counts to get the total
}


/// Get the index in the coverage array for a given position.
fn get_array_index(
    genome_info: &Vec<(String, usize, usize)>,
    chrom: usize,
    position: usize,
    bin_width: usize,
) -> Option<usize> {
    if position >= genome_info[chrom].1 {
        None
    }else {
        Some(genome_info[chrom].2 + position / bin_width)
    }
    
}


/// Write the coverage data to a BEDGraph file.
fn write_bedgraph(
    file_path: &str,                         // Output file path
    genome_info: &Vec<(String, usize, usize)>, // Genome information (name, length, bin offset)
    coverage: &Vec<u32>,                   // Coverage data array (assuming coverage for each bin)
    bin_width: usize,                        // Bin width
) -> std::io::Result<()> {
    // Open the output file
    let mut file = File::create(file_path)?;

    // Iterate over the genome info and coverage data
    let mut bin_start: usize = 0; // Start position of each bin
    for (chrom_index, (chrom_name, chrom_length, bin_offset)) in genome_info.iter().enumerate() {
        let mut bin_end = 0; // End position of each bin
        
        // Iterate through the coverage data for the current chromosome
        for bin_index in 0..(chrom_length / bin_width) {
            // Calculate start and end positions for each bin
            bin_start = bin_index * bin_width;
            bin_end = (bin_index + 1) * bin_width;

            if bin_end > *chrom_length {
                bin_end = *chrom_length; // Adjust for last bin that may not fill up to full width
            }

            // Get the coverage value for this bin (you'll need to map the correct coverage)
            let coverage_value = *coverage.get(bin_offset + bin_index).unwrap_or(&0);

            // Write to the BedGraph file
            if coverage_value > 0 {
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