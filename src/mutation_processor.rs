//use regex::Regex;
//use rustody::mapping_info::MappingInfo;
//use rustody::genes_mapper::cigar::{Cigar, CigarEnum};
use needletail::parse_fastx_file;
use std::collections::HashMap;

use rust_htslib::bam::{ Read, Reader, Header };

use crate::read_data::ReadData;

use rustody::singlecelldata::IndexedGenes;
use rustody::singlecelldata::SingleCellData;
use rustody::mapping_info::MappingInfo;
use rustody::singlecelldata::cell_data::GeneUmiHash;


pub struct MutationProcessor {
    /// the offsets for each chr
    genome_offsets: Vec<usize>,
    // the chr names (for matching to the fasta names)
    chr_names: HashMap<String, usize>,
    /// the real genome vector
    genome: Vec<u8>,
} 

impl MutationProcessor {


    /// Initializes the MutationProcessor from a BAM header and a FASTA file
    pub fn new( bam_file: &str, fasta: &str ) -> Result<Self, String> {
        
        let reader = match Reader::from_path(bam_file) {
            Ok(r) => r,
            Err(e) => panic!("Error opening BAM file: {}", e),
        };


        let header = Header::from_template(reader.header());

        let header_map = header.to_hashmap();
        let mut total_length = 0;

        let mut chr_names = HashMap::new();
        let mut genome_offsets = Vec::new();
        let mut i = 0;

        // Parse the BAM header and extract chromosome names and lengths
        if let Some(reference_info) = header_map.get("SQ") {
            for record in reference_info {
                if let Some(length) = record.get("LN") {
                    // The "LN" tag holds the length of the chromosome
                    if let Ok(length_int) = length.parse::<usize>() {
                        if let Some(name) = record.get("SN") {
                            chr_names.insert(name.to_string(), i);
                            genome_offsets.push(total_length);
                            total_length += length_int as usize;
                            i +=1;
                        }
                    }
                }
            }
        }

        // Initialize the genome vector
        let mut genome = vec![0; total_length]; // Initialize with zeros

        // Load the FASTA file and populate the genome sequence
        let mut fasta_file = parse_fastx_file(fasta).map_err(|_| "Invalid FASTA file path".to_string())?;
        let mut found_chrs = 0;

        while let Some(e_record) = fasta_file.next() {
            let seqrec = e_record.map_err(|_| "Error parsing FASTA record".to_string())?;
            match std::str::from_utf8(seqrec.id()) {
                Ok(st) => {
                    // Extract chromosome name from the FASTA header
                    if let Some(orig_gname) = st.split(' ').next() {
                        if let Some(&chr_id) = chr_names.get(orig_gname) {
                            // We need that sequence
                            let offset = genome_offsets[chr_id];
                            let length = if chr_id + 1 < genome_offsets.len() {
                                genome_offsets[chr_id + 1] - offset
                            } else {
                                total_length - offset
                            };

                            // Check length consistency
                            if length != seqrec.seq().len() {
                                return Err(format!(
                                    "chr {} in fasta is not the expected size from the bam header bam: {} != fasta: {}",
                                    orig_gname,
                                    length,
                                    seqrec.seq().len()
                                ));
                            } else {
                                // Copy sequence into the genome vector
                                genome[offset..(offset + length)]
                                    .copy_from_slice(&seqrec.seq());
                                found_chrs += 1;
                            }
                        }
                    }
                }
                Err(_) => return Err("Invalid FASTA sequence name".to_string()),
            }
        }

        // Ensure all chromosomes from BAM header are found in the FASTA file
        if found_chrs != chr_names.len() {
            return Err(format!(
                "Missing {} chromosomes expected from the bam file header: {:?}",
                chr_names.len() - found_chrs,
                chr_names
            ));
        }

        // Return the initialized MutationProcessor
        Ok(MutationProcessor {
            genome_offsets,
            chr_names,
            genome,
        })
    }

    pub fn handle_mutations(
        &self,
        data: &ReadData,
        gene_id: &str,
        mut_idx: &mut IndexedGenes,
        mut_gex: &mut SingleCellData,
        mapping_info: &mut MappingInfo,
        cell_id: &u64,
    ) {
        
        for mutation_name in self.get_all_mutations(
            &data.chromosome,
            data.start.try_into().unwrap(),
            &data.cigar,
        ) {
            let snip = format!("{gene_id}/{}/{}", data.chromosome, mutation_name);
            let mut_id = mut_idx.get_gene_id(&snip);
            let ghum = GeneUmiHash(mut_id, data.umi);

            let _ = mut_gex.try_insert(cell_id, ghum, mapping_info);
        }

    }

    pub fn get_all_mutations(&self, chr_name: &str, start: usize, md_tag: &str) -> Vec<String> {
        let mut mutations = Vec::new();
        
        // Retrieve chromosome name
        let chr_id = self.chr_names.get(chr_name).unwrap_or_else(|| panic!("I copuld not find chr {} in my data!", chr_name));
        
        let mut current_pos = start -1;
        let mut i = 0;
        
        // Iterate over the MD tag, which contains the match/mismatch information
        //let segments = md_tag.split(|c| c == 'A' || c == 'T' || c == 'C' || c == 'G');
        let mut match_len = 0;

        for c in md_tag.chars() {
            match c {
                // Match or mismatch (number indicates match length)
                '0'..='9' => {
                    match_len = match_len * 10 + c.to_digit(10).unwrap() as usize;
                }
                // Mismatch - the character in the MD string
                'A' | 'T' | 'C' | 'G' => {
                    // First, move the current position based on the match length
                    if match_len > 0 {
                        current_pos += match_len;
                        match_len = 0;  // Reset match length after moving
                    }
                    // When we encounter a mismatch, we construct the mutation string
                    //println!("I am getting the nucleotide at position {}: {}", self.genome_offsets[*chr_id] + current_pos, char::from(self.genome[self.genome_offsets[*chr_id] + current_pos]));
                    let ref_base = self.genome[self.genome_offsets[*chr_id] + current_pos ];// Get the reference base at current position
                    let alt_base = c;
                    mutations.push(format!(
                        "{}:g.{}{}>{}",
                        chr_name,
                        current_pos,
                        ref_base as char,
                        alt_base
                    ));

                    // Update the current position
                    current_pos += 1;
                }
                // Handle deletions, where ^ indicates deleted bases in the reference
                '^' => {
                    i += 1; // Skip the '^'
                    let mut deleted_bases = String::new();
                    while i < md_tag.len() && md_tag.chars().nth(i).unwrap().is_alphabetic() {
                        deleted_bases.push(md_tag.chars().nth(i).unwrap());
                        i += 1;
                    }

                    // We are not changing the current position since the bases were deleted in the reference
                    mutations.push(format!(
                        "{}:g.{}del{}",
                        chr_name,
                        current_pos + 1, // Position in 1-based coordinate
                        deleted_bases
                    ));
                }
                // Skip anything else - should not be anything in fact
                _ => {}
            }
        }

        mutations
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_mutation_processor_new() {
        let fasta_path = "testData/mutation_test.fa.gz";
        let bam_path = "testData/mutation_test.bam"; // Or use a dummy path if unused

        assert!(Path::new(fasta_path).exists(), "FASTA test file not found.");
        assert!(Path::new(bam_path).exists(), "BAM test file not found.");

        match MutationProcessor::new(bam_path, fasta_path) {
            Ok(processor) => {
                assert!(!processor.genome.is_empty(), "Genome should not be empty.");
                assert!(
                    !processor.chr_names.is_empty(),
                    "Chromosome names should be loaded."
                );
                println!("Loaded genome with {} bases.", processor.genome.len());
            }
            Err(e) => panic!("MutationProcessor failed to load: {}", e),
        }
    }
    

    #[test]
    fn test_handle_mutations() {
        // Set up dummy genome
        let fasta_path = "testData/mutation_test.fa.gz";
        /*
>chr1
00000000001111111111222222222233
01234567890123456789012345678901
ACTGACTGACTGACTGACTGACTGACTGACTG


        */
        let bam_path = "testData/mutation_test.bam"; // Not used, can be dummy
        /*
@HD VN:1.6  SO:coordinate
@SQ SN:chr1 LN:32
@PG ID:samtools PN:samtools VN:1.19.2   CL:samtools view -b mutation_test.sam
@PG ID:samtools.1   PN:samtools PP:samtools VN:1.19.2   CL:samtools view -h mutation_test.bam
@PG ID:samtools.2   PN:samtools PP:samtools.1   VN:1.19.2   CL:samtools view -h mutation_test.sam
@PG ID:samtools.3   PN:samtools PP:samtools.2   VN:1.19.2   CL:samtools view -h mutation_test.bam
read1:AGCTGCTGAGATCGATACAGTAATTGGTCA    0   chr1    1   60  10M *   0   0   ACTGACTTAC  *   MD:Z:6T3
read2:AGCTGCTGAGATCGATACAGTTTTTGGTCA    0   chr1    4   60  10M *   0   0   ACTGACTTAC  *   MD:Z:2T7
read3:AGCTGCTGAGATCGATACAGTTAAAGGTCA    0   chr1    4   60  10M *   0   0   ACTGACTTAC  *   MD:Z:4C5
        */
        let processor = MutationProcessor::new(bam_path, fasta_path)
            .expect("Failed to load MutationProcessor");

        let mut reader = match Reader::from_path(bam_path) {
            Ok(r) => r,
            Err(e) => panic!("Error opening BAM file: {}", e),
        };

        let mut mapping_info = MappingInfo::new(None, 3.0, 0, None);
        let mut gex = SingleCellData::new(1);
        let mut idx = IndexedGenes::empty(Some(0));
        let cell_id = 1_u64;


        for r in reader.records() {
            // Read a record from BAM file
            let record = match r {
                Ok(r) => r,
                Err(e) => panic!("I could not collect a read: {e:?}"),
            };
            let mut ref_id_to_name = HashMap::<i32, String>::new();
            ref_id_to_name.insert( 0, "chr1".to_string() );
            match ReadData::from_singlecell_bowtie2(&record, &ref_id_to_name) {
                Ok((_id, read)) => {
                    // Process the read if successful
                    processor.handle_mutations(&read, "unknown", &mut idx, &mut gex, &mut mapping_info, &cell_id);
                },
                Err(e) => {
                    // Handle the error, log it, and continue processing
                    panic!("Error processing BAM record: {:?}", e);

                }
            }
            
        }

        assert_eq!( idx.get_all_gene_names().len(), 2, "detected only one mutation {:?}", idx.get_all_gene_names() );
        assert_eq!( idx.get_all_gene_names()[0], "unknown/chr1/chr1:g.7G>T" );
        //assert_eq!( idx.get_all_gene_names()[1], "unknown/chr1/chr1:g.7G>T" );
        assert_eq!( idx.get_all_gene_names()[1], "unknown/chr1/chr1:g.8A>C" );

    }

}