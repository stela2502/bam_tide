//use regex::Regex;
//use rustody::mapping_info::MappingInfo;
//use rustody::genes_mapper::cigar::{Cigar, CigarEnum};
use needletail::parse_fastx_file;
use std::collections::HashMap;

use rust_htslib::bam::Header;

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
    pub fn new( header: &Header, fasta: &str ) -> Result<Self, String> {
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
        
        for mutation_name in processor.get_all_mutations(
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
                    let ref_base = self.genome[self.genome_offsets[*chr_id] + current_pos];// Get the reference base at current position
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
