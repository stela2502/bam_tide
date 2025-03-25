use regex::Regex;
use rustody::mapping_info::MappingInfo;
use rustody::genes_mapper::cigar::{Cigar, CigarEnum};


pub struct MutationProcessor {
    pub quality_cutoff: usize, // Mean quality score threshold
}

impl MutationProcessor {

    /// the main entry into the mutations.
    /// Give me the bam entries start, the Cigar,  sequence and quality scores and I give you a vector detected mutations
    pub fn get_all_mutations(&self, bam_start:usize, cigar:&str, sequence:&[u8], qual:&[u8], mapping_info: &mut MappingInfo ) -> Vec<String>{
        let cigar_obj = Cigar::new("");

        let mut current_position = bam_start; // Start from the provided position
        let mut current_loc_pos = 0;
        let mut ret = Vec::new();


        for cap in cigar_obj.str_to_tuple_vec(cigar, false) {

            if cap.option != CigarEnum::Nothing {
                if current_loc_pos + cap.len() > qual.len() {
                    panic!("MutationProcessor::get_all_mutations - cigar {cigar} lead to a out of range issue: {current_loc_pos} + {} ({})> {}", 
                        cap.len() ,current_loc_pos + cap.len(),qual.len() );
                }
                if let Some(name) = self.mutation_name( &cap.option, current_position, current_loc_pos, cap.len(), sequence, qual, mapping_info) {
                    ret.push(name)
                }else {
                    mapping_info.report("not mutated");
                }
                if cap.option.adds_to_read(true) {
                    current_position += cap.len();
                    current_loc_pos += cap.len();
                }
            }
        }
        ret
    }

    fn mutation_name(
        &self,
        mutation_type: &CigarEnum,
        start: usize,
        local_start: usize,
        length: usize,
        seq: &[u8], // Accept slice of u8 directly
        qual: &[u8], // Accept slice of u8 directly
        mapping_info: &mut MappingInfo,
    ) -> Option<String> {
        // Check the quality of nucleotides involved in the mutation
        if start+length > qual.len() {
            return None
        }
        let quality_check = self.check_quality(local_start, length, qual);

        // If the mutation quality is below the cutoff, skip it
        if quality_check < self.quality_cutoff {
            mapping_info.report(&format!("mutation quality {} - failed",quality_check ) );
            return None;
        }

        // Handle each mutation type
        match mutation_type {
            CigarEnum::Mismatch  => self.mismatch_name(start, local_start,length, seq),
            CigarEnum::Deletion  => self.deletion_name(start, length),
            CigarEnum::Insertion => self.insertion_name(start, local_start, length, seq),
            _ => None, // If the mutation type is unhandled, return None
        }
    }

    fn check_quality(&self, start: usize, length: usize, qual: &[u8]) -> usize {
        // Substring the quality scores involved in the mutation
        
        let mutation_qual = &qual[start..start + length];

        // Calculate the mean quality score for this substring
        let total_quality: usize = mutation_qual.iter().map(|&c| (c as usize) - 33).sum();
        total_quality / mutation_qual.len()
    }

    fn mismatch_name(&self, start: usize, local_start:usize, length: usize, seq: &[u8]) -> Option<String> {
        let mismatch_seq = &seq[local_start..local_start + length];
        // Convert the mutated sequence slice into a String when needed
        Some(format!(
            "snp/{}/{}",
            start + 1, // Position is 1-based in mutation naming
            String::from_utf8_lossy(mismatch_seq) // Convert slice to String only for the mutated part
        ))
    }

    fn deletion_name(&self, start: usize, length: usize) -> Option<String> {
        Some(format!(
            "ins/{}/{}",
            start + 1, // Position is 1-based in mutation naming,
            length,
        ))
    }

    fn insertion_name(&self, start: usize, local_start:usize, length: usize, seq: &[u8]) -> Option<String> {
        let inserted_seq = &seq[local_start..local_start + length];
        // Convert the inserted sequence slice into a String when needed
        Some(format!(
            "del/{}/{}",
            start + 1, // Position is 1-based in mutation naming
            String::from_utf8_lossy(inserted_seq), // Convert slice to String only for the mutated part
        ))
    }
}
