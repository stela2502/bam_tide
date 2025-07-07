use crate::gtf::exon_iterator::{ExonIterator,ReadResult};
use crate::read_data::ReadData;
use crate::mutation_processor::MutationProcessor;
use crate::gtf_logics::MatchType;


use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::IndexedGenes;
use rustody::mapping_info::MappingInfo;
use rustody::int_to_str::IntToStr;


#[derive(Debug, PartialEq )]
pub enum QueryErrors {
    OutOfGenes,
    ChrNotFound,
    NoMatch,
}

pub trait FeatureMatcher: Sync + Send  + std::fmt::Display {
    /// process feature gets one or two BAM reads and needs to figure out how to add them to the exp_gex and mut_gex single cells data obejcts.

    fn init_search(
        &self,
        chr: &str,
        start: usize,
        iterator: &mut ExonIterator,
    ) -> Result<(), QueryErrors>;

    /*fn process_buffer<T>(
        &self,
        data: &(ReadData, Option<ReadData>),
        mutations: &Option<MutationProcessor>,
        iterator: &mut ExonIterator,
        exp_gex: &mut SingleCellData,
        exp_idx: &mut IndexedGenes,
        mut_gex: &mut SingleCellData,
        mut_idx: &mut IndexedGenes,
        mapping_info: &mut MappingInfo,
        match_type: &MatchType,
    )
    where
        T: FeatureMatcher + Send + Sync;
    */
    fn extract_gene_ids(
        &self,
        read_result: &Option<Vec<ReadResult>>,
        data: &ReadData,
        mapping_info: &mut MappingInfo,
    ) -> Vec<String>;

    fn parse_cell_id(cell_id_str: &str) -> Result<u64, String> {
        cell_id_str.parse::<u64>().or_else(|_| {

            let ret = match IntToStr::new(cell_id_str.as_bytes().to_vec(), 32){
                Ok(obj) => Ok(obj.into_u64()),
                Err(e) => Err( format!("cell_name could not be parsed to u64: {e}") )
            };
            //println!("I am trying to encode the cell ID {} as 2bit u64 {:?}", cell_id_str, ret);
            ret
        })
    }

    fn process_feature(
        &self,
        data: &(ReadData, Option<ReadData>),
        mutations: &Option<MutationProcessor>,
        iterator: &mut ExonIterator,
        exp_gex: &mut SingleCellData,
        exp_idx: &mut IndexedGenes,
        mut_gex: &mut SingleCellData,
        mut_idx: &mut IndexedGenes,
        mapping_info: &mut MappingInfo,
        match_type: &MatchType,
    );


}