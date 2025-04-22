use crate::gtf::exon_iterator::{ExonIterator,ReadResult};

#[derive(Debug, PartialEq)]
pub enum QueryErrors {
    OutOfGenes,
    ChrNotFound,
    NoMatch,
}


pub trait FeatureMatcher {
    /// process feature gets one or two BAM reads and needs to figure out how to add them to the exp_gex and mut_gex single cells data obejcts.

    fn init_search(
        &self,
        chr: &str,
        start: usize,
        iterator: &mut ExonIterator,
    ) -> Result<(), QueryErrors>;

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