use crate::ont_normalizer::fastq_record::FastqRecord;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BarcodeMatcherSpec {
    Any,
    Whitelist {
        path: PathBuf,
        max_mismatches: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrimerPart {
    Fixed {
        name: String,
        seq: Vec<u8>,
        max_mismatches: usize,
    },
    Random {
        name: String,
        min_len: usize,
        max_len: usize,
    },
    CellId {
        name: String,
        len: usize,
        matcher: BarcodeMatcherSpec,
    },
    Umi {
        name: String,
        len: usize,
    },
    PolyT {
        min_len: usize,
        max_non_t: usize,
    },
    Insert,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerStructure {
    pub name: String,
    pub adapter: Vec<u8>,
    pub adapter_suffix_min_len: usize,
    pub max_adapter_mismatches: usize,
    pub parts: Vec<PrimerPart>,
}

#[derive(Debug, Clone)]
pub struct PrimerHit {
    pub structure_name: String,
    pub orientation: Orientation,

    /// Start of the primer structure.
    /// This closes the previous insert.
    pub primer_start: usize,

    /// End of the primer structure.
    pub primer_end: usize,

    /// Start of the biological insert.
    /// This opens the next emitted molecule.
    pub insert_start: usize,

    pub cell_id: Option<FastqRecord>,
    pub umi: Option<FastqRecord>,
}

#[derive(Debug, Clone)]
pub struct PrimerSplit {
    pub original_read_id: String,
    pub insert: FastqRecord,
    pub cell_id: Option<FastqRecord>,
    pub umi: Option<FastqRecord>,
    pub orientation: Orientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AdapterHit {
    pub primer_start: usize,
    pub adapter_end: usize,
}

#[derive(Debug, Default)]
pub(crate) struct MatchState {
    pub cell_seq: Vec<u8>,
    pub cell_qual: Vec<u8>,
    pub umi_seq: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub insert_start: Option<usize>,
}