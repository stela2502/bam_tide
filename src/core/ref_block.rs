//ref_block.rs

use rust_htslib::bam::Record;
use rust_htslib::bam::record::Cigar;

/// A contiguous reference interval covered by an alignment.
/// Coordinates are 0-based, half-open: [start, end).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RefBlock {
    pub start: u32,
    pub end: u32,
}

impl RefBlock {
    /// Create a new block, assuming start < end.
    pub fn new(start: u32, end: u32) -> Self {
        Self { start, end }
    }
}

/// Convert a BAM record into reference-covered blocks using its CIGAR string.
///
/// - Only reference-consuming *and* coverage-producing ops (M, =, X) create blocks.
/// - Ops that consume reference but do not represent coverage (D, N) advance the reference
///   without emitting blocks.
/// - Other ops (I, S, H, P) do not consume reference and do not emit blocks.
/// - Output blocks are 0-based, half-open and adjacent blocks are merged.
///
/// This function assumes the record is mapped; callers should filter unmapped reads beforehand.
pub fn record_to_blocks(rec: &Record) -> Vec<RefBlock> {
    // BAM position is 0-based; for unmapped reads it's typically -1.
    // We assume filtering already removed unmapped records, but guard anyway.
    let pos0 = rec.pos();
    if pos0 < 0 {
        return Vec::new();
    }

    let mut ref_pos: u64 = pos0 as u64;
    let mut blocks: Vec<RefBlock> = Vec::new();

    for op in rec.cigar().iter() {
        match *op {
            // Coverage on reference: emit a block
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_u64 = len as u64;
                if len_u64 == 0 {
                    continue;
                }

                let start = ref_pos;
                let end = ref_pos + len_u64;

                // Merge with previous block if adjacent.
                if let Some(last) = blocks.last_mut() {
                    if last.end as u64 == start {
                        // Safe because end fits in u64; we clamp to u32 below.
                        last.end = end.min(u32::MAX as u64) as u32;
                    } else {
                        blocks.push(RefBlock::new(
                            start.min(u32::MAX as u64) as u32,
                            end.min(u32::MAX as u64) as u32,
                        ));
                    }
                } else {
                    blocks.push(RefBlock::new(
                        start.min(u32::MAX as u64) as u32,
                        end.min(u32::MAX as u64) as u32,
                    ));
                }

                // Advance reference position
                ref_pos = end;
            }

            // Consumes reference but is not "covered bases"
            // (deletions and splices/introns)
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += len as u64;
            }

            // Does not consume reference; no coverage emitted
            Cigar::Ins(_len)
            | Cigar::SoftClip(_len)
            | Cigar::HardClip(_len)
            | Cigar::Pad(_len) => {
                // no-op for reference position
            }

            // Some CIGAR implementations include Back (rare). If present, treat conservatively.
            #[allow(unreachable_patterns)]
            Cigar::Back(len) => {
                // Move backwards on reference (rare/edge-case). Avoid underflow.
                let l = len as u64;
                ref_pos = ref_pos.saturating_sub(l);
            }
        }
    }

    blocks
}



#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn fake_record(pos: i64, cigar: CigarString) -> Record {
        let mut rec = Record::new();
        rec.set_pos(pos);
        rec.set_cigar(&cigar);
        rec
    }

    #[test]
    fn test_simple_match() {
        let rec = fake_record(100, CigarString(vec![Cigar::Match(10)]));
        let blocks = record_to_blocks(&rec);
        assert_eq!(blocks, vec![RefBlock { start: 100, end: 110 }]);
    }

    #[test]
    fn test_splice() {
        let rec = fake_record(
            100,
            CigarString(vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)]),
        );
        let blocks = record_to_blocks(&rec);
        assert_eq!(
            blocks,
            vec![
                RefBlock { start: 100, end: 105 },
                RefBlock { start: 115, end: 120 }
            ]
        );
    }
}

