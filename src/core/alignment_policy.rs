//core.rs
use rust_htslib::bam::Record;
use crate::cli::CoverageCli;


/// SAM flag bits (u16)
const FLAG_PAIRED: u16        = 0x1;
const FLAG_UNMAPPED: u16      = 0x4;
const FLAG_SECONDARY: u16     = 0x100;
const FLAG_QCFAIL: u16        = 0x200;
const FLAG_DUPLICATE: u16     = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const FLAG_READ2: u16         = 0x80; // "second in pair"


#[derive(Clone, Copy, Debug)]
pub struct AlignmentPolicy {
    pub min_mapq: u8,
    pub include_secondary: bool,
    pub include_supplementary: bool,
    pub include_duplicates: bool,
    pub only_r1: bool,

    /// Final, effective exclude mask used by `passes_filter`.
    /// This is computed from booleans unless overridden by CLI.
    sam_flag_exclude:u16,
}

impl AlignmentPolicy {

    /// CLI constructor:
    /// - builds default mask from booleans
    /// - if `cli.sam_flag_exclude != 0`, that overrides the computed mask
    pub fn from_cli(cli: &crate::cli::CoverageCli) -> Self {
        let computed = Self::computed_exclude_mask_bamcoverage(cli);

        let effective = cli.sam_flag_exclude.unwrap_or(computed);

        Self {
            min_mapq: cli.min_mapping_quality,
            include_secondary: cli.include_secondary,
            include_supplementary: cli.include_supplementary,
            include_duplicates: cli.include_duplicates,
            only_r1: cli.only_r1,
            sam_flag_exclude: effective,
        }
    }


    fn computed_exclude_mask_bamcoverage(cli: &CoverageCli) -> u16 {
        let mut mask =
              FLAG_UNMAPPED
            | FLAG_SECONDARY
            | FLAG_QCFAIL
            | FLAG_DUPLICATE
            | FLAG_SUPPLEMENTARY;

        // Allow overrides
        if cli.include_secondary {
            mask &= !FLAG_SECONDARY;
        }
        if cli.include_supplementary {
            mask &= !FLAG_SUPPLEMENTARY;
        }
        if cli.include_duplicates {
            mask &= !FLAG_DUPLICATE;
        }
        if cli.only_r1 {
            mask |= FLAG_READ2;
        }

        mask
    }

    
    /// Filters on:
    /// - unmapped (always excluded here)
    /// - effective SAM exclude mask (computed or overridden)
    /// - MAPQ
    #[inline]
    pub fn passes_filter(&self, rec: &Record) -> bool {

        let flag: u16 = rec.flags();

        // Single bitmask check replaces all the redundant per-flag checks.
        if (flag & self.sam_flag_exclude) != 0 {
            return false;
        }

        if rec.mapq() < self.min_mapq {
            return false;
        }

        true
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Record;

    /// Create a minimal BAM record with configurable flags and MAPQ.
    /// This is only for testing filter logic, not alignment semantics.
    fn fake_record(flags: u16, mapq: u8) -> Record {
        let mut rec = Record::new();
        rec.set_flags(flags);
        rec.set_mapq(mapq);
        rec
    }

    #[test]
    fn test_passes_basic_record() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: false,
            include_duplicates: false,
            only_r1: false,
        };

        let rec = fake_record(0, 30);
        assert!(policy.passes_filter(&rec));
    }

    #[test]
    fn test_filters_unmapped() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: true,
            include_supplementary: true,
            include_duplicates: true,
            only_r1: false,
        };

        // 0x4 = unmapped
        let rec = fake_record(0x4, 30);
        assert!(!policy.passes_filter(&rec));
    }

    #[test]
    fn test_filters_secondary() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: true,
            include_duplicates: true,
            only_r1: false,
        };

        // 0x100 = secondary
        let rec = fake_record(0x100, 30);
        assert!(!policy.passes_filter(&rec));
    }

    #[test]
    fn test_allows_secondary_if_enabled() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: true,
            include_supplementary: true,
            include_duplicates: true,
            only_r1: false,
        };

        let rec = fake_record(0x100, 30);
        assert!(policy.passes_filter(&rec));
    }

    #[test]
    fn test_filters_supplementary() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: true,
            include_supplementary: false,
            include_duplicates: true,
            only_r1: false,
        };

        // 0x800 = supplementary
        let rec = fake_record(0x800, 30);
        assert!(!policy.passes_filter(&rec));
    }

    #[test]
    fn test_filters_duplicates() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: true,
            include_supplementary: true,
            include_duplicates: false,
            only_r1: false,
        };

        // 0x400 = duplicate
        let rec = fake_record(0x400, 30);
        assert!(!policy.passes_filter(&rec));
    }

    #[test]
    fn test_only_r1_filters_r2() {
        let policy = AlignmentPolicy {
            min_mapq: 0,
            include_secondary: true,
            include_supplementary: true,
            include_duplicates: true,
            only_r1: true,
        };

        // 0x80 = second in pair
        let rec = fake_record(0x80, 30);
        assert!(!policy.passes_filter(&rec));
    }

    #[test]
    fn test_min_mapq() {
        let policy = AlignmentPolicy {
            min_mapq: 20,
            include_secondary: true,
            include_supplementary: true,
            include_duplicates: true,
            only_r1: false,
        };

        let good = fake_record(0, 30);
        let bad = fake_record(0, 10);

        assert!(policy.passes_filter(&good));
        assert!(!policy.passes_filter(&bad));
    }
}
