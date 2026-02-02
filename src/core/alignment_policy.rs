//core.rs
use rust_htslib::bam::Record;

#[derive(Clone, Copy, Debug)]
pub struct AlignmentPolicy {
    pub min_mapq: u8,
    pub include_secondary: bool,
    pub include_supplementary: bool,
    pub include_duplicates: bool,
    pub only_r1: bool,
}

impl AlignmentPolicy {
    pub fn from_cli(cli: &crate::cli::CoverageCli) -> Self {
        Self {
            min_mapq: cli.min_mapping_quality,
            include_secondary: cli.include_secondary,
            include_supplementary: cli.include_supplementary,
            include_duplicates: cli.include_duplicates,
            only_r1: cli.only_r1,
        }
    }
    /// Check whether a BAM record should be included given this policy.
    /// Filters on unmapped/secondary/supplementary/duplicate flags, MAPQ, and optional "read1 only".
    pub fn passes_filter(&self, rec: &Record) -> bool {
        if rec.is_unmapped() {
            return false;
        }
        if !self.include_secondary && rec.is_secondary() {
            return false;
        }
        if !self.include_supplementary && rec.is_supplementary() {
            return false;
        }
        if !self.include_duplicates && rec.is_duplicate() {
            return false;
        }
        if self.only_r1 && !rec.is_first_in_template() {
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
