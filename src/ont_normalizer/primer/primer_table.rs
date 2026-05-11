use std::collections::HashMap;

use int_to_str::IntToStr;

#[derive(Debug, Clone)]
pub struct PrimerKmerInfo {
    pub primer_start_offset: i32,
    pub kmer_start: usize,
}

#[derive(Debug, Clone)]
pub struct PrimerStartCandidate {
    pub primer_start: i32,
    pub primer_end: usize,
    pub support: usize,
    pub first_kmer_start: usize,
    pub last_kmer_start: usize,
}

#[derive(Debug, Clone)]
pub struct PrimerTable {
    pub primer: Vec<u8>,
    pub kmer_size: usize,
    pub infos: Vec<PrimerKmerInfo>,
    pub table: HashMap<Vec<u8>, Vec<usize>>,
}

impl PrimerTable {
    pub fn new<T: AsRef<[u8]>>(primer: T, kmer_size: usize) -> Result<Self, String> {
        let primer = primer.as_ref().to_vec();

        if kmer_size == 0 {
            return Err("kmer_size must be > 0".to_string());
        }

        if kmer_size % 4 != 0 {
            return Err("kmer_size must currently be a multiple of 4".to_string());
        }

        let mut this = Self {
            primer,
            kmer_size,
            infos: Vec::new(),
            table: HashMap::new(),
        };

        this.build()?;
        Ok(this)
    }

    fn build(&mut self) -> Result<(), String> {
        if self.primer.len() < self.kmer_size {
            return Ok(());
        }

        for kmer_start in 0..=(self.primer.len() - self.kmer_size) {
            let word = IntToStr::enc_bytes(
                &self.primer[kmer_start..kmer_start + self.kmer_size],
            )?;

            let info_id = self.infos.len();

            self.infos.push(PrimerKmerInfo {
                primer_start_offset: -(kmer_start as i32),
                kmer_start,
            });

            self.table.entry(word).or_default().push(info_id);
        }

        Ok(())
    }

    pub fn detect_primer_start(
        &self,
        read: &IntToStr,
    ) -> Option<PrimerStartCandidate> {
        self.detect_candidates(read).into_iter().max_by_key(|c| {
            (
                c.support,
                c.last_kmer_start.saturating_sub(c.first_kmer_start),
            )
        })
    }

    pub fn detect_candidates(
        &self,
        read: &IntToStr,
    ) -> Vec<PrimerStartCandidate> {
        let read_seq = read.to_string(read.size);
        let read_seq = read_seq.as_bytes();

        if read_seq.len() < self.kmer_size {
            return Vec::new();
        }

        let mut starts: HashMap<i32, PrimerStartAccumulator> = HashMap::new();

        for read_kmer_start in 0..=(read_seq.len() - self.kmer_size) {
            let Ok(word) = IntToStr::enc_bytes(
                &read_seq[read_kmer_start..read_kmer_start + self.kmer_size],
            ) else {
                continue;
            };

            let Some(info_ids) = self.table.get(&word) else {
                continue;
            };

            for &info_id in info_ids {
                let info = &self.infos[info_id];

                let primer_start =
                    read_kmer_start as i32 + info.primer_start_offset;

                let acc = starts
                    .entry(primer_start)
                    .or_insert_with(PrimerStartAccumulator::new);

                acc.add(info.kmer_start);
            }
        }

        let mut out = Vec::new();

        for (primer_start, acc) in starts {
            let aligned_primer_end = primer_start + self.primer.len() as i32;

            if aligned_primer_end <= 0 {
                continue;
            }

            let primer_end = aligned_primer_end as usize;

            out.push(PrimerStartCandidate {
                primer_start,
                primer_end,
                support: acc.support,
                first_kmer_start: acc.first_kmer_start,
                last_kmer_start: acc.last_kmer_start,
            });
        }

        out.sort_by_key(|c| c.primer_start);
        out
    }
}

#[derive(Debug, Clone)]
struct PrimerStartAccumulator {
    support: usize,
    first_kmer_start: usize,
    last_kmer_start: usize,
}

impl PrimerStartAccumulator {
    fn new() -> Self {
        Self {
            support: 0,
            first_kmer_start: usize::MAX,
            last_kmer_start: 0,
        }
    }

    fn add(&mut self, kmer_start: usize) {
        self.support += 1;
        self.first_kmer_start = self.first_kmer_start.min(kmer_start);
        self.last_kmer_start = self.last_kmer_start.max(kmer_start);
    }
}


#[cfg(test)]
mod primer_table_tests {
    use super::*;
    use int_to_str::IntToStr;

    const PRIMER: &[u8] = b"ACGTGACTGACCTGATCGTA"; // 20 bp

    fn detect(read: &[u8]) -> Vec<PrimerStartCandidate> {
        let table = PrimerTable::new(PRIMER, 4).unwrap();
        let enc = IntToStr::new(read);
        table.detect_candidates(&enc)
    }

    fn best(read: &[u8]) -> PrimerStartCandidate {
        let table = PrimerTable::new(PRIMER, 4).unwrap();
        let enc = IntToStr::new(read);
        table.detect_primer_start(&enc)
            .expect("expected primer candidate")
    }

    #[test]
    fn test_clean_single_primer_hit() {
        let read = b"TTTTTTTTTTACGTGACTGACCTGATCGTACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

        let hit = best(read);

        assert_eq!(hit.primer_start, 10);
        assert_eq!(hit.primer_end, 30);
        assert!(
            hit.support >= 14,
            "unexpectedly low support for clean primer: {hit:?}"
        );
    }

    #[test]
    fn test_partial_primer_hit_still_predicts_full_primer_end() {
        // Primer starts at 10, but the first 6 bp are missing/damaged.
        // Remaining visible part starts at read pos 16.
        let read = b"TTTTTTTTTTNNNNNNCTGACCTGATCGTACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

        let hit = best(read);

        assert_eq!(hit.primer_start, 10);
        assert_eq!(hit.primer_end, 30);
        assert!(
            hit.support >= 5,
            "partial primer should still have support: {hit:?}"
        );
    }

    #[test]
    fn test_single_mismatch_keeps_candidate_alive() {
        // One mismatch inside the primer: ACGTGACTGACC TGATCGTA
        //                          mutated: ACGTGACTGACC AGATCGTA
        let read = b"GGGGGGGGGGACGTGACTGACCAGATCGTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        let hit = best(read);

        assert_eq!(hit.primer_start, 10);
        assert_eq!(hit.primer_end, 30);
        assert!(
            hit.support >= 10,
            "single mismatch should leave enough exact 4bp words: {hit:?}"
        );
    }

    #[test]
    fn test_one_n_in_primer_region_keeps_candidate_alive() {
        // N is encoded as A by IntToStr. This is not true wildcard matching,
        // but one N should only damage local 4bp words.
        let read = b"GGGGGGGGGGACGTGACTGANCNTGATCGTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let read = read.iter().copied().filter(|b| *b != b' ').collect::<Vec<_>>();

        let hit = best(&read);

        assert_eq!(hit.primer_start, 10);
        assert_eq!(hit.primer_end, 30);
        assert!(
            hit.support >= 5,
            "single N should leave enough exact 4bp words: {hit:?}"
        );
    }

    #[test]
    fn test_two_binding_sites_are_detected() {
        let read = b"AAAAAACGTGACTGACCTGATCGTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTGACTGACCTGATCGTAGGGGGGGGGG";

        let candidates = detect(read);

        let starts: Vec<i32> = candidates
            .iter()
            .filter(|c| c.support >= 10)
            .map(|c| c.primer_start)
            .collect();

        assert!(
            starts.contains(&6),
            "missing first primer site; candidates={candidates:#?}"
        );

        assert!(
            starts.contains(&66),
            "missing second primer site; candidates={candidates:#?}"
        );
    }
    fn assert_detects_primer(
        name: &str,
        read: &[u8],
        expected_start: i32,
        expected_end: usize,
        min_support: usize,
    ) {
        let candidates = detect(read);

        let found = candidates.iter().any(|hit| {
            hit.primer_start == expected_start
                && hit.primer_end == expected_end
                && hit.support >= min_support
        });

        assert!(
            found,
            "{name}: expected primer candidate start={expected_start}, end={expected_end}, min_support={min_support}; candidates={candidates:#?}"
        );
    }

    #[test]
    fn test_three_increasingly_degenerate_10x_like_reads_still_detect_primer() {
        let reads: Vec<(&str, &[u8], i32, usize, usize)> = vec![
            (
                "clean",
                b"AAAAAACGTGACTGACCTGATCGTAAACCGGTTAACCGGTTTTTTTTTTTTTTTTTTTTTTTGGGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGG",
                6,
                26,
                10,
            ),
            (
                "one_mismatch",
                b"AAAAAACGTGACTGACCAGATCGTAAACCGGTTAACCGGTTTTTTTTTTTTTTTTTTTTTTTGGGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGG",
                6,
                26,
                5,
            ),
            (
                "partial_plus_n",
                b"AAAAAANNNNNNCTGACCTGATCGTAAACCGGTTAACCGGTTTTTTTTTTTTTTTTTTTTTTTGGGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGG",
                6,
                26,
                4,
            ),
        ];

        for (name, read, expected_start, expected_end, min_support) in reads {
            assert_detects_primer(
                name,
                read,
                expected_start,
                expected_end,
                min_support,
            );
        }
    }
}