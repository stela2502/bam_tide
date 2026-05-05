const DNA_COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];

    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }

    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b't';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'a';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'c';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'g';
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';

    table
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

impl Orientation {
    pub fn as_str(self) -> &'static str {
        match self {
            Orientation::Forward => "forward",
            Orientation::ReverseComplement => "reverse_complement",
        }
    }
}

#[derive(Debug, Clone)]
pub struct AdapterHit {
    pub start: usize,
    pub end: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone)]
pub struct Cassette {
    pub orientation: Orientation,
    pub adapter_start: usize,
    pub adapter_end: usize,
    pub cb: Vec<u8>,
    pub cb_qual: Vec<u8>,
    pub umi: Vec<u8>,
    pub umi_qual: Vec<u8>,
    pub poly_t_start: usize,
    pub poly_t_len: usize,
    pub segment_start: usize,
    pub segment_end: usize,
}

#[derive(Debug, Clone)]
pub struct CassetteExtractor {
    adapter: Vec<u8>,
    min_adapter_match: usize,
    cb_len: usize,
    umi_len: usize,
    poly_t_min: usize,
    poly_t_window: usize,
    min_transcript_len: usize,
    max_adapter_mismatches: usize,
}

impl CassetteExtractor {
    pub fn new(
        adapter: Vec<u8>,
        min_adapter_match: usize,
        cb_len: usize,
        umi_len: usize,
        poly_t_min: usize,
        poly_t_window: usize,
        min_transcript_len: usize,
        max_adapter_mismatches: usize,
    ) -> Self {
        Self {
            adapter,
            min_adapter_match,
            cb_len,
            umi_len,
            poly_t_min,
            poly_t_window,
            min_transcript_len,
            max_adapter_mismatches,
        }
    }

    pub fn extract_both_orientations(
        &self,
        seq: &[u8],
        qual: &[u8],
        too_short_after_adapter: &mut u64,
        failed_poly_t: &mut u64,
    ) -> Vec<Cassette> {
        let mut out = Vec::new();

        out.extend(self.extract_orientation(
            seq,
            qual,
            Orientation::Forward,
            too_short_after_adapter,
            failed_poly_t,
        ));

        let (rc_seq, rc_qual) = self.revcomp_with_qual(seq, qual);

        out.extend(self.extract_orientation(
            &rc_seq,
            &rc_qual,
            Orientation::ReverseComplement,
            too_short_after_adapter,
            failed_poly_t,
        ));

        out
    }

    fn find_adapter_hits(&self, seq: &[u8]) -> Vec<AdapterHit> {
        let mut hits = Vec::new();

        if self.adapter.is_empty() || self.min_adapter_match == 0 {
            return hits;
        }

        let min_len = self.min_adapter_match.min(self.adapter.len());
        let adapter_len = self.adapter.len();

        for pos in 0..seq.len() {
            let mut best_hit: Option<AdapterHit> = None;

            for len in (min_len..=adapter_len).rev() {
                if pos + len > seq.len() {
                    continue;
                }

                let adapter_suffix = &self.adapter[adapter_len - len..];

                let Some(mismatches) =
                    self.adapter_window_matches(&seq[pos..pos + len], adapter_suffix)
                else {
                    continue;
                };

                let suffix_offset = adapter_len - len;

                let adapter_start = pos.saturating_sub(suffix_offset);

                // Important:
                // the adapter end is fixed by the matched suffix end.
                // If we matched the last `len` bases of the adapter,
                // then pos + len is the expected full adapter end.
                let adapter_end = pos + len;

                best_hit = Some(AdapterHit {
                    start: adapter_start,
                    end: adapter_end,
                    mismatches,
                });

                break;
            }

            if let Some(hit) = best_hit {
                hits.push(hit);
            }
        }

        self.remove_overlapping_hits(hits)
    }

    fn adapter_window_matches(&self, read_window: &[u8], adapter_window: &[u8]) -> Option<usize> {
        if read_window.len() != adapter_window.len() {
            return None;
        }

        let mut mismatches = 0usize;

        for (read_base, adapter_base) in read_window.iter().zip(adapter_window.iter()) {
            if read_base.eq_ignore_ascii_case(adapter_base) {
                continue;
            }

            // ONT ambiguity: treat N as unknown, not as evidence against the hit.
            if read_base.eq_ignore_ascii_case(&b'N') || adapter_base.eq_ignore_ascii_case(&b'N') {
                continue;
            }

            mismatches += 1;

            if mismatches > self.max_adapter_mismatches {
                return None;
            }
        }

        Some(mismatches)
    }

    pub fn revcomp_with_qual(&self, seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut rc_seq = Vec::with_capacity(seq.len());
        let mut rc_qual = Vec::with_capacity(qual.len());

        for base in seq.iter().rev() {
            rc_seq.push(Self::complement(*base));
        }

        for q in qual.iter().rev() {
            rc_qual.push(*q);
        }

        (rc_seq, rc_qual)
    }

    fn extract_orientation(
        &self,
        seq: &[u8],
        qual: &[u8],
        orientation: Orientation,
        too_short_after_adapter: &mut u64,
        failed_poly_t: &mut u64,
    ) -> Vec<Cassette> {
        let hits = self.find_adapter_hits(seq);
        let mut out = Vec::new();

        for (i, hit) in hits.iter().enumerate() {
            let cb_start = hit.end;
            let cb_end = cb_start + self.cb_len;
            let umi_start = cb_end;
            let umi_end = umi_start + self.umi_len;
            let poly_t_start = umi_end;

            if poly_t_start + self.poly_t_window > seq.len() {
                *too_short_after_adapter += 1;
                continue;
            }

            let Some(poly_t_start) = self.find_poly_t_start(seq, poly_t_start) else {
                *failed_poly_t += 1;
                continue;
            };

            let t_count = seq[poly_t_start..poly_t_start + self.poly_t_window]
                .iter()
                .filter(|base| base.eq_ignore_ascii_case(&b'T'))
                .count();

            if t_count < self.poly_t_min {
                *failed_poly_t += 1;
                continue;
            }

            let poly_t_len = seq[poly_t_start..]
                .iter()
                .take_while(|base| base.eq_ignore_ascii_case(&b'T'))
                .count();

            let after_polyt = poly_t_start + poly_t_len;
            let transcript_len = seq.len().saturating_sub(after_polyt);

            if transcript_len < self.min_transcript_len {
                *too_short_after_adapter += 1;
                continue;
            }

            let segment_end = hits.get(i + 1).map(|next| next.start).unwrap_or(seq.len());

            if segment_end <= hit.start {
                *too_short_after_adapter += 1;
                continue;
            }

            out.push(Cassette {
                orientation,
                adapter_start: hit.start,
                adapter_end: hit.end,
                cb: seq[cb_start..cb_end].to_vec(),
                cb_qual: qual[cb_start..cb_end].to_vec(),
                umi: seq[umi_start..umi_end].to_vec(),
                umi_qual: qual[umi_start..umi_end].to_vec(),
                poly_t_start,
                poly_t_len,
                segment_start: hit.start,
                segment_end,
            });
        }

        out
    }

    fn remove_overlapping_hits(&self, hits: Vec<AdapterHit>) -> Vec<AdapterHit> {
        let mut filtered: Vec<AdapterHit> = Vec::new();

        for hit in hits {
            if let Some(last) = filtered.last_mut()
                && hit.start < last.end
            {
                let hit_len = hit.end - hit.start;
                let last_len = last.end - last.start;

                if hit_len > last_len {
                    *last = hit;
                }

                continue;
            }

            filtered.push(hit);
        }

        filtered
    }

    #[inline(always)]
    fn complement(base: u8) -> u8 {
        DNA_COMPLEMENT[base as usize]
    }

    fn find_poly_t_start(&self, seq: &[u8], expected_start: usize) -> Option<usize> {
        if expected_start + self.poly_t_window > seq.len() {
            return None;
        }

        let t_count = seq[expected_start..expected_start + self.poly_t_window]
            .iter()
            .filter(|base| base.eq_ignore_ascii_case(&b'T'))
            .count();

        if t_count >= self.poly_t_min {
            Some(expected_start)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn extractor() -> CassetteExtractor {
        CassetteExtractor::new(b"CTACACGACGCTCTTCCGATCT".to_vec(), 13, 16, 12, 10, 14, 4, 2)
    }

    fn q(len: usize) -> Vec<u8> {
        vec![30; len]
    }

    #[test]
    fn extracts_full_adapter_example() {
        let ex = extractor();

        let seq = b"NNNNCTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTTTTTTTTTTTACGTACGTACGT";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            1,
            "expected one cassette, got {cassettes:#?}"
        );
        assert_eq!(too_short, 0);
        assert_eq!(failed_poly_t, 0);

        let c = &cassettes[0];

        assert_eq!(c.orientation, Orientation::Forward);
        assert_eq!(c.adapter_start, 4);
        assert_eq!(c.adapter_end, 26);
        assert_eq!(&c.cb, b"GCATTAACAATAGACC");
        assert_eq!(&c.umi, b"TGTTGGAGACGC");
        assert_eq!(c.poly_t_start, 54);
        assert!(c.poly_t_len >= 24);
        assert_eq!(c.segment_start, c.adapter_start);
        assert_eq!(c.segment_end, seq.len());
    }

    #[test]
    fn extracts_adapter_suffix_example() {
        let ex = extractor();

        let seq =
            b"AAAAGCTCTTCCGATCTGGGTTTGCAGCCTAAAGGTAAGGGCTCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGCCCC";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            1,
            "expected one suffix cassette, got {cassettes:#?}"
        );

        let c = &cassettes[0];

        assert_eq!(c.orientation, Orientation::Forward);
        assert_eq!(c.adapter_start, 0);
        assert_eq!(c.adapter_end, 17);
        assert_eq!(&c.cb, b"GGGTTTGCAGCCTAAA");
        assert_eq!(&c.umi, b"GGTAAGGGCTCA");
        assert_eq!(c.poly_t_start, 45);
        assert!(c.poly_t_len >= 30);
    }

    #[test]
    fn extracts_reverse_complement_cassette() {
        let ex = extractor();

        let forward_molecule =
            b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTTTTTTTTTTTACGTACGT";
        let forward_qual = q(forward_molecule.len());

        let (reverse_seq, reverse_qual) = ex.revcomp_with_qual(forward_molecule, &forward_qual);

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes = ex.extract_both_orientations(
            &reverse_seq,
            &reverse_qual,
            &mut too_short,
            &mut failed_poly_t,
        );

        assert_eq!(
            cassettes.len(),
            1,
            "expected one reverse-complement cassette, got {cassettes:#?}"
        );

        let c = &cassettes[0];

        assert_eq!(c.orientation, Orientation::ReverseComplement);
        assert_eq!(c.adapter_start, 0);
        assert_eq!(c.adapter_end, 22);
        assert_eq!(&c.cb, b"GCATTAACAATAGACC");
        assert_eq!(&c.umi, b"TGTTGGAGACGC");
    }

    #[test]
    fn splits_two_cassettes_from_one_read() {
        let ex = extractor();

        let mol1 = b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTAAAAAA";
        let mol2 = b"GCTCTTCCGATCTGGGTTTGCAGCCTAAAGGTAAGGGCTCATTTTTTTTTTTTTTCCCCCC";
        let seq = [mol1.as_slice(), mol2.as_slice()].concat();
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(&seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            2,
            "expected two cassettes, got {cassettes:#?}"
        );

        assert_eq!(cassettes[0].orientation, Orientation::Forward);
        assert_eq!(cassettes[1].orientation, Orientation::Forward);

        assert_eq!(&cassettes[0].cb, b"GCATTAACAATAGACC");
        assert_eq!(&cassettes[0].umi, b"TGTTGGAGACGC");

        assert_eq!(&cassettes[1].cb, b"GGGTTTGCAGCCTAAA");
        assert_eq!(&cassettes[1].umi, b"GGTAAGGGCTCA");

        assert_eq!(cassettes[0].segment_end, cassettes[1].adapter_start);
        assert_eq!(cassettes[1].segment_end, seq.len());
    }

    #[test]
    fn counts_failed_polyt() {
        let ex = extractor();

        let seq = b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCACGTACGTACGTACGT";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(cassettes.len(), 0);
        assert_eq!(too_short, 0);
        assert_eq!(failed_poly_t, 1);
    }

    #[test]
    fn counts_too_short_after_adapter() {
        let ex = extractor();

        let seq = b"CTACACGACGCTCTTCCGATCTGCATTAAC";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(cassettes.len(), 0);
        assert_eq!(too_short, 1);
        assert_eq!(failed_poly_t, 0);
    }

    #[test]
    fn accepts_ns_inside_full_adapter() {
        let ex = extractor();

        // Full adapter:
        // CTACACGACGCTCTTCCGATCT
        //
        // Corrupted with Ns:
        // CTACACGACGNTCTTCCGATCT
        let seq = b"AAAAC T A C A C G A C G N T C T T C C G A T C TGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTACGTACGT"
            .iter()
            .copied()
            .filter(|b| *b != b' ')
            .collect::<Vec<u8>>();

        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(&seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            1,
            "N inside adapter should be tolerated, got {cassettes:#?}"
        );

        let c = &cassettes[0];
        assert_eq!(c.orientation, Orientation::Forward);
        assert_eq!(&c.cb, b"GCATTAACAATAGACC");
        assert_eq!(&c.umi, b"TGTTGGAGACGC");
    }

    #[test]
    fn accepts_ns_inside_suffix_adapter() {
        let ex = extractor();

        // Suffix adapter:
        // GCTCTTCCGATCT
        //
        // Corrupted with N:
        // GCTCTTCNGATCT
        let seq = b"AAAAGCTCTTCNGATCTGGGTTTGCAGCCTAAAGGTAAGGGCTCATTTTTTTTTTTTTTGGGGCCCC";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            1,
            "N inside suffix adapter should be tolerated, got {cassettes:#?}"
        );

        let c = &cassettes[0];
        assert_eq!(c.orientation, Orientation::Forward);
        assert_eq!(&c.cb, b"GGGTTTGCAGCCTAAA");
        assert_eq!(&c.umi, b"GGTAAGGGCTCA");
    }

    #[test]
    fn rejects_polyt_too_soon_after_adapter() {
        let ex = extractor();

        // Adapter followed almost immediately by polyT.
        // This should NOT be accepted as cassette, because there is no room for CB+UMI.
        let seq = b"CTACACGACGCTCTTCCGATCTTTTTTTTTTTTTTACGTACGTACGT";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(cassettes.len(), 0);
        assert!(
            too_short > 0 || failed_poly_t > 0,
            "polyT too soon must be counted as a failure"
        );
    }

    #[test]
    fn rejects_no_gene_sequence_after_polyt() {
        let ex = extractor();

        // Valid adapter + CB + UMI + polyT, but no transcript bases after polyT.
        let seq = b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTT";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(
            cassettes.len(),
            0,
            "cassette without transcript sequence after polyT must be rejected"
        );

        assert!(
            too_short > 0,
            "missing transcript sequence after polyT should increment too_short_after_adapter"
        );
    }

    #[test]
    fn rejects_missing_polyt_after_cb_umi() {
        let ex = extractor();

        // Valid adapter + CB + UMI, but then transcript-like sequence instead of polyT.
        let seq = b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCACGTACGTACGTACGT";
        let qual = q(seq.len());

        let mut too_short = 0;
        let mut failed_poly_t = 0;

        let cassettes =
            ex.extract_both_orientations(seq, &qual, &mut too_short, &mut failed_poly_t);

        assert_eq!(cassettes.len(), 0);
        assert_eq!(too_short, 0);
        assert_eq!(
            failed_poly_t, 1,
            "missing polyT after CB+UMI should increment failed_poly_t"
        );
    }
}
