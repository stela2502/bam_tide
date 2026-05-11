use crate::ont_normalizer::fastq_record::FastqRecord;
use crate::ont_normalizer::primer::model::{
    AdapterHit, BarcodeMatcherSpec, MatchState, Orientation, PrimerHit, PrimerPart,
    PrimerSplit, PrimerStructure,
};
use anyhow::{bail, Result};
use mapping_info::MappingInfo;

#[derive(Debug, Clone)]
pub struct PrimerDetector {
    pub structures: Vec<PrimerStructure>,
}

impl PrimerDetector {
    pub fn new(structures: Vec<PrimerStructure>) -> Self {
        Self { structures }
    }

    pub fn tenx_v3() -> Self {
        Self::new(vec![PrimerStructure::tenx_v3()])
    }

    pub fn split_read(
        &self,
        read: &FastqRecord,
        stats: &mut MappingInfo,
    ) -> Result<Vec<PrimerSplit>> {
        match self.split_oriented_read(read, Orientation::Forward, stats) {
            Ok(splits) => Ok(splits),
            Err(_) => {
                let rc = read.revcomp(format!("{}/rc", read.id));
                self.split_oriented_read(&rc, Orientation::ReverseComplement, stats)
            }
        }
    }

    fn split_oriented_read(
        &self,
        read: &FastqRecord,
        orientation: Orientation,
        stats: &mut MappingInfo,
    ) -> Result<Vec<PrimerSplit>> {

        #[derive(Debug)]
        struct PendingPrimer {
            insert_start: usize,
            cell_id: Option<FastqRecord>,
            umi: Option<FastqRecord>,
        }

        let mut out = Vec::new();
        let mut cursor = 0usize;
        let mut pending: Option<PendingPrimer> = None;

        while let Some(hit) =
            self.detect_next_primer(read, cursor, orientation, stats)
        {
            if let Some(prev) = pending.take() {
                if hit.primer_start > prev.insert_start {
                    out.push(PrimerSplit {
                        original_read_id: read.id.clone(),

                        insert: read.clipped(
                            format!("{}/mol{}", read.id, out.len()),
                            prev.insert_start,
                            hit.primer_start,
                        ),

                        cell_id: prev.cell_id,
                        umi: prev.umi,

                        orientation,
                    });
                } else {
                    stats.report("empty_insert_before_next_primer");
                }
            } else if hit.primer_start > cursor {
                stats.report("prefix_before_first_primer");
            }

            pending = Some(PendingPrimer {
                insert_start: hit.insert_start,
                cell_id: hit.cell_id,
                umi: hit.umi,
            });

            cursor = hit.insert_start;

            if cursor >= read.len() {
                break;
            }
        }
        

        if let Some(prev) = pending.take() {
            if prev.insert_start < read.len() {
                out.push(PrimerSplit {
                    original_read_id: read.id.clone(),

                    insert: read.clipped(
                        format!("{}/mol{}", read.id, out.len()),
                        prev.insert_start,
                        read.len(),
                    ),

                    cell_id: prev.cell_id,
                    umi: prev.umi,

                    orientation,
                });
            } else {
                stats.report("empty_final_insert");
            }
        }

        let orient = match orientation {
            Orientation::Forward => "forward_splits",
            Orientation::ReverseComplement => "reverse_splits",
        };
        if out.is_empty() {
            bail!("no primer-compatible insert detected in read {}", read.id);
        }

        *stats
            .error_counts
            .entry("emitted_splits".to_string())
            .or_insert(0) +=  out.len();
        *stats
            .error_counts
            .entry(orient.to_string())
            .or_insert(0) +=  out.len();

        if out.len() > 1 {
            stats.report("multi_split_reads");
        }

        Ok(out)
    }

    fn detect_next_primer(
        &self,
        read: &FastqRecord,
        cursor: usize,
        orientation: Orientation,
        stats: &mut MappingInfo,
    ) -> Option<PrimerHit> {
        let mut best: Option<PrimerHit> = None;

        for structure in &self.structures {
            if let Some(hit) = structure.detect_next(read, cursor, orientation, stats) {
                let replace = best
                    .as_ref()
                    .map(|old| hit.primer_start < old.primer_start)
                    .unwrap_or(true);

                if replace {
                    best = Some(hit);
                }
            }
        }

        best
    }
}

impl PrimerStructure {
    pub fn detect_next(
        &self,
        read: &FastqRecord,
        cursor: usize,
        orientation: Orientation,
        stats: &mut MappingInfo,
    ) -> Option<PrimerHit> {
        let adapter_hit = self.find_next_adapter(read, cursor)?;

        self.match_after_adapter(read, adapter_hit, orientation, stats)
    }

    fn find_next_adapter(&self, read: &FastqRecord, cursor: usize) -> Option<AdapterHit> {
        if self.adapter.is_empty() || read.is_empty() || cursor >= read.len() {
            return None;
        }

        let min_len = self.adapter_suffix_min_len.min(self.adapter.len());
        let mut best: Option<AdapterHit> = None;

        for pos in cursor..read.len() {
            for suffix_len in (min_len..=self.adapter.len()).rev() {
                if pos + suffix_len > read.len() {
                    continue;
                }

                let suffix = &self.adapter[self.adapter.len() - suffix_len..];
                let candidate = &read.seq[pos..pos + suffix_len];

                if mismatches_n_tolerant(candidate, suffix) <= self.max_adapter_mismatches {
                    let hit = AdapterHit {
                        primer_start: pos,
                        adapter_end: pos + suffix_len,
                    };

                    let replace = best
                        .as_ref()
                        .map(|old| hit.primer_start < old.primer_start)
                        .unwrap_or(true);

                    if replace {
                        best = Some(hit);
                    }
                }
            }

            if best.is_some() {
                break;
            }
        }

        best
    }

    fn match_after_adapter(
        &self,
        read: &FastqRecord,
        adapter_hit: AdapterHit,
        orientation: Orientation,
        stats: &mut MappingInfo,
    ) -> Option<PrimerHit> {
        let mut pos = adapter_hit.adapter_end;
        let mut state = MatchState::default();

        for part in &self.parts {
            pos = part.match_at(read, pos, &mut state, stats)?;

            if matches!(part, PrimerPart::Insert) {
                break;
            }
        }

        let insert_start = state.insert_start.unwrap_or(pos);

        Some(PrimerHit {
            structure_name: self.name.clone(),
            orientation,
            primer_start: adapter_hit.primer_start,
            primer_end: insert_start,
            insert_start,
            cell_id: state.cell_id(),
            umi: state.umi(),
        })
    }
}

impl PrimerPart {
    fn match_at(
        &self,
        read: &FastqRecord,
        pos: usize,
        state: &mut MatchState,
        stats: &mut MappingInfo,
    ) -> Option<usize> {
        match self {
            PrimerPart::Fixed {
                name,
                seq,
                max_mismatches,
            } => {
                if pos + seq.len() > read.len() {
                    stats.report(&format!("too_short_before_fixed_{name}"));
                    return None;
                }

                let candidate = &read.seq[pos..pos + seq.len()];

                if mismatches_n_tolerant(candidate, seq) > *max_mismatches {
                    stats.report(&format!("failed_fixed_{name}"));
                    return None;
                }

                Some(pos + seq.len())
            }

            PrimerPart::Random { .. } => {
                panic!("PrimerPart::Random is not implemented yet");
            }

            PrimerPart::CellId { name, len, matcher } => {
                if pos + len > read.len() {
                    stats.report(&format!("too_short_cell_id_{name}"));
                    return None;
                }

                match matcher {
                    BarcodeMatcherSpec::Any => {}
                    BarcodeMatcherSpec::Whitelist { .. } => {
                        panic!("Barcode whitelist matching is not implemented yet");
                    }
                }

                state.cell_seq.extend_from_slice(&read.seq[pos..pos + len]);
                state.cell_qual.extend_from_slice(&read.qual[pos..pos + len]);

                Some(pos + len)
            }

            PrimerPart::Umi { name, len } => {
                if pos + len > read.len() {
                    stats.report(&format!("too_short_umi_{name}"));
                    return None;
                }

                state.umi_seq.extend_from_slice(&read.seq[pos..pos + len]);
                state.umi_qual.extend_from_slice(&read.qual[pos..pos + len]);

                Some(pos + len)
            }

            PrimerPart::PolyT { min_len, max_non_t } => {
                let Some(poly_len) = Self::match_polyt(read, pos, *min_len, *max_non_t) else {
                    stats.report("failed_poly_t");
                    return None;
                };

                Some(pos + poly_len)
            }

            PrimerPart::Insert => {
                state.insert_start = Some(pos);
                Some(pos)
            }
        }
    }

    fn match_polyt(
        read: &FastqRecord,
        start: usize,
        min_len: usize,
        _max_non_t: usize,
    ) -> Option<usize> {
        if start >= read.len() {
            return None;
        }

        let mut len = 0usize;

        for &b in &read.seq[start..] {
            let u = b.to_ascii_uppercase();

            if u == b'T' || u == b'N' {
                len += 1;
            } else {
                break;
            }
        }

        (len >= min_len).then_some(len)
    }
}

impl MatchState {
    fn cell_id(&self) -> Option<FastqRecord> {
        (!self.cell_seq.is_empty()).then(|| {
            FastqRecord::new("cell_id", &self.cell_seq, &self.cell_qual)
        })
    }

    fn umi(&self) -> Option<FastqRecord> {
        (!self.umi_seq.is_empty()).then(|| {
            FastqRecord::new("umi", &self.umi_seq, &self.umi_qual)
        })
    }
}

fn mismatches_n_tolerant(a: &[u8], b: &[u8]) -> usize {
    assert_eq!(a.len(), b.len());

    a.iter()
        .zip(b.iter())
        .filter(|(x, y)| {
            let x = x.to_ascii_uppercase();
            let y = y.to_ascii_uppercase();
            x != b'N' && y != b'N' && x != y
        })
        .count()
}


#[cfg(test)]
mod tests {
    use super::*;

    fn q(len: usize) -> Vec<u8> {
        vec![30; len]
    }

    fn read(id: &str, seq: &[u8]) -> FastqRecord {
        FastqRecord::new(id, seq, &q(seq.len()))
    }

    fn det() -> PrimerDetector {
        PrimerDetector::tenx_v3()
    }

    #[test]
    fn splits_single_10x_forward() {
        let rec = read(
            "read1",
            b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTACGTACGT",
        );

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rec, &mut stats).unwrap();

        assert_eq!(splits.len(), 1);
        assert_eq!(splits[0].original_read_id, "read1");
        assert_eq!(splits[0].orientation, Orientation::Forward);
        assert_eq!(splits[0].insert.seq, b"ACGTACGT");
        assert_eq!(splits[0].cell_id.as_ref().unwrap().seq, b"GCATTAACAATAGACC");
        assert_eq!(splits[0].umi.as_ref().unwrap().seq, b"TGTTGGAGACGC");

        assert_eq!(stats.get_issue_count("emitted_splits"), 1);
        assert_eq!(stats.get_issue_count("forward_splits"), 1);
    }

    #[test]
    fn splits_two_concatenated_molecules() {
        let mol1 = b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mol2 = b"GCTCTTCCGATCTGGGTTTGCAGCCTAAAGGTAAGGGCTCATTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let seq = [mol1.as_slice(), mol2.as_slice()].concat();

        let rec = read("chimera1", &seq);

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rec, &mut stats).unwrap();

        assert_eq!(splits.len(), 2);

        assert_eq!(splits[0].insert.seq, b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); // shorter as the tnext tests clips the full theoretical primer length here..
        assert_eq!(splits[0].cell_id.as_ref().unwrap().seq, b"GCATTAACAATAGACC");
        assert_eq!(splits[0].umi.as_ref().unwrap().seq, b"TGTTGGAGACGC");

        assert_eq!(splits[1].insert.seq, b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
        assert_eq!(splits[1].cell_id.as_ref().unwrap().seq, b"GGGTTTGCAGCCTAAA");
        assert_eq!(splits[1].umi.as_ref().unwrap().seq, b"GGTAAGGGCTCA");

        assert_eq!(stats.get_issue_count("emitted_splits"), 2, "{splits:?}\n{stats}");
        assert_eq!(stats.get_issue_count("multi_split_reads"), 1);
    }

    #[test]
    fn reverse_complement_is_normalized() {
        let forward = read(
            "read1",
            b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTACGTACGT",
        );

        let rc = forward.revcomp("read1_rc");

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rc, &mut stats).unwrap();

        assert_eq!(splits.len(), 1);
        assert_eq!(splits[0].orientation, Orientation::ReverseComplement);
        assert_eq!(splits[0].insert.seq, b"ACGTACGT");
        assert_eq!(splits[0].cell_id.as_ref().unwrap().seq, b"GCATTAACAATAGACC");
        assert_eq!(splits[0].umi.as_ref().unwrap().seq, b"TGTTGGAGACGC");

        assert_eq!(stats.get_issue_count("reverse_splits"), 1);
    }

    #[test]
    fn suffix_adapter_detection_works() {
        let rec = read(
            "suffix",
            b"GCTCTTCCGATCTGGGTTTGCAGCCTAAAGGTAAGGGCTCATTTTTTTTTTTTTTGGGGCCCC",
        );

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rec, &mut stats).unwrap();

        assert_eq!(splits.len(), 1);
        assert_eq!(splits[0].insert.seq, b"GGGGCCCC");
        assert_eq!(splits[0].cell_id.as_ref().unwrap().seq, b"GGGTTTGCAGCCTAAA");
        assert_eq!(splits[0].umi.as_ref().unwrap().seq, b"GGTAAGGGCTCA");
    }

    #[test]
    fn reports_prefix_before_first_primer() {
        let rec = read(
            "prefix",
            b"NNNNCTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTACGT",
        );

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rec, &mut stats).unwrap();

        assert_eq!(splits.len(), 1);
        assert_eq!(splits[0].insert.seq, b"ACGT");
        assert_eq!(stats.get_issue_count("prefix_before_first_primer"), 1);
    }

    #[test]
    fn fails_without_polyt() {
        let rec = read(
            "no_polyt",
            b"CTACACGACGCTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCACGTACGTACGT",
        );

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let result = det().split_read(&rec, &mut stats);

        assert!(result.is_err());
        assert!(stats.get_issue_count("failed_poly_t") > 0);
    }

    #[test]
    fn fails_without_primer() {
        let rec = read("junk", b"ACGTACGTACGTACGT");

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let result = det().split_read(&rec, &mut stats);

        assert!(result.is_err());
    }

    #[test]
    fn accepts_n_inside_adapter() {
        let rec = read(
            "adapter_n",
            b"CTACACGACGNTCTTCCGATCTGCATTAACAATAGACCTGTTGGAGACGCTTTTTTTTTTTTTTACGT",
        );

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let splits = det().split_read(&rec, &mut stats).unwrap();

        assert_eq!(splits.len(), 1);
        assert_eq!(splits[0].insert.seq, b"ACGT");
        assert_eq!(splits[0].cell_id.as_ref().unwrap().seq, b"GCATTAACAATAGACC");
    }

    #[test]
    #[should_panic(expected = "PrimerPart::Random is not implemented yet")]
    fn random_part_panics_for_now() {
        let structure = PrimerStructure {
            name: "random-test".to_string(),
            adapter: b"AAA".to_vec(),
            adapter_suffix_min_len: 3,
            max_adapter_mismatches: 0,
            parts: vec![
                PrimerPart::Random {
                    name: "stagger".to_string(),
                    min_len: 1,
                    max_len: 3,
                },
                PrimerPart::CellId {
                    name: "CB".to_string(),
                    len: 4,
                    matcher: BarcodeMatcherSpec::Any,
                },
                PrimerPart::Umi {
                    name: "UMI".to_string(),
                    len: 4,
                },
                PrimerPart::PolyT {
                    min_len: 4,
                    max_non_t: 0,
                },
                PrimerPart::Insert,
            ],
        };

        let detector = PrimerDetector::new(vec![structure]);
        let rec = read("random", b"AAANACGTGGCCTTTTAAAA");

        let mut stats = MappingInfo::new(None, 0.0, 0);
        let _ = detector.split_read(&rec, &mut stats);
    }
}