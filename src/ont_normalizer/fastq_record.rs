//use crate::ont_normalizer::primer::{PrimerDetector, PrimerHit, Orientation};
//use mapping_info::MappingInfo;



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


#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}


impl FastqRecord {

    pub fn new(id: impl Into<String>, seq: &[u8], qual: &[u8]) -> Self {
        assert_eq!(
            seq.len(),
            qual.len(),
            "FASTQ sequence and quality lengths differ"
        );

        Self {
            id: id.into(),
            seq: seq.to_vec(),
            qual: qual.to_vec(),
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn is_empty(&self) -> bool{
        self.seq.is_empty()
    }

/*
    pub fn split_by_primer_detector(
        &self,
        detector: &PrimerDetector,
        stats: &mut MappingInfo,
    ) -> Vec<FastqRecord> {
        let hits = detector.detect_both_orientations(self, stats);

        hits.iter()
            .enumerate()
            .map(|(i, hit)| self.to_normalized_primer_read(i, hit))
            .collect()
    }

    fn to_normalized_primer_read(&self, i: usize, hit: &PrimerHit) -> FastqRecord {
        let oriented = match hit.orientation {
            Orientation::Forward => self.clone(),
            Orientation::ReverseComplement => self.revcomp(format!("{}/rc", self.id)),
        };

        let cell = String::from_utf8_lossy(&hit.cell_id.unwrap_or(FastqRecord::default()).seq);
        let umi = String::from_utf8_lossy(&hit.umi.unwrap_or(FastqRecord::default()).seq);

        oriented.clipped(
            format!("{}_{} CB:{} UMI:{}", self.id, i, cell, umi),
            hit.read_start,
            hit.segment_end,
        )
    }
*/    

    pub fn from_bam_record(rec: &rust_htslib::bam::Record) -> Self {
        let seq = rec.seq().as_bytes();
        Self::new(
            String::from_utf8_lossy(rec.qname()).to_string(),
            &seq,
            rec.qual(),
        )
    }

    pub fn clipped(&self, id: impl Into<String>, start: usize, end: usize) -> Self {
        assert!(start <= end, "invalid FASTQ clip range: {start}..{end}");
        assert!(end <= self.seq.len(), "FASTQ clip end outside sequence");

        Self::new(id, &self.seq[start..end], &self.qual[start..end])
    }

    pub fn qual_string(&self) -> String {
        self.qual
            .iter()
            .map(|q| q.saturating_add(33) as char)
            .collect()
    }

    pub fn revcomp(&self, id: impl Into<String>) -> Self {
        let mut rc_seq = Vec::with_capacity(self.seq.len());
        let mut rc_qual = Vec::with_capacity(self.qual.len());

        for base in self.seq.iter().rev() {
            rc_seq.push(Self::complement(*base));
        }

        for q in self.qual.iter().rev() {
            rc_qual.push(*q);
        }

        Self::new(id, &rc_seq, &rc_qual)
    }

    #[inline(always)]
    fn complement(base: u8) -> u8 {
        DNA_COMPLEMENT[base as usize]
    }

    pub fn write_to<W: std::io::Write>(
        &self,
        writer: &mut W,
        qual_buf: &mut Vec<u8>,
    ) -> std::io::Result<()> {
        writer.write_all(b"@")?;
        writer.write_all(self.id.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&self.seq)?;
        writer.write_all(b"\n+\n")?;

        qual_buf.clear();
        qual_buf.extend(self.qual.iter().map(|q| q.saturating_add(33)));

        writer.write_all(qual_buf)?;
        writer.write_all(b"\n")?;

        Ok(())
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_record_basic() {
        let seq = b"ACGT";
        let qual = &[30, 31, 32, 33]; // Phred scores

        let rec = FastqRecord::new("read1", seq, qual);

        assert_eq!(rec.id, "read1");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, qual);
    }

    #[test]
    fn test_fastq_record_clipping() {
        let rec = FastqRecord::new(
            "read1",
            b"ACGTACGT",
            &[10, 11, 12, 13, 14, 15, 16, 17],
        );

        let clipped = rec.clipped("read1/mol1", 2, 6);

        assert_eq!(clipped.id, "read1/mol1");
        assert_eq!(clipped.seq, b"GTAC");
        assert_eq!(clipped.qual, vec![12, 13, 14, 15]);
    }

    #[test]
    fn test_fastq_record_write() {
        let rec = FastqRecord::new(
            "read1",
            b"ACGT",
            &[0, 1, 2, 3], // will become ! " # $
        );

        let mut out = Vec::new();
        let mut qual_buf = Vec::new();

        rec.write_to(&mut out, &mut qual_buf).unwrap();

        let result = String::from_utf8(out).unwrap();

        assert_eq!(
            result,
            "@read1\nACGT\n+\n!\"#$\n"
        );
    }

    #[test]
    #[should_panic(expected = "FASTQ sequence and quality lengths differ")]
    fn test_fastq_len_mismatch_panics() {
        FastqRecord::new("bad", b"ACGT", &[1, 2, 3]);
    }
    #[test]
    fn test_fastq_record_revcomp() {
        let rec = FastqRecord::new(
            "read1",
            b"ACGTNacgtn",
            &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        );

        let rc = rec.revcomp("read1/rc");

        assert_eq!(rc.id, "read1/rc");
        assert_eq!(rc.seq, b"nacgtNACGT");
        assert_eq!(rc.qual, vec![10, 9, 8, 7, 6, 5, 4, 3, 2, 1]);
    }
}