use std::fmt;
use std::fmt::Display;


#[derive(Debug, PartialEq, Clone, PartialOrd, Hash)]
pub struct ChrArea{
    pub start:usize,
    pub end:usize,
    pub chr:String
}

impl Display for ChrArea{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end )?;
        Ok(())
    }
}

impl Ord for ChrArea {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First, compare by chr
        let chr_cmp = self.chr.cmp(&other.chr);
        if chr_cmp != std::cmp::Ordering::Equal {
            return chr_cmp;
        }
        // If chr is the same, compare by start position
        self.start.cmp(&other.start)
    }
}

impl Eq for ChrArea {}  // Implementing Eq is necessary for Ord to work

impl ChrArea {
    pub fn match_pos( &self, pos:usize ) -> bool {
        self.start <= pos && pos < self.end
    }
}