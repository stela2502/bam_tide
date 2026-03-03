//subset_bam.rs
// bam_tide/src/subset_bam.rs

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rust_htslib::bam::record::{Aux, Record};

/// Maps tag values -> group id (outfile index), plus holds the derived outfile names.
#[derive(Debug, Default)]
pub struct Subsetter {
    tags: HashMap<String, usize>,
    group_count: usize,
    pub ofile_names: Vec<String>,
    strip_suffix: Option<String>,
}

impl Subsetter {
    pub fn new() -> Self {
        Self {
            tags: HashMap::new(),
            group_count: 0,
            ofile_names: Vec::with_capacity(128),
            strip_suffix: None,
        }
    }

    /// Optional convenience: strip a suffix from tag values before matching.
    /// Example: `.with_strip_suffix("-1")`
    pub fn with_strip_suffix(mut self, suffix: impl Into<String>) -> Self {
        self.strip_suffix = Some(suffix.into());
        self
    }

    /// Read one text file with one tag value per line.
    /// This file becomes ONE group/output BAM (group id = current group_count).
    ///
    /// `prefix` is the outfile prefix (can include a directory + stem prefix).
    /// Outfile name is: `{prefix}{file_stem}.bam`
    pub fn read_simple_list<P: AsRef<Path>>(
        &mut self,
        bc_file: P,
        prefix: &str,
    ) -> std::io::Result<()> {
        let bc_file = bc_file.as_ref();

        let file = File::open(bc_file)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let mut tag_value = line?;
            tag_value = tag_value.trim().to_string();
            if tag_value.is_empty() {
                continue;
            }

            if let Some(suf) = &self.strip_suffix {
                if let Some(stripped) = tag_value.strip_suffix(suf) {
                    tag_value = stripped.to_string();
                }
            }

            // If duplicates occur, last one would win. Better to fail fast:
            if self.tags.contains_key(&tag_value) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Duplicate tag value '{}' encountered (already assigned to a group)",
                        tag_value
                    ),
                ));
            }

            self.tags.insert(tag_value, self.group_count);
        }

        let stem = bc_file
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("Could not determine file_stem for {}", bc_file.display()),
                )
            })?;

        self.ofile_names.push(format!("{prefix}{stem}.bam"));
        self.group_count += 1;
        Ok(())
    }

    #[inline]
    pub fn num_groups(&self) -> usize {
        self.group_count
    }

    /// Returns the group id if the record tag matches.
    #[inline]
    pub fn process_record(&self, record: &Record, tag: &[u8; 2]) -> Option<usize> {
        let tag_value = get_tag_value(record, tag)?;
        self.tags.get(&tag_value).copied()
    }
}

#[inline]
fn get_tag_value(record: &Record, tag: &[u8; 2]) -> Option<String> {
    let aux = record.aux(tag).ok()?;
    match aux {
        Aux::String(s) => Some(s.to_string()),
        _ => None,
    }
}