use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam;
//use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{Read, Writer};

//use std::collections::BTreeMap;
use std::fs::{self};
//use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use bam_tide::multi_subset_bam::Subsetter;

/// Split a BAM file into multiple BAMs based on tag values.
///
/// Each input list file defines one output BAM.
/// Each line in the list file is one accepted tag value.
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    bam: PathBuf,

    /// BAM tag to match (must be two characters, e.g. CB, CR, UB)
    #[arg(short, long, default_value = "CR")]
    tag: String,

    /// One or more files containing tag values (one per line)
    #[arg(short = 'l', long = "list", required = true)]
    lists: Vec<PathBuf>,

    /// Output prefix (directory + filename prefix)
    #[arg(short, long)]
    prefix: String,

    /// Also write unmatched reads
    #[arg(long)]
    write_unmatched: bool,
}

fn parse_tag(tag: &str) -> Result<[u8; 2]> {
    let bytes = tag.as_bytes();
    if bytes.len() != 2 {
        anyhow::bail!("Tag must be exactly 2 characters (got '{}')", tag);
    }
    Ok([bytes[0], bytes[1]])
}

fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("creating directory {}", parent.display()))?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let now = SystemTime::now();
    let cli = Cli::parse();

    let tag = parse_tag(&cli.tag)?;

    let mut reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("opening BAM {}", cli.bam.display()))?;

    let header = bam::Header::from_template(reader.header());

    let mut subsetter = Subsetter::new();

    for list in &cli.lists {
        subsetter.read_simple_list(list, &cli.prefix)?;
    }

    // Prepare output writers
    let mut writers: Vec<Writer> = Vec::new();

    for path in &subsetter.ofile_names {
        ensure_parent_dir(path)?;
        writers.push(
            Writer::from_path(path, &header, bam::Format::Bam)
                .with_context(|| format!("creating {}", path.display()))?,
        );
    }

    // Optional unmatched writer
    let mut unmatched_writer = if cli.write_unmatched {
        let path = PathBuf::from(format!("{}unmatched.bam", cli.prefix));
        ensure_parent_dir(&path)?;
        Some(
            Writer::from_path(&path, &header, bam::Format::Bam)
                .context("creating unmatched BAM")?,
        )
    } else {
        None
    };

    let mut total = 0usize;
    let mut matched = 0usize;
    let mut unmatched = 0usize;

    for rec in reader.records() {
        let record = rec?;
        total += 1;

        if let Some(id) = subsetter.process_record(&record, &tag) {
            writers[id].write(&record)?;
            matched += 1;
        } else {
            unmatched += 1;
            if let Some(w) = unmatched_writer.as_mut() {
                w.write(&record)?;
            }
        }
    }

    // Timing
    if let Ok(elapsed) = now.elapsed() {
        let ms = elapsed.as_millis();
        let sec = ms / 1000;
        let min = sec / 60;
        let hr = min / 60;

        eprintln!(
            "\nProcessed {total} reads\nMatched: {matched}\nUnmatched: {unmatched}\nTime: {hr}h {}m {}s\n",
            min % 60,
            sec % 60
        );
    }

    for path in &subsetter.ofile_names {
        eprintln!("Wrote {}", path.display());
    }

    Ok(())
}
