//bam-subset-tag.rs
// bam_tide/src/bin/bam-subset-tag.rs

use clap::Parser;
use rust_htslib::{bam, bam::Read};

use std::fs;
use std::path::Path;
use std::time::SystemTime;

use bam_tide::subset_bam::Subsetter;
use bam_tide::compute_io_threads;

/// Split one BAM into multiple BAMs by a 2-char BAM tag (e.g. CR/CB/BC).
/// Each values-file corresponds to one output BAM and contains one tag value per line.
#[derive(Parser, Debug)]
#[command(disable_version_flag = true)]
#[command(version = "0.2.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// Input BAM
    #[arg(short, long)]
    bam: String,

    /// BAM tag to look for (must be exactly 2 chars)
    #[arg(default_value = "CR", short, long)]
    tag: String,

    /// threads for bam read process (default 4)
    #[arg(default_value_t = 4, short, long)]
    threads: usize,

    /// Text files: each file is one group, one tag value per line
    #[arg(short, long, num_args(1..), value_delimiter = ' ')]
    values: Vec<String>,

    /// Output prefix (can include directory). Actual names become: `<ofile><values_file_stem>.bam`
    #[arg(short, long)]
    ofile: String,

    /// Strip this suffix from tag values before matching (e.g. "-1")
    #[arg(long)]
    strip_suffix: Option<String>,
}

fn main() -> anyhow::Result<()> {

    let now = SystemTime::now();
    let opts: Opts = Opts::parse();

    if opts.tag.len() != 2 {
        anyhow::bail!("The tag needs to be exactly two chars long - not {}", opts.tag);
    }
    let tag: [u8; 2] = opts.tag.as_bytes().try_into().unwrap();

    // Reader + header template
    let mut reader = bam::Reader::from_path(&opts.bam)?;
    let hts_threads = compute_io_threads(opts.threads);
    if let Err(e) = reader.set_threads(hts_threads) {
        eprintln!(
            "Warning: failed to enable HTSlib threading ({}). Continuing single-threaded.",
            e
        )
    };
    let header = bam::Header::from_template(reader.header());

    // Build subsetter
    let mut subsetter = Subsetter::new();
    if let Some(s) = &opts.strip_suffix {
        subsetter = subsetter.with_strip_suffix(s.clone());
    }

    for fname in &opts.values {
        subsetter.read_simple_list(fname, &opts.ofile)?;
    }

    // Ensure output dir exists
    let outpath = Path::new(&opts.ofile)
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Could not get parent path of outfile prefix {}", opts.ofile))?;
    if fs::metadata(outpath).is_err() {
        fs::create_dir_all(outpath)?;
    }

    // Create writers
    let mut ofiles = Vec::with_capacity(subsetter.num_groups());
    for fname in &subsetter.ofile_names {
        ofiles.push(bam::Writer::from_path(fname, &header, bam::Format::Bam)?);
    }

    let mut reads: u64 = 0;
    let mut lines: u64 = 0;

    let mut per_group_selected: Vec<u64> = vec![0; subsetter.num_groups()];

    for rec in reader.records() {
        let record = rec?;
        lines += 1;

        if let Some(ofile_id) = subsetter.process_record(&record, &tag) {
            reads += 1;
            per_group_selected[ofile_id] += 1;
            ofiles[ofile_id].write(&record)?;
        }
    }
    
    if let Ok(elapsed) = now.elapsed() {
        let mut milli = elapsed.as_millis();
        let mil = milli % 1000;
        milli = (milli - mil) / 1000;

        let sec = milli % 60;
        milli = (milli - sec) / 60;

        let min = milli % 60;
        milli = (milli - min) / 60;

        println!(
            "\nSelected {reads}/{lines} reads in {milli}h {min}min {sec}s {mil}ms\n"
        );
    }

    for (i, n) in per_group_selected.iter().enumerate() {
        println!("{:>12}  {}", n, subsetter.ofile_names[i]);
    }

    Ok(())
}
