use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::Parser;

use bam_tide::compare_report::CompareReport;

// bigtools
use bigtools::{BigWigRead, ChromInfo};
use bigtools::utils::reopen::ReopenableFile;

type BwReader = BigWigRead<ReopenableFile>;

/// Verify Rust bam-coverage against python bamCoverage BigWig
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Python-generated BigWig
    #[arg(long)]
    python_bw: PathBuf,

    /// Rust-generated BigWig
    #[arg(long)]
    rust_bw: PathBuf,

    /// Bin width (must match coverage generation)
    #[arg(long, default_value_t = 50)]
    bin_width: u32,

    /// Epsilon threshold
    #[arg(long, default_value_t = 1e-5)]
    eps: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut py = BigWigRead::open_file(&args.python_bw)
        .with_context(|| format!("open python BW {}", args.python_bw.display()))?;

    let mut rs = BigWigRead::open_file(&args.rust_bw)
        .with_context(|| format!("open rust BW {}", args.rust_bw.display()))?;

    // Read chrom lists
    let py_chroms: Vec<ChromInfo> = py.chroms().iter().map(|s| s.clone()).collect();
    let rs_chroms: Vec<ChromInfo> = rs.chroms().iter().map(|s| s.clone()).collect();

    // Build rust chrom lookup (borrowed, no clones)
    let mut rs_map: HashMap<String, u32> = HashMap::new();
    for c in &rs_chroms {
        rs_map.insert(c.name.clone(), c.length);
    }

    // Validate chrom sets + lengths
    for c in &py_chroms {
        match rs_map.get(c.name.as_str()) {
            None => bail!("chrom {} missing in rust BigWig", c.name),
            Some(&len2) if len2 != c.length => bail!(
                "chrom {} length mismatch: python={} rust={}",
                c.name,
                c.length,
                len2
            ),
            _ => {}
        }
    }

    let mut total = CompareReport::default();

    println!("CHR\tmax_abs\tmean_abs\trmse\tfrac_over\tpearson rho");

    for c in &py_chroms {
        let chr = c.name.as_str();
        let chr_len = c.length;

        let py_bins = bins_from_bigwig(&mut py, chr, chr_len, args.bin_width)?;
        let rs_bins = bins_from_bigwig(&mut rs, chr, chr_len, args.bin_width)?;

        let mut chr_rep = CompareReport::default();
        chr_rep.update_from_bins(&py_bins, &rs_bins, args.eps);

        let (frac, mean, rmse) = chr_rep.finish();

        println!(
            "{}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.4}\t{:.4}",
            chr, chr_rep.max_abs, mean, rmse, frac, chr_rep.pearson()
        );

        total.merge(&chr_rep);
    }

    let (frac, mean, rmse) = total.finish();
    println!(
        "TOTAL\t{:.3e}\t{:.3e}\t{:.3e}\t{:.4}\t{:.4}",
        total.max_abs, mean, rmse, frac, total.pearson()
    );

    Ok(())
}

/// Create a binned representation of a BigWig chromosome
fn bins_from_bigwig(
    bw: &mut BwReader,
    chr: &str,
    chr_len: u32,
    bin_w: u32,
) -> Result<Vec<f64>> {
    let nbins = ((chr_len as u64 + bin_w as u64 - 1) / bin_w as u64) as usize;
    let mut bins = vec![0.0_f64; nbins];

    let mut it = bw
        .get_interval(chr, 0, chr_len)
        .with_context(|| format!("intervals {chr}:0-{chr_len}"))?;

    while let Some(iv) = it.next() {
        let iv = iv.unwrap();
        let mut s = iv.start;
        let e = iv.end.min(chr_len);
        let v = iv.value as f64;

        if e <= s {
            continue;
        }

        while s < e {
            let bin_id = (s / bin_w) as usize;
            let bin_start = (bin_id as u32) * bin_w;
            let bin_end = (bin_start + bin_w).min(chr_len);

            let seg_end = e.min(bin_end);
            let overlap = (seg_end - s) as f64;

            bins[bin_id] += v * overlap;
            s = seg_end;
        }
    }

    // normalize to mean per base in bin
    for (i, x) in bins.iter_mut().enumerate() {
        let bin_start = (i as u32) * bin_w;
        let bin_end = (bin_start + bin_w).min(chr_len);
        let denom = (bin_end - bin_start) as f64;
        if denom > 0.0 {
            *x /= denom;
        }
    }

    Ok(bins)
}
