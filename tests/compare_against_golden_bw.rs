use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::tempdir;

// bigtools read API
use bigtools::BigWigRead;

fn exe_path(name: &str) -> PathBuf {
    let base = if cfg!(debug_assertions) {
        PathBuf::from("./target/debug")
    } else {
        PathBuf::from("./target/release")
    };
    base.join(name)
}

fn run_ok(cmd: &mut Command) -> Result<(), String> {
    let out = cmd.output().map_err(|e| format!("failed to spawn: {e:?}"))?;
    if !out.status.success() {
        return Err(format!(
            "command failed (exit={:?})\ncmd: {:?}\nstdout:\n{}\nstderr:\n{}",
            out.status.code(),
            cmd,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    Ok(())
}

/// Compute weighted mean value in [start,end) from bigWig.
/// If no intervals overlap, returns 0.0 (bin is empty).
fn bw_weighted_mean(
    bw: &mut BigWigRead<&File>,
    chr: &str,
    start: u32,
    end: u32,
) -> Result<f32, String> {
    let mut it = bw
        .get_interval(chr, start, end)
        .map_err(|e| format!("bigwig get_interval failed for {chr}:{start}-{end}: {e:?}"))?;

    let mut sum = 0.0_f64;
    let mut covered = 0u64;

    while let Some(item) = it.next() {
        let v = item.map_err(|e| format!("bigwig interval item error: {e:?}"))?;
        // overlap with requested window
        let s = start.max(v.start);
        let e = end.min(v.end);
        if e > s {
            let w = (e - s) as u64;
            covered += w;
            sum += (v.value as f64) * (w as f64);
        }
    }

    if covered == 0 {
        Ok(0.0)
    } else {
        Ok((sum / covered as f64) as f32)
    }
}

/// Pearson correlation over paired vectors (ignores all-zero case gracefully).
fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;

    let mut num = 0.0;
    let mut den_x = 0.0;
    let mut den_y = 0.0;
    for i in 0..n {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    if den_x == 0.0 || den_y == 0.0 {
        return f64::NAN;
    }
    num / (den_x.sqrt() * den_y.sqrt())
}

#[test]
fn rust_bigwig_matches_golden_deeptools() -> Result<(), String> {
    // --- config ---
    let bam = Path::new("legacy/testData/subset.bam");
    let golden = Path::new("legacy/testData/subset.deepTools.bw");
    let bin_width: u32 = 50;

    // Thresholds (tune once you have a stable baseline)
    let eps_abs: f64 = 1e-3;
    let max_frac_bad: f64 = 1e-4;
    let min_corr: f64 = 0.99;

    if !bam.exists() {
        return Err(format!("missing test BAM: {}", bam.display()));
    }
    if !golden.exists() {
        return Err(format!(
            "missing golden bigWig: {}\n\
             Create it once using deepTools and commit it.\n\
             Example:\n\
             bamCoverage -b {} -o {} --outFileFormat bigwig --binSize {} --normalizeUsing None --samFlagExclude 3328 --minMappingQuality 0",
            golden.display(),
            bam.display(),
            golden.display(),
            bin_width
        ));
    }

    // --- run your tool to produce out.bw ---
    let tmp = tempdir().map_err(|e| format!("tempdir: {e:?}"))?;
    let out_bw = tmp.path().join("rust.bw");

    let exe = exe_path("bam-coverage");
    if !exe.exists() {
        return Err(format!(
            "missing binary: {} (build first: cargo build --release)",
            exe.display()
        ));
    }

    // Adjust flags to match your CLI.
    // IMPORTANT: choose settings matching the golden bamCoverage command:
    // - same bin width
    // - same normalize meaning (ideally mean coverage == deepTools None)
    // - same filtering policy (secondary/supp/dups) if your tool supports it
    run_ok(Command::new(&exe)
        .args([
            "-b", bam.to_str().unwrap(),
            "-o", out_bw.to_str().unwrap(),
            "-w", &bin_width.to_string(),
            // make sure your output is bigWig in this binary / mode
            // and normalization corresponds to deepTools "None"
            "-n", "not",
            "--min-mapping-quality", "0",
            "--include-secondary", "false",
            "--include-supplementary", "false",
            "--include-duplicates", "false",
        ])
    )?;

    // --- open bigWigs ---
    let gf = File::open(golden).map_err(|e| format!("open golden bw: {e:?}"))?;
    let of = File::open(&out_bw).map_err(|e| format!("open rust bw: {e:?}"))?;

    let mut bw_g = BigWigRead::open(&gf).map_err(|e| format!("BigWigRead open golden: {e:?}"))?;
    let mut bw_o = BigWigRead::open(&of).map_err(|e| format!("BigWigRead open rust: {e:?}"))?;

    // Iterate chromosomes from golden (source of truth)
    let chroms = bw_g.chroms().to_vec();

    let mut x = Vec::<f64>::new();
    let mut y = Vec::<f64>::new();

    let mut bins_with_signal = 0u64;
    let mut bad_bins = 0u64;
    let mut mean_abs_sum = 0.0_f64;
    let mut max_abs = 0.0_f64;

    for c in chroms {
        let chr = c.name.as_str();
        let len = c.length;

        let mut start = 0u32;
        while start < len {
            let end = (start + bin_width).min(len);

            let vg = bw_weighted_mean(&mut bw_g, chr, start, end)? as f64;
            let vo = bw_weighted_mean(&mut bw_o, chr, start, end)? as f64;

            if vg != 0.0 || vo != 0.0 {
                bins_with_signal += 1;
                let abs = (vo - vg).abs();
                mean_abs_sum += abs;
                if abs > max_abs {
                    max_abs = abs;
                }
                if abs > eps_abs {
                    bad_bins += 1;
                }

                x.push(vg);
                y.push(vo);
            }

            start = end;
        }
    }

    let mean_abs = if bins_with_signal > 0 {
        mean_abs_sum / bins_with_signal as f64
    } else {
        0.0
    };
    let frac_bad = if bins_with_signal > 0 {
        bad_bins as f64 / bins_with_signal as f64
    } else {
        0.0
    };
    let corr = pearson_corr(&x, &y);

    eprintln!(
        "golden bigWig comparison:\n  bins_with_signal={}\n  bad_bins={}\n  frac_bad={:.3e}\n  mean_abs={:.6}\n  max_abs={:.6}\n  corr={:.6}",
        bins_with_signal, bad_bins, frac_bad, mean_abs, max_abs, corr
    );

    if bins_with_signal == 0 {
        return Err("comparison found zero signal bins; test BAM/golden likely wrong".into());
    }
    if frac_bad > max_frac_bad {
        return Err(format!("frac_bad too high: {frac_bad:.3e} > {max_frac_bad:.3e}"));
    }
    if corr.is_nan() || corr < min_corr {
        return Err(format!("corr too low: {corr} < {min_corr}"));
    }

    Ok(())
}

