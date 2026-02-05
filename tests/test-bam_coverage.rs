use std::process::Command;
use std::path::{Path,PathBuf};

#[test]
fn test_bam_coverage_matches_deeptools_bigwig() {
    let bam = "legacy/testData/subset.bam";

    let out_dir = PathBuf::from("legacy/testData");
    let py_bw = out_dir.join("_ref_deeptools.bw");
    let rs_bw = out_dir.join("_test_bam_coverage.bw");

    // ------------------------------------------------------------------
    // 1) Create reference BigWig with deeptools bamCoverage
    // ------------------------------------------------------------------
    if !Path::new(&py_bw).exists() {
        let status = Command::new("bamCoverage")
            .args([
                "-b", bam,
                "-o", py_bw.to_str().unwrap(),
                "--binSize", "50",
                "--normalizeUsing", "None",
                "--minMappingQuality", "0",
                //"--ignoreDuplicates",
                "--extendReads", "0",
                "--numberOfProcessors", "4",
            ])
            .status()
            .expect("failed to run deeptools bamCoverage");

        assert!(status.success(), "deeptools bamCoverage failed");
    }

    // ------------------------------------------------------------------
    // 2) Create BigWig with our bam-coverage
    // ------------------------------------------------------------------
    let exe = env!("CARGO_BIN_EXE_bam-coverage");

    let status = Command::new(exe)
        .args([
            "-b", bam,
            "-o", rs_bw.to_str().unwrap(),
            "-w", "50",
            "--min-mapping-quality", "0",
        ])
        .status()
        .expect("failed to run bam-coverage");

    assert!(status.success(), "bam-coverage failed");

    // ------------------------------------------------------------------
    // 3) Compare using bw-compare
    // ------------------------------------------------------------------
    let cmp = env!("CARGO_BIN_EXE_bw-compare");

    let output = Command::new(cmp)
        .args([
            "--python-bw", py_bw.to_str().unwrap(),
            "--rust-bw", rs_bw.to_str().unwrap(),
            "--bin-width", "50",
            "--eps", "1e-4",
        ])
        .output()
        .expect("failed to run bw-compare");

    if !output.status.success() {
        panic!(
            "bw-compare failed\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }
    let stderr = String::from_utf8_lossy(&output.stderr);

    parse_total_line( &stderr );

}

fn parse_total_line(stdout: &str) -> (f64, f64, f64, f64) {
    for line in stdout.lines() {
        if line.starts_with("TOTAL") {
            let parts: Vec<&str> = line.split('\t').collect();
            assert!(
                parts.len() >= 5,
                "Malformed TOTAL line: {line}"
            );

            let max_abs: f64 = parts[1].parse().unwrap();
            let mean_abs: f64 = parts[2].parse().unwrap();
            let rmse: f64 = parts[3].parse().unwrap();
            let frac_over: f64 = parts[4].parse().unwrap();

            return (max_abs, mean_abs, rmse, frac_over);
        }
    }
    panic!("No TOTAL line found in bw-compare output:\n{stdout}");
}