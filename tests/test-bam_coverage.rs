//test-bam_coverage.rs
use std::process::Command;

fn parse_bedgraph_lines(s: &str) -> Vec<(String, u32, u32, i64)> {
    s.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with("track"))
        .map(|line| {
            let mut it = line.split('\t');
            let chr = it.next().unwrap().to_string();
            let start: u32 = it.next().unwrap().parse().unwrap();
            let end: u32 = it.next().unwrap().parse().unwrap();

            // your output appears integer-like (e.g. 54), so parse as i64.
            // If you later emit floats, switch to f32 and epsilon compare.
            let val: i64 = it.next().unwrap().parse().unwrap();

            (chr, start, end, val)
        })
        .collect()
}

#[test]
fn test_bam_coverage_bedgraph_runs_and_outputs_lines() {
    // Cargo provides the built binary path.
    let exe = env!("CARGO_BIN_EXE_bam-coverage");

    // Use a temp output file
    //let out = tempfile::NamedTempFile::new().expect("tempfile");
    //let out_path = out.path().to_path_buf();
    let out_path = "legacy/testData/_out_bam_coverage.bedgraph";

    let output = Command::new(exe)
        .args([
            "-b", "legacy/testData/bowtie_chrm_mutated.bam",
            "-o", out_path,
            "-w", "10",
            "--min-mapping-quality", "0",
        ])
        .output()
        .expect("Failed to run bam-coverage");

    if !output.status.success() {
        panic!(
            "bam-coverage failed.\nstatus: {:?}\nstdout:\n{}\nstderr:\n{}",
            output.status.code(),
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }

    // Read output file and do basic sanity checks
    let s = std::fs::read_to_string(&out_path).expect("read output");
    let rows = parse_bedgraph_lines(&s);
    // Compare a small stable prefix of non-zero bins
    // 2) Exact first 10 lines (your `head`)
    let expected_head: Vec<(String, u32, u32, i64)> = vec![
        ("chr1".to_string(), 0, 248_956_422, 0),
        ("chr2".to_string(), 0, 242_193_529, 0),
        ("chr3".to_string(), 0, 198_295_559, 0),
        ("chr4".to_string(), 0, 190_214_555, 0),
        ("chr5".to_string(), 0, 181_538_259, 0),
        ("chr6".to_string(), 0, 170_805_979, 0),
        ("chr7".to_string(), 0, 159_345_973, 0),
        ("chr8".to_string(), 0, 145_138_636, 0),
        ("chr9".to_string(), 0, 138_394_717, 0),
        ("chr10".to_string(), 0, 133_797_422, 0),
    ];

    for (i, exp) in expected_head.iter().enumerate() {
        assert_eq!(&rows[i], exp, "Mismatch in head line {i}");
    }

    assert!(
        rows.len() >= 1335,
        "Too few rows: got {}, expected at least {}",
        rows.len(),
        1493
    );

    // 3) First 10 non-zero lines (your grep -v "\s0$" | head)
    let nonzero: Vec<(String, u32, u32, i64)> = rows
        .iter()
        .cloned()
        .filter(|(_, _, _, v)| *v != 0)
        .collect();

    let expected_chrM_prefix: Vec<(String, u32, u32, i64)> = vec![
        ("chrM".to_string(), 0, 10, 7),
        ("chrM".to_string(), 10, 20, 35),
        ("chrM".to_string(), 20, 30, 77),
        ("chrM".to_string(), 30, 40, 107),
        ("chrM".to_string(), 40, 50, 156),
        ("chrM".to_string(), 50, 60, 192),
        ("chrM".to_string(), 60, 70, 224),
        ("chrM".to_string(), 70, 80, 243),
        ("chrM".to_string(), 80, 90, 217),
        ("chrM".to_string(), 90, 100, 178),
    ];

    assert!(
        nonzero.len() >= 1122,
        "Too few non-zero rows: got {}, expected at least {}",
        nonzero.len(),
        1280
    );

    for (i, exp) in expected_chrM_prefix.iter().enumerate() {
        assert_eq!(&nonzero[i], exp, "Mismatch in first non-zero line {i}");
    }

}