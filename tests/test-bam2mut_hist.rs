use std::process::Command;

#[test]
fn test_bam2mut_hist_output() {
    let output = Command::new("target/release/bam2mut_hist")
        .args(&["--input", "testData/bowtie_chrm_mutated.bam", "-b", "10", "-h"])
        .output()
        .expect("Failed to run bam2mut_hist");

    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);

    let expected = vec![1904, 1416, 1422, 1504, 1402, 1452, 1399, 1395, 1671, 1603];

    for (i, &val) in expected.iter().enumerate() {
        let expected_line = format!("{:2}: {}", i, val);
        assert!(
            stdout.contains(&expected_line),
            "Missing line in output: '{}'",
            expected_line
        );
    }
}