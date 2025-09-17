
use std::process::Command;
use std::path::Path;
use std::fs;
#[allow(unused_imports)]
use std::process::exit;

#[test]
fn test_bam2bigwig_creates_output() {
    let input_bam = "testData/test.bam";
    let output_bw = "testData/output/test.bw";

    // Remove output file if it exists (cleanup from prior runs)
    if Path::new(output_bw).exists() {
        fs::remove_file(output_bw).expect("Failed to remove old test output file");
    }

    let is_release_mode = !cfg!(debug_assertions);

    if ! is_release_mode {
        eprintln!("Test should be re-run in release mode (speed!)");
        exit(0);
    }

    let command = if is_release_mode {
        "./target/release/bam2bedgraph"
    } else {
        "./target/debug/bam2bedgraph"
    };

    // Run the compiled binary (assumes it's in target/debug or target/release)
    let status = Command::new( command )
        .args([
            "-b", input_bam,
            "-o", output_bw,
            "--umi-tag", "No"
        ])
        .status()
        .expect("Failed to execute bam2bigwig binary");

    assert!(status.success(), "bam2bigwig did not exit successfully");

    // Check if the output file was created
    assert!(
        Path::new(output_bw).exists(),
        "Output file was not created: {}",
        output_bw
    );

    // Clean up the output file after test
    fs::remove_file(output_bw).expect("Failed to remove test output file after check");
}
