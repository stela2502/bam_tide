// Import necessary modules
use std::process::Command;
use std::fs::File;
use std::fs;
use std::io::BufReader;
use std::io::BufRead;
use std::collections::HashMap;
use std::process::exit;
use std::path::PathBuf;

// bigwig should work as well - but it creates enormouse files as it adds the zeros in, too.
#[test]
fn test_bam2bedgraph() {

    let is_release_mode = !cfg!(debug_assertions);

    if ! is_release_mode {
        eprintln!("Test should be re-run in release mode (speed!)");
        exit(0);
    }

    let command = if is_release_mode {
        "./target/release/bowtie2mtx"
    } else {
        "./target/debug/bowtie2mtx"
    };

    let args = &[
        "-b", "testData/mutation_test.bam",
        "-o", "testData/output/mutations",
        "-f", "testData/mutation_test.fa.gz", 
        "--min-umi", "0",
    ]; 

    // check if the cell names make sense
    let file_path = PathBuf::from("testData/output/mutations");
    if file_path.exists() {
        fs::remove_dir_all(&file_path).expect("Failed to remove existing directory");
    }

    // Execute the command with the provided arguments
    let output = Command::new( command ).args( args )
        .output()
        .map_err(|e| {
            eprintln!("Failed to execute command: {}", e);
            e
        }).unwrap();

    let cmd = format!("{} {}", command, args.join(" "));
    if !output.status.success() {
        eprintln!("Command failed: {}", cmd);
        // Handle failure accordingly
    }else {
        println!("{}", cmd );
    }

    // Check if the command was successful (exit code 0)
    assert!(output.status.success(), "{:?}", String::from_utf8_lossy(&output.stderr));
    // Convert output to string
    let output_str = String::from_utf8_lossy(&output.stdout);

    let file_path = PathBuf::from("testData/output/mutations/bowtie2mtx/mutations");

    let expected_files = vec![
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
    ];

    for file in expected_files {
        let file = file_path.join(file);
        assert!(file.exists(), "Expected file {:?} was not created", file);
    }

}