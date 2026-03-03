
// Import necessary modules
use std::fs::{self};
use std::path::Path;
use std::collections::BTreeMap;
use assert_cmd::cargo::cargo_bin_cmd;


/// Parse:
///   "Selected 37/8292 reads"
/// and subsequent lines like:
///   "11  testData/output/two_clusters_barcodes.bam"
///   "26  testData/output/two_clusters_barcodes2.bam"
fn parse_cli_report(stdout: &str) -> (usize, usize, BTreeMap<String, usize>) {
    let mut selected: Option<usize> = None;
    let mut total: Option<usize> = None;
    let mut per_file: BTreeMap<String, usize> = BTreeMap::new();

    for raw in stdout.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }

        // "Selected 37/8292 reads in 0h 0min 0s 8ms"
        if let Some(rest) = line.strip_prefix("Selected ") {
            // find the first token that looks like "X/Y"
            let xy_token = rest
                .split_whitespace()
                .find(|tok| tok.contains('/') )
                .unwrap_or_else(|| panic!("Failed to find X/Y token in Selected line: {line:?}"));

            let (a, b) = xy_token
                .split_once('/')
                .unwrap_or_else(|| panic!("Failed to split X/Y token in Selected line: {line:?}"));

            let a = a.trim().parse::<usize>()
                .unwrap_or_else(|_| panic!("Failed to parse selected from: {line:?}"));
            let b = b.trim().parse::<usize>()
                .unwrap_or_else(|_| panic!("Failed to parse total from: {line:?}"));

            selected = Some(a);
            total = Some(b);
            continue;
        }

        // "<count> <path>"
        // (your output has variable whitespace indentation, so split_whitespace is ideal)
        let mut it = line.split_whitespace();
        let first = it.next();
        let second = it.next();

        if let (Some(cnt_str), Some(path_str)) = (first, second) {
            if let Ok(cnt) = cnt_str.parse::<usize>() {
                // If there are additional whitespace-separated pieces, they belong to the path only if
                // you ever print paths containing spaces. If not, this is enough.
                // If you DO want to support spaces in paths, see note below.
                per_file.insert(path_str.to_string(), cnt);
            }
        }
    }

    let selected = selected.expect("Did not find a 'Selected X/Y reads' line in stdout");
    let total = total.expect("Did not find a 'Selected X/Y reads' line in stdout");

    (selected, total, per_file)
}

fn assert_cli_numerics(
    stdout: &str,
    expected_selected: usize,
    expected_total: usize,
    expected_files: &[(&str, usize)],
) {
    let (selected, total, per_file) = parse_cli_report(stdout);

    assert_eq!(selected, expected_selected, "Selected count mismatch.\nstdout:\n{stdout}");
    assert_eq!(total, expected_total, "Total count mismatch.\nstdout:\n{stdout}");

    // Check file counts (both presence + value)
    for (path, exp_cnt) in expected_files {
        let got = per_file.get(*path).copied();
        assert_eq!(
            got,
            Some(*exp_cnt),
            "Per-file count mismatch for {path:?}.\nParsed: {per_file:#?}\nstdout:\n{stdout}"
        );
    }

    // Optional: ensure sum of per-file counts equals selected
    let sum: usize = per_file.values().sum();
    assert_eq!(
        sum, selected,
        "Sum(per-file counts) != Selected.\nParsed: {per_file:#?}\nstdout:\n{stdout}"
    );
}

#[test]
fn test_multi_subset_bam() {

    let path = Path::new("testData/output");
    if path.exists() {
        fs::remove_dir_all(path).expect("Failed to remove directory");
    }


    let args = &[
        "-b", "testData/bam_subset_test.bam",
        "-v", "testData/barcodes.txt","testData/barcodes2.txt",
        "-o", "testData/output/two_clusters_",
    ];

    // Execute the command with the provided arguments
    let mut cmd = cargo_bin_cmd!( "bam-subset-tag" );
    cmd.args( args );
    let exe = cmd.get_program().to_string_lossy();

    // build printable command
    let cmdline = format!(
        "{} {}",
        exe,
        args.iter().map(|a| format!("{:?}", a)).collect::<Vec<_>>().join(" ")
    );

    let output = cmd.output()
        .map_err(|e| {
            eprintln!("Failed to execute command\n{cmdline}: {}", e);
            e
        }).unwrap();


    // Check if the command was successful (exit code 0)
    if !output.status.success() {
        panic!(
            "Command failed!\n{}\n\nstdout:\n{}\n\nstderr:\n{}",
            cmdline,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }

    let expect = vec!["two_clusters_barcodes.bam","two_clusters_barcodes2.bam"];

    for file in expect{
        let ofile= "testData/output/".to_string() +file;
        assert!(Path::new( &ofile ).exists(), "Expected outfile {} does not exist!", ofile);
    }

    // Now verify numerics:
    assert_cli_numerics(
        &String::from_utf8_lossy(&output.stdout),
        37,     // expected selected
        8292,   // expected total
        &[
            ("testData/output/two_clusters_barcodes.bam", 11),
            ("testData/output/two_clusters_barcodes2.bam", 26),
        ],
    );

}
