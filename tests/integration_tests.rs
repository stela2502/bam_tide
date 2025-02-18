use std::fs;
use std::collections::HashMap;
use std::path::Path;

use bam_tide::bed_data::BedData;
use bam_tide::gtf_logics::{AnalysisType};

#[test]
fn test_bam_to_bigwig() {
    let bam_file = "testData/test.bam"; // Replace with a small BAM test file
    let bigwig_file = "testData/output/test.wig";

    // Ensure output folder exists
    fs::create_dir_all("testData/output/").unwrap();

    let path = Path::new(bigwig_file);
    if path.exists() {
        if let Err(e) = fs::remove_file(path) {
            eprintln!("Failed to remove file {}: {}", bigwig_file, e);
        }
    }

    // Run the function
    let data =  BedData::new( bam_file, 50, 1, &AnalysisType::Bulk, &b"CR", &b"No" );
    let result = BedData::write_bedgraph( &data, bigwig_file );

    // Assert success
    assert!(result.is_ok());

    // Optionally, validate the BigWig file
    assert!(fs::metadata(bigwig_file).is_ok());
}