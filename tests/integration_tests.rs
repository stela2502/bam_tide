use std::fs;
use std::collections::HashMap;

use bam_tide::bed_data::BedData;

#[test]
fn test_bam_to_bigwig() {
    let bam_file = "tests/test.bam"; // Replace with a small BAM test file
    let bigwig_file = "test_output/test.bigwig";

    // Ensure output folder exists
    fs::create_dir_all("test_output").unwrap();

    // Run the function
    let data =  BedData::new( bam_file, 50, 1 );
    let result = BedData::write_bedgraph( &data, bigwig_file );

    // Assert success
    assert!(result.is_ok());

    // Optionally, validate the BigWig file
    assert!(fs::metadata(bigwig_file).is_ok());
}