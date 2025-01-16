use bam_tide;
use std::fs;
use std::collections::HashMap;


#[test]
fn test_calculate_total_bins(){
    #[test]
    fn test_calculate_total_bins() {
        // Create a HashMap representing chromosome sizes
        let mut genome_info = HashMap::new();
        genome_info.insert("chr1".to_string(), 100, 0  );
        genome_info.insert("chr2".to_string(), 1444, (100 + 50 - 1) /50 );
        genome_info.insert("chr3".to_string(), 3456, (1444 + 50 - 1) / 50 );
        
        let bin_width = 50;
        
        // Test case: Calculate total bins
        let total_bins = bam_tide::calculate_total_bins(&genome_info, bin_width);
        
        // Expected number of bins:
        let expected_bins = (100 + 50 - 1) / 50 + (1444 + 50 - 1) / 50 + (3456 + 50 - 1) / 50;
        assert_eq!(total_bins, expected_bins, "The total bins calculated do not match the expected value.");
    }
}

#[test]
fn test_bam_to_bigwig() {
    let bam_file = "test_data/test.bam"; // Replace with a small BAM test file
    let bigwig_file = "test_output/test.bigwig";

    // Ensure output folder exists
    fs::create_dir_all("test_output").unwrap();

    // Run the function
    let result = bam_tide::bam_to_bedgraph(bam_file, bigwig_file, 50);

    // Assert success
    assert!(result.is_ok());

    // Optionally, validate the BigWig file
    assert!(fs::metadata(bigwig_file).is_ok());
}