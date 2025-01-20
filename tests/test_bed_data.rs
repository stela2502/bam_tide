#[cfg(test)]
mod tests {
    use super::*;

    use bam_tide::bed_data::{BedData, Value}; // Adjust the path if needed
    use bam_tide::data_iter::DataIter;


    #[test]
    fn run_data_iter_test() {
        // Prepare mock data for testing
        let genome_info = vec![
            ("chr1".to_string(), 990, 0),
            ("chr2".to_string(), 2000, 5),
        ];

        let search = BedData::genome_info_to_search( &genome_info );

        let coverage_data = vec![10, 20, 30, 40, 50, 1,2,3,4,5,6,7,8,9,10 ]; // Example coverage data
        let bin_width = 200;

        // Initialize the mock BedData struct
        let bed_data = BedData{
            genome_info,
            search,
            coverage_data,
            bin_width,
            threads:1,
        };

        // Create an instance of DataIter
        let mut data_iter = DataIter::new(&bed_data);

        // Run the test and collect results
        let mut results = Vec::new();
        while let Some((chrom, value)) = data_iter.next() {
            results.push((chrom, value));
        }

        // Print results for inspection
        println!("Results: {:?}", results);

        // Add assertions to check expected results
        assert_eq!(results.len(), 15); // We expect 15 results
        assert_eq!(results[0].0, "chr1");
        assert_eq!(&format!("{:?}",results[0].1), "Value { start: 0, end: 200, value: 10.0 }" );

        assert_eq!(results[4].0, "chr1");
        assert_eq!(&format!("{:?}",results[4].1), "Value { start: 800, end: 990, value: 50.0 }" );

        assert_eq!(results[5].0, "chr2");
        assert_eq!(&format!("{:?}",results[5].1), "Value { start: 0, end: 200, value: 1.0 }" );

    }
}