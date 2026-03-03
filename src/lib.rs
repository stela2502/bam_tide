pub mod cli;
pub mod core;
pub mod bed_data;
pub mod data_iter;
pub mod compare_report;
pub mod subset_bam;


pub fn compute_io_threads(user_threads: usize) -> usize {
    user_threads
    /*match user_threads {
        0 | 1 => 1,
        2 => 1,
        _ => user_threads - 1,
    }*/
}