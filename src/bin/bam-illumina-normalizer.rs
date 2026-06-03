use anyhow::Result;
use bam_tide::illumina_normalizer::{cli::Cli, IlluminaNormalizer};

fn main() -> Result<()> {
    let cli = Cli::parse_args();
    let mut normalizer = IlluminaNormalizer::from_cli(cli)?;
    normalizer.run()?;

    eprintln!("{}", normalizer.stats_report());

    let log_path = normalizer.config().out.with_extension("log.tsv");
    normalizer.stats().report_to_csv(log_path.to_str().unwrap());

    Ok(())
}
