use anyhow::{Context, Result};
use bam_tide::illumina_normalizer::{cli::Cli, IlluminaNormalizer};

fn main() -> Result<()> {
    let cli = Cli::parse_args();

    let mut normalizer = IlluminaNormalizer::from_cli(cli)
        .context("failed to initialize Illumina normalizer")?;

    normalizer
        .run()
        .context("Illumina normalizer failed while processing input FASTQ files")?;

    eprintln!("{}", normalizer.stats_report());

    let log_path = normalizer.config().out.with_extension("log.tsv");
    let log_path = log_path
        .to_str()
        .context("stats log path is not valid UTF-8")?;

    normalizer.stats().report_to_csv(log_path);

    Ok(())
}
