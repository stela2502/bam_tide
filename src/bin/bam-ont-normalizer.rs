use anyhow::Result;
use bam_tide::ont_normalizer::{OntNormalizer, cli::Cli};

fn main() -> Result<()> {
    let cli = Cli::parse_args();

    let mut normalizer = OntNormalizer::from_cli(cli);
    normalizer.run()?;

    eprintln!("{}", normalizer.stats());

    Ok(())
}
