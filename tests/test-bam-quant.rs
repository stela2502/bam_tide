use assert_cmd::Command;
use std::path::{Path, PathBuf};
use std::process::Command as StdCommand;

const PBMC_URL: &str = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam";

// GENCODE Human Release 49, “Basic gene annotation CHR” GTF on EMBL-EBI FTP. :contentReference[oaicite:1]{index=1}
const GENCODE_BASIC_CHR_GTF_GZ: &str = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz";

fn run(cmd: &mut StdCommand) {
    let status = cmd.status().expect("failed to run command");
    assert!(status.success(), "command failed: {:?}", cmd);
}

fn download_if_needed(url: &str, path: &Path) {
    if path.exists() {
        eprintln!("✓ cached: {}", path.display());
        return;
    }
    eprintln!("⬇ downloading: {}", url);
    run(StdCommand::new("wget").args(["-O", path.to_str().unwrap(), url]));
}

fn require_exe(name: &str) {
    let exists = std::process::Command::new("which")
        .arg(name)
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    assert!(exists, "required executable not found in PATH: {}", name);
}

fn require_tools() {
    require_exe("gtf_splice_index");
    require_exe("samtools");
    require_exe("wget");
    require_exe("bash");
    require_exe("gzip");
    require_exe("awk");
}

/// Create chr1-only GTF (keeps all header lines, filters feature lines by seqname == "chr1").
/// Works for typical GENCODE GTFs where seqname is column 1.
fn subset_gtf_chr1_if_needed(gtf_gz: &Path, out_gtf: &Path) {
    if out_gtf.exists() {
        eprintln!("✓ cached: {}", out_gtf.display());
        return;
    }

    eprintln!("⚙ subsetting chr1 GTF: {}", out_gtf.display());

    // gzip -dc <in> | awk '(^#) || ($1=="chr1")' > out
    let script = format!(
        r#"
set -euo pipefail
gzip -dc "{in_gz}" \
  | awk 'BEGIN{{FS="\t"}} /^#/ {{print; next}} ($1=="chr1") {{print}}' \
  > "{out_gtf}"
"#,
        in_gz = gtf_gz.display(),
        out_gtf = out_gtf.display(),
    );

    run(StdCommand::new("bash").arg("-c").arg(script));
}

/// Build your splice index if missing.
/// Assumes you have a workspace binary named `gtf_splice_index` with a `build` subcommand:
///   gtf_splice_index build -a <annotation.gtf> -i <index.dat>
fn build_index_if_needed(gtf_chr1: &Path, index_path: &Path) {
    if index_path.exists() {
        eprintln!("✓ cached: {}", index_path.display());
        return;
    }

    eprintln!("⚙ building index: {}", index_path.display());

    let mut cmd = assert_cmd::Command::new("gtf_splice_index");

    cmd.args([
        "build",
        "-a",
        gtf_chr1.to_str().unwrap(),
        "-i",
        index_path.to_str().unwrap(),
    ]);

    cmd.assert().success();
}

fn make_subset_bam_if_needed(input_bam: &Path, output_bam: &Path, n: u64) {
    let bai = output_bam.with_extension("bam.bai");

    if output_bam.exists() && bai.exists() {
        eprintln!("✓ cached: {}", output_bam.display());
        if !bai.exists() {
            let script = format!(
                r#"
samtools index {bam_out}
"#,
                bam_out = output_bam.display()
            );
            run(StdCommand::new("bash").arg("-lc").arg(script));
            eprintln!(" create missing bam index")
        }
        return;
    }

    eprintln!(
        "⚙ creating subset BAM (first {n} alignments): {}",
        output_bam.display()
    );

    // samtools view -h | awk keep headers + first N alignments | samtools view -b | samtools index
    let script = format!(
        r#"
set -euo pipefail

samtools view -h "{in_bam}" \
  | awk -v n={n} 'BEGIN{{c=0}} /^@/ {{print; next}} {{c++; if (c<=n) print; else exit}}' \
  | samtools view -b -o "{out_bam}" -
samtools index "{out_bam}"
"#,
        in_bam = input_bam.display(),
        out_bam = output_bam.display(),
        n = n
    );

    run(StdCommand::new("bash").arg("-lc").arg(script));
}

fn output_nonempty(outdir: &Path) -> bool {
    if !outdir.exists() {
        return false;
    }
    std::fs::read_dir(outdir)
        .ok()
        .and_then(|mut it| it.next())
        .is_some()
}

#[test]
#[ignore = "Downloads data + requires wget/gzip/awk/samtools. Run manually."]
fn bam_quant_full_pipeline_cached_with_auto_index() {
    // ---- knobs ----
    let n_reads: u64 = std::env::var("BAM_QUANT_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2_000_000);

    require_tools(); // check for all external tools

    // ---- persistent cache ----
    let cache_dir = PathBuf::from("/tmp/bam_quant_cache");
    std::fs::create_dir_all(&cache_dir).unwrap();

    // Inputs
    let pbmc_bam = cache_dir.join("pbmc_1k.bam");
    let gtf_gz = cache_dir.join("gencode.basic.chr.gtf.gz");
    let gtf_chr1 = cache_dir.join("gencode.basic.chr1.gtf");

    // Derived artifacts
    let subset_bam = cache_dir.join(format!("pbmc_subset_{n_reads}.bam"));
    let index_path = cache_dir.join("splice_index_chr1.dat");

    // Output
    let outdir = cache_dir.join(format!("out_{n_reads}_chr1"));
    std::fs::create_dir_all(&outdir).unwrap();

    // ---- step 1: download inputs ----
    download_if_needed(PBMC_URL, &pbmc_bam);
    download_if_needed(GENCODE_BASIC_CHR_GTF_GZ, &gtf_gz);

    // ---- step 2: subset GTF to chr1 ----
    subset_gtf_chr1_if_needed(&gtf_gz, &gtf_chr1);

    // ---- step 3: build index (via gtf_splice_index) ----
    build_index_if_needed(&gtf_chr1, &index_path);

    // ---- step 4: subset BAM + index it ----
    make_subset_bam_if_needed(&pbmc_bam, &subset_bam, n_reads);

    // ---- step 5: run bam-quant only if outputs missing ----
    if output_nonempty(&outdir) {
        eprintln!(
            "✓ outputs already exist, skipping bam-quant: {}",
            outdir.display()
        );
    } else {
        eprintln!("🚀 running bam-quant…");

        let mut cmd = Command::cargo_bin("bam-quant").unwrap();
        cmd.args([
            "--bam",
            subset_bam.to_str().unwrap(),
            "--index",
            index_path.to_str().unwrap(),
            "--outpath",
            outdir.to_str().unwrap(),
        ]);

        cmd.assert().success();
    }

    // ---- sanity ----
    assert!(output_nonempty(&outdir), "bam-quant produced no output");
}
