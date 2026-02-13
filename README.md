[![Rust](https://github.com/stela2502/bam_tide/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/bam_tide/actions/workflows/rust.yml)

# bam_tide

Fast, reproducible BAM → binned coverage exporters and validation utilities written in Rust.

This repository currently ships two main command-line tools:

- **`bam-coverage`** — compute binned coverage from a position-sorted BAM and write **bedGraph** or **BigWig** (depending on the build/features of your binary).
- **`bw-compare`** — validate Rust-generated BigWigs against a reference (typically **deeptools `bamCoverage`**) and report per-chromosome and global statistics, including a final **`TOTAL`** line.

> Goal: provide **bamCoverage-compatible filtering semantics** (include/exclude SAM flags, MAPQ, duplicates, secondary/supplementary) with **much higher performance** and simple, scriptable reproducibility.

---

## Installation

### Option 1: Build from source (recommended)

Prerequisites:
- Rust toolchain (stable) via `rustup`
- A C toolchain (for some dependencies) as needed on your system

```bash
git clone https://github.com/stela2502/bam_tide.git
cd bam_tide
cargo build --release
```

Binaries will be available here:

```bash
./target/release/bam-coverage
./target/release/bw-compare
```

### Option 2: Use inside a container (optional)

If you run on HPC and prefer containers, you can build an Apptainer/Singularity image around a release build of these binaries. (A recipe is not included yet—see **Future developments**.)

---

## Quickstart

### 1) Generate coverage with `bam-coverage`

**Minimal example:**

```bash
./target/release/bam-coverage \
  --bam input.sorted.bam \
  --outfile sample.bw
```

**Common options:**
- `--width 50` to change bin size (default: 50)
- `--min-mapping-quality 30` to require MAPQ ≥ 30
- `--sam-flag-exclude 256` to exclude secondary alignments (deeptools-compatible)
- `--sam-flag-include 64` to include only Read1

Example (exclude secondary + supplementary alignments):

```bash
./target/release/bam-coverage \
  --bam input.sorted.bam \
  --outfile sample.bw \
  --sam-flag-exclude 2304 \
  --width 50
```

> `--sam-flag-exclude` is a **bitmask**: any read with `(flag & mask) != 0` is discarded.

---

### 2) Validate results with `bw-compare`

Compare a Python (deeptools) BigWig to a Rust BigWig:

```bash
./target/release/bw-compare \
  --python-bw python_sample.bw \
  --rust-bw rust_sample.bw
```

This prints a per-chromosome report and finishes with a **`TOTAL`** line that summarizes all chromosomes.

**Write the report to a file:**

```bash
./target/release/bw-compare \
  --python-bw python_sample.bw \
  --rust-bw rust_sample.bw \
  --outfile compare_report.txt
```

See the full benchmark results here:

[Benchmark results](Benchmark.md)

---

## Command reference

### `bam-coverage --help`

Shared CLI options for coverage exporters (bedGraph / bigWig)

```
Usage: bam-coverage [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>
          Input BAM file (sorted by chromosome position)

  -o, --outfile <OUTFILE>
          Output file (bedGraph or BigWig depending on binary)

  -n, --normalize <NORMALIZE>
          Normalize the data somehow

          [default: not]
          [possible values: not, rpkm, cpm, bpm, rpgc]

  -w, --width <WIDTH>
          Bin width for coverage calculation

          [default: 50]

      --only-r1
          Collect only R1 areas

      --min-mapping-quality <MIN_MAPPING_QUALITY>
          Minimum mapping quality to include a read

          [default: 0]

      --include-secondary
          Include secondary alignments

      --include-supplementary
          Include supplementary alignments

      --include-duplicates
          Include duplicate-marked reads

      --sam-flag-exclude <SAM_FLAG_EXCLUDE>
          Exclude reads with ANY of these SAM flag bits set (equivalent to deeptools --samFlagExclude).

          The value is a bitmask of SAM flags. Any read with (read_flag & mask) != 0 will be discarded.

          Examples:
            256  -> exclude secondary alignments
            512  -> exclude QC-failed reads
            1024 -> exclude PCR/optical duplicates
            2048 -> exclude supplementary alignments
            2816 -> exclude secondary + QC-fail + supplementary

          Default: None (no flag-based exclusion, matches bamCoverage defaults)

      --sam-flag-include <SAM_FLAG_INCLUDE>
          Include only reads that have ALL of these SAM flag bits set. Applied after the exclusion test (equivalent to deeptools --samFlagInclude).

          The value is a bitmask of SAM flags. A read is kept only if: (read_flag & mask) == mask

          Examples:
            64  -> include only read1
            128 -> include only read2
            2   -> include only properly paired reads

          Default: None (no include constraint, matches bamCoverage defaults)

  -h, --help
          Print help (see a summary with '-h')
```

---

### `bw-compare --help`

Compare two BigWig files (typically python bamCoverage vs. rust bam-coverage) and report per-chromosome and global differences.

The tool bins both BigWigs with the same bin width and compares the values position by position. It reports several statistics describing how different the signals are.

Output columns:
- `n_over_eps`       Number of bins where `|python - rust| > eps`
- `frac_n_over_eps`  Fraction of bins over eps
- `mean_abs`         Mean absolute difference
- `var_abs`          Variance of absolute differences
- `rmse`             Root mean squared error
- `max_abs`          Maximum absolute difference
- `pearson_rho`      Pearson correlation between tracks

A final **TOTAL** line summarizes all chromosomes.

If `--outfile` is not given, a report file will be created automatically:
`bw_compare_<rust_basename>_w<bin_width>.txt`

```
Usage: bw-compare [OPTIONS] --python-bw <FILE> --rust-bw <FILE>

Options:
      --python-bw <FILE>
          Python-generated BigWig (reference)

      --rust-bw <FILE>
          Rust-generated BigWig (to be validated)

      --bin-width <INT>
          Bin width used during coverage generation (must match both files)

          [default: 50]

      --eps <FLOAT>
          Epsilon threshold for counting a bin as different (|python - rust| > eps)

          [default: 0.00001]

      --outfile <FILE>
          Optional output file for the comparison report.

          If not provided, the report is written to:
          bw_compare_<rust_basename>_w<bin_width>.txt

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

---

## Reproducibility notes

To get reproducible and comparable coverage tracks:

1. **Use a position-sorted BAM**
   - `bam-coverage` expects the BAM to be sorted by chromosome position.
2. **Match bin width**
   - `bam-coverage --width` and `bw-compare --bin-width` must match (and must also match your `bamCoverage --binSize` if you compare to deeptools).
3. **Match filtering semantics**
   - MAPQ threshold: `--min-mapping-quality`
   - Include/exclude: `--sam-flag-exclude` / `--sam-flag-include`
   - Secondary/supplementary/duplicates: `--include-secondary`, `--include-supplementary`, `--include-duplicates`
4. **Normalization**
   - Ensure both tools use the same normalization scheme (e.g. CPM/RPKM/BPM/RPGC) and genome size settings (where applicable).
5. **Floating point**
   - Minor floating point differences can occur; use `bw-compare --eps` to set a meaningful tolerance and rely on correlation (`pearson_rho`) plus RMSE/mean absolute difference for validation.

### Example validation workflow

```bash
# Python reference (deeptools)
bamCoverage input.sorted.bam --samFlagExclude 256 --binSize 50 -o python_input_flag256.bw

# Rust candidate
./target/release/bam-coverage -b input.sorted.bam --sam-flag-exclude 256 --width 50 -o rust_input_flag256.bw

# Compare
./target/release/bw-compare --python-bw python_input_flag256.bw --rust-bw rust_input_flag256.bw --bin-width 50
```

---

## Benchmarks (example)

In a small validation experiment across multiple SAM flag exclusion masks, the Rust implementation produced numerically indistinguishable results (Pearson r ≈ 1.0) while being **~10× faster** and using **~10× less peak memory** than deeptools `bamCoverage` (exact numbers depend on dataset, IO, and chosen flags).

---

## Future developments

Planned improvements to increase parity with deeptools `bamCoverage` and to improve reproducibility/UX:

- **Feature parity** (selected examples)
  - more read extension / fragment handling options (paired-end fragment coverage semantics)
  - region-restricted output (chromosome/interval subsets)
  - smoothing / rolling aggregation options
  - advanced normalization configurations (explicit effective genome size, RPGC parameters)
  - blacklist / region exclusion
- **Better packaging**
  - published releases with prebuilt binaries for common Linux targets
  - container recipes for Apptainer/Singularity and Docker
- **Performance and QA**
  - expanded test suite (known-answer tests and randomized property tests)
  - continuous benchmarking and regression checks
  - improved bw-compare reporting (optional CSV/TSV output, plots)

If you rely on a particular `bamCoverage` option that is currently missing, please open an issue describing:
- the exact deeptools CLI you use
- a small test BAM (or synthetic minimal example)
- the expected behavior

---

## License

See `LICENSE` in this repository.

---

## Acknowledgements

- Inspired by the deeptools `bamCoverage` interface and semantics.
- Uses the Rust ecosystem for HTS parsing and BigWig writing (see `Cargo.toml` for dependencies).
