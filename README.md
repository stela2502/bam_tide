# bam_tide

Fast, memory-efficient BAM coverage calculation in Rust.

`bam_tide` provides a modern, CIGAR-aware implementation of genome coverage binning and exports results as **bedGraph** or **bigWig**.  
It is designed as a high-performance alternative to tools like **deepTools `bamCoverage`**, with explicit and reproducible coverage semantics.

---

## ✨ Features

- 🚀 **Very fast** BAM → bedGraph / bigWig conversion
- 🧬 **CIGAR-correct** reference coverage (M/= /X only)
- 📦 Low memory footprint
- 📏 Fixed-width genomic binning
- 🔬 Numerically comparable to `deepTools bamCoverage`
- 🧪 Built-in benchmarking and regression tests

---

## 📦 Installation

### From source (recommended)

Requires:
- Rust ≥ 1.70
- `cargo`
- `samtools` (for test data / benchmarking)
- `deepTools` (optional, for benchmarking comparison)

```bash
git clone https://github.com/stela2502/bam_tide.git
cd bam_tide
cargo build --release
```

Binary will be at:

```bash
./target/release/bam-coverage
```

---

## ▶️ Usage

### Basic usage

```bash
bam-coverage \
  -b input.bam \
  -o output.bedgraph \
  -w 50
```

Arguments:

| Option | Description |
|--------|-------------|
| `-b, --bam` | Input BAM file (sorted by coordinate) |
| `-o, --outfile` | Output bedGraph or bigWig |
| `-w, --width` | Bin width (default: 50 bp) |
| `-n, --normalize` | Normalization mode (`not`, `cpm`, `rpkm`, `bpm`) |
| `--min-mapping-quality` | Filter reads by MAPQ |
| `--include-secondary` | Include secondary alignments |
| `--include-supplementary` | Include supplementary alignments |
| `--include-duplicates` | Include duplicate-marked reads |

---

## 🧮 Coverage semantics

`bam_tide` computes **binwise coverage as the number of alignment blocks overlapping each bin**:

- CIGAR operations counted: `M`, `=`, `X`
- CIGAR operations ignored: `I`, `S`, `H`, `P`
- Reference skips (`N`) and deletions (`D`) advance reference position but do not contribute coverage
- Each block contributes **+1** to every bin it overlaps

This matches the effective behavior of `deepTools bamCoverage` when run without smoothing, normalization, or fragment extension.

---

## 🧪 Benchmark vs deepTools

Benchmark command:

```bash
bash legacy/benchmark.sh legacy/testData/subset.bam bench_out
```

Results:

```
Rust:
Elapsed: 0.23 s
RSS:     5.5 MB

deepTools:
Elapsed: 38.53 s
RSS:     73 MB
```

### Numeric comparison

```bash
python3 legacy/compare_bedgraph_bins.py bench_out/rust.bedgraph bench_out/python.bedgraph 50
```

Output:

```
== bedGraph binwise comparison ==
bin_width       50
total_bins      54564613
bins_with_signal 854
bad_bins        39
mean_abs        5.57e-06
corr            1.0
```

This shows:
- Perfect correlation with deepTools output
- Only 39 bins differ slightly out of >54 million
- Numerical differences are negligible

---

## 🏃 Running your own benchmark

Requirements:
- `deepTools` installed (`bamCoverage`)
- `/usr/bin/time`

```bash
bash legacy/benchmark.sh your_data.bam my_benchmark 50 4
python3 legacy/bench_table.py my_benchmark
```

Compare outputs:

```bash
python3 legacy/compare_bedgraph_bins.py \
  my_benchmark/rust.bedgraph \
  my_benchmark/python.bedgraph \
  50
```

---

## 🧪 Testing

Run unit and integration tests:

```bash
cargo test
```

This includes:
- CIGAR → block conversion tests
- binning logic tests
- regression tests for bedGraph and bigWig output

---

## 📊 Output formats

Supported outputs:
- `bedGraph`
- `bigWig`

Format is inferred from output filename extension.

---

## 🧠 Design

Pipeline:

```
BAM → AlignmentPolicy filter
    → record_to_blocks()  (CIGAR → RefBlock[])
    → BedData::add_ref_blocks()
    → DataIter
    → bedGraph / bigWig
```

Core structure:

```rust
struct RefBlock {
    start: u32,
    end: u32,   // half-open [start, end)
}
```

Bins accumulate **per-block** overlap (not per-base).

---

## 📜 License

MIT License.

---

## 👤 Author

Stefan Lang  
Division of Molecular Hematology  
Lund University  

---

## 🧬 Citation

If you use this tool in a paper, please cite the repository:

```
https://github.com/stela2502/bam_tide
```

(Manuscript forthcoming.)

---

## 💬 Why “bam_tide”?

Because it’s fast, memory-efficient, and turns BAMs into coverage tracks in one clean sweep 🌊
