# 🌊 bam_tide

**bam_tide** is a fast, memory-efficient Rust tool for converting BAM files into binned genomic coverage tracks (e.g. bedGraph / BigWig-ready data).

This repository contains:

* ✅ **Refactored implementation** at the repository root
* 🗂 **Legacy implementation** preserved in [`legacy/`](./legacy) for reference and reproducibility

The refactor focuses on:

* clearer internal structure
* simpler APIs
* improved performance
* maintainability for future features

---

## 🚧 Status

⚠️ **Active refactor / development**

The new implementation is NOT functional and still evolving. Interfaces and CLI options may change until the first stable refactored release.

---

## 📁 Repository Layout

```text
.
├── src/            # Refactored bam_tide implementation
├── tests/          # Tests for the new implementation
├── legacy/         # Original implementation (frozen)
├── Cargo.toml
└── README.md
```

---

## ✨ Features

* Fast streaming BAM processing
* Genome binning with configurable bin width
* Low memory footprint (no full-genome arrays)
* Designed for large BAM files (scRNA-seq, bulk RNA-seq, ChIP-seq, ATAC-seq)
* Suitable for downstream BigWig generation
* Written in Rust for performance and safety

---

## 📦 Installation

Clone and build from source:

```bash
git clone https://github.com/stela2502/bam_tide.git
cd bam_tide
cargo build --release
```

The binary will be located at:

```text
target/release/bam_tide
```

---

## ▶️ Usage

Example:

```bash
bam_tide \
  --bam input.bam \
  --bin-width 100 \
  --out coverage.bedgraph
```

(Exact CLI options may change during refactoring; run `bam_tide --help` for current usage.)

---

## 🧬 Legacy Implementation

The original implementation is preserved under:

```text
legacy/
```

This allows:

* reproducibility of previous results
* comparison with the refactored implementation
* gradual migration of features

The legacy code is frozen and no longer actively developed.

---

## 🛣 Roadmap

Planned improvements:

* Stable CLI interface
* Direct BigWig output
* Improved multi-threading
* Region / chromosome filtering
* Better test coverage

---

## 📜 License

MIT

---

## 🙌 Acknowledgements

Developed as part of ongoing work on efficient genomics tooling in Rust.

