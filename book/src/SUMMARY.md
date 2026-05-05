# bam_tide

`bam_tide` is a Rust toolkit for fast and flexible processing of BAM files across bulk, single-cell, and long-read sequencing workflows.

It is built around a simple idea:

> Work directly on BAM files, make assumptions explicit, and keep the logic transparent.

---

## 🧠 Why bam_tide exists

Modern sequencing workflows often rely on complex pipelines that:

- convert data between formats (BAM → FASTQ → BAM)
- enforce strict reference models
- hide how reads are interpreted
- are difficult to debug or extend

This becomes especially problematic for:

- single-cell data (CB / UB handling)
- long-read data (ONT)
- custom annotations and transcript models

`bam_tide` takes a different approach.

---

## 🔬 Core ideas

### BAM-first

BAM is treated as the central data format.

- no unnecessary conversion to FASTQ
- all processing happens directly on alignments
- metadata (tags, flags, structure) is preserved

---

### Explicit models

Instead of hiding logic inside pipelines, `bam_tide` makes it visible:

- how reads are matched to genes
- how splice structures are interpreted
- how SNPs are detected
- how tags (CB, UB, etc.) are used

---

### Composable tools

Each tool has a clear role:

- prepare data
- transform reads
- quantify signal
- validate results

Tools are designed to be combined:

```text
ONT BAM
  ↓
bam-ont-normalizer
  ↓
external mapping
  ↓
bam-subset-tag
  ↓
bam-quant
  ↓
analysis
```

---

### Performance

`bam_tide` is written in Rust:

- streaming BAM processing
- multi-threaded execution
- low memory overhead

This enables fast processing of large datasets.

---

## 🔧 Overview of tools

`bam_tide` currently includes:

### Core tools

- **bam-coverage**  
  Generate coverage tracks from BAM files.

- **bam-quant**  
  Perform single-cell-aware quantification against a splice index.

- **bam-ont-normalizer**  
  Normalize ONT reads into a consistent structure.

---

### Infrastructure

- **gtf-splice-index**  
  Build a splice index from annotation files.

---

### Utilities

- **bam-subset-tag**  
  Split BAM files based on tag values.

---

### Validation / development

- **bw-compare**  
  Compare coverage tracks.

- **bam-quant-testdata**  
  Generate synthetic test data with known ground truth.

---

## 🧬 Typical workflows

### Coverage

```text
BAM → bam-coverage → BigWig
```

---

### Single-cell quantification

```text
BAM → gtf-splice-index → bam-quant → matrix
```

---

### ONT integration

```text
ONT BAM → bam-ont-normalizer → bam-quant
```

---

### Debugging / filtering

```text
BAM → bam-subset-tag → smaller BAM → analysis
```

---

## ⚠️ Status

`bam_tide` is under active development.

- `bam-coverage` is the most stable component
- `bam-quant` is evolving
- ONT and SNP workflows are experimental

Expect changes and iteration.

---

## 🧪 Testing philosophy

`bam_tide` includes synthetic data generation to enable:

- end-to-end validation
- reproducible debugging
- regression testing

This is critical for verifying correctness in complex workflows.

---

## 🧠 Design philosophy

- **BAM-first** instead of pipeline-first  
- **explicit logic** instead of hidden heuristics  
- **small tools** instead of monolithic frameworks  
- **performance-aware** but not opaque  
- **debuggable at every step**

---

## 👤 Author

Stefan Lang  
Division of Molecular Hematology, Lund University

---

## Where to go next

- [Installation](installation.md)
- [Tools](SUMMARY.md)
- Individual tool documentation
