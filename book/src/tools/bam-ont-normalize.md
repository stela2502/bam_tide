# bam-ont-normalizer

`bam-ont-normalizer` processes ONT-derived BAM files and converts reads into a simplified, consistent structure suitable for downstream single-cell-style workflows.

It is designed to clean up long-read sequencing data so that it behaves more like the assumptions made by typical 10x / short-read pipelines.

## Status

Early-stage / experimental.

The core idea is stable, but heuristics for detecting structure (barcodes, UMI, transcript boundaries) are still evolving.

This tool should be used for exploration and method development, not yet for production pipelines without validation.

## What it does

`bam-ont-normalizer`:

- reads ONT BAM files (e.g. from Dorado)
- identifies 10x-style structures within reads:
  - cell barcode (CB)
  - UMI (UB)
  - polyT / adapter regions
- normalizes reads into a consistent layout
- ensures each read represents a single structured molecule
- writes a cleaned BAM with standardized tags

The goal is to turn heterogeneous ONT reads into something that downstream tools can interpret reliably.

## Why use it?

ONT BAM files often:

- contain complex alignments (indels, large skips, multi-block structures)
- include multiple logical “segments” in a single read
- vary in orientation and structure
- lack consistent tagging (CB / UB)

This makes direct use in single-cell pipelines difficult.

`bam-ont-normalizer` attempts to:

- simplify reads
- enforce a consistent structure
- recover useful metadata (CB / UB)
- reduce ambiguity before quantification

## Conceptual model

The tool assumes that valid reads should look like:

```
[adapter] [cell barcode] [UMI] [polyT] [transcript sequence]
```

The normalization process:

1. detects candidate regions within the read
2. aligns / matches against expected patterns
3. extracts CB and UMI
4. enforces a single logical read per molecule
5. removes or trims inconsistent parts
6. outputs a normalized representation

## Output

The output is a BAM file where:

- each read is simplified
- CB / UB tags are set where possible
- orientation is consistent
- ambiguous or malformed reads are filtered or trimmed

This BAM is intended for direct use with:

- `bam-quant`
- downstream single-cell workflows
- exploratory analysis

## Basic usage

```bash
bam-ont-normalizer \
  --bam input.bam \
  --outfile normalized.bam
```

## Typical workflow

```bash
bam-ont-normalizer \
  --bam dorado_unmapped.bam \
  --outfile normalized.bam

bam-quant \
  --bam normalized.bam \
  --index splice.index \
  --outpath quant
```

## Design goals

- operate directly on BAM (no FASTQ roundtrip)
- recover structure from noisy reads
- produce deterministic, reproducible output
- remain transparent and debuggable

## Limitations

- relies on heuristic pattern detection
- barcode/UMI detection may fail in low-quality reads
- complex structural artifacts may not be fully resolved
- assumptions are currently tuned for 10x-like chemistry

## Philosophy

Rather than forcing ONT data into existing pipelines,  
`bam-ont-normalizer` adapts the data to make its structure explicit.

This makes it possible to:

- test new workflows
- integrate long-read data into single-cell pipelines
- debug assumptions about read structure

## Future directions

- improved adapter / barcode detection
- better handling of multi-segment reads
- tighter integration with `bam-quant`
- optional QC reporting

## See also

- [`bam-quant`](./bam-quant.md)
- [`bam-coverage`](./bam-coverage.md)