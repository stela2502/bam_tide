# bam-quant

`bam-quant` performs single-cell-aware quantification of BAM files against a GTF-derived splice index.

It is designed as a flexible alternative to fixed pipelines such as Cell Ranger, allowing direct control over how reads are interpreted, matched, and counted.

## Status

`bam-quant` is experimental.

The core ideas and data model are stable, but details of matching, filtering, and output may evolve.

Results should always be validated against known datasets.

## What it does

`bam-quant`:

1. Reads a BAM file (typically 10x-style)
2. Extracts cell barcodes (CB) and UMIs (UB)
3. Matches reads against a splice index derived from a GTF
4. Classifies reads (exonic, intronic, junction-aware)
5. Optionally matches SNPs (if VCF + genome are provided)
6. Produces sparse matrices (scdata format)

It supports:

- gene-level quantification
- transcript-level quantification
- optional intronic separation
- optional SNP ref/alt counting

## Why use it?

Standard pipelines like Cell Ranger:

- require strict reference formats
- hide matching logic
- limit flexibility

`bam-quant` is built for:

- custom GTFs
- alternative transcript models
- exploratory single-cell analysis
- integration with SNP-aware workflows
- explicit control over read classification

## Conceptual model

`bam-quant` works on two key abstractions:

### SplicedRead

A read represented as genomic blocks (exons/introns) derived from BAM alignment.

This is used for:

- gene matching
- transcript matching
- splice-aware classification

### AlignedRead (optional)

A refined, genome-aware representation of the read.

Used for:

- SNP matching
- base-level inspection
- quality filtering

### Splice Index

Built from a GTF file.

Contains:

- genes
- transcripts
- exon structures
- genomic bins for fast lookup

### SNP Index (optional)

Built from a VCF file.

Used to:

- identify SNP positions
- classify observed bases as ref/alt/other

## Pipeline

The internal pipeline is:

1. Load splice index (GTF-derived)
2. Optionally load genome FASTA
3. Optionally load SNP index (VCF)
4. Stream BAM records:
   - extract CB / UB
   - build `SplicedRead`
   - optionally build `AlignedRead`
5. Chunk processing (parallel):
   - match reads to genes/transcripts
   - classify matches
   - optionally match SNPs
6. Accumulate counts per cell
7. Write sparse matrices

## Basic usage

```bash
bam-quant \
  --bam input.bam \
  --index splice.index \
  --outpath quant
```

## Gene vs transcript mode

```bash
--quant-mode gene
--quant-mode transcript
```

Gene mode collapses transcript structure.  
Transcript mode retains isoform-level resolution.

## Intronic reads

```bash
--split-intronic
```

Separates intronic counts into a separate matrix.

Note:

> Intronic classification is currently strict and may undercount in noisy data.

## SNP-aware mode

Requires both genome and VCF:

```bash
bam-quant \
  --bam input.bam \
  --index splice.index \
  --genome genome.fa \
  --vcf variants.vcf \
  --outpath quant
```

Produces:

- main matrix
- SNP reference matrix
- SNP alternate matrix

## Matching behaviour

Matching is controlled by parameters such as:

- strand requirement
- junction chain matching
- overhang tolerance
- allowed intronic gaps

These affect how strictly reads must match transcript structures.

## Output

`bam-quant` writes sparse matrices in `scdata` format.

Typical outputs:

- gene / transcript matrix
- optional intronic matrix
- optional SNP ref matrix
- optional SNP alt matrix

These can be imported into:

- Python (AnnData)
- R (Seurat / SingleCellExperiment)

## Performance

Designed for:

- large BAM files
- multi-core processing
- streaming execution
- minimal memory overhead

Parallelization is done in chunks using Rayon.

## Limitations

- matching rules are still evolving
- intronic classification is strict
- SNP matching depends on read quality and alignment correctness
- requires consistent CB/UB tagging in BAM

## Philosophy

`bam-quant` does not try to hide complexity.

Instead, it exposes:

- how reads are interpreted
- how matches are decided
- how counts are generated

This makes it suitable for:

- debugging pipelines
- custom analysis workflows
- research use cases

## See also

- [`bam-coverage`](./bam-coverage.md)
- [`bam-ont-normalizer`](./bam-ont-normalizer.md)