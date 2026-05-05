# gtf-splice-index

`gtf-splice-index` builds and inspects splice indexes from genome annotation files.

The resulting index is used by `bam-quant` to match aligned reads against genes, transcripts, exons, introns, and splice junction structures.

## Status

Experimental, but central to `bam-quant`.

The index format and matching logic may still evolve.

## What it does

`gtf-splice-index` reads an annotation file and builds a serialized splice index.

Supported annotation-style inputs include:

- GTF
- GFF
- GFF3

The index stores:

- genes
- transcripts
- exon blocks
- chromosome / contig structure
- binned genomic lookup tables

This makes read matching much faster during quantification.

## Why it exists

`bam-quant` should not parse a full GTF file for every run.

Instead, the annotation is preprocessed once:

```bash
gtf-splice-index build \
  --annotation genes.gtf \
  --index genes.splice.idx
```

Then reused:

```bash
bam-quant \
  --bam sample.bam \
  --index genes.splice.idx \
  --outpath quant
```

## Build an index

```bash
gtf-splice-index build \
  --annotation genes.gtf \
  --index genes.splice.idx
```

By default, the bin width is:

```text
1,000,000 bp
```

You can change it:

```bash
gtf-splice-index build \
  --annotation genes.gtf \
  --index genes.splice.idx \
  --bin-width 1000000
```

## Inspect an index

```bash
gtf-splice-index stats \
  --index genes.splice.idx
```

This prints a summary of the serialized index, for example:

```text
SpliceIndex: 78691 genes, 507365 transcripts, 25 chromosomes, bin_width=1000000 bp
```

## Attribute keys

Different annotation sources use different attribute names.

For standard GTF files, the defaults are:

| Meaning | Default key |
|---|---|
| Gene ID | `gene_id` |
| Gene name | `gene_name` |
| Transcript ID | `transcript_id` |
| Transcript name | `transcript_name` |
| Exon feature type | `exon` |

For GFF3-style annotations, exon-to-transcript linkage usually uses:

| Meaning | Default key |
|---|---|
| Parent transcript | `Parent` |

## Custom attribute keys

You can override attribute keys if your annotation uses non-standard names.

Example:

```bash
gtf-splice-index build \
  --annotation annotation.gff3 \
  --index annotation.splice.idx \
  --gene-id-key ID \
  --gene-name-key Name \
  --transcript-id-key transcript_id \
  --parent-key Parent
```

Multiple keys can be provided.

Example:

```bash
gtf-splice-index build \
  --annotation annotation.gtf \
  --index annotation.splice.idx \
  --gene-id-key gene_id geneID GeneID \
  --gene-name-key gene_name geneName Name
```

The index builder will try the provided keys when extracting identifiers and names.

## Exon feature types

By default, only features with type `exon` are treated as exon blocks.

You can provide additional feature types:

```bash
gtf-splice-index build \
  --annotation annotation.gff3 \
  --index annotation.splice.idx \
  --exon-feature-type exon CDS
```

Use this carefully. Including `CDS` changes the biological meaning of the index and may affect quantification.

## Relationship to bam-quant

`bam-quant` depends on this index.

The splice index allows `bam-quant` to:

- find candidate genes quickly
- match reads to transcript structures
- classify exon-compatible reads
- detect exact junction chains
- distinguish intronic from exonic signal

## Notes and limitations

The quality of the index depends on the quality of the annotation.

Common problems include:

- missing transcript IDs
- inconsistent gene/transcript relationships
- mixed GTF/GFF3 conventions
- contig names that do not match the BAM
- annotations with unexpected feature types

If read matching gives surprising results, inspect both:

1. the annotation used to build the index
2. the chromosome names in the BAM header

## See also

- [`bam-quant`](./bam-quant.md)