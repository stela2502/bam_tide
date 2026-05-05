# bam-coverage

`bam-coverage` summarizes aligned BAM reads into genomic coverage tracks.

It is intended as a fast Rust-based alternative to tools such as `deepTools bamCoverage`, with a focus on performance, simple command-line usage, and direct BAM-to-track conversion.

## Status

`bam-coverage` is currently the most mature tool in `bam_tide`.

It is suitable for regular use, but output should still be validated when changing parameters, reference genomes, or downstream workflows.

## What it does

`bam-coverage` reads an input BAM file and produces binned coverage output.

Typical use cases include:

- creating genome-wide coverage tracks
- exporting coverage for visualization
- comparing BAM-level signal between samples
- producing `bedGraph` or `BigWig`-style output for downstream tools

## Why use it?

`bam-coverage` was written because many existing coverage tools are flexible but slow, especially when used repeatedly in larger workflows.

The design goals are:

- fast BAM processing
- low overhead
- simple CLI usage
- useful defaults
- direct support for common coverage-track workflows

## Basic usage

```bash
bam-coverage \
  --bam sample.bam \
  --outfile sample.bw
```

The output format is inferred from the output filename where supported.

For example:

```bash
bam-coverage \
  --bam sample.bam \
  --outfile sample.bedgraph
```

## Unsorted BAM input

One important feature of `bam-coverage` is that it can process unsorted BAM files.

```bash
bam-coverage \
  --bam input.unsorted.bam \
  --outfile sample.bw
```

This can be useful when working with intermediate BAM files where sorting would be expensive or unnecessary.

## Typical workflow

A common workflow is:

1. align reads to a reference genome
2. produce a BAM file
3. run `bam-coverage`
4. visualize or compare the resulting coverage track

Example:

```bash
bam-coverage \
  --bam aligned.bam \
  --outfile aligned.coverage.bw
```

## Relationship to deepTools bamCoverage

`bam-coverage` targets a similar problem space as `deepTools bamCoverage`:

- read BAM
- summarize coverage
- export a genome-wide signal track

The main difference is implementation philosophy.

`bam-coverage` is written in Rust and is designed to be fast and lightweight. It does not aim to reimplement every option from `deepTools`, but instead focuses on the operations that are most useful in high-throughput BAM processing workflows.

## Performance

Performance depends on:

- BAM size
- compression level
- storage speed
- number of threads
- output format
- whether the BAM is sorted
- reference/genome size

In practical use, `bam-coverage` is expected to be much faster than Python-based coverage workflows for common BAM summarization tasks.

## Output

Depending on the selected output path and enabled features, `bam-coverage` can write coverage-style track files such as:

- `bedGraph`
- `BigWig`

These files can be used in genome browsers or downstream analysis workflows.

## Notes and limitations

`bam-coverage` is intended to be a practical high-performance coverage tool, not a full clone of every `bamCoverage` option.

When reproducing published analyses, always verify that binning, filtering, normalization, and output settings match the original workflow.

## See also

- [`bw-compare`](./bw-compare.md), if available in your build
- [`bam-quant`](./bam-quant.md)
- [`bam-ont-normalizer`](./bam-ont-normalizer.md)