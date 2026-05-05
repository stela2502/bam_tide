# bw-compare

`bw-compare` compares coverage tracks (e.g. BigWig or bedGraph) to evaluate differences between methods or tools.

It is primarily intended to validate `bam-coverage` against established tools such as deepTools.

## Status

Utility tool.

Stable enough for internal validation and benchmarking, but not intended as a general-purpose track comparison framework.

## What it does

`bw-compare` takes two coverage tracks and computes differences.

Typical use cases:

- compare `bam-coverage` output vs `deepTools bamCoverage`
- detect systematic biases
- verify correctness after code changes
- benchmark parameter effects

## Why it exists

When reimplementing core functionality (like coverage calculation), you need a way to verify:

- numerical agreement
- systematic differences
- edge-case behaviour

`bw-compare` provides a simple way to do that.

## Basic usage

```bash
bw-compare \
  --bw1 sample_rust.bw \
  --bw2 sample_python.bw
```

## Typical workflow

```bash
bam-coverage \
  --bam input.bam \
  --outfile rust.bw

bamCoverage \
  -b input.bam \
  -o python.bw

bw-compare \
  --bw1 rust.bw \
  --bw2 python.bw
```

## Output

The tool reports:

- overall similarity metrics
- per-region differences
- potential outliers

(Exact output depends on implementation details.)

## Interpretation

Small differences are expected due to:

- floating point rounding
- binning strategies
- handling of edge cases

Large differences usually indicate:

- parameter mismatch
- filtering differences
- bugs

## Limitations

- assumes comparable input tracks
- does not normalize automatically
- not intended for large-scale statistical comparison
- primarily a developer / validation tool

## Role in bam_tide

`bw-compare` is not a primary analysis tool.

It is part of the development workflow:

- validate new features
- ensure reproducibility
- guard against regressions

## See also

- [`bam-coverage`](./bam-coverage.md)