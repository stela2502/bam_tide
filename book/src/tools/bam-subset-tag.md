# bam-subset-tag

`bam-subset-tag` splits a BAM file into multiple BAM files based on the value of a given tag (e.g. CB, CR, UB).

Each input list file defines one output BAM.  
Each line in a list file represents one accepted tag value.

## Status

Utility tool.

Stable and intended for regular use in debugging, filtering, and workflow preparation.

## What it does

- reads a BAM file
- extracts a specified tag from each read
- matches the tag value against one or more input lists
- writes reads into separate BAM files per list

Optionally:

- writes unmatched reads into a separate BAM

## Basic usage

```bash
bam-subset-tag \
  --bam input.bam \
  --tag CR \
  --list cluster_a.txt cluster_b.txt \
  --prefix out/subset_
```

## Input format

Each list file contains one tag value per line:

```text
AAACCCAAGGAGAGTA
TTTGGGCCCAAATTTG
...
```

Empty lines and lines starting with `#` are ignored.

## Output

For each list file:

```text
out/subset_<listname>.bam
```

Example:

```bash
--list cluster_a.txt cluster_b.txt
--prefix out/subset_
```

Produces:

```text
out/subset_cluster_a.bam
out/subset_cluster_b.bam
```

If enabled:

```text
out/subset_unmatched.bam
```

## Options

```text
-b, --bam <BAM>            Input BAM file

-t, --tag <TAG>            BAM tag to match (default: CR)

-l, --list <LISTS>         One or more files containing tag values

-p, --prefix <PREFIX>      Output prefix (directory + filename prefix)

    --write-unmatched      Also write unmatched reads
```

## Typical use cases

### Extract specific cells

```bash
bam-subset-tag \
  --bam full.bam \
  --tag CB \
  --list selected_cells.txt \
  --prefix subset/
```

### Split clusters into separate BAMs

```bash
bam-subset-tag \
  --bam full.bam \
  --tag CR \
  --list cluster_a.txt cluster_b.txt cluster_c.txt \
  --prefix clusters/
```

### Debug problematic reads

```bash
bam-subset-tag \
  --bam full.bam \
  --tag UB \
  --list suspicious_umis.txt \
  --prefix debug/ \
  --write-unmatched
```

## Supported tag types

The tool supports multiple BAM tag types:

- string tags (e.g. CB, CR, UB)
- numeric tags (converted to string)
- character tags

Matching is exact.

## Notes

- tag must be exactly two characters
- input BAM must contain the specified tag
- chromosome naming / alignment is not modified
- output BAMs retain original alignments and flags

## Performance

- streaming BAM processing
- minimal memory overhead
- suitable for large BAM files
- performance scales with disk IO

## Relationship to other tools

`bam-subset-tag` is often used together with:

- [`bam-quant`](./bam-quant.md)  
  → subset cells before quantification

- [`bam-ont-normalizer`](./bam-ont-normalizer.md)  
  → filter normalized ONT reads

- [`bam-coverage`](./bam-coverage.md)  
  → compute coverage on subsets

## Design philosophy

`bam-subset-tag` is intentionally simple:

- no fuzzy matching
- no implicit grouping
- no hidden assumptions

Each output BAM corresponds exactly to a user-defined set of tag values.

This makes it reliable for debugging and reproducible filtering.

## See also

- [`bam-quant`](./bam-quant.md)
- [`bam-coverage`](./bam-coverage.md)
- [`bam-ont-normalizer`](./bam-ont-normalizer.md)