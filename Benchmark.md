# bam_tide Benchmark Summary

## Overview

We compared **bam_tide (Rust)** against **deepTools bamCoverage
(Python)** using two datasets:

1.  GL000220.1 (very small chromosome fragment)
2.  pbmc_1k (realistic medium-sized dataset)

For each dataset, multiple SAM flag exclusion scenarios were tested:

-   0
-   256
-   512
-   1024
-   2048

Metrics collected:

-   Absolute differences between output BigWig tracks
-   Pearson correlation
-   Runtime
-   Peak memory usage

------------------------------------------------------------------------

## How to run the benchmark

### Requirements

-   bam_tide (bam-coverage in PATH)
-   deepTools (bamCoverage in PATH)
-   time (GNU /usr/bin/time recommended)

### Example run

``` bash
/usr/bin/time -v bam-coverage -b /tmp/pbmc_1k.bam -o /tmp/rust_pbmc_1k.bw -w 50 -t 3
/usr/bin/time -v bamCoverage -b /tmp/pbmc_1k.bam -o /tmp/python_pbmc_1k.bw --binSize 50
```

------------------------------------------------------------------------

## Multi-threading performance (Rust)

bam_tide supports multi-core processing via the `-t` parameter.

### pbmc_1k results (bin width = 50)

| Threads | Runtime (s) | Eff. cores | Speedup vs t=1 |
|---------|------------|------------|----------------|
| 1       | 76.1       | 1.57       | 1.00×          |
| 3       | 47.8       | 2.78       | 1.59×          |
| 5       | 54.3       | 2.60       | 1.40×          |

### Interpretation

-   Multi-threading provides a clear speedup.
-   Optimal performance for this dataset is around **3 threads**.
-   Increasing threads beyond this introduces overhead and reduces
    performance.
-   Effective CPU usage reaches \~280%, showing strong parallel
    execution.

------------------------------------------------------------------------

## 1. GL000220.1 (Small reference)

### Key results

| Metric                  | Python    | Rust (single processor) |
|:------------------------|:----------|:------------------------|
| Runtime                 | ~46–50 s  | ~4.2–4.4 s              |
| Peak RAM                | ~73 MB    | ~7.3–7.8 MB             |
| Pearson correlation     | 1.0000    | 1.0000                  |
| Max absolute difference | 0         | 0                       |

### Interpretation

-   Rust is \~10× faster
-   Rust uses \~10× less memory
-   Outputs are identical

------------------------------------------------------------------------

## 2. pbmc_1k (Realistic dataset)

### Key results

| Metric                  | Python   | Rust   |
|:------------------------|:---------|:-------|
| Runtime                 | 721.8 s  | 91.3 s |
| Peak RAM                | 190 MB   | 232 MB |
| Pearson correlation     | 1.0      | 1.0    |
| Max absolute difference | 0        | 0      |

### Interpretation

-   Rust is \~7--8× faster
-   Slightly higher memory usage due to in-memory processing
-   Perfect agreement with reference

------------------------------------------------------------------------

## Cross-dataset Summary

| Property    | GL000220.1  | pbmc_1k           |
|:------------|:------------|:------------------|
| Speedup     | ~10×        | ~7–8×             |
| Memory      | ~10× lower  | ~20–30% higher    |
| Correlation | 1.0         | 1.0               |
| Differences | 0           | 0                 |

------------------------------------------------------------------------

## Main Conclusions

### Speed

-   bam_tide is consistently much faster than deepTools

### Parallelism

-   Multi-core scaling is effective up to a dataset-dependent limit
-   Optimal thread count should be chosen empirically

### Accuracy

-   Outputs are identical (correlation 1.0)

------------------------------------------------------------------------

## Overall Assessment

bam_tide provides:

-   Major speed improvements
-   Effective multi-core utilization
-   Identical output

This makes it a strong high-performance alternative to bamCoverage.
