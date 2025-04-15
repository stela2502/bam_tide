# bam_tide

The two programs bam2bigwig and bam2bedgraph do exactly what they claim to do - they parse a bam file and count the reads.
At the moment no normalization is implemented and there is also no way to influence how the reads are collected:

Each read will be counted for the whole region starting at feature.start to feature.start + read.sequence.length.

This project has been started as bam to bigwig conversion took quite a lot of time on our system.
It is in a very early state and not tested a lot.

The other incuded program is bam2mtx and this quantifies a variety of bam files and put the resulting counts into the widely accepted [Matrix Market MTX format](https://math.nist.gov/MatrixMarket/formats.html) that is also created by e.g. CellRanger and is accepted as input from e.g. both the Python/Scanpy and R/Seurat packges for single cell analysis.

# Installation

You need to [install the Rust compiler](https://www.rust-lang.org/tools/install) to compile this software.

The easiest way to compile and install this package is:

```
cargo install --git https://github.com/stela2502/bam_tide.git
```

# Usage

## bam2bedgraph

This program collects all reads from the bam file and exports only those regions that did match to a read in the bedgraph.

```
bam2bedgraph -h
A command-line tool for converting BAM coverage data to BigWig format

Usage: bam2bedgraph [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>                      Input BAM file (sorted by chromosome position)
  -o, --outfile <OUTFILE>              Output BigWig file
  -c, --cell-tag <CELL_TAG>            tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
  -u, --umi-tag <UMI_TAG>              tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
  -a, --analysis-type <ANALYSIS_TYPE>  Collect single cell info or bulk [default: bulk] [possible values: bulk, single-cell]
  -n, --normalize <NORMALIZE>          Normalize the data somehow [default: not] [possible values: not]
  -w, --width <WIDTH>                  Bin width for coverage calculation (default: 50bp) [default: 50]
  -h, --help                           Print help
  -V, --version                        Print version

```

## bam2bigwig

This program collects all reads from the bam file and exports all possible regions to the bigwig out file.

```
bam2bigwig -h
A command-line tool for converting BAM coverage data to BigWig format

Usage: bam2bigwig [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>                      Input BAM file (sorted by chromosome position)
  -o, --outfile <OUTFILE>              Output BigWig file
  -c, --cell-tag <CELL_TAG>            tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
  -u, --umi-tag <UMI_TAG>              tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
  -a, --analysis-type <ANALYSIS_TYPE>  Collect single cell info or bulk [default: bulk] [possible values: bulk, single-cell]
  -n, --normalize <NORMALIZE>          Normalize the data somehow [default: not] [possible values: not]
  -w, --width <WIDTH>                  Bin width for coverage calculation (default: 50bp) [default: 50]
  -h, --help                           Print help
  -V, --version                        Print version

```

## bam2mtx

```
bam2mtx -h
This tool can quantify a bam file from both single cell data as well as bulk data. The default is to quantify single cell data from using the velocyto input format with these default settings: --cell_tag CB --umi_tag UB --analysis_type single-cell --match_type exact --gtf_type genes. Switching to --analysis_type "bulk" will create the mtx out files for one single cell id: "1" Switching to --match_type "overlap" will quantify reads in a sticky way - any overlap will be called a match - even paired ones!! Switching to --gtf_type "exons" will quantify exons instead of genes - make sure our exons are uniquely named!

Usage: bam2mtx [OPTIONS] --bam <BAM> --gtf <GTF> --outpath <OUTPATH> --min-umi <MIN_UMI>

Options:
  -b, --bam <BAM>                      the bam file to quantify
  -g, --gtf <GTF>                      the gtf file fitting to the Bam file (text or gzipped)
  -o, --outpath <OUTPATH>              the outpath
  -m, --min-umi <MIN_UMI>              the minimum (UMI) reads per cell (sample + genes + antibody combined)
  -n, --num-proc <NUM_PROC>            used processor cores (default all)
  -c, --cell-tag <CELL_TAG>            tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
  -u, --umi-tag <UMI_TAG>              tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
  -q, --qual <QUAL>                    For mutation collection please give me a quality cutoff for accepting a nucl as a valid mutation (20? 30?)
  -a, --analysis-type <ANALYSIS_TYPE>  Collect single cell info or bulk [default: single-cell] [possible values: bulk, single-cell]
      --match-type <MATCH_TYPE>        Match only inside exons or overlapping? [default: exact] [possible values: overlap, exact]
      --gtf-type <GTF_TYPE>            Group gtf info into genes or e.g. singel cell quantifications (genes) or use each exon on it's own for e.g. TE analyses (exon) [default: genes] [possible values: genes, exons]
  -h, --help                           Print help
  -V, --version                        Print version

```

## bowtie2mtx

```
bowtie2mtx -h
This tool can quantify a bam file from both single cell data as well as bulk data (untested). The main aim of this is to quantify mutations in the bam file. with these default settings: --match_type exact --gtf_type genes. Switching to --analysis_type "bulk" will create the mtx out files for one single cell id: "1" (untested) Switching to --match_type "overlap" will quantify reads in a sticky way - any overlap will be called a match - even paired ones!! Switching to --gtf_type "exons" will quantify exons instead of genes - make sure our exons are uniquely named!

Usage: bowtie2mtx [OPTIONS] --bam <BAM> --gtf <GTF> --outpath <OUTPATH> --min-umi <MIN_UMI>

Options:
  -b, --bam <BAM>                      the bam file to quantify
  -f, --fasta <FASTA>                  the fasta file needed for mutation anaylsis
  -g, --gtf <GTF>                      the gtf file fitting to the Bam file (text or gzipped)
  -o, --outpath <OUTPATH>              the outpath
  -m, --min-umi <MIN_UMI>              the minimum (UMI) reads per cell (sample + genes + antibody combined)
  -n, --num-proc <NUM_PROC>            used processor cores (default all)
  -q, --qual <QUAL>                    For mutation collection please give me a quality cutoff for accepting a nucl as a valid mutation (20? 30?)
  -a, --analysis-type <ANALYSIS_TYPE>  Collect single cell info or bulk [default: single-cell] [possible values: bulk, single-cell]
      --match-type <MATCH_TYPE>        Match only inside exons or overlapping? [default: exact] [possible values: overlap, exact]
      --gtf-type <GTF_TYPE>            Group gtf info into genes or e.g. singel cell quantifications (genes) or use each exon on it's own for e.g. TE analyses (exon) [default: genes] [possible values: genes, exons]
  -h, --help                           Print help
  -V, --version                        Print version
```

# Thanks

Jack Huey for help on how to use his bigtools library ( https://github.com/jackh726/bigtools ).