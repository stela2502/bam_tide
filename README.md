# bam_tide

The two programs bam2bigwig and bam2bedgraph do exactly what they claim to do - they parse a bam file and count the reads.
At the moment no normalization is implemented and there is also no way to influence how the reads are collected:

Each read will be counted for the whole region starting at feature.start to feature.start + read.sequence.length.


This project has been started as bam to bigwig conversion took quite a lot of time on our system.
It is in a very early state and not tested a lot.

# Insttallation

You need the Rust compiler to compile this software.
The easiest way to install it is to state

```
cargo install --git https://github.com/stela2502/bam_tide.git
```

# Usage

There are two programs included that do the same thing, but create two different outputs:

## bam2bedgraph

This program collects all reads from the bam file and exports only those regions that did match to a read in the bedgraph.

```
bam2bedgraph  -h
A command-line tool for converting BAM coverage data to BigWig format

Usage: bam2bedgraph [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>          Input BAM file (sorted by chromosome position)
  -o, --outfile <OUTFILE>  Output BigWig file
  -w, --width <WIDTH>      Bin width for coverage calculation (default: 50bp) [default: 50]
  -h, --help               Print help
  -V, --version            Print version
```

## bam2bigwig

This program collects all reads from the bam file and exports all possible regions to the bigwig out file.

```
bam2bigwig  -h
A command-line tool for converting BAM coverage data to BigWig format

Usage: bam2bigwig [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>          Input BAM file (sorted by chromosome position)
  -o, --outfile <OUTFILE>  Output BigWig file
  -w, --width <WIDTH>      Bin width for coverage calculation (default: 50bp) [default: 50]
  -h, --help               Print help
  -V, --version            Print version
```


# Thanks

Jack Huey for help on how to use his bigtools library ( https://github.com/jackh726/bigtools ).