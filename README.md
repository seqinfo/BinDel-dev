# PPDx

## Prerequisites
1. System with Slurm Workload Manager
2. Coordinates of regions to analyze (in .bed format)

## Dependencies
1. [BSgenome.Hsapiens.UCSC.hg38](https://anaconda.org/bioconda/bioconductor-bsgenome.hsapiens.ucsc.hg38)
2. [Cromwell](https://anaconda.org/bioconda/cromwell)
4. [GenomicAlignments](https://anaconda.org/bioconda/bioconductor-genomicalignments)
5. [Rsamtools](https://anaconda.org/bioconda/bioconductor-rsamtools)
6. [Tidyverse](https://anaconda.org/r/r-tidyverse)

```./create_env.sh```under the workflows can be used to install required dependencies (requires [Anaconda](https://anaconda.org/)).

## Usage
1. Pre-process the samples to .bam format (not part of this package)
2. Create reference
3. Analyze the sample of interest


## Scripts
1. ```count_reads.R``` - Script for binning and counting reference group reads.
2. ```analyse_bam.R``` - Script for analyzing .bam file with reference group.
