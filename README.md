# PPDx

## Prerequisites
1. System with Slurm Workload Manager
2. Coordinates of regions to analyze (in .bed format)


## Dependencies
1. [BSgenome.Hsapiens.UCSC.hg38](https://anaconda.org/bioconda/bioconductor-bsgenome.hsapiens.ucsc.hg38)
2. [Cromwell](https://anaconda.org/bioconda/cromwell)
3. [GenomicAlignments](https://anaconda.org/bioconda/bioconductor-genomicalignments)
4. [here](https://anaconda.org/conda-forge/r-here)
5. [Rsamtools](https://anaconda.org/bioconda/bioconductor-rsamtools)
6. [Scales](https://anaconda.org/r/r-scales)
7. [Tidyverse](https://anaconda.org/r/r-tidyverse)

```./create_env.sh```under the workflows can be used to install required dependencies (requires [Anaconda](https://anaconda.org/)).

## Usage
1. Pre-process the samples to .bam format (not part of this package)
2. Create reference
3. Analyze the sample of interest


## Scripts
1. ```count_reads.R``` - Script for binning and counting reads.
2. ```analyse_bam.R``` - Script for analyzing .bam file with reference group.
3. ```util.R``` - Central script for providing tools to bin, count and find GC%.
4. ```test.R``` - Convenience script for locally testing.
5. ```archived.R``` - Possible future additions. Currently has chi-squared variation transformation, which is still in the development phase.

