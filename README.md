# SyndromeDetector

## Prerequisites
1. System with Slurm Workload Manager
2. Coordinates of regions to analyze (in .bed format)

## Dependencies
1. [BSgenome.Hsapiens.UCSC.hg38](https://anaconda.org/bioconda/bioconductor-bsgenome.hsapiens.ucsc.hg38)
2. [Cromwell](https://anaconda.org/bioconda/cromwell)
3. [Fuzzyjoin](https://anaconda.org/conda-forge/r-fuzzyjoin)
4. [GenomicAlignments](https://anaconda.org/bioconda/bioconductor-genomicalignments)
5. [here](https://anaconda.org/conda-forge/r-here)
6. [Rsamtools](https://anaconda.org/bioconda/bioconductor-rsamtools)
7. [Tidyverse](https://anaconda.org/r/r-tidyverse)
8  [Txdb.Hsapiens.UCSC.hg38.Knowngene](https://anaconda.org/bioconda/bioconductor-txdb.hsapiens.ucsc.hg38.knowngene)
9. [Gridextra](https://anaconda.org/r/r-gridextra)
10. [DNAcopy](https://anaconda.org/bioconda/bioconductor-dnacopy)

```./create_env.sh```under the workflows can be used to install required dependencies (requires [Anaconda](https://anaconda.org/)).

## Usage
1. Pre-process the samples to .bam format (not part of this package)
2. Create reference
3. Analyze the sample of interest

#### Reference creation
1. Fill ```inputs_ref.json``` with sample paths to use as a reference group and coordinates with the key ```Main.bin.regions``` over which the analysis will be conducted.
2. Update ```create_ref.sh``` slurm comments, ```cromwell.conf``` and ```options_ref.json``` (if needed, not mandatory).
3. Run ```./create_ref.sh```

#### Analyze
1. Fill ```inputs.json``` with sample paths to analyze, provide the coordinates used in the reference creation with the key ```Main.coordinates``` and the coordinates of interest ```Main.coordinates_of_interest```, which must overlap the ```Main.coordinates```. ```Main.coordinates_of_interest``` will be used for visualization.
2. Update ```analyze_bam.sh``` slurm comments, ```cromwell.conf``` and ```options.json``` (if needed, not mandatory).
3. Run ```./analyze_bam.sh```

## Scripts
1. ```count_reads.R``` - Script for binning and counting reads.
2. ```analyse_bam.R``` - Script for analyzing .bam file with reference group. Currently supports local Z-score, reference based Z-score and ZZ-score. Before calculating these scores, reads are GC-corrected, normalized by sample and normalized by bin length.
3. ```visualize.R``` - Script for visualization purposes.
4. ```util.R``` - Central script for providing tools to bin, count and find GC%.
5. ```test.R``` - Convenience script for locally testing.
6. ```archived.R``` - Possible future additions. Currently has chi-squared variation transformation, which is still in the development phase.

## Extra information
1. Currently provided bed regions are from A534N.
2. D129N	T15
3. C075N	AS?
4. B869N	AS?
5. E121N	low risk
6. E149N	low risk

