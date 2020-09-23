#!/bin/bash

# Exit when any command fails
set -e

conda create -n detector bioconductor-rsamtools r-tidyverse r-fs bioconductor-genomicalignments bioconductor-bsgenome.hsapiens.ucsc.hg38 cromwell r-here r-gridextra r-purrr bioconductor-dnacopy bioconductor-txdb.hsapiens.ucsc.hg38.knowngene r-fuzzyjoin r-scales

conda install -c hcc r-depmixs4
