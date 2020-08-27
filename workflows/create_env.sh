#!/bin/bash

# Exit when any command fails
set -e

conda create -n detector bioconductor-rsamtools r-tidyverse r-fs bioconductor-genomicalignments bioconductor-bsgenome.hsapiens.ucsc.hg38 cromwell r-here r-gridextra r-purrr

