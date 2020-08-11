#!/bin/bash

#SBATCH --job-name=AnalyzeBAM
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --partition=amd

# Exit when any command fails
set -e

module load python-3.6.3
source activate detector

bam_location=$1
sample_name=$2


Rscript ../R/analyse_bam.R $bam_location "coordinates/chr15.bed" "reference.tsv"

Rscript ../R/visualize.R results.$sample_name.tsv "coordinates/as_pws.bed"

