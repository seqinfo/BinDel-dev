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
analyze_locations_bed="coordinates/chr15.bed"
reference_location="reference.tsv"


Rscript analyze_bam.R $bam_location $analyze_locations_bed $reference_location

Rscript visualize.R results.$sample_name.tsv $sample_name "coordinates/as_pws.bed" "AS/PWS"

