#!/bin/bash

#SBATCH --job-name=AnalyzeBAM
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=22:10:00
#SBATCH --partition=amd

# Exit when any command fails
set -e

module load python-3.6.3
source activate detector

export analyser=$(realpath ../../R/analyse_bam.R)

cromwell -Dconfig.file=cromwell.conf -Xmx12g run -i inputs.json -o options.json \
  analyze.wdl

