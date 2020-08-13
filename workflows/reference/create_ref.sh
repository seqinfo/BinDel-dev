#!/bin/bash

#SBATCH --job-name=SyndromeDetectorRef
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=5-00:00:00
#SBATCH --partition=amd

# Exit when any command fails
set -e

module load python-3.6.3
source activate detector

export count_reads=$(realpath ../../R/count_reads.R)

cromwell -Dconfig.file=cromwell.conf -Xmx12g run -i inputs_ref.json -o options_ref.json \
  reference.wdl

