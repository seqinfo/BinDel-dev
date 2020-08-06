library(Rsamtools)
library(tidyverse)
source("util.R")

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and output locations.", call.=FALSE)
  }

# TODO: rm
# bam_location <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/A749N.bam"# args[1]
# bed_location <- "coordinates/chr15.bed" #args[2]
# out_location <- "test/" #args[3]

bam_location <- args[1]
bed_location <- args[2]
out_location <- args[3]

bam_name <- basename(bam_location)

sample <- readGAlignments(bam_location)
bed <- read_tsv(bed_location)
reads <- count(bed, sample, bam_name)

write_tsv(path = paste0(out_location, "count.", basename(bam_name), ".tsv"), x = reads)
