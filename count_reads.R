library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and output locations.", call.=FALSE)
  }

bam_location <- args[1]
bed_location <- args[2]
out_location <- args[3]

bam_name <- basename(bam_location)

sample <- readGAlignments(bam_location)
bed <- read_tsv(bed_location)

reads <- bed %>% 
  mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
  mutate(normalized_by_sample = reads/sum(reads)) %>% 
  mutate(sample = bam_name)

write_tsv(path = paste0(out_location, "count.", basename(bam_name), ".tsv"), x = reads)

