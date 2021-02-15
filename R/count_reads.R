library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)
library(Rsamtools)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)


if (length(args) != 2) {
  stop("Please provide (1) binnable BAM file (.bam) and (2) regions to bin (.bed)",
       call. = FALSE)
}


read_bam_counts <- function(bam_location) {
  param <-
    ScanBamParam(
      flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
      what = c("pos")
    )
  
  return(readGAlignments(bam_location, param = param))
}



bin_counts <- function(bam_location, bed) {
  bam <- read_bam_counts(bam_location)
  binned_counts <- bed %>%
    mutate(reads = assay(
      summarizeOverlaps(makeGRangesFromDataFrame(.), bam, mode = "IntersectionStrict")
    )) %>%
    mutate(sample = basename(bam_location))
  
  return(binned_counts)
}



bam_location <- args[1]
bed_location <- args[2]

bed <- read_tsv(bed_location)

reads_per_bin <- bin_counts(bam_location, bed) %>%
  group_by(sample) %>% 
  mutate(sample = group_indices(., sample)) %>%
  mutate(sample = paste0("ref.", sample)) %>% 
  ungroup()


write_tsv(path = paste0(basename(bam_location), ".tsv"), x = reads_per_bin)
