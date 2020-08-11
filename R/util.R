library(Rsamtools)
library(dplyr)
library(readr)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg38)

bin_counts <- function(bam_location, bed_location) {
  bed <- read_tsv(bed_location)
  bam <- read_bam_counts(bam_location)
  
  binned_counts <- bed %>%
    mutate(reads = assay(
      summarizeOverlaps(makeGRangesFromDataFrame(.), bam, mode = "IntersectionStrict")
    )) %>%
    mutate(sample = basename(bam_location))
  
  return(binned_counts)
}


find_gc <- function(bed_location) {
  reads <- read_tsv(bed_location) %>%
    mutate(gc = letterFrequency(
      getSeq(BSgenome.Hsapiens.UCSC.hg38,
             GRanges(
               chromosome, IRanges(start = start, end = end)
             )),
      "GC",
      as.prob = T
    )) %>%
    mutate(gc = round(gc, 1))
}


read_bam_counts <- function(bam_location) {
  param <-
    ScanBamParam(
      flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
      what = c("pos")
    )
  
  return(readGAlignments(bam_location, param = param))
}

