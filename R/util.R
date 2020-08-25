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


#' Sum of Z-scores of a n^th long window. 
#' 
#' For example, region with length 4 with bin Z-scores of "-5, -5, -9, -10" 
#' would yield Z-score of -29.
#' 
#' With Z-scores of "-5, -5, 0, -10" it would be 0 as 0 denotes that region is not 
#' sequential or has other issues/limitations/filters.
#' 
#' Takes an n number of dplyr columns as a parameters.
#'
windowed_sum <- function(...) {
  p <- 1
  s <- 0
  for (i in list(...)) {
    p <- i * p
    s <- i + s
  }
  return(if_else(p == 0, 0, s, missing = 0))
}

#' Check if a region is eligible to count as an region.
#' 
#' For example, region with length 4 with bin Z-scores of "-5, -5, -9, -10" 
#' would yield 1.
#' 
#' With Z-scores of "-5, -5, 0, -10" it would yield 0 as 0 denotes that region is not 
#' eligible to be considered as an region.
#' 
#' #' Takes an n number of dplyr columns as a parameters.
windowed_count <- function(...) {
  p <- 1
  for (i in list(...)) {
    p <- i * p
  }
  return(if_else(p == 0, 0, 1, missing = 0))
}

