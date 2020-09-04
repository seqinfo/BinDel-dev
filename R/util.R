library(BSgenome.Hsapiens.UCSC.hg38)
library(DNAcopy)
library(dplyr)
library(GenomicAlignments)
library(purrr)
library(readr)
library(Rsamtools)
library(tidyr)

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


# https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
# https://stackoverflow.com/questions/61425521/challenge-replacing-soft-deprecated-funs-within-mutate-at
lagTransformation <- function(ds, n) {
  # this function creats lag transformation of dataframe
  # args:
  # ds : Dataset
  # n : number of lags
  lags <- seq(n)
  lag_names <-
    paste("lag", formatC(lags, width = nchar(max(lags)), flag = "0"), sep = "")
  lag_functions <-
    purrr::map(lags, ~ function(x)
      dplyr::lag(x, .x)) %>%
    setNames(lag_names)
  ds <- ds %>% mutate_at(vars(names(ds)), lag_functions)
  return(ds)
}