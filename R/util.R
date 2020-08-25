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


#' Calculate consecutive negative Z-score bin stats (sum and number of occurrences of length 1...9)
#'  using windowed_sum and windowed_count up to bin length of 9
calculate_neg_bin_stat <- function(data){
  group_z_aggerated <- data %>%
    group_by(sample, focus, start) %>%
    mutate(score = abs(min(z_score_ref, 0, na.rm = T))) %>%
    ungroup() %>%
    group_by(sample, focus) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(
      n1 = lead(x = score, n = 1, default = 0),
      n2 = lead(x = score, n = 2, default = 0),
      n3 = lead(x = score, n = 3, default = 0),
      n4 = lead(x = score, n = 4, default = 0),
      n5 = lead(x = score, n = 5, default = 0),
      n6 = lead(x = score, n = 6, default = 0),
      n7 = lead(x = score, n = 7, default = 0),
      n8 = lead(x = score, n = 8, default = 0),
      n9 = lead(x = score, n = 9, default = 0)
    ) %>%
    ungroup() %>% 
    group_by(sample, focus) %>%
    summarise(
      sum1 = sum(windowed_sum(n1)),
      count1 = sum(windowed_count(n1)),
      
      sum2 = sum(
        windowed_sum(n1, n2),
        windowed_sum(n2, n3),
        windowed_sum(n3, n4),
        windowed_sum(n4, n5),
        windowed_sum(n5, n6),
        windowed_sum(n6, n7),
        windowed_sum(n7, n8),
        windowed_sum(n8, n9)
      ),
      count2 = sum(
        windowed_count(n1, n2),
        windowed_count(n2, n3),
        windowed_count(n3, n4),
        windowed_count(n4, n5),
        windowed_count(n5, n6),
        windowed_count(n6, n7),
        windowed_count(n7, n8),
        windowed_count(n8, n9)
      ),
      
      sum3 = sum(
        windowed_sum(n1, n2, n3),
        windowed_sum(n2, n3, n4),
        windowed_sum(n3, n4, n5),
        windowed_sum(n4, n5, n6),
        windowed_sum(n5, n6, n7),
        windowed_sum(n6, n7, n8),
        windowed_sum(n7, n8, n9)
      ),
      count3 = sum(
        windowed_count(n1, n2, n3),
        windowed_count(n2, n3, n4),
        windowed_count(n3, n4, n5),
        windowed_count(n4, n5, n6),
        windowed_count(n5, n6, n7),
        windowed_count(n6, n7, n8),
        windowed_count(n7, n8, n9)
      ),
      
      sum4 = sum(
        windowed_sum(n1, n2, n3, n4),
        windowed_sum(n2, n3, n4, n5),
        windowed_sum(n3, n4, n5, n6),
        windowed_sum(n4, n5, n6, n7),
        windowed_sum(n5, n6, n7, n8),
        windowed_sum(n6, n7, n8, n9)
      ),
      count4 = sum(
        windowed_count(n1, n2, n3, n4),
        windowed_count(n2, n3, n4, n5),
        windowed_count(n3, n4, n5, n6),
        windowed_count(n4, n5, n6, n7),
        windowed_count(n5, n6, n7, n8),
        windowed_count(n6, n7, n8, n9)
      ),
      sum5 = sum(
        windowed_sum(n1, n2, n3, n4, n5),
        windowed_sum(n2, n3, n4, n5, n6),
        windowed_sum(n3, n4, n5, n6, n7),
        windowed_sum(n4, n5, n6, n7, n8),
        windowed_sum(n5, n6, n7, n8, n9)
      ),
      count5 = sum(
        windowed_count(n1, n2, n3, n4, n5),
        windowed_count(n2, n3, n4, n5, n6),
        windowed_count(n3, n4, n5, n6, n7),
        windowed_count(n4, n5, n6, n7, n8),
        windowed_count(n5, n6, n7, n8, n9)
      ),
      sum6 = sum(
        windowed_sum(n1, n2, n3, n4, n5, n6),
        windowed_sum(n2, n3, n4, n5, n6, n7),
        windowed_sum(n3, n4, n5, n6, n7, n8),
        windowed_sum(n4, n5, n6, n7, n8, n9)
      ),
      count6 = sum(
        windowed_count(n1, n2, n3, n4, n5, n6),
        windowed_count(n2, n3, n4, n5, n6, n7),
        windowed_count(n3, n4, n5, n6, n7, n8),
        windowed_count(n4, n5, n6, n7, n8, n9)
      ),
      sum7 = sum(
        windowed_sum(n1, n2, n3, n4, n5, n6, n7),
        windowed_sum(n2, n3, n4, n5, n6, n7, n8),
        windowed_sum(n3, n4, n5, n6, n7, n8, n9)
      ),
      count7 = sum(
        windowed_count(n1, n2, n3, n4, n5, n6, n7),
        windowed_count(n2, n3, n4, n5, n6, n7, n8),
        windowed_count(n3, n4, n5, n6, n7, n8, n9)
      ),
      sum8 = sum(
        windowed_sum(n1, n2, n3, n4, n5, n6, n7, n8),
        windowed_sum(n2, n3, n4, n5, n6, n7, n8, n9)
      ),
      count8 = sum(
        windowed_count(n1, n2, n3, n4, n5, n6, n7, n8),
        windowed_count(n2, n3, n4, n5, n6, n7, n8, n9)
      ),
      sum9 = sum(windowed_sum(n1, n2, n3, n4, n5, n6, n7, n8, n9)),
      count9 = sum(windowed_count(n1, n2, n3, n4, n5, n6, n7, n8, n9))
    ) %>%
    ungroup()
  
  return(group_z_aggerated)
}
