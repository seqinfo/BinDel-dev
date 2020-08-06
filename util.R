library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)
require(BSgenome.Hsapiens.UCSC.hg38)

count <- function(bed, sample, bam_name){
  reads <- bed %>% 
    mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
    mutate(sample = bam_name) %>% 
    group_by(chromosome, start, end) %>% # Calculate GC%
    mutate(gc = letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                       GRanges(chromosome, 
                                               IRanges(start=start, end=end), strand="+", as.character=T)),
                                "GC", 
                                as.prob = T)) %>% 
    ungroup() %>% # GC correct DOI: 10.1038/s41598-017-02031-5
    mutate(gc = round(gc, 1)) %>% 
    group_by(gc) %>% 
    mutate(avg_gc_bin = mean(reads)) %>% 
    ungroup() %>% 
    mutate(factor = avg_gc_bin / mean(reads)) %>% 
    mutate(gc_corrected = round(reads * factor, 0)) %>% 
    mutate(normalized_by_sample = gc_corrected / sum(gc_corrected)) %>% 
    select(-avg_gc_bin, -factor)
  return(reads)
}

