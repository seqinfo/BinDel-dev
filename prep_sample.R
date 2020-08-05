library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)


bam_location <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/60569250_S11.bam"
bed_location <- "coordinates/chr15.bed"
out_location <- "test/"

bam_name <- basename(bam_location)

sample <- readGAlignments(bam_location)
bed <- read_tsv(bed_location)

reads <- bed %>% 
  mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
  mutate(normalized_by_sample = reads/sum(reads)) %>% 
  mutate(sample = bam_name)

write_tsv(path = paste0(out_location, "count.",basename(bam_name),".tsv"), x = reads)

