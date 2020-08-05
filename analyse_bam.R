library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)

reference_location <- "reference.tsv"
bam_location <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/A749N.bam"
syndrome_bed <- "coordinates/as_pws.bed"
analyze_locations_bed <- "coordinates/chr15.bed"


sample <- readGAlignments(bam_location)
bam_name <- basename(bam_location)


syndrome_coordinates <- read_tsv(syndrome_bed) %>% 
  mutate(region = "syndrome")

analyze_coordinates <- read_tsv(analyze_locations_bed) %>% 
  mutate(region = "reference")

locations <- analyze_coordinates %>% 
  anti_join(syndrome_coordinates, by = c("chromosome", "start", "end")) %>% 
  bind_rows(syndrome_coordinates)


counts <- locations %>% 
  mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
  mutate(normalized_by_sample = reads/sum(reads)) %>% 
  mutate(sample = bam_name) %>% 
  mutate(type = "analysis")
  

reference <- read_tsv(reference_location) 
  
  
results <- reference %>% 
  bind_rows(counts) %>% 
  group_by(chromosome, start, end) %>% 
  mutate(z_score = (normalized_by_sample - mean(normalized_by_sample))/sd(normalized_by_sample)) %>% 
  ungroup() %>% 
  filter(!is.na(type)) %>% 
  select(-type)

