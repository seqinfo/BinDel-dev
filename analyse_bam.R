library(Rsamtools)
library(tidyverse)
source("util.R")

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and reference locations.", call.=FALSE)
}

# TODO: rm
# bam_location <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/B869N.bqsr.bam"# args[1]
# analyze_locations_bed <- "coordinates/chr15.bed" #args[2]
# reference_location <- "reference.tsv" #args[3]

bam_location <-  args[1]
analyze_locations_bed <- args[2]
reference_location <- args[3]

bam_name <- basename(bam_location)


sample <- readGAlignments(bam_location)
bed <- read_tsv(analyze_locations_bed) 
reads <- count(bed, sample, bam_name)
reference <- read_tsv(reference_location) 

ref_size <- reference %>% 
  select(sample) %>% 
  distinct() %>% 
  nrow(.) + 1 # add also the sample under analysis for chi squared calculation

total_bins <- reference %>% 
  select(chromosome, start, end) %>% 
  group_by(chromosome, start, end) %>% 
  distinct() %>% 
  nrow(.) 
  
# Degrees of freedom
df <- ref_size  - 1


results <- reference %>% 
  # Z-score calculation with reference.
  bind_rows(reads) %>% 
  group_by(chromosome, start, end) %>% 
  mutate(z_score_ref = (normalized_by_sample - mean(normalized_by_sample)) / sd(normalized_by_sample)) %>%
  ungroup() %>% 
  # Chi squared variation reduction: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431782/
  group_by(chromosome, start, end) %>%  
  mutate(numerator = sum(gc_corrected)) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(denominator = sum(gc_corrected)) %>% 
  ungroup() %>% 
  mutate(numerator = numerator / (ref_size * total_bins)) %>% 
  mutate(denominator = denominator / total_bins) %>% 
  mutate(factor = numerator / denominator) %>% 
  mutate(on = factor * gc_corrected) %>% # Normalized counts
  group_by(chromosome, start, end) %>% 
  mutate(expected_on = mean(on)) %>% # Expected normalized count value
  ungroup() %>% 
  group_by(chromosome, start, end) %>%  
  # Calculate chi-squared
  mutate(chi_squared = ((expected_on - on)^2) / expected_on) %>% 
  # transform to a standard normal distribution N(0, 1)
  mutate(chi_z_score = (chi_squared - df) / sqrt(2 * df)) %>% 
  ungroup() %>%
  filter(sample == bam_name) %>% 
  # Calculate local Z-score
  mutate(local_z_score = (normalized_by_sample - mean(normalized_by_sample)) / sd(normalized_by_sample)) %>% 
  # Filter columns
  select(chromosome, start, end, sample, z_score_ref, local_z_score, chi_z_score)

write_tsv(results, paste0("results.", bam_name, ".tsv"))

