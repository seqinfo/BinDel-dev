library(Rsamtools)
library(tidyverse)
source("util.R")

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and reference locations.", call.=FALSE)
}

# TODO: rm, used only for local testing.
# bam_location <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/60569250_S11.bam"# args[1]
# analyze_locations_bed <- "coordinates/chr15.bed" #args[2]
# reference_location <- "reference.tsv" #args[3]

bam_location <-  args[1]
analyze_locations_bed <- args[2]
reference_location <- args[3]

bam_name <- basename(bam_location)

# BAM file
sample <- readGAlignments(bam_location)

# Bins to analyze
bed <- read_tsv(analyze_locations_bed) 

# Analyzable sample read counts for bins
reads <- count(bed, sample, bam_name)

# Reference sample to be used to calculate z-scores
reference <- read_tsv(reference_location) 


results <- reference %>% 
  # Add sample under analysis
  bind_rows(reads) %>% 
  
  # GC correct (sample wise) (PMID: 28500333 and PMID: 20454671)
  group_by(sample, gc) %>% 
  mutate(avg_reads_gc_interval = mean(reads)) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(weights = mean(reads) / avg_reads_gc_interval) %>% 
  mutate(gc_corrected = reads * weights) %>% 
  
  # Normalize by sample (important to make samples comparable)
  mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>% 
  ungroup() %>% 
  
  # Z-score calculation with reference (bin wise)
  group_by(chromosome, start, end) %>% 
  mutate(z_score_ref = (gc_corrected - mean(gc_corrected)) / sd(gc_corrected)) %>%
  ungroup() %>% 
  
  # Calculate local Z-score (sample wise)
  group_by(sample) %>% 
  mutate(local_z_score = (gc_corrected - mean(gc_corrected)) / sd(gc_corrected)) %>% 
  ungroup() %>% 
  
  # Keep in the output only the analyzable sample
  filter(sample == bam_name) %>% 
  
  # Filter columns
  select(chromosome, start, end, reads, gc, sample, z_score_ref, local_z_score)

write_tsv(results, paste0("results.", bam_name, ".tsv"))

