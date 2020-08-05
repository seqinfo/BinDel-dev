library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and reference locations.", call.=FALSE)
}

bam_location <- args[1]
analyze_locations_bed <- args[2]
reference_location <- args[3]


sample <- readGAlignments(bam_location)
bam_name <- basename(bam_location)


# Coordinates for use
locations <- read_tsv(analyze_locations_bed) %>% 
  mutate(region = "reference")

# Count reads from the sample under analysis.
counts <- locations %>% 
  mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
  mutate(normalized_by_sample = reads/sum(reads)) %>% 
  mutate(sample = bam_name) %>% 
  mutate(type = "analysis")
  

reference <- read_tsv(reference_location) 

ref_size <- reference %>% 
  select(sample) %>% 
  distinct() %>% 
  nrow(.)

total_bins <- reference %>% 
  nrow(.) / ref_size
  
# Degrees of freedom
df <- ref_size - 1

# We will add sample under analysis.
df <- df + 1

# Z-score calculation with reference.
results <- reference %>% 
  bind_rows(counts) %>% 
  group_by(chromosome, start, end) %>% 
  mutate(z_score_over_ref = (normalized_by_sample - mean(normalized_by_sample)) / sd(normalized_by_sample)) %>% 
  ungroup() %>% # Chi squared: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431782/
  group_by(chromosome, start, end) %>%  # the observed read counts are first normalized by multiplying them with a normalization factor. 
  mutate(factor = sum(reads) / (ref_size * total_bins)) %>% # Factor numerator (sum of specific bin divided by ...)
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(factor = factor / (sum(reads) / total_bins)) %>% # Numerator / (the denominator: sum of sample reads divided by...)
  ungroup() %>% 
  mutate(on = factor * reads) %>% # Normalized counts
  select(-factor) %>%  # Factor is not needed
  group_by(chromosome, start, end) %>% 
  mutate(expected_on = mean(on)) %>% # Expected value is the average of each bin
  ungroup() %>% 
  mutate(chi_squared = ((expected_on - on)^2) / expected_on) %>% # Calculate chi-squared
  mutate(chi_squared_standard_normal_z_score = (chi_squared - df) / sqrt(2 * df)) %>% # transform to a standard normal distribution N(0, 1) %>% 
  filter(!is.na(type)) %>% 
  select(-type) %>% 
  mutate(z_score_sample_within = (normalized_by_sample - mean(normalized_by_sample)) / sd(normalized_by_sample)) %>% 
  write_tsv(paste0("results.", bam_name,".tsv"))

