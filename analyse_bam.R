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
  mutate(type = "analysis") %>% 
  group_by(chromosome, start, end) %>% # Calculate GC%
  mutate(gc = letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                     GRanges(chromosome, 
                                             IRanges(start=start, end=end), strand="+", as.character=T)),
                              "GC", 
                              as.prob = T)) %>% 
  ungroup() %>% # LOESS GC correct PMC3130771
  mutate(P = predict(loess(gc ~ reads, .))) %>% 
  mutate(M = median(reads)) %>% 
  mutate(factor = M/P) %>% 
  mutate(gc_corrected = reads * factor)
  

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
  mutate(z_score_over_gc = (gc_corrected - mean(gc_corrected)) / sd(gc_corrected)) %>% 
  ungroup() %>% 
  # Chi squared: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431782/
  group_by(chromosome, start, end) %>%  
  # the observed read counts are first normalized by multiplying them with a normalization factor. 
  # Factor numerator (sum of specific bin divided by ...)
  mutate(factor = sum(gc_corrected) / (ref_size * total_bins)) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  # Numerator / (the denominator: sum of sample reads divided by...)
  mutate(factor = factor / (sum(gc_corrected) / total_bins)) %>% 
  ungroup() %>% 
  mutate(on = factor * gc_corrected) %>% # Normalized counts
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

