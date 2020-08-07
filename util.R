library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg38)

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
    mutate(gc = round(gc, 1))
  return(reads)
}

# ref_size <- reference %>% 
#   select(sample) %>% 
#   distinct() %>% 
#   nrow(.) + 1 # add also the sample under analysis for chi squared calculation

# total_bins <- reference %>% 
#   select(chromosome, start, end) %>% 
#   group_by(chromosome, start, end) %>% 
#   distinct() %>% 
#   nrow(.) 

# Degrees of freedom
# df <- ref_size - 1

chi_squared_z_score <- function(dataframe, degrees_of_freedom, total_bins, ref_size){
  # Chi squared variation reduction: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431782/
  
  # Sum of bin i over all the samples
  dataframe %>% group_by(chromosome, start, end) %>%  
    mutate(numerator = sum(gc_corrected)) %>% 
    ungroup() %>% 
    mutate(numerator = numerator / (ref_size * total_bins)) %>% 
    group_by(sample) %>% 
    
    # Average of sample normalized count.
    mutate(denominator = mean(gc_corrected)) %>%
    ungroup() %>%
    
    # Multiply reads by normalization factor
    mutate(on = gc_corrected * numerator / denominator) %>%
    group_by(chromosome, start, end) %>% 
    mutate(expected_on = mean(on)) %>% # Expected normalized count value
    ungroup() %>% 
    
    # Calculate chi-squared
    mutate(dif = (expected_on - on)^2) %>% 
    mutate(chi_squared = ((expected_on - on)^2) / expected_on) %>% 
    
    # chi correction:
    # mutate(chi_corrected_count = gc_corrected * chi_squared / degrees_of_freedom) %>% 
    
    # transform to a standard normal distribution N(0, 1)
    ungroup() %>% 
    mutate(chi_z_score = (chi_squared - degrees_of_freedom) / sqrt(2 * degrees_of_freedom)) %>% 
    ungroup()
  
  return(dataframe)
}

