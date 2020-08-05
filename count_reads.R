library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)
require(BSgenome.Hsapiens.UCSC.hg38)

args = commandArgs(trailingOnly=TRUE)

# Check that at least three arguments are supplied.
if (length(args) != 3) {
  stop("Please provide .bam, .bed and output locations.", call.=FALSE)
  }

bam_location <- args[1]
bed_location <- args[2]
out_location <- args[3]

bam_name <- basename(bam_location)

sample <- readGAlignments(bam_location)
bed <- read_tsv(bed_location)


# LOESS: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3130771/

reads <- bed %>% 
  mutate(reads = assay(summarizeOverlaps(makeGRangesFromDataFrame(.), sample, mode="IntersectionStrict"))) %>% 
  mutate(normalized_by_sample = reads/sum(reads)) %>% 
  mutate(sample = bam_name) %>% 
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
  

write_tsv(path = paste0(out_location, "count.", basename(bam_name), ".tsv"), x = reads)
