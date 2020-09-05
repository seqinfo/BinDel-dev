source(here::here("R/util.R"))

args = commandArgs(trailingOnly = TRUE)


if (length(args) != 3) {
  stop(
    "Please provide (1) analyzable BAM (.bam), (2) analyzable regions (.bed) and (3) reference file (.tsv).",
    call. = FALSE
  )
}


bam_location <-  args[1]
bed_location <- args[2]
reference_location <- args[3]


binned_reads <- bin_counts(bam_location, bed_location) # Bin BAM
queried_gc <- find_gc(bed_location) # Find GC for the locations
reference <-
  read_tsv(reference_location) # Reference samples to be used to calculate z-scores


# Merge: BAM + reference + GC
merged <- reference %>%
  mutate(sample = paste0("ref.set.", sample)) %>% # reference samples names must not overlap with the analyzable sample.
  bind_rows(binned_reads) %>%
  left_join(queried_gc)


# GC-correct (sample wise) (PMID: 28500333 and PMID: 20454671)
gc_corrected <- merged %>%
  group_by(sample, gc) %>%
  mutate(avg_reads_gc_interval = mean(reads)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(weights = mean(reads) / avg_reads_gc_interval) %>%
  mutate(gc_corrected = reads * weights) %>%
  ungroup()


# Normalize by sample (important to make samples comparable)
gc_corrected <- gc_corrected %>%
  group_by(sample) %>%
  mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
  ungroup()


# Normalize by bin length
bin_length_normalized <- gc_corrected %>%
  mutate(gc_corrected = gc_corrected / (end - start))


# Calculate reference group statistics
sample_only <- bin_length_normalized %>%
  filter(sample == basename(bam_location))


# Calculate ref set (each) bin SD and mean
without_sample <- bin_length_normalized %>%
  filter(sample != basename(bam_location)) %>%
  ungroup()

# Calculate each reference bin i mean
ref_bins <- without_sample %>%
  group_by(chromosome, start) %>%
  summarise(expected = mean(gc_corrected)) %>%
  ungroup()


reference_bin_info <- without_sample %>%
  group_by(focus, start, end) %>%
  mutate(mean_ref_bin = mean(gc_corrected)) %>%
  mutate(mean_ref_sd = sd(gc_corrected)) %>%
  ungroup() %>%
  select(focus, start, end, mean_ref_bin, mean_ref_sd) %>%
  distinct()


bin_length_normalized <- without_sample %>%
  bind_rows(sample_only) %>%
  left_join(reference_bin_info)


# Z-score calculation with reference (bin wise)
results <- bin_length_normalized %>%
  mutate(z_score_ref = (gc_corrected - mean_ref_bin) / mean_ref_sd)


# Calculate expected value and actual value ratio
results <- ref_bins %>%
  right_join(results, by = c("chromosome", "start")) %>%
  mutate(ratio = log(gc_corrected / expected, base = 2))


# Clean the output
results <- results %>%
  filter(sample == basename(bam_location)) %>%  # Keep in the output only the analyzable sample
  select(chromosome,
         start,
         end,
         reads,
         gc,
         sample,
         z_score_ref,
         ratio)


# Calculate aberrations with circular binary segmentation
# PMID: 15475419
CNA.object <- CNA(
  genomdat = results$ratio,
  chrom = results$chromosome,
  maploc = results$start,
  data.type = "logratio",
  sampleid = basename(bam_location)
)

smoothed <- smooth.CNA(
  CNA.object,
  smooth.region = 10,
  outlier.SD.scale = 4,
  smooth.SD.scale = 2,
  trim = 0.025
)

segments <- segment(smoothed, verbose = 3, nperm = 1000)


write_tsv(results, paste0(basename(bam_location), ".results.tsv"))
write_tsv(segments$output, paste0(basename(bam_location), ".segments.tsv"))
