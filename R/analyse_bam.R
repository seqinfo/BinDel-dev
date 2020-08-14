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


# Remove extreme reference bins compared to other bins
sample_only <- bin_length_normalized %>%
  filter(sample == basename(bam_location))


# Remove bin with extreme values from the reference group.
cut_off <- 3
without_sample <- bin_length_normalized %>%
  filter(sample != basename(bam_location)) %>%
  group_by(chromosome, start, end) %>%
  mutate(mean = mean(gc_corrected)) %>%
  mutate(sd = sd(gc_corrected)) %>%
  filter((mean - cut_off * sd < gc_corrected) &
           (gc_corrected < mean + cut_off * sd)) %>%
  ungroup()

bin_length_normalized <- without_sample %>%
  bind_rows(sample_only)


# Z-score calculation with reference (bin wise)
results <- bin_length_normalized %>%
  group_by(chromosome, start, end) %>%
  mutate(z_score_ref = (gc_corrected - mean(gc_corrected)) / sd(gc_corrected)) %>%
  ungroup()


# Calculate local Z-score (sample wise)
results <- results %>%
  group_by(sample) %>%
  mutate(local_z_score = (gc_corrected - mean(gc_corrected)) / sd(gc_corrected)) %>%
  ungroup()


# Calculate local ZZ-score (sample wise)
results <- results %>%
  group_by(sample) %>%
  mutate(zz_score = (z_score_ref - mean(z_score_ref)) / sd(z_score_ref)) %>%
  ungroup()


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
         local_z_score,
         zz_score)


write_tsv(results, paste0("results.", basename(bam_location), ".tsv"))

