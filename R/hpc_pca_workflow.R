library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)

find_gc <- function(bed) {
  reads <- bed %>%
    mutate(gc = letterFrequency(
      getSeq(BSgenome.Hsapiens.UCSC.hg38,
             GRanges(chr, IRanges(
               start = start, end = end
             ))),
      "GC",
      as.prob = T
    )) %>%
    mutate(gc = round(gc, 1))
}

ref_location <-
  "/gpfs/hpc/samba/CCHT/service_lab/science/CNV/SyndromeDetector/workflows/reference/ref_outputs/reference.tsv"
ref_samples <- "reference_669.txt"
out_name <- "output_pca_669_pca160.txt"
nComp <- 160
reference <- read_tsv(ref_samples)

samples <-
  read_tsv(ref_location) %>%
  left_join(find_gc(
    select(.data = ., chr, start, end, focus) %>%
      distinct(chr, start, end, focus)
  )) %>%
  mutate(sample = str_remove(sample, ".bam")) %>%
  # GC correct
  group_by(sample, gc) %>%
  mutate(avg_reads_gc_interval = mean(reads)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(weights = mean(reads) / avg_reads_gc_interval) %>%
  mutate(gc_corrected = reads * weights) %>%
  filter(!is.na(gc_corrected)) %>%
  ungroup() %>%
  # Sample read count correct
  group_by(sample) %>%
  mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
  ungroup() %>%
  # Sample bin length correct:
  mutate(gc_corrected = gc_corrected / (end - start)) %>%
  group_by(sample, chr) %>%
  mutate(reads_chr = sum(reads)) %>%
  ungroup() %>%
  # Optimize memory
  select(-end,-gc,-reads,-avg_reads_gc_interval,-weights) %>%
  # mutate(sample = str_replace(sample, "_S[0-9]+$", "")) %>%
  left_join(reference) %>%
  arrange(desc(reference))


# Pivot wide for PCA normalization
wider <- samples %>%
  select(focus, start, sample, reference, gc_corrected) %>%
  pivot_wider(
    names_from = c(focus, start),
    id_cols = c(sample, reference),
    values_from = gc_corrected,
    names_sep = ":"
  )

# https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
# Train PCA
ref <- wider %>%
  filter(reference) %>%
  select(-reference, -sample)

mu <- colMeans(ref, na.rm = T)
refPca <- prcomp(ref)


Xhat <- refPca$x[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
Xhat <- scale(Xhat, center = -mu, scale = FALSE)

# Use trained PCA on other samples
pred <- wider %>%
  filter(is.na(reference)) %>%
  select(-reference, -sample)

Yhat <-
  predict(refPca, pred)[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
Yhat <- scale(Yhat, center = -mu, scale = FALSE)

# Actual PCA normalization and conversion back to long:
normalized <-
  bind_rows(as.data.frame(as.matrix(ref) / as.matrix(Xhat)),
            as.data.frame(as.matrix(pred) / as.matrix(Yhat))) %>%
  pivot_longer(
    names_sep = ":",
    names_to = c("focus", "start"),
    cols = everything(),
    values_to = "gc_corrected"
  )

normalized$sample <- samples$sample
normalized$reference <- samples$reference
normalized$chr <- samples$chr

# Clean memory footprint
rm(wider)
rm(pred)
rm(ref)
rm(Yhat)
rm(Xhat)

# Calculate each reference bin i mean and filter out high variance and low mean
reference <- normalized %>%
  filter(reference) %>%
  group_by(chr, start) %>%
  summarise(mean_ref_bin = mean(gc_corrected),
            mean_ref_sd = sd(gc_corrected)) %>%
  ungroup()

filtered <- reference # %>%
#group_by(chr) %>%
#filter(mean_ref_bin > mean(mean_ref_bin)) %>%
#filter(mean_ref_sd < mean(mean_ref_sd))


samples <- normalized %>%
  right_join(filtered) %>%
  mutate(z_score = (gc_corrected - mean_ref_bin) / mean_ref_sd) %>%
  mutate(under_median = as.integer(gc_corrected < mean_ref_bin)) %>%
  mutate(over_median = as.integer(gc_corrected >= mean_ref_bin)) %>%
  filter(!is.na(z_score)) %>%
  mutate(reference = !is.na(reference)) %>%
  group_by(sample, focus) %>%
  summarise(
    z_score_neg = sum(case_when(z_score > 0 ~ NA_real_,
                                z_score <= 0 ~ z_score), na.rm = T) / sqrt(n()),
    
    z_score_pos = sum(case_when(z_score > 0 ~ z_score,
                                z_score <= 0 ~ NA_real_), na.rm = T) / sqrt(n()),
    z_score_PPDX = sum(z_score) / sqrt(n()),
    under_median = sum(under_median),
    over_median = sum(over_median),
    mid_range = (max(z_score) + min(z_score)) / 2
  ) %>%
  mutate(
    z_score_PPDX_norm =  (z_score_PPDX + 1) / (under_median + 2),
    z_score_PPDX_norm2 =  (z_score_PPDX + 1) / (over_median + 2),
    dup_score = (z_score_pos + 1) / (-1 * z_score_neg + 2)
  ) %>%
  select(-z_score_neg,-z_score_pos)



write_tsv(samples, out_name)
