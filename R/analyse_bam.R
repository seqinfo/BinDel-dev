library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(GenomicAlignments)
library(Rsamtools)

args <- commandArgs(trailingOnly = TRUE)


if (length(args) != 2) {
  stop("Please provide (1) analyzable BAM (.bam) and (2) reference file (.tsv).",
       call. = FALSE)
}


bam_location <- args[1]
sample_name <- basename(bam_location)

reference_location <- args[2]


# Useful for 45,X detection (and aneuploidy detection), but can hide CNVs.
# Not recommended for micro-deletion detection.
bin_filter_on <- FALSE

# Has important effect for micro-deletions detection
use_pca <- TRUE
# Output probabilities of not being reference.
output_mahalanobis <- TRUE
# Use GC-correct (recommended to use)
do_gc_correct <- TRUE
# If PCA normalization is used, how many components?
# Tested to have effect between 20-160.
nComp <- 160

# Save plot
create_plot <- TRUE
# Include reference in the output?
include_reference <- TRUE

# Optimize RAM usage
clean_env <- TRUE

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

bin_counts <- function(bam_location, bed) {
  bam <- readGAlignments(bam_location, param =  ScanBamParam(
    flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
    what = c("pos")
  ))
  
  binned_counts <- bed %>%
    mutate(reads = assay(
      summarizeOverlaps(makeGRangesFromDataFrame(.), bam, mode = "IntersectionStrict")
    )) %>%
    mutate(sample = basename(bam_location))
  
  return(binned_counts)
}

# Reference samples to be used to calculate z-scores
reference <-
  read_tsv(reference_location) %>%
  mutate(reference = TRUE)


# Bin BAM under investigation
binned_reads <- bin_counts(
  bam_location,
  reference %>%
    select(chr, start, end, focus) %>%
    distinct(chr, start, end, focus)
) %>%
  mutate(reference = FALSE)


# Merge: BAM + reference
samples <- reference %>%
  # Note, reference samples names must not overlap with the analyzable sample.
  dplyr::bind_rows(binned_reads)


if (clean_env) {
  rm(binned_reads)
}

# GC-correct (sample wise) (PMID: 28500333 and PMID: 20454671)
if (do_gc_correct) {
  # Find GC% of the genome (HG38)
  samples <- samples %>%
    left_join(find_gc(
      select(.data = ., chr, start, end, focus) %>%
        distinct(chr, start, end, focus)
    )) %>%
    # Do GC correct
    group_by(sample, gc) %>%
    mutate(avg_reads_gc_interval = mean(reads)) %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(weights = mean(reads) / avg_reads_gc_interval) %>%
    mutate(gc_corrected = reads * weights) %>%
    filter(!is.na(gc_corrected)) %>%
    ungroup()
  
} else{
  samples <- samples %>%
    mutate(gc_corrected = reads) %>%
    filter(!is.na(gc_corrected))
}

samples <- samples %>%
  # Sample read count correct
  group_by(sample) %>%
  mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
  ungroup() %>%
  # Sample bin length correct
  mutate(gc_corrected = gc_corrected / (end - start)) %>%
  # Optimize memory
  select(chr, focus, start, sample, reference, gc_corrected)


if (use_pca) {
  # For PCA sort ()
  samples <- samples %>%
    arrange(reference)
  
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
    select(-reference,-sample)
  
  mu <- colMeans(ref, na.rm = T)
  refPca <- prcomp(ref)
  
  
  Xhat <- refPca$x[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
  Xhat <- scale(Xhat, center = -mu, scale = FALSE)
  
  # Use trained PCA on other samples
  pred <- wider %>%
    filter(!reference) %>%
    select(-reference,-sample)
  
  Yhat <-
    predict(refPca, pred)[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
  Yhat <- scale(Yhat, center = -mu, scale = FALSE)
  
  # Actual PCA normalization and conversion back to long:
  normalized <-
    bind_rows(as.data.frame(as.matrix(pred) / as.matrix(Yhat)),
              as.data.frame(as.matrix(ref) / as.matrix(Xhat))) %>%
    pivot_longer(
      names_sep = ":",
      names_to = c("focus", "start"),
      cols = everything(),
      values_to = "gc_corrected"
    )
  
  normalized$sample <- samples$sample
  normalized$reference <- samples$reference
  normalized$chr <- samples$chr
  samples <- normalized
  
  # Clean memory footprint
  if (clean_env) {
    rm(wider)
    rm(pred)
    rm(ref)
    rm(Yhat)
    rm(Xhat)
    rm(normalized)
  }
  
}

# Calculate each reference bin i mean and filter out high variance and low mean
reference <- samples %>%
  filter(reference) %>%
  group_by(chr, start) %>%
  summarise(mean_ref_bin = mean(gc_corrected),
            mean_ref_sd = sd(gc_corrected)) %>%
  ungroup()


if (bin_filter_on) {
  filtered <- reference %>%
    group_by(chr) %>%
    filter(mean_ref_bin > mean(mean_ref_bin)) %>%
    filter(mean_ref_sd < mean(mean_ref_sd)) %>%
    ungroup()
} else{
  filtered <- reference
}


samples <- samples %>%
  right_join(filtered) %>%
  mutate(
    z_score = (gc_corrected - mean_ref_bin) / mean_ref_sd,
    over_median = as.integer(gc_corrected >= mean_ref_bin)
  ) %>%
  filter(!is.na(z_score)) %>%
  group_by(sample, focus, reference) %>%
  summarise(z_score_PPDX = sum(z_score) / sqrt(n()),
            over_median = sum(over_median)) %>%
  mutate(z_score_PPDX_norm =  (z_score_PPDX + 1) / (over_median + 2)) %>%
  group_by(reference, focus) %>%
  mutate(mean_x = mean(z_score_PPDX_norm),
         mean_y = mean(z_score_PPDX)) %>%
  ungroup()

if (clean_env) {
  rm(filtered)
}

if (output_mahalanobis) {
  samples <- samples %>%
    group_by(focus) %>%
    group_split() %>%
    map_dfr( ~ {
      cov <-
        cov(.x %>%
              filter(reference) %>%
              select(z_score_PPDX_norm, z_score_PPDX))
      
      center <- .x %>%
        filter(reference) %>%
        select(mean_x, mean_y) %>%
        distinct()
      
      distances <- mahalanobis(
        .x %>% select(z_score_PPDX_norm, z_score_PPDX),
        c(center$mean_x, center$mean_y),
        cov
      )
      
      bind_cols(
        sample = .x$sample,
        focus = .x$focus,
        reference = .x$reference,
        z_score_PPDX_norm = .x$z_score_PPDX_norm,
        z_score_PPDX = .x$z_score_PPDX,
        p = -log10(pchisq(
          distances, df = 2, lower.tail = FALSE
        ) + 1e-100)
      )
    })
  
}

main_theme <- theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


if (create_plot) {
  last_plot <- ggplot(samples, aes(x = focus,
                                   y = p)) +
    geom_point(
      data = samples %>% filter(reference),
      aes(alpha = 0.5),
      color = "grey",
      shape = 16
    ) +
    geom_point(
      data = samples %>% filter(!reference),
      color = "black",
      shape = 17
    ) +
    ylab("Probability of interest") +
    ggtitle(sample_name) +
    main_theme +
    scale_color_grey()
  
  
  ggsave(paste0(sample_name, ".png"), last_plot)
  
}

if (include_reference) {
  write_tsv(samples, paste0(sample_name, ".tsv"))
} else{
  write_tsv(samples %>% filter(!reference), paste0(sample_name, ".tsv"))
}
