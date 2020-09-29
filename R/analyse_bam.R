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


sample_name = basename(bam_location)
binned_reads <- bin_counts(bam_location, bed_location) # Bin BAM
queried_gc <- find_gc(bed_location) # Find GC for the locations
reference <-
  readr::read_tsv(reference_location) # Reference samples to be used to calculate z-scores


# Merge: BAM + reference + GC
merged <- reference %>%
  dplyr::mutate(sample = paste0("ref.set.", sample)) %>% # reference samples names must not overlap with the analyzable sample.
  dplyr::bind_rows(binned_reads) %>%
  dplyr::left_join(queried_gc)


# GC-correct (sample wise) (PMID: 28500333 and PMID: 20454671)
gc_corrected <- merged %>%
  dplyr::group_by(sample, gc) %>%
  dplyr::mutate(avg_reads_gc_interval = mean(reads)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(weights = mean(reads) / avg_reads_gc_interval) %>%
  dplyr::mutate(gc_corrected = reads * weights) %>%
  dplyr::filter(!is.na(gc_corrected)) %>%
  dplyr::ungroup()


# Normalize by sample (important to make samples comparable)
gc_corrected <- gc_corrected %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
  dplyr::ungroup()


# Normalize by bin length
bin_length_normalized <- gc_corrected %>%
  dplyr::mutate(gc_corrected = gc_corrected / (end - start))


# Calculate reference group statistics
sample_only <- bin_length_normalized %>%
  dplyr::filter(sample == sample_name) %>%
  dplyr::mutate(reference = FALSE)


# Calculate ref set (each) bin SD and mean
without_sample <- bin_length_normalized %>%
  dplyr::filter(sample != sample_name) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(reference = TRUE)


# Calculate each reference bin i mean
ref_bins <- without_sample %>%
  dplyr::group_by(chromosome, start) %>%
  dplyr::summarise(expected = mean(gc_corrected)) %>%
  dplyr::ungroup()


reference_bin_info <- without_sample %>%
  dplyr::group_by(focus, start, end) %>%
  dplyr::mutate(mean_ref_bin = mean(gc_corrected)) %>%
  dplyr::mutate(mean_ref_sd = sd(gc_corrected)) %>%
  dplyr::ungroup() %>%
  dplyr::select(focus, start, end, mean_ref_bin, mean_ref_sd) %>%
  dplyr::distinct()


bin_length_normalized <- without_sample %>%
  dplyr::bind_rows(sample_only) %>%
  dplyr::left_join(reference_bin_info)


# Z-score calculation with reference (bin wise)
results <- bin_length_normalized %>%
  dplyr::mutate(z_score_ref = (gc_corrected - mean_ref_bin) / mean_ref_sd)


# Calculate expected value and actual value ratio
results <- ref_bins %>%
  dplyr::right_join(results, by = c("chromosome", "start")) %>%
  dplyr::mutate(ratio = log(gc_corrected / expected, base = 2))


# HMM
results <- results %>%
  filter(chromosome != focus) %>% 
  tidyr::drop_na() %>%
  dplyr::arrange(desc(reference, sample, focus, start), .by_group = TRUE)


# Mann–Whitney U test
results <- results %>%
  dplyr::group_by(chromosome, start) %>%
  dplyr::mutate(Mann_Whitney = wilcox.test(gc_corrected ~ reference, exact = FALSE)$p.value) %>%
  dplyr::ungroup()


results <- results %>%
  dplyr::group_by(focus) %>%
  tidyr::nest() %>%
  dplyr::mutate(HMM = purrr::map(data, function(df)
    depmixS4::posterior(
      depmixS4::fit(
        depmixS4::depmix(
          list(ratio ~ 1, Mann_Whitney ~ 1),
          family = list(gaussian(), gaussian()),
          nstates = 2,
          data = df
        )
        ,
        verbose = 1
      )
    ))) %>%
  tidyr::unnest(cols = c(data, HMM)) %>%
  dplyr::mutate(HMM = state) %>%
  dplyr::ungroup()


MW_count <- results %>%
  dplyr::mutate(MW = round(-log10(Mann_Whitney), 3), sign = sign(ratio)) %>%
  dplyr::group_by(sample, focus, sign, HMM) %>%
  dplyr::summarise(sum = sum(MW) / n(), count = n()) %>%
  dplyr::ungroup()


MW_stats <- MW_count %>%
  filter(sample != sample_name) %>%
  dplyr::group_by(focus, sign, HMM) %>%
  summarise(mean_sum = mean(sum), sd_sum = sd(sum)) %>%
  ungroup() %>%
  right_join(MW_count) %>%
  filter(sample == sample_name) %>%
  mutate((sum - mean_sum) / sd_sum)


# Output sample info only
results <- results %>%
  dplyr::filter(sample == sample_name)  # Keep in the output only the analyzable sample


# Clean the output
results <- results %>%
  dplyr::select(
    chromosome,
    start,
    end,
    reads,
    gc_corrected,
    gc,
    sample,
    z_score_ref,
    ratio,
    Mann_Whitney,
    HMM
  )


# Calculate aberrations with circular binary segmentation
# PMID: 15475419
CNA.object <- DNAcopy::CNA(
  genomdat = -log10(results$Mann_Whitney),
  chrom = results$chromosome,
  maploc = results$start,
  data.type = "logratio",
  sampleid = sample_name
)


smoothed <- DNAcopy::smooth.CNA(
  CNA.object,
  smooth.region = 10,
  outlier.SD.scale = 4,
  smooth.SD.scale = 2,
  trim = 0.025
)


segments <-
  DNAcopy::segment(smoothed, verbose = 1, nperm = 10000)$output %>%
  dplyr::filter(loc.start != loc.end)


readr::write_tsv(results, paste0(sample_name, ".results.tsv"))
readr::write_tsv(MW_stats, paste0(sample_name, ".metrics.tsv"))
readr::write_tsv(segments, paste0(sample_name, ".segments.tsv"))
