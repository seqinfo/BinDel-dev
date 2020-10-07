source(here::here("R/util.R"))
hmm_script_location <- here::here("R/hmm.py")
stats_script_location <- here::here("R/stats.py")

args <- commandArgs(trailingOnly = TRUE)


if (length(args) != 2) {
  stop("Please provide (1) analyzable BAM (.bam) and (2) reference file (.tsv).",
       call. = FALSE)
}


bam_location <- args[1]
reference_location <- args[2]

reference <-
  readr::read_tsv(reference_location) # Reference samples to be used to calculate z-scores

bed <- reference %>%
  select(chromosome, start, end, focus) %>%
  distinct(chromosome, start, end, focus)


sample_name <- basename(bam_location)
binned_reads <- bin_counts(bam_location, bed) # Bin BAM
queried_gc <- find_gc(bed) # Find GC for the locations


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


# Normalize by sample
gc_corrected <- gc_corrected %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
  dplyr::ungroup() %>% # And by bin length:
  dplyr::mutate(gc_corrected = gc_corrected / (end - start))


# Calculate reference group statistics
sample_only <- gc_corrected %>%
  dplyr::filter(sample == sample_name) %>%
  dplyr::mutate(reference = FALSE)


# Calculate ref set (each) bin SD and mean
without_sample <- gc_corrected %>%
  dplyr::filter(sample != sample_name) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(reference = TRUE)


# Calculate each reference bin i mean
ref_bins <- without_sample %>%
  dplyr::group_by(chromosome, start) %>%
  dplyr::summarise(expected = mean(gc_corrected),
                   sd = sd(gc_corrected)) %>%
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
  dplyr::right_join(reference_bin_info) %>% 
  filter(sd <= 1.5e-08)


# Z-score calculation with reference (bin wise)
results <- bin_length_normalized %>%
  dplyr::mutate(z_score_ref = (gc_corrected - mean_ref_bin) / mean_ref_sd)


# Calculate expected value and actual value ratio
results <- ref_bins %>%
  dplyr::right_join(results, by = c("chromosome", "start")) %>%
  dplyr::mutate(ratio = log(gc_corrected / expected, base = 2))


# Mann–Whitney U test
results <- results %>%
  dplyr::group_by(chromosome, start) %>%
  dplyr::mutate(Mann_Whitney = -log10(wilcox.test(gc_corrected ~ reference, exact = FALSE)$p.value)) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na() %>%
  dplyr::arrange(desc(sample, focus, start))


# Continue with sample only
results <- results %>%
  dplyr::filter(sample == sample_name)


# HMM
temp_tsv <- paste0(sample_name, ".temp")
hmm_temp <- paste0(temp_tsv, ".hmm")


readr::write_tsv(results, temp_tsv)
command <-
  paste("python", hmm_script_location, "-i", temp_tsv, "-o", hmm_temp)


if (system(command) == 0) {
  results <- readr::read_tsv(hmm_temp)
  
} else {
  stop("hmm.py did not finish with expected exit code!")
}


# Output, plots and additional statistics


# Add state names for plotting
results <- results %>%
  dplyr::mutate(HMM = paste0("S", HMM))


# Statistics:
stats <- results %>%
  tidyr::pivot_wider(id_cols = c(focus, HMM, start)) %>%
  dplyr::arrange(desc(focus, start))


temp_tsv <- paste0(sample_name, ".temp")
stats_temp <- paste0(temp_tsv, ".stats")
readr::write_tsv(stats, temp_tsv)


command <-
  paste("python",
        stats_script_location,
        "-i",
        temp_tsv,
        "-o",
        stats_temp)

if (system(command) == 0) {
  stats <- readr::read_tsv(stats_temp) %>%
    dplyr::group_by(focus, HMM, length) %>%
    dplyr::summarise(n = n())
  
  readr::write_tsv(stats, paste0(sample_name, ".stats.tsv"))
  
} else {
  stop("stats.py did not finish with expected exit code!")
}


# Plots
pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10,
  height = 20
)

ggplot(results, aes(x = start, y = ratio)) +
  geom_line(size = 0.001,
            alpha = 0.5,
            color = "grey") +
  geom_point(aes(x = start, y = ratio, color = HMM),
             size = 1,
             alpha = 1) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(facets = ~ focus,
             scales = "free",
             ncol = 3) +
  scale_color_manual(values = c("S0" = "red", "S1" = "grey", "S2" = "purple")) +
  main_theme +
  get_line(1) +
  get_line(0) +
  get_line(-1)

ggplot(results, aes(x = start, y = sd)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ focus, scales = "free") +
  main_theme +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6))

ggplot(stats %>% filter(HMM == "S0"), aes(x = length, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ focus + HMM, scales = "free") +
  scale_x_continuous(breaks = pretty_breaks()) +
  main_theme

ggplot(stats %>% filter(HMM == "S1"), aes(x = length, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ focus + HMM, scales = "free") +
  scale_x_continuous(breaks = pretty_breaks()) +
  main_theme

ggplot(stats %>% filter(HMM == "S2"), aes(x = length, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ focus + HMM, scales = "free") +
  scale_x_continuous(breaks = pretty_breaks()) +
  main_theme

dev.off()


# Clean the output
results <- results %>%
  dplyr::select(
    chromosome,
    focus,
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

readr::write_tsv(results, paste0(sample_name, ".results.tsv"))
