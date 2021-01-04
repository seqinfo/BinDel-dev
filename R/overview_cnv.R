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

main_theme <- theme_bw() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5
    )
  )

cut_off_ff <- 7

fetal_fractions <- read_tsv("seqff.niptify1.tsv") %>%
  select(sample, SeqFF1) %>%
  mutate(SeqFF1 = SeqFF1 * 100)

reference <- readr::read_tsv("reference_est.txt")

samples <- readr::read_tsv("reference.tsv") %>%
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
  select(-end, -gc, -reads, -avg_reads_gc_interval, -weights) %>%
  left_join(fetal_fractions) %>%
  mutate(sample = str_replace(sample, "_S[0-9]+$", "")) %>%
  left_join(reference) %>%
  # SeqFF filter
  filter(SeqFF1 > cut_off_ff | reference)


# Calculate each reference bin i mean and filter out high variance and low mean
reference <- samples %>%
  filter(reference) %>%
  group_by(chr, start) %>%
  summarise(mean_ref_bin = mean(gc_corrected),
            mean_ref_sd = sd(gc_corrected)) %>%
  ungroup()

filtered <- reference %>%
  group_by(chr) %>%
  filter(mean_ref_bin > mean(mean_ref_bin)) %>%
  filter(mean_ref_sd < mean(mean_ref_sd))

samples <- samples %>%
  right_join(filtered) %>%
  mutate(z_score = (gc_corrected - mean_ref_bin) / mean_ref_sd) %>%
  mutate(ratio = log(gc_corrected / mean_ref_bin, base = 2)) %>%
  mutate(under_median = as.integer(gc_corrected < mean_ref_bin)) %>%
  filter(!is.na(z_score)) %>%
  mutate(reference = !is.na(reference)) %>%
  ungroup() %>%
  group_by(sample, focus) %>%
  mutate(sd = sd(gc_corrected)) %>%
  ungroup() %>%
  group_by(sample, focus, reads_chr, sd) %>%
  summarise(
    z_score_PPDX = sum(z_score) / sqrt(n()),
    under_median = sum(under_median)
  ) %>%
  mutate(z_score_PPDX_norm =  z_score_PPDX / under_median) %>%
  left_join(fetal_fractions)


write_tsv(samples, "output.tsv")
