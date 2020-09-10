library(tidyverse)
library(gridExtra)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fuzzyjoin)

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Please provide (1) results file and (2) segments file in .tsv and (3) regions of interest in .bed format.",
    call. = FALSE
  )
}

results_location <- args[1]
segments_location <- args[2]
regions_of_interest_location <- args[3]


# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text = l)
}


results <- read_tsv(results_location) %>%
  dplyr::left_join((read_tsv(regions_of_interest_location))
                   , by =  c("chromosome", "start", "end")) %>%
  separate(
    chromosome,
    remove = F,
    sep = "chr",
    into = c("temp", "chr_number"),
    convert = T
  ) %>%
  dplyr::select(-temp) %>%
  dplyr::arrange(chr_number) %>%
  dplyr::filter(!is.na(z_score_ref))


focuses <- results %>%
  dplyr::select(chromosome, start, focus)


segments <- read_tsv(segments_location) %>%
  dplyr::mutate(chromosome = chrom) %>%
  dplyr::mutate(start = loc.start) %>%
  dplyr::left_join(focuses) %>%
  separate(
    chromosome,
    remove = F,
    sep = "chr",
    into = c("temp", "chr_number"),
    convert = T
  ) %>%
  dplyr::select(-temp) %>%
  dplyr::arrange(chr_number)


theme <- theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 10)
  )


zero_line <-
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey",
    size = 1
  )


one_line <-
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey",
    size = 1
  )


minus_line <-
  geom_hline(
    yintercept = -1,
    linetype = "dashed",
    color = "grey",
    size = 1
  )


box.plot.chr <-
  ggplot(results %>% filter(chromosome == focus),
         aes(fct_reorder(focus, chr_number), ratio)) +
  geom_jitter(
    shape = 16,
    position = position_jitter(0.2),
    aes(color = chromosome),
    alpha = 0.3
  ) +
  geom_boxplot(aes(fill = chromosome), alpha = 0.5) +
  theme +
  zero_line +
  minus_line +
  one_line


box.plot.target <-
  ggplot(results %>% filter(chromosome != focus),
         aes(fct_reorder(focus, chr_number), ratio)) +
  geom_jitter(
    shape = 16,
    position = position_jitter(0.2),
    aes(color = chromosome),
    alpha = 0.3
  ) +
  geom_boxplot(aes(fill = chromosome), alpha = 0.5) +
  theme +
  zero_line +
  minus_line +
  one_line


theme <- theme_bw() +
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


overall <- ggplot(results, aes(x = start, y = ratio)) +
  geom_point(aes(color = ifelse(
    focus != chromosome, ifelse(ratio < -0.5, 'red', "grey"), 'grey'
  )), size = 1, alpha = 1) +
  scale_x_continuous(n.breaks = 10, labels = fancy_scientific) +
  facet_wrap(facets = vars(chr_number),
             scales = "free",
             ncol = 2) +
  scale_color_identity() +
  geom_segment(
    data = segments %>%
      filter(abs(seg.mean) > 0.1),
    aes(
      x = loc.start,
      y = seg.mean,
      xend = loc.end,
      yend = seg.mean
    ),
    size = 1
  ) +
  theme +
  zero_line +
  minus_line +
  one_line


target_results <- results %>%
  dplyr::filter(chromosome != focus)

temp <- target_results %>%
  dplyr::select(chromosome) %>%
  dplyr::distinct()

genes <-
  as.data.frame(transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)) %>%
  dplyr::mutate(chromosome = seqnames) %>%
  dplyr::right_join(temp, by = "chromosome") %>%
  fuzzy_left_join(
    target_results,
    .,
    by = c(
      "start" = "start",
      "end" = "end",
      "chromosome" = "chromosome"
    ),
    match_fun = list(`>`, `<`, `==`)
  )



filtered_genes <- genes %>%
  dplyr::select(chromosome.y, start.y, end.y, tx_name, focus) %>%
  dplyr::filter(!is.na(tx_name)) %>%
  mutate(chromosome = chromosome.y,
         start = start.y,
         end = end.y) %>%
  dplyr::select(chromosome, start, end, tx_name, focus) %>%
  distinct()


targets <-
  ggplot(target_results, aes(x = start, y = ratio)) +
  zero_line +
  minus_line +
  one_line +
  geom_point(aes(color = ifelse(ratio < -0.5, 'red', "grey")), size = 1, alpha = 1) +
  geom_line(aes(color = "grey"), size = 0.001, alpha = 0.5) +
  scale_x_continuous(n.breaks = 10, labels = fancy_scientific) +
  geom_segment(
    data = filtered_genes,
    aes(
      x = start,
      y = -1,
      xend = end,
      yend = -1,
      color = "blue"
    ),
    size = 1
  ) +
  facet_wrap(facets = vars(focus),
             scales = "free",
             ncol = 2) +
  scale_color_identity() +
  geom_segment(
    data = segments %>%
      filter(chromosome != focus),
    aes(
      x = loc.start,
      y = seg.mean,
      xend = loc.end,
      yend = seg.mean
    ),
    size = 1
  ) +
  geom_label(
    data = segments %>% 
      filter(chromosome != focus),
    aes(
      x = loc.start + ((loc.end - loc.start) / 2),
      y = 0.05,
      label = paste0(num.mark, "/", seg.mean)
    )
  ) +
  theme


sample_name <- basename(results_location)

write_tsv(genes, paste0(sample_name, ".genes.tsv"))

pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10
)


grid.arrange(box.plot.chr, box.plot.target, ncol = 1)
overall
targets

dev.off()
