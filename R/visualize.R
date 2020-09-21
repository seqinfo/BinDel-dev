library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(scales)

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




get_line <- function(y) {
  return(geom_hline(
    yintercept = y,
    linetype = "dashed",
    color = "grey",
    size = 1
  ))
}


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
  get_line(1) +
  get_line(0) +
  get_line(-1)


overall <-
  ggplot(results, aes(
    x = start,
    y = ratio,
    color = focus != chromosome
  )) +
  geom_point(size = 1, alpha = 1) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(facets = vars(chr_number),
             scales = "free",
             ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  theme +
  get_line(1) +
  get_line(0) +
  get_line(-1)



target_results <- results %>%
  dplyr::filter(chromosome != focus)

temp <- target_results %>%
  dplyr::select(chromosome) %>%
  dplyr::distinct()

genes <-
  as.data.frame(transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)) %>%
  dplyr::mutate(chromosome = seqnames) %>%
  dplyr::right_join(temp, by = "chromosome")


targets <-
  ggplot(target_results, aes(x = start, y = ratio)) +
  geom_point(size = 1, alpha = 1) +
  geom_line(size = 0.001, alpha = 0.5) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(facets = vars(focus),
             scales = "free",
             ncol = 2) +
  scale_color_identity() +
  theme +
  get_line(1) +
  get_line(0) +
  get_line(-1)



pvalues <-
  ggplot(results, aes(
    x = start,
    y = -log10(Mann_Whitney),
    color = focus != chromosome
  )) +
  geom_point(size = 1, alpha = 1) +
  geom_segment(
    data = segments %>%
      filter(focus != chromosome),
    aes(
      x = loc.start,
      y = -log10(seg.mean),
      xend = loc.end,
      yend = -log10(seg.mean),
      
    )
  ) +
  geom_segment(
    data = genes,
    aes(
      x = start,
      y = 0,
      xend = end,
      yend = 0,
      
    ),
    color = "red",
    size = 1
  ) +
  facet_wrap(facets = vars(chromosome),
             scales = "free",
             ncol = 2) +
  theme +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_color_brewer(palette = "Paired")


sample_name <- basename(results_location)

write_tsv(genes, paste0(sample_name, ".genes.tsv"))

pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10
)


box.plot.chr
overall
targets
pvalues

dev.off()
