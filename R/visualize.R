library(tidyverse)
library(gridExtra)


args = commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Please provide (1) results file in .tsv and (2) regions of interest in .bed format.",
    call. = FALSE
  )
}

results_location <- args[1]
regions_of_interest_location = args[2]


# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text = l)
}


results <- read_tsv(results_location) %>%
  left_join((read_tsv(regions_of_interest_location))
            , by =  c("chromosome", "start", "end")) %>%
  separate(
    chromosome,
    remove = F,
    sep = "chr",
    into = c("temp", "chr_number"),
    convert = T
  ) %>%
  select(-temp) %>%
  arrange(chr_number)


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
    color = "red",
    size = 1
  )



box.plot.chr <-
  ggplot(results %>% filter(chromosome == focus),
         aes(fct_reorder(focus, chr_number), z_score_ref)) +
  geom_jitter(
    shape = 16,
    position = position_jitter(0.2),
    aes(color = chromosome),
    alpha = 0.3
  ) +
  geom_boxplot(aes(fill = chromosome), alpha = 0.5) +
  ggtitle("Z-score ref chromosomes") +
  theme +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 0.5,
    vjust = 0.5,
  )) +
  zero_line +
  scale_y_continuous(n.breaks = 8, limits = c(-5, 5))


box.plot.target <-
  ggplot(results %>% filter(chromosome != focus),
         aes(fct_reorder(focus, chr_number), z_score_ref)) +
  geom_jitter(
    shape = 16,
    position = position_jitter(0.2),
    aes(color = chromosome),
    alpha = 0.3
  ) +
  geom_boxplot(aes(fill = chromosome), alpha = 0.5) +
  ggtitle("Z-score ref targets") +
  theme +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 0.5,
    vjust = 0.5,
  )) +
  zero_line


overall <- ggplot(results, aes(x = start, y = z_score_ref)) +
  geom_point(aes(color = focus), size = 0.001, alpha = 0.5) +
  geom_line(aes(color = focus), size = 0.001, alpha = 0.5) +
  geom_boxplot(
    aes(color = focus),
    alpha = 0.5,
    outlier.shape = NA,
    position = "identity"
  ) +
  scale_x_continuous(n.breaks = 10, labels = fancy_scientific) +
  facet_wrap(facets = vars(chr_number),
             scales = "free",
             ncol = 2) +
  ggtitle("Combined") +
  theme +
  theme(axis.text.x = element_text(size = 6)) +
  zero_line 



targets <-
  ggplot(results %>% filter(chromosome != focus),
         aes(x = start, y = z_score_ref)) +
  geom_point(aes(color = chromosome), size = 0.001, alpha = 0.5) +
  geom_line(aes(color = chromosome), size = 0.001, alpha = 0.5) +
  geom_boxplot(
    aes(color = chromosome),
    alpha = 0.5,
    outlier.shape = NA,
    position = "identity"
  ) +
  scale_x_continuous(n.breaks = 10, labels = fancy_scientific) +
  facet_wrap(facets = vars(focus),
             scales = "free",
             ncol = 2) +
  ggtitle("Targets") +
  theme +
  theme(axis.text.x = element_text(size = 6)) +
  zero_line



sample_name <- basename(results_location)


pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10
)


grid.arrange(box.plot.chr, box.plot.target, ncol = 1)
overall
targets

dev.off()

