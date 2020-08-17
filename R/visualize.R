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


results <- read_tsv(results_location) %>%
  left_join((read_tsv(regions_of_interest_location))
            , by =  c("chromosome", "start", "end"))


theme <- theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 10)
  )


# Box plots
###############
zero_line <-
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "red")

z.score.ref.plot <-
  ggplot(results, aes(x = focus, y = z_score_ref)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("Z-score ref") +
  theme +
  zero_line


z.score.local.plot <-
  ggplot(results, aes(x = focus, y = local_z_score)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("Z-score local") +
  theme +
  zero_line


zz.score.plot <-
  ggplot(results, aes(x = focus, y = zz_score)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("ZZ-score") +
  theme +
  zero_line


# Density plots
###############

# Z-score
z.score.ref.dens <-
  ggplot(results, aes(x = z_score_ref)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("Z-score ref") +
  theme

# Z-score local scatter
z.score.local.dens <-
  ggplot(results, aes(x = local_z_score)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("Z-score local") +
  theme

# ZZ-score
zz.score.dens <-
  ggplot(results, aes(x = zz_score)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("ZZ score") +
  theme

sample_name <- basename(results_location)

pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10
)
grid.arrange(
  z.score.local.plot,
  z.score.local.dens,
  
  z.score.ref.plot,
  z.score.ref.dens,
  
  zz.score.plot,
  zz.score.dens,
  nrow = 3
)

dev.off()
