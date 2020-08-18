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


zero_line <-
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "red")


create_box_plot <- function(data, x, y, title) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_jitter(
      shape = 16,
      position = position_jitter(0.2),
      aes(color = x),
      alpha = 0.3
    ) +
    geom_boxplot(aes(color = x), alpha = 0.5) +
    ggtitle(title) +
    theme +
    zero_line
  return(p)
}


create_density <- function(data, x, title) {
  p <- ggplot(results, aes(x = x)) +
    geom_density(aes(color = focus), size = 0.5) +
    ggtitle(title) +
    theme
  
  return(p)
}


create_scatter_plot <- function(data, x, y, title) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = focus), alpha = 0.5) +
    facet_grid(cols = vars(focus), scales = "free") +
    ggtitle(title) +
    theme +
    zero_line
  return(p)
}


# Box plots
box.local <-
  create_box_plot(results, results$focus, results$local_z_score, "Z-score local")
box.ref <-
  create_box_plot(results, results$focus, results$z_score_ref, "Z-score ref")
box.zz <-
  create_box_plot(results, results$focus, results$zz_score, "ZZ-score")


# Density plots
dens.local <-
  create_density(results, results$local_z_score, "Z-score local")
dens.ref <-
  create_density(results, results$z_score_ref, "Z-score ref")
dens.zz <-
  create_density(results, results$zz_score, "ZZ score")


# Scatter plots
scatter.local <-
  create_scatter_plot(results, results$start, results$local_z_score, "Z-score local")
scatter.ref <-
  create_scatter_plot(results,
                      results$start,
                      results$z_score_ref,
                      "Z-score reference")
scatter.zz <-
  create_scatter_plot(results, results$start, results$zz_score, "ZZ-score reference")


# Scatter plots |abs|
scatter.local.abs <-
  create_scatter_plot(results,
                      results$start,
                      abs(results$local_z_score),
                      "|Z|-score local")
scatter.ref.abs <-
  create_scatter_plot(results,
                      results$start,
                      abs(results$z_score_ref),
                      "|Z|-score reference")
scatter.zz.abs <-
  create_scatter_plot(results,
                      results$start,
                      abs(results$zz_score),
                      "|ZZ|-score reference")



sample_name <- basename(results_location)

pdf(
  file = paste0(sample_name, ".pdf"),
  title = sample_name,
  width = 10
)
grid.arrange(box.local,
             dens.local,
             
             box.ref,
             dens.ref,
             
             box.zz,
             dens.zz,
             nrow = 3)

grid.arrange(scatter.local,
             scatter.ref,
             scatter.zz,
             nrow = 3)

grid.arrange(scatter.local.abs,
             scatter.ref.abs,
             scatter.zz.abs,
             nrow = 3)

dev.off()
