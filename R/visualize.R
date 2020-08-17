library(tidyverse)
library(grid)


start_x <- -15
end_x <- 15

start_y <- -15
end_y <- 15


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
  left_join((
    read_tsv(regions_of_interest_location))
  , by =  c("chromosome", "start", "end"))


# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., cols = 1) {
  plots <- c(list(...))
  
  numPlots = length(plots)
  
  layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                   ncol = cols,
                   nrow = ceiling(numPlots / cols))
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  
  for (i in 1:numPlots) {
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    print(plots[[i]],
          vp = viewport(
            layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col
          ))
  }
}


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
    plot.title = element_text(size = 10)
  )


# Box plots
###############
z.score.ref.plot <-
  ggplot(results, aes(x = focus, y = z_score_ref)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("Z-score ref") +
  theme +
  scale_y_continuous(limits = c(start_y, end_y, 1))


z.score.local.plot <-
  ggplot(results, aes(x = focus, y = local_z_score)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("Z-score local") +
  theme +
  scale_y_continuous(limits = c(start_y, end_y, 1))


zz.score.plot <-
  ggplot(results, aes(x = focus, y = zz_score)) +
  geom_boxplot(aes(color = focus)) +
  ggtitle("ZZ-score") +
  theme +
  scale_y_continuous(limits = c(start_y, end_y, 1))


# Density plots
###############

# Z-score
z.score.ref.dens <-
  ggplot(results, aes(x = z_score_ref)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("Z-score ref") +
  theme +
  scale_x_continuous(limits = c(start_x, end_x, 1))

# Z-score local scatter
z.score.local.dens <-
  ggplot(results, aes(x = local_z_score)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("Z-score local") +
  theme +
  scale_x_continuous(limits = c(start_x, end_x, 1))

# ZZ-score
zz.score.dens <-
  ggplot(results, aes(x = zz_score)) +
  geom_density(aes(color = focus), size = 0.5) +
  ggtitle("ZZ score") +
  theme +
  scale_x_continuous(limits = c(start_x, end_x, 1))

sample_name <- basename(results_location)


pdf(file = paste0(sample_name, ".pdf"), title = sample_name)
multiplot(
  z.score.local.plot,
  z.score.local.dens,
  
  z.score.ref.plot,
  z.score.ref.dens,
  
  zz.score.plot,
  zz.score.dens,
  
  cols = 3
)
dev.off()
