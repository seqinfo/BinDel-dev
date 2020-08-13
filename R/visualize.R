library(tidyverse)
library(grid)

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
    read_tsv(regions_of_interest_location) %>%
      mutate(region = "Focus")
  )
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
  ggplot(results, aes(x = region, y = z_score_ref)) +
  geom_boxplot(aes(color = region)) +
  ggtitle("Z-score ref") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))


z.score.local.plot <-
  ggplot(results, aes(x = region, y = local_z_score)) +
  geom_boxplot(aes(color = region)) +
  ggtitle("Z-score local") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))


zz.score.plot <-
  ggplot(results, aes(x = region, y = zz_score)) +
  geom_boxplot(aes(color = region)) +
  ggtitle("ZZ-score") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))


# Scatter plots
###############

# Z-score ref scatter
z_score_ref.na_mean <- results %>%
  filter(is.na(region)) %>%
  select(z_score_ref) %>%
  summarise(mean(z_score_ref))

z_score_ref.syndrome_mean <- results %>%
  filter(!is.na(region)) %>%
  select(z_score_ref) %>%
  summarise(mean(z_score_ref))


z.score.ref.scatter <-
  ggplot(results, aes(x = start, y = z_score_ref)) +
  geom_point(aes(color = region), size = 0.5) +
  geom_hline(yintercept = z_score_ref.na_mean$`mean(z_score_ref)`, color =
               "blue") +
  geom_hline(yintercept = z_score_ref.syndrome_mean$`mean(z_score_ref)`, color =
               "red") +
  ggtitle("Z-score ref") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))

# Z-score local scatter
z.score.local.na_mean <- results %>%
  filter(is.na(region)) %>%
  select(local_z_score) %>%
  summarise(mean(local_z_score))

z.score.local.syndrome_mean <- results %>%
  filter(!is.na(region)) %>%
  select(local_z_score) %>%
  summarise(mean(local_z_score))

z.score.local.scatter <-
  ggplot(results, aes(x = start, y = local_z_score)) +
  geom_point(aes(color = region), size = 0.5) +
  geom_hline(yintercept = z.score.local.na_mean$`mean(local_z_score)`,
             color = "blue") +
  geom_hline(yintercept = z.score.local.syndrome_mean$`mean(local_z_score)`,
             color = "red") +
  ggtitle("Z-score local") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))

# ZZ-score  scatter
zz.score.na_mean <- results %>%
  filter(is.na(region)) %>%
  select(zz_score) %>%
  summarise(mean(zz_score))

zz.score.syndrome_mean <- results %>%
  filter(!is.na(region)) %>%
  select(zz_score) %>%
  summarise(mean(zz_score))

zz.score.scatter <-
  ggplot(results, aes(x = start, y = zz_score)) +
  geom_point(aes(color = region), size = 0.5) +
  geom_hline(yintercept = zz.score.na_mean$`mean(zz_score)`,
             color = "blue") +
  geom_hline(yintercept = zz.score.syndrome_mean$`mean(zz_score)`,
             color = "red") +
  ggtitle("ZZ score") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))

sample_name <- basename(results_location)

pdf(file = paste0(sample_name, ".pdf"), title = sample_name)
multiplot(
  z.score.local.plot,
  z.score.local.scatter,
  
  z.score.ref.plot,
  z.score.ref.scatter,
  
  zz.score.plot,
  zz.score.scatter,
  
  cols = 3
)
dev.off()
