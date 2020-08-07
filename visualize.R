library(tidyverse)
library(grid)

if (length(args) != 4) {
  stop("Please provide 'results.tsv', 'sample name', 'syndrome regions' and their 'name'", call.=FALSE)
}

results <- read_tsv(read_tsv(args[1])) 
sample_name = args[2]
syndrome_regions = args[3]
syndrome_region_name = args[4]

# results <- read_tsv("results.A749N.bam.tsv")# read_tsv(args[1])
# sample_name = "A749" #args[2]
# syndrome_regions = "coordinates/as_pws.bed"#args[3]
# syndrome_region_name = "AS/PWS"#args[4]

coords_syndrome <- read_tsv(syndrome_regions) %>% 
  mutate(region = syndrome_region_name)


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
    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
  }
}


theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        plot.title = element_text(size=10))


results <- results %>% 
  left_join(coords_syndrome, by =  c("chromosome", "start", "end"))

# Box plots
###############
z.score.ref.plot <- ggplot(results, aes(x=region, y=z_score_ref)) + 
  geom_boxplot(aes(color = region)) +
  ggtitle("Z-score ref") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))
  

z.score.local.plot <- ggplot(results, aes(x=region, y=local_z_score)) + 
  geom_boxplot(aes(color = region)) +
  ggtitle("Z-score local") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))
  

#chi.z.score.plot <- ggplot(results, aes(x=region, y=chi_z_score)) + 
#  geom_boxplot(aes(color = region)) +
#  ggtitle("Z-score chi") +
#  theme +
#  scale_y_continuous(limits = c(-10, 10, 1))

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
  

z.score.ref.scatter <- ggplot(results, aes(x=start, y=z_score_ref)) + 
  geom_point(aes(color = region), size = 0.5) +
  geom_hline(yintercept = z_score_ref.na_mean$`mean(z_score_ref)`, color="blue") +
  geom_hline(yintercept = z_score_ref.syndrome_mean$`mean(z_score_ref)`, color="red") +
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

z.score.local.scatter <- ggplot(results, aes(x=start, y=local_z_score)) + 
  geom_point(aes(color = region), size = 0.5) +
  geom_hline(yintercept = z.score.local.na_mean$`mean(local_z_score)`, color="blue") +
  geom_hline(yintercept = z.score.local.syndrome_mean$`mean(local_z_score)`, color="red") +
  ggtitle("Z-score local") +
  theme +
  scale_y_continuous(limits = c(-10, 10, 1))

# Z-score chi scatter
#z.score.chi.na_mean <- results %>% 
#  filter(is.na(region)) %>% 
#  select(chi_z_score) %>% 
#  summarise(mean(chi_z_score))

#z.score.chi.syndrome_mean <- results %>% 
#  filter(!is.na(region)) %>% 
#  select(chi_z_score) %>% 
#  summarise(mean(chi_z_score))

#chi.z.score.scatter <- ggplot(results, aes(x=start, y=chi_z_score)) + 
#  geom_point(aes(color = region), size = 0.5) +
#  geom_hline(yintercept = z.score.chi.na_mean$`mean(chi_z_score)`, color="blue") +
#  geom_hline(yintercept = z.score.chi.syndrome_mean$`mean(chi_z_score)`, color="red") +
#  ggtitle("Z-score chi") +
#  theme +
#  scale_y_continuous(limits = c(-10, 10, 1))

pdf(file = paste0(sample_name, ".pdf"), title = sample_name)
multiplot(z.score.local.plot, 
          z.score.local.scatter,          
          
          z.score.ref.plot,
          z.score.ref.scatter, 
          
          #chi.z.score.plot, 
          #chi.z.score.scatter, 
          cols = 2)
dev.off()

