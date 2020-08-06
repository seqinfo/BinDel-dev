library(tidyverse)

results <- read_tsv("results.B869N.bqsr.bam.tsv")


coords_syndrome <- read_tsv("coordinates/as_pws.bed") %>% 
  mutate(region = "AS/PWS")

results <- results %>% 
  left_join(coords_syndrome, by =  c("chromosome", "start", "end"))


z.score.plot <- ggplot(results, aes(x=region, y=z_score_sample_within)) + 
  geom_boxplot(aes(color = region)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(-10, 10, 1))
z.score.plot
