# Stouffer's_Z-score_method
# sum((results$z_score_ref))/sqrt(nrow(results))

# ref_size <- reference %>%
#   select(sample) %>%
#   distinct() %>%
#   nrow(.) + 1 # add also the sample under analysis for chi squared calculation

# total_bins <- reference %>%
#   select(chromosome, start, end) %>%
#   group_by(chromosome, start, end) %>%
#   distinct() %>%
#   nrow(.)

# Degrees of freedom
# df <- ref_size - 1

# chi_squared_z_score <- function(dataframe, degrees_of_freedom, total_bins, ref_size){
#  # Chi squared variation reduction: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5431782/
#
#  # Sum of bin i over all the samples
#  dataframe %>% group_by(chromosome, start, end) %>%
#    mutate(numerator = sum(gc_corrected)) %>%
#    ungroup() %>%
#    mutate(numerator = numerator / (ref_size * total_bins)) %>%
#    group_by(sample) %>%
#
#    # Average of sample normalized count.
#    mutate(denominator = mean(gc_corrected)) %>%
#    ungroup() %>%
#
#    # Multiply reads by normalization factor
#    mutate(on = gc_corrected * numerator / denominator) %>%
#    group_by(chromosome, start, end) %>%
#    mutate(expected_on = mean(on)) %>% # Expected normalized count value
#    ungroup() %>%
#
#    # Calculate chi-squared
#    mutate(dif = (expected_on - on)^2) %>%
#    mutate(chi_squared = ((expected_on - on)^2) / expected_on) %>%
#
#    # chi correction:
#    # mutate(chi_corrected_count = gc_corrected * chi_squared / degrees_of_freedom) %>%
#
#    # transform to a standard normal distribution N(0, 1)
#    ungroup() %>%
#    mutate(chi_z_score = (chi_squared - degrees_of_freedom) / sqrt(2 * degrees_of_freedom)) %>%
#    ungroup()
#
#  return(dataframe)
#}



# chi.z.score.plot <- ggplot(results, aes(x=region, y=chi_z_score)) +
#  geom_boxplot(aes(color = region)) +
#  ggtitle("Z-score chi") +
#  theme +
#  scale_y_continuous(limits = c(-10, 10, 1))


# Z-score chi scatter
# z.score.chi.na_mean <- results %>%
#  filter(is.na(region)) %>%
#  select(chi_z_score) %>%
#  summarise(mean(chi_z_score))

# z.score.chi.syndrome_mean <- results %>%
#  filter(!is.na(region)) %>%
#  select(chi_z_score) %>%
#  summarise(mean(chi_z_score))

# chi.z.score.scatter <- ggplot(results, aes(x=start, y=chi_z_score)) +
#  geom_point(aes(color = region), size = 0.5) +
#  geom_hline(yintercept = z.score.chi.na_mean$`mean(chi_z_score)`, color="blue") +
#  geom_hline(yintercept = z.score.chi.syndrome_mean$`mean(chi_z_score)`, color="red") +
#  ggtitle("Z-score chi") +
#  theme +
#  scale_y_continuous(limits = c(-10, 10, 1))


# lags <- seq(9)
# lead_names <- paste("lead", formatC(lags, width = nchar(max(lags)), flag = "0"), sep = "")
# lead_functions <- purrr::map(lags, ~ function(x)
# dplyr::lead(x, .x)) %>%
#  setNames(lead_names)


# test <- results %>%
#  filter(focus == "AS/PWS") %>%  #TODO: rm filter
#  select(sample, start, focus, z_score_ref) %>%
#  group_by(sample, focus) %>%
#  arrange(start, .by_group = TRUE) %>%
#  mutate_at(vars(z_score_ref), lead_functions) %>%
#  ungroup() %>%
#  pivot_longer(cols = all_of(lead_names), names_to = "lead_name", values_to = "lead_value")


# test3 <- test %>%
#  group_by(sample, start, focus) %>%
#  filter(!is.na(lead_value)) %>%
#  filter(lead_value < 0) %>%
#  summarise(bins_sum = sum(lead_value))

# temp <- test3 %>%
#  filter(sample != basename(bam_location))

# mean_bin = mean(temp$bins_sum)
# sd_bin = sd(temp$bins_sum)


# Z-score calculation with reference (bin wise)
# test4 <- test3 %>%
# mutate(z_score = (bins_sum - mean_bin) / sd_bin)


# library(ggplot2)
# to.plot <- k %>%
  # filter(sample != basename(bam_location)) %>%
  # filter(focus == "AS/PWS") %>%
  # select(start, ratio, sample)


# theme <- theme_bw() +
#  theme(
#    legend.position = "none",
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    axis.title.x = element_blank(),
#    axis.title.y = element_blank(),
#    plot.title = element_text(size = 10)
#  )

# line3 <-
#  geom_hline(
#    yintercept = 2,
#    linetype = "dashed",
#    color = "red",
#    size = 2
#  )

# line33 <-
#  geom_hline(
#    yintercept = -2,
#    linetype = "dashed",
#    color = "red",
#    size = 2
#  )

# ggplot(to.plot, aes(x = start, y = ratio)) +
#  geom_line() +
#  #line3 +
#  #line33 +
#  theme

