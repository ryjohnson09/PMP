# negctrl_abun_prev.R
# Ryan Johnson
# 14 June 2020
# Compare the abundance and prevalence of neg-ctrl ASVs between various groups

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)
library(RColorBrewer)
library(cowplot)

# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds")
ps2_tibble <- psmelt(ps2)

# Select neg-ctrl data ----------------------------------
bad_asvs <- ps2_tibble %>% 
  filter(Group == "Neg-Ctrl") %>% 
  filter(Abundance > 0)

bad_asvs_relabun <- bad_asvs %>% 
  group_by(Sample) %>% 
  mutate(rel_abun_neg_ctrl = (Abundance / sum(Abundance)) * 100) %>% 
  ungroup() %>% 
  select(OTU, rel_abun_neg_ctrl) %>% 
  group_by(OTU) %>% 
  summarise(average_relabun_negctrl = mean(rel_abun_neg_ctrl))

# Select same ASVs in non-neg-ctrl data -----------------------------
bad_asvs_good_data <- ps2_tibble %>% 
  filter(Group != "Neg-Ctrl") %>% 
  group_by(Sample) %>% 
  mutate(rel_abun_good = (Abundance / sum(Abundance)) * 100) %>% 
  filter(OTU %in% bad_asvs$OTU) %>% 
  ungroup()

bad_asvs_good_data_relabun <- bad_asvs_good_data %>% 
  select(OTU, rel_abun_good) %>% 
  mutate(rel_abun_good = if_else(is.nan(rel_abun_good), 0, rel_abun_good)) %>% 
  group_by(OTU) %>% 
  summarise(average_relabun_good = mean(rel_abun_good))

# Relative Abundance Plot ----------------------
rel_abun_good_bad <- bad_asvs_relabun %>% 
  full_join(bad_asvs_good_data_relabun)

rel_abun_good_bad_plot <- rel_abun_good_bad %>% 
  ggplot(aes(x = average_relabun_negctrl, y = average_relabun_good)) +
  geom_point(size = 2.5, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Negative Controls",
       y = "Real Samples",
       title = "ASV Average Relative Abundance") +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )

rel_abun_good_bad_plot

#########################################################################

# Bad ASVs prevalence ----------------------------------
bad_asvs_prev <- bad_asvs %>% 
  select(OTU, Sample, Abundance) %>% 
  mutate(sample_presence = if_else(Abundance > 0, 1, 0)) %>% 
  group_by(OTU) %>% 
  summarise(bad_prevalence = sum(sample_presence) / length(unique(bad_asvs$Sample)) * 100)

# Bad ASVs good data prevalence ----------------------------------
bad_asvs_good_data_prev <- bad_asvs_good_data %>% 
  select(OTU, Sample, Abundance) %>% 
  mutate(sample_presence = if_else(Abundance > 0, 1, 0)) %>% 
  group_by(OTU) %>% 
  summarise(good_prevalence = sum(sample_presence) / length(unique(bad_asvs_good_data$Sample)) * 100)

# Plot Prevalence -----------------------
prev_good_bad <- bad_asvs_prev %>% 
  full_join(bad_asvs_good_data_prev)

prev_good_bad_plot <- prev_good_bad %>% 
  ggplot(aes(x = bad_prevalence, y = good_prevalence)) +
  geom_jitter(width = 0, height = 0, size = 2.5, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Negative Controls",
       y = "Real Samples",
       title = "ASV Prevalence") +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )

prev_good_bad_plot

#########################################################################

# Bad ASVs abundance per sample ----------------------------------
bad_asvs_abund <- bad_asvs %>% 
  select(OTU, Sample, Abundance) %>% 
  group_by(OTU) %>% 
  summarise(bad_abundance = sum(Abundance) / length(unique(bad_asvs$Sample)))

# Bad ASVs good data abundance ----------------------------------
bad_asvs_good_data_abund <- bad_asvs_good_data %>% 
  select(OTU, Sample, Abundance) %>% 
  group_by(OTU) %>% 
  summarise(good_abundance = sum(Abundance) / length(unique(bad_asvs_good_data$Sample)))

# Plot Abundance -----------------------
abund_good_bad <- bad_asvs_abund %>% 
  full_join(bad_asvs_good_data_abund)

abund_good_bad_plot <- abund_good_bad %>% 
  ggplot(aes(x = bad_abundance, y = good_abundance)) +
  geom_point(size = 2.5, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Negative Controls",
       y = "Real Samples",
       title = "Average ASV Abundance / Sample") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks = seq(0, 800, 100), limits = c(0, 800)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks = seq(0, 800, 100), limits = c(0,800)) +
    theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

abund_good_bad_plot

# Merge Plots -------------------------------
my_plot <- plot_grid(rel_abun_good_bad_plot, 
                     prev_good_bad_plot,
                     abund_good_bad_plot)

my_plot

ggsave2(filename = "results/figs/contam_abun_prev.png",
        plot = my_plot)



