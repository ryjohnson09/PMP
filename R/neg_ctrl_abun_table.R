# neg_ctrl_abun_table.R
# Ryan Johnson
# 07 August 2020
# Create table showing abundance of neg ctrl ASVs in neg ctrl samples
#   and in good sample

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)

# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds")
ps2_tibble <- psmelt(ps2)

# Select neg-ctrl data ----------------------------------
bad_asvs <- ps2_tibble %>% 
  filter(Group == "Neg-Ctrl") %>% 
  filter(Abundance > 0)

# Bad ASVs abundance per sample ----------------------------------
bad_asvs_abund <- bad_asvs %>% 
  select(OTU, Sample, Abundance) %>% 
  group_by(OTU) %>% 
  summarise(bad_abundance = sum(Abundance) / length(unique(bad_asvs$Sample)))

# Add in each negctrl sample
bad_asvs_abund_sample <- bad_asvs %>% 
  select(OTU, Sample, Abundance) %>% 
  pivot_wider(names_from = Sample, 
              values_from = Abundance,
              values_fill = 0) %>% 
  full_join(bad_asvs_abund) %>% 
  select(OTU, "Avg NegCtrl Abundance" = bad_abundance, everything())


# Select same ASVs in non-neg-ctrl data -----------------------------
bad_asvs_good_data <- ps2_tibble %>% 
  filter(Group != "Neg-Ctrl") %>% 
  group_by(Sample) %>% 
  mutate(rel_abun_good = (Abundance / sum(Abundance)) * 100) %>% 
  filter(OTU %in% bad_asvs$OTU) %>% 
  ungroup()

# Bad ASVs good data abundance ----------------------------------
bad_asvs_good_data_abund <- bad_asvs_good_data %>% 
  select(OTU, Sample, Abundance) %>% 
  group_by(OTU) %>% 
  summarise(good_abundance = sum(Abundance) / length(unique(bad_asvs_good_data$Sample)))

# Add in each good sample
bad_asvs_good_abund_sample <- bad_asvs_good_data %>% 
  select(OTU, Sample, Abundance) %>% 
  pivot_wider(names_from = Sample, 
              values_from = Abundance,
              values_fill = 0) %>% 
  full_join(bad_asvs_good_data_abund) %>% 
  select(OTU, "Avg Good Abundance" = good_abundance, everything())


# Final table data ----------------------------
asv_abund <- bad_asvs_abund_sample %>% 
  full_join(bad_asvs_good_abund_sample) %>% 
  left_join(unique(select(bad_asvs, OTU, Genus))) %>% 
  select(OTU, Genus, everything())

write_csv(asv_abund, "results/tables/ASV_abun.csv")
