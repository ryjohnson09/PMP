# abundance_table.R
# Ryan Johnson
# 16 Sept 2020
# Generate table of ASV rel abun
# "a table of the relative abundances of 
# sequence features (columns) in each sample (rows)."

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)

# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds") %>% 
  psmelt()

# Get OTU and bacterial name key -------------------------------------------
asv_name <- ps2 %>% 
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  unique()

write_csv(asv_name, "results/tables/ASV_Taxa_key.csv")


# Convert to Relative abundance ------------------------------------------
ps2_relabun <- ps2 %>% 
  group_by(Sample) %>% 
  mutate(sample_total = sum(Abundance)) %>% 
  mutate(rel_abun = (Abundance / sample_total) * 100) %>% 
  select(OTU, Sample, Group, Diagnosis, rel_abun) %>% 
  ungroup()

# Convert columns to rel_abun values
ps2_relabun_wide <- ps2_relabun %>% 
  pivot_wider(names_from = OTU, values_from = rel_abun, values_fill = 0)

# Write Table
write_csv(ps2_relabun_wide, "results/tables/ASV_relabun.csv")
