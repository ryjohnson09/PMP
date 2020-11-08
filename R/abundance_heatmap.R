# abundance_heatmap.R
# Ryan Johnson
# 30 May 2020
# Analyze polished microbiome data to produce 
#  relative abundance plots at phylum, genus, and ASV
#  level

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)
library(ggthemes)
library(RColorBrewer)
library(scales)



# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds") %>% 
  psmelt()


# Convert to Relative abundance ------------------------------------------
ps2_relabun <- ps2 %>% 
  group_by(Sample) %>% 
  mutate(sample_total = sum(Abundance)) %>% 
  mutate(rel_abun = (Abundance / sample_total) * 100) %>% 
  select(OTU, Sample, Abundance, sample_total, rel_abun, everything())


# Top ASVs ---------------------------------------------------------------
top_N_asv <- ps2_relabun %>% 
  group_by(OTU) %>% 
  summarise(asv_mean_abun = mean(rel_abun)) %>% 
  arrange(desc(asv_mean_abun)) %>% 
  head(100)


# Filter for top ASVs -----------------------------------------------------
ps2_relabun_top_ASV <- ps2_relabun %>% 
  filter(OTU %in% top_N_asv$OTU)


# Plot --------------------------------------------------------------------
my_breaks <- c(0, 0.01, 0.1, 1, 10, 100)

ps2_topASV_heatmap <- ps2_relabun_top_ASV %>% 
  ggplot(aes(x = Sample, y = OTU, fill = rel_abun)) +
  geom_tile() +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, 
                       trans = scales::pseudo_log_trans(sigma = 0.001)) +
  facet_wrap(~ Group, scales = "free") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ps2_topASV_heatmap


