# abundance_heatmap_negctrl.R
# Ryan Johnson
# 30 May 2020
# Analyze polished microbiome data to produce 
#  relative abundance plots at phylum, genus, and ASV
#  level. Only looking at ASV present in negctrl.

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


# negctrl ASVs ---------------------------------------------------------------
neg_ctrl_asv <- ps2_relabun %>% 
  filter(Gender == "Neg-Ctrl" & Abundance > 0)


# Filter for top ASVs -----------------------------------------------------
ps2_relabun_negctrl_ASVs <- ps2_relabun %>% 
  filter(OTU %in% neg_ctrl_asv$OTU)


# cluster -----------------------------------------------------------------
asv_matrix <- ps2_relabun_negctrl_ASVs %>% 
  select(OTU, rel_abun, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = rel_abun) %>% 
  column_to_rownames("OTU") %>% 
  as.matrix()

asv_cluster <- hclust(dist(asv_matrix))

OTU_order <- asv_cluster$labels[asv_cluster$order]


# Plot --------------------------------------------------------------------
my_breaks <- c(0, 0.01, 0.1, 1, 10, 100)
ps2_relabun_negctrl_ASVs$OTU <- factor(ps2_relabun_negctrl_ASVs$OTU, levels = OTU_order)

ps2_negctrlASV_heatmap <- ps2_relabun_negctrl_ASVs %>% 
  ggplot(aes(x = Sample, y = OTU, fill = rel_abun)) +
  geom_tile() +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, 
                       trans = scales::pseudo_log_trans(sigma = 0.001)) +
  #facet_wrap(~ Gender, scales = "free") +
  labs(title = "All ASV present in Neg-Ctrls") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ps2_negctrlASV_heatmap


