# abundance_barplots.R
# Ryan Johnson
# 14 March 2020
# Analyze polished microbiome data to produce 
#  relative abundance plots at phylum and genus
#  level

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)
library(ggthemes)
library(RColorBrewer)



# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds")

# Variables
num_top_genus <- 15 # How many most abundant genera to consider in plots



# Covert to Relative Abundance ----------------------------------------------
# Combine data to Phylum Level
ps_phylum <- tax_glom(ps2, "Phylum", NArm = TRUE)

# Convert phylum data to Relative Abundance
ps_phylum_relabun <- transform_sample_counts(ps_phylum, 
                                             function(OTU) OTU/sum(OTU) * 100)
# Convert phylum relative abundance to data frame
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)

# Combine data to Genus Level
ps_genus <- tax_glom(ps2, "Genus", NArm = TRUE)

# Convert genus data to relative abundance
ps_genus_relabun <- transform_sample_counts(ps_genus, 
                                            function(OTU) OTU/sum(OTU) * 100)

# Get top N genera
topN_genus <- names(sort(taxa_sums(ps_genus), decreasing=TRUE))[1:num_top_genus]

# Extract the top N genera
ps_genus_topN <- prune_taxa(topN_genus, ps_genus_relabun)

# Convert top genus relative abundance to data frame
taxa_abundance_table_genus <- psmelt(ps_genus_topN)


# Phylum Level Bar Plots -----------------------------------------------------
# Get fill colors
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
# Phylum Level By Patient
barplot_phylum_patient <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Sample Name",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level") +
  ylim(0,100) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = getPalette(length(unique(taxa_abundance_table_phylum$Phylum)))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = here("results", "figs", "barplot_phylum_patient.png"), 
       plot = barplot_phylum_patient, 
       height = 8,
       width = 15)




# Genus Level Bar Plots -----------------------------------------------------
# Genus level by patient
barplot_genus_patient <- taxa_abundance_table_genus %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Sample",
       y = "Relative Abundance",
       title = paste0("Genus (top ", num_top_genus, ") Relative Abundance")) +
  ylim(0,100) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = getPalette(num_top_genus)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = here("results", "figs", "barplot_genus_patient.png"), 
       plot = barplot_genus_patient, 
       height = 8,
       width = 15)

