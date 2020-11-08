# remove_contaminants.R
# Ryan Johnson
# 3 May 2020
# Remove all ASV associated with PCR neg and Extract Neg

## Load Packages -------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(here)
library(RColorBrewer)

# Load Data and set Variables -----------------------------------------------
ps2 <- readRDS("data/processed/ps2.rds")
ps2_tibble <- psmelt(ps2)

# Select neg-ctrl data ----------------------------------
bad_data <- ps2_tibble %>% 
  filter(Group == "Neg-Ctrl") %>% 
  filter(Abundance > 0) %>% 
  pull(OTU)

# Filter out bad otus from ps2 -------------------------
ps2_good <- ps2_tibble #%>% 
  # filter(!OTU %in% bad_data)

# Create table of number of reads per sample after removing contaminants -------
ps2_before <- ps2_tibble %>% 
  group_by(Sample) %>% 
  summarise(before_total_reads = sum(Abundance))

ps2_after <- ps2_good %>% 
  group_by(Sample) %>% 
  summarise(after_total_reads = sum(Abundance))

ps2_before_after <- ps2_before %>% 
  full_join(ps2_after) %>% 
  mutate(percent_remaining = (after_total_reads / before_total_reads) * 100)

# Write to results
write_csv(x = ps2_before_after,
          path = here("results", "tables", "contaminant_removal.csv"))

# Convert to Relative Abundance ---------------------------
ps2_good_relabun <- ps2_good %>% 
  group_by(Sample) %>% 
  mutate(rel_abun = (Abundance / sum(Abundance)) * 100) %>% 
  select(OTU, Sample, rel_abun, everything()) %>% 
  ungroup()

# Try plotting phylum ------------------
# Get fill colors
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

phylum_abun_plot <- ps2_good_relabun %>% 
  group_by(Sample, Phylum, Group) %>% 
  summarise(phylum_abun = sum(rel_abun, na.rm = T), 
            phylum_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(total_reads = sum(phylum_reads)) %>% 
  ggplot(aes(x = Sample, y = phylum_abun, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Group, scales = "free") +
  scale_fill_manual(values = getPalette(length(unique(ps2_good_relabun$Phylum)))) +
  ylim(0,100) +
  labs(x = "Sample Name",
       y = "Relative Abundance",
       title = "Relative Abundance at Phylum Level") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

phylum_abun_plot

ggsave(filename = here("results", "figs", "barplot_phylum_patient_clean.png"), 
       plot = phylum_abun_plot, 
       height = 8,
       width = 15)


# Try plotting genus ------------------
# Get top N Genera
top_genera_number <- 15

top_N_genera <- ps2_good_relabun %>% 
  filter(!is.na(Genus)) %>% 
  filter(!is.na(rel_abun)) %>% 
  group_by(Genus) %>% 
  summarise(genera_mean_abun = mean(rel_abun)) %>% 
  arrange(desc(genera_mean_abun)) %>% 
  head(top_genera_number)

# Top N general
genus_abun_plot <- ps2_good_relabun %>% 
  filter(Genus %in% top_N_genera$Genus) %>% 
  group_by(Sample, Genus, Group) %>% 
  summarise(genus_abun = sum(rel_abun, na.rm = T), 
            genus_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(total_reads = sum(genus_reads)) %>% 
  ggplot(aes(x = Sample, y = genus_abun, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Group, scales = "free") +
  labs(title = paste0("Top ", top_genera_number, " Genera")) +
  scale_fill_manual(values = getPalette(top_genera_number)) +
  labs(x = "Sample",
       y = "Relative Abundance",
       title = paste0("Genus (top ", top_genera_number, ") Relative Abundance")) +
  ylim(0,100) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

genus_abun_plot

ggsave(filename = here("results", "figs", "barplot_genus_patient_clean.png"), 
       plot = genus_abun_plot, 
       height = 8,
       width = 15)
