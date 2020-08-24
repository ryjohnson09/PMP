# Load Data and set Variables -----------------------------------------------

library(tidyverse)
library(shiny)
library(phyloseq)
library(RColorBrewer)
library(plotly)

#ps2 <- readRDS("data/processed/ps2.rds")
ps2_tibble <- read_csv("ASV_Filtering_Shiny/ps2_tibble.csv")

prev_scale_Good <- c(0, 100)
prev_scale_Bad <- c(0, 100)

abun_scale_Good <- c(0, 650)
abun_scale_Bad <- c(5, 650)


bad_asvs <- ps2_tibble %>% 
  filter(Gender == "Neg-Ctrl") %>% 
  filter(Abundance > 0)

bla <- ps2_tibble %>% 
  filter(OTU %in% c(bad_asvs$OTU)) %>% 
  select(OTU, Gender, Sample, Abundance, Genus) %>% 
  mutate(Sample_Type = case_when(Gender == "Neg-Ctrl" ~ "Neg-Ctrl",
                                 Gender != "Neg-Ctrl" | is.na(Gender) ~ "Good")) %>% 
  mutate(Presence = ifelse(Abundance > 0, 1, 0)) %>% 
  group_by(Sample_Type, OTU, Genus) %>%
  summarise(mean_abundance = mean(Abundance), 
            sum_prevalence = sum(Presence), 
            n = n(),
            perc_prevalence = sum(Presence) / n() * 100) %>% 
  ungroup() %>% 
  select(Sample_Type, OTU, mean_abundance, perc_prevalence, Genus) %>% 
  pivot_wider(names_from = Sample_Type, values_from = c(mean_abundance, perc_prevalence)) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Other", levels(Genus)[Genus])) %>%
  mutate(keep = ifelse(
    between(mean_abundance_Good, abun_scale_Good[1], abun_scale_Good[2]) &
    between(`mean_abundance_Neg-Ctrl`, abun_scale_Bad[1], abun_scale_Bad[2]) &
    between(perc_prevalence_Good, prev_scale_Good[1], prev_scale_Good[2]) &
    between(`perc_prevalence_Neg-Ctrl`, prev_scale_Bad[1], prev_scale_Bad[2]), "Yes", "No")          
  )


prev_plot <- bla %>% 
  ggplot(aes(x = `perc_prevalence_Neg-Ctrl`, y = perc_prevalence_Good, color = keep, tooltip = Genus)) +
  geom_jitter_interactive(height = 0.3, width = 1)

abun_plot <- bla %>% 
  ggplot(aes(x = `mean_abundance_Neg-Ctrl`, y = mean_abundance_Good, color = keep, tooltip = Genus)) +
  geom_point()

girafe(code = print(prev_plot))
abun_plot







getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

# Filter out bad otus from ps2
ps2_good_relabun <- ps2_tibble %>% 
  # Remove contaminant ASVs
  filter(!OTU %in% filter(bla, keep == "No")$OTU) %>% 
  # Convert to Relative Abundance
  select(-OTU) %>% 
  group_by(Sample) %>% 
  mutate(ASV_rel_abun = (Abundance / sum(Abundance)) * 100) %>% 
  ungroup() %>% 
  # If no Genus taxon, convert to "Other"
  mutate(Genus = ifelse(is.na(Genus), "Other", levels(Genus)[Genus])) %>% 
  group_by(Sample, Genus, Gender) %>% 
  summarise(Genus_rel_abun = sum(ASV_rel_abun)) %>% 
  ungroup()

# Get top N genera
top_N_genera <- ps2_good_relabun %>% 
  group_by(Genus) %>% 
  summarise(genera_mean_abun = mean(Genus_rel_abun, na.rm = TRUE)) %>% 
  arrange(desc(genera_mean_abun)) %>% 
  head(15)

# Top N genera plot
genus_abun_plot <- ps2_good_relabun %>% 
  filter(Genus %in% top_N_genera$Genus) %>% 
  group_by(Sample, Genus, Gender) %>% 
  summarise(genus_abun = sum(Genus_rel_abun, na.rm = T)) %>%  
            # genus_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  #group_by(Sample) %>% 
  #mutate(total_reads = sum(genus_reads)) %>% 
  ggplot(aes(x = Sample, y = genus_abun, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Gender, scales = "free") +
  labs(title = paste0("Top ", 25, " Genera")) +
  scale_fill_manual(values = getPalette(25)) +
  labs(x = "Sample",
       y = "Relative Abundance",
       title = paste0("Genus (top ", 25, ") Relative Abundance")) +
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
