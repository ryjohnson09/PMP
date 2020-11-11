# Load Data and set Variables -----------------------------------------------

library(tidyverse)
library(shiny)
library(phyloseq)
library(RColorBrewer)
library(plotly)
library(ggiraph)
library(Cairo)
library(decontam)

# Set prev and abun scale values
prev_scale_Good <- c(0, 100)
prev_scale_Bad <- c(0, 100)

abun_scale_Good <- c(5, 650)
abun_scale_Bad <- c(0, 650)

# Read in Data ----------------------------------------
ps2 <- readRDS("data/processed/ps2.rds")

# Determine Contaminants (neg-ctrl vs PMP/control samples) ----------------------
# Remove Mock
ps2_noMock <- subset_samples(ps2, Group != "Mock")
ps2_noMock_tibble <- as_tibble(psmelt(ps2_noMock))

# Call contaminants
sample_data(ps2_noMock)$is_neg <- sample_data(ps2_noMock)$True_Sample != TRUE
ps2_noMock_contam <- as_tibble(isContaminant(ps2_noMock, method="prevalence", neg="is_neg"), rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  rename(Contams_negctrl_vs_pmp = contaminant)

# Merge
ps2_contam_temp <- ps2_noMock_tibble %>% 
  left_join(select(ps2_noMock_contam, OTU, Contams_negctrl_vs_pmp), by = "OTU")



# Determine Contaminants (neg-ctrl vs mock samples) ----------------------
# Select Mock & Neg-Ctrl
ps2_Mock <- subset_samples(ps2, Group %in% c("Mock", "Neg-Ctrl"))

# Call contaminants
sample_data(ps2_Mock)$is_neg <- sample_data(ps2_Mock)$True_Sample != TRUE
ps2_Mock_contam <- as_tibble(isContaminant(ps2_Mock, method="prevalence", neg="is_neg"), rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  rename(Contams_negctrl_vs_mock = contaminant)

# Merge
ps2_contam <- ps2_contam_temp %>% 
  left_join(select(ps2_Mock_contam, OTU, Contams_negctrl_vs_mock), by = "OTU") %>% 
  # Final contam call
  mutate(contaminant = case_when(
    Contams_negctrl_vs_mock | Contams_negctrl_vs_pmp ~ "Yes",
    !Contams_negctrl_vs_mock & !Contams_negctrl_vs_pmp ~ "No"))

# Note: ps2_contam has the mock community removed and indicates if 
#  and ASV is a contaminant either by comparing neg-ctrl to PMP/control
#  and comparing mock vs neg-ctrl


# Prev and Abundance filtering ----------------------------------------
# Select contaminant ASVs in Neg-Ctrl Samples
bad_asvs <- ps2_contam %>% 
  filter(Group == "Neg-Ctrl") %>% 
  filter(Abundance > 0)

# Get mean abundance and percent prevalence of contaminant ASVs
contam_ASVs <- ps2_contam %>% 
  filter(OTU %in% c(bad_asvs$OTU)) %>% 
  select(OTU, Group, Sample, Abundance, Genus, contaminant) %>% 
  mutate(Sample_Type = case_when(Group == "Neg-Ctrl" ~ "Neg-Ctrl",
                                 Group != "Neg-Ctrl" | is.na(Group) ~ "Good")) %>% 
  mutate(Presence = ifelse(Abundance > 0, 1, 0)) %>% 
  group_by(Sample_Type, OTU, Genus, contaminant) %>%
  summarise(mean_abundance = mean(Abundance), 
            sum_prevalence = sum(Presence), 
            n = n(),
            perc_prevalence = sum(Presence) / n() * 100) %>% 
  ungroup() %>% 
  select(Sample_Type, OTU, mean_abundance, perc_prevalence, Genus, contaminant) %>% 
  pivot_wider(names_from = Sample_Type, values_from = c(mean_abundance, perc_prevalence)) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
  mutate(keep = ifelse(
    between(mean_abundance_Good, abun_scale_Good[1], abun_scale_Good[2]) &
    between(`mean_abundance_Neg-Ctrl`, abun_scale_Bad[1], abun_scale_Bad[2]) &
    between(perc_prevalence_Good, prev_scale_Good[1], prev_scale_Good[2]) &
    between(`perc_prevalence_Neg-Ctrl`, prev_scale_Bad[1], prev_scale_Bad[2]), "Yes", "No")          
  )

abun_plot <- contam_ASVs %>% 
  ggplot(aes(x = `mean_abundance_Neg-Ctrl`,
             y = mean_abundance_Good,
             fill = keep)) +
  # Lightly fill contaminants
  geom_point(data = filter(contam_ASVs, contaminant == "Yes"), 
             shape = 21, 
             size = 2, 
             alpha = 0.2) +
  # Dark fill non-contaminants
  geom_point(data = filter(contam_ASVs, contaminant == "No"), 
             shape = 21, 
             size = 2, 
             alpha = 0.7) +
  theme_minimal() +
  labs(x = "Negative Controls",
       y = "Real Samples") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 400, 100), limits = c(0, 400)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 400, 100), limits = c(0,400)) +
  scale_fill_manual(values=c("#d73027", "#4575b4"),
                     name="Include or Remove\nContaminant ASV",
                     breaks=c("No", "Yes"),
                     labels=c("Remove", "Include")) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

abun_plot

# Set seed to make jitter reproducible
set.seed(1234)
prev_plot <- contam_ASVs %>% 
  ggplot(aes(x = `perc_prevalence_Neg-Ctrl`,
             y = perc_prevalence_Good,
             fill = keep)) +
  # Lightly fill contaminants
  geom_jitter(data = filter(contam_ASVs, contaminant == "Yes"), 
              width = 2, height = 0.5,
              shape = 21,
              size = 2,
              alpha = 0.2) +
  # Dark fill non-contaminants
  geom_jitter(data = filter(contam_ASVs, contaminant == "No"), 
              width = 2, height = 0.5,
              shape = 21,
              size = 2,
              alpha = 0.7) +
  theme_minimal() +
  labs(x = "Negative Controls",
       y = "Real Samples") +
  scale_fill_manual(values=c("#d73027", "#4575b4"),
                     name="Include or Remove\nContaminant ASV",
                     breaks=c("No", "Yes"),
                     labels=c("Remove", "Include")) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

prev_plot





# Abundance Bar Plots --------------------------------------
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

# Filter out bad otus from ps2
ps2_good_relabun <- psmelt(ps2) %>% 
  # Remove contaminant ASVs
  filter(!OTU %in% filter(contam_ASVs, keep == "No" & contaminant == "No")$OTU) %>% 
  # Convert to Relative Abundance
  select(-OTU) %>% 
  group_by(Sample) %>% 
  mutate(ASV_rel_abun = (Abundance / sum(Abundance)) * 100) %>% 
  ungroup() %>% 
  # If no Genus taxon, convert to "Other"
  mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>% 
  group_by(Sample, Genus, Group) %>% 
  summarise(Genus_rel_abun = sum(ASV_rel_abun)) %>% 
  ungroup()

# Get top N genera
top_N_genera <- ps2_good_relabun %>% 
  group_by(Genus) %>% 
  summarise(genera_mean_abun = mean(Genus_rel_abun, na.rm = TRUE)) %>% 
  arrange(desc(genera_mean_abun)) %>% 
  head(10)

# Top N genera plot
genus_abun_plot <- ps2_good_relabun %>% 
  filter(Genus %in% top_N_genera$Genus) %>% 
  group_by(Sample, Genus, Group) %>% 
  summarise(genus_abun = sum(Genus_rel_abun, na.rm = T)) %>%  
            # genus_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  #group_by(Sample) %>% 
  #mutate(total_reads = sum(genus_reads)) %>% 
  ggplot(aes(x = Sample, y = genus_abun, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Group, scales = "free") +
  labs(title = paste0("Top ", 25, " Genera")) +
  scale_fill_manual(values = getPalette(10)) +
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
