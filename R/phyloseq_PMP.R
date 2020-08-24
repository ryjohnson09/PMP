# Phyloseq PMP Analysis
# 14 March 2020
# Ryan Johnson
# Process the DADA2 analyzed PMP sequence dat

library(tidyverse)
library(phyloseq)
library(here)
library(readxl)

# Read in DADA2 Data ------------------------------------
seqtab <- readRDS(here("data", "processed", "seqtab_PMP.rds"))
tax <- readRDS(here("data", "processed", "tax_PMP.rds"))


# Create Sample Metadata --------------------------------
metadata <- read_excel(here("data", "raw", "PMP samples in plate.xlsx")) %>% 
  select(Patient...6, Gender, Diagnosis)

colnames(metadata) <- c("Sample_Name", "Gender", "Diagnosis")

# Fix sample names to match rownames(seqtab)
metadata_final <- metadata %>% 
  filter(!Sample_Name == "Nothing") %>% 
  mutate(Sample_Name = str_replace_all(Sample_Name, " ", "-")) %>% 
  mutate(Sample_Name = str_replace_all(Sample_Name, "PMP-CON", "CON")) %>% 
  mutate(Sample_Name = str_replace_all(Sample_Name, "Mock-10\\^", "Mock-E")) %>%
  mutate(Gender = ifelse(str_detect(Sample_Name, "^Mock"), "Mock", Gender)) %>% 
  mutate(Gender = ifelse(str_detect(Sample_Name, "^CON"), "Control", Gender)) %>% 
  mutate(Gender = ifelse(str_detect(Sample_Name, "^OV"), "OV", Gender)) %>% 
  full_join(tibble(Sample_Name = rownames(seqtab))) %>% 
  filter(Sample_Name %in% rownames(seqtab)) %>% 
  column_to_rownames("Sample_Name")

# Convert to data frame
metadata_final <- data.frame(metadata_final)

# Save metada file
write_csv(rownames_to_column(metadata_final, "Sample Name"), 
          path = "results/tables/PMP_metadata.csv")


# Make Phyloseq Object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(metadata_final), 
               tax_table(tax))
ps


# Filter out Phyla with Low Prevalence -----------------------------------
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Remove samples with 0 reads from phyloseq object
ps <- prune_samples(sample_sums(ps) > 0, ps)

# Filter out any ambiguous Phyla
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Generate Prevalence Functions
prevdf <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})

# Filter out phyla based on above findings
## This can be adjusted
filterPhyla <- c()#"Chloroflexi",
                 #"Euryarchaeota") # Remove Archaea
                 #"Fusobacteria",
                 #"Spirochaetes",
                 #"Synergistetes",
                 #"Tenericutes")

ps1 <- subset_taxa(ps, !Phylum %in% filterPhyla)
ps1


# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

# Plot (x is total # of reads in phyla among ALL samples,
#  y is prevalence among all samples divided by total sample #)
prev_plot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps), color=Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  #geom_hline(yintercept = 0.00, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + 
  theme(legend.position="none")

ggsave(filename = here("results", "figs", "Phylum_prevalence_filtering.png"),
       plot = prev_plot,
       height = 8,
       width = 10)


# Filter out low prevalence taxa (present in less than N% of samples)
prevalenceThreshold <- 0.00 * nsamples(ps)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 <- prune_taxa(keepTaxa, ps)


# Save new ps2 object
saveRDS(ps2, here("data", "processed", "ps2.rds"))

