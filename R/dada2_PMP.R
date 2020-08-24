# dada2_PMP.R
# Author: Ryan Johnson
# date: 2 March 2020
# Purpose: Analyze the raw reads from the PMP study
#  using DADA2. Primer sequences have been removed from 
#  raw reads.

# Set up Environment -----------------------------------------------
library(tidyverse)
library(dada2); packageVersion("dada2")
library(here); here()
library(ShortRead)

# Get path to raw reads
raw_reads <- here("data", "raw", "raw_reads")

# Get list of forward/reverse reads
f_reads <- sort(list.files(raw_reads, pattern="_R1_001.fastq", full.names = TRUE))
r_reads <- sort(list.files(raw_reads, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names
sample_names <- basename(f_reads) %>% 
  str_split(pattern = "_") %>% 
  map(function(x) x[[1]]) %>% 
  unlist()

# Get directories and file names for filtered seqs
if(!dir.exists(here("data", "processed", "filtered_reads"))){
  dir.create(here("data", "processed", "filtered_reads"))
}

f_reads_filt <- here("data", "processed", "filtered_reads", paste0(sample_names, "_F_filt.fastq.gz"))
r_reads_filt <- here("data", "processed", "filtered_reads", paste0(sample_names, "_R_filt.fastq.gz"))

names(f_reads_filt) <- sample_names
names(r_reads_filt) <- sample_names

# ensure that filtered reads dirs is empty
if(length(dir(here("data", "processed", "filtered_reads"))) > 0){
  unlink(dir(here("data", "processed", "filtered_reads"), full.names = T))
}

# Count reads per sample to start -------------------------------------
init_read_count <- tibble(sample_id = character(0), num_reads = numeric(0))

for (i in seq_along(f_reads)){
  seq_row <- tibble(sample_id = sample_names[i],
                    num_reads = length(readFastq(f_reads[i])))
  init_read_count <- init_read_count %>% 
    full_join(seq_row, by = c("sample_id", "num_reads"))
  rm(seq_row)
}

# Write to tables directory
if(!dir.exists(here("results", "tables"))){
  dir.create(here("results", "tables"))
}

write_csv(x = init_read_count,
          path = here("results", "tables", "initial_read_counts.csv"))



# Remove samples with low # of reads --------------------------------
## This can be adjusted
samples_to_remove <- init_read_count %>% 
  filter(num_reads <= 0) %>% # remove samples with less than 500 reads
  pull(sample_id)

message(paste0("Sample(s) removed (<500 reads): ", paste(samples_to_remove, collapse = ", ")))

f_reads <- f_reads[!sample_names %in% samples_to_remove]
f_reads_filt <- f_reads_filt[!sample_names %in% samples_to_remove]
r_reads <- r_reads[!sample_names %in% samples_to_remove]
r_reads_filt <- r_reads_filt[!sample_names %in% samples_to_remove]



# Quality plots for raw reads -------------------------------------
if(!dir.exists(here("results", "figs"))){
  dir.create(here("results", "figs"))
}

f_reads_qualplot <- plotQualityProfile(f_reads[c(1,5,10,20,30,40,50,60,70)])
r_reads_qualplot <- plotQualityProfile(r_reads[c(1,5,10,20,30,40,50,60,70)])

# Save plots in results
ggsave(filename = here("results", "figs", "f_reads_qualplot.png"),
       plot = f_reads_qualplot,
       device = "png")
ggsave(filename = here("results", "figs", "r_reads_qualplot.png"), 
       plot = r_reads_qualplot,
       device = "png")


# Quality Filtering ------------------------------------------
filt_out <- filterAndTrim(f_reads, 
                          f_reads_filt, 
                          r_reads, 
                          r_reads_filt, 
                          truncLen = c(200, 190),
                          maxN = 0, 
                          maxEE = c(2, 2), 
                          truncQ = 2, 
                          rm.phix = TRUE,
                          compress = TRUE, 
                          multithread = TRUE)

# convert to tibble and write results
filt_out <- as_tibble(filt_out) %>% 
  mutate(sample_name = names(f_reads_filt)) %>% 
  select(sample_name, reads.in, reads.out)

write_csv(x = filt_out, 
          path = here("results", "tables", "intial_filtering_results.csv"))


# Remove any samples with less than N reads after filtering
## This can be adjusted
low_read_samps <- filt_out %>% 
  #filter(reads.out < 100) %>%
  filter(reads.out == 0) %>% 
  pull(sample_name)

low_read_samps_path_f <- f_reads_filt[
  purrr::map_lgl(f_reads_filt, ~any(str_detect(.x, low_read_samps)))]

low_read_samps_path_r <- r_reads_filt[
  purrr::map_lgl(r_reads_filt, ~any(str_detect(.x, low_read_samps)))]

low_reads_samps_path <- c(low_read_samps_path_f, low_read_samps_path_r)

file.remove(low_reads_samps_path)

# update filt paths
f_reads_filt <- f_reads_filt[f_reads_filt %in% dir(here("data",
                                                        "processed",
                                                        "filtered_reads"),
                                                       full.names = T)]

r_reads_filt <- r_reads_filt[r_reads_filt %in% dir(here("data",
                                                        "processed",
                                                        "filtered_reads"),
                                                   full.names = T)]


# Learn Error rates ----------------------------------------------
errF <- learnErrors(f_reads_filt, multithread=TRUE)
errR <- learnErrors(r_reads_filt, multithread=TRUE)

f_error_plot <- plotErrors(errF, nominalQ=TRUE)
r_error_plot <- plotErrors(errR, nominalQ = TRUE)

ggsave(plot = f_error_plot, 
       filename = here("results", "figs", "f_error_plot.png"))
ggsave(plot = r_error_plot, 
       filename = here("results", "figs", "r_error_plot.png"))


# Sample Inference ------------------------------------------
dadaFs <- dada(f_reads_filt, err=errF, multithread=TRUE)
dadaRs <- dada(r_reads_filt, err=errR, multithread=TRUE)


# Merge reads ----------------------------------------------
mergers <- mergePairs(dadaFs, f_reads_filt, 
                      dadaRs, r_reads_filt, 
                      verbose=TRUE)


# Construct sequence table --------------------------------
seq_tab <- makeSequenceTable(mergers)
dim(seq_tab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seq_tab)))

# Trim seq_tab to correct length
# Our amplicons should be around 253bp. We see a large chunk at 
# 221bp. I think this is likely chimeras or contaminants. We 
# need to remove these.
seq_tab2 <- seq_tab[,nchar(colnames(seq_tab)) %in% 251:256]



# Remove chimeras --------------------------------------------
seq_tab_nochim <- removeBimeraDenovo(seq_tab2, 
                                     method="consensus", 
                                     multithread=TRUE, 
                                     verbose = TRUE)

dim(seq_tab2)
dim(seq_tab_nochim)

sum(seq_tab_nochim)/sum(seq_tab2)


# Track Reads through pipline --------------------------------
# update filt_out
filt_out2 <- filt_out %>% 
  filter(sample_name %in% rownames(seq_tab_nochim))

getN <- function(x) sum(getUniques(x))

track <- cbind(filt_out2, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seq_tab_nochim))

# Change Column names
colnames(track) <- c("Sample Name", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

# Convert to tibble
track_tibble <- as_tibble(track)

# Write to results
write_csv(x = track_tibble,
          path = here("results", "tables", "track_reads.csv"))



# Assign taxonomy --------------------------------------------
tax <- assignTaxonomy(seqs = seq_tab_nochim, 
                      refFasta = "data/raw/silva_nr_v132_train_set.fa.gz", 
                      multithread = TRUE)

# Preview the taxonomy assignments
tax.print <- tax
rownames(tax.print) <- NULL
head(tax.print)


# Save Results to disk ----------------------------------
saveRDS(seq_tab_nochim, "data/processed/seqtab_PMP.rds")
saveRDS(tax, "data/processed/tax_PMP.rds")
