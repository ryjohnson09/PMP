# sample_material_remaining.R
# Ryan Johnson
# June 18, 2020
# Samples with material left over, how many reads left?

library(tidyverse)

final_reads <- read_csv("results/tables/track_reads.csv")

material_remaining <- c("PMP-367", "PMP-211", "PMP-478", "PMP-357", 
                        "PMP-227", "PMP-254", "PMP-231", "PMP-150", 
                        "PMP-423", "PMP-193", "CON-29", "PMP-501", 
                        "PMP-325", "USC-01", "PMP-359", "PMP-485", 
                        "CON-17", "PMP-433", "PMP-420", "PMP-136", 
                        "CON-18", "OV-02", "PMP-123", "PMP-127", 
                        "CON-24", "OV-01", "PMP-140", "CON-16", 
                        "CON-28", "PMP-185", "PMP-354", "PMP-492", 
                        "MUC-OV", "PMP-361", "PMP-187", "PMP-465", 
                        "PMP-186", "PMP-460", "PMP-407", "PMP-362", 
                        "PMP-482", "PMP-245", "PMP-436", "PMP-475", 
                        "OV-05")

final_reads_materials <- final_reads %>% 
  filter(`Sample Name` %in% material_remaining)

write_csv(final_reads_materials, "results/tables/samples_remaining.csv")
