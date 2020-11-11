# Server ------------------------------------------------------------------
server <- function(input, output) {
  
  ## Read in Data ----------------------------------------
  ps2 <- readRDS("ps2.rds")
  
  
  ## Determine Contaminants (neg-ctrl vs PMP/control samples) ----------------------
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
  
  
  ## Determine Contaminants (neg-ctrl vs mock samples) ------------------------------
  # Select Mock & Neg-Ctrl
  ps2_Mock <- subset_samples(ps2, Group %in% c("Mock", "Neg-Ctrl"))
  
  # Call contaminants
  sample_data(ps2_Mock)$is_neg <- sample_data(ps2_Mock)$True_Sample != TRUE
  ps2_Mock_contam <- as_tibble(isContaminant(ps2_Mock, method="prevalence", neg="is_neg"), rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    rename(Contams_negctrl_vs_mock = contaminant)
  
  ## Merge (final tibble)
  ps2_contam <- ps2_contam_temp %>% 
    left_join(select(ps2_Mock_contam, OTU, Contams_negctrl_vs_mock), by = "OTU") %>% 
    # Final contam call
    mutate(contaminant = case_when(
      Contams_negctrl_vs_mock | Contams_negctrl_vs_pmp ~ "Yes",
      !Contams_negctrl_vs_mock & !Contams_negctrl_vs_pmp ~ "No"))
  
  # Prev and Abundance filtering ----------------------------------------
  # Select contaminant ASVs in Neg-Ctrl Samples
  bad_asvs <- ps2_contam %>% 
    filter(Group == "Neg-Ctrl") %>% 
    filter(Abundance > 0)
  
  # Get mean abundance and percent prevalence of contaminant ASVs
  contam_ASVs <- reactive({
    ps2_contam %>% 
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
        between(mean_abundance_Good, input$g_abun[1], input$g_abun[2]) &
          between(`mean_abundance_Neg-Ctrl`, input$n_abun[1], input$n_abun[2]) &
          between(perc_prevalence_Good, input$g_prev[1], input$g_prev[2]) &
          between(`perc_prevalence_Neg-Ctrl`, input$n_prev[1], input$n_prev[2]), "Yes", "No"))
  })
  
  
  ## Abundance and Prevalence Scatter Plots -----------------------------------
  # Abundance plot
  output$abun_plot <- renderGirafe({
    abun_plot1 <- contam_ASVs() %>% 
      ggplot(aes(x = `mean_abundance_Neg-Ctrl`,
                 y = mean_abundance_Good,
                 fill = keep,
                 tooltip = Genus)) +
      # Lightly fill contaminants
      geom_point_interactive(data = filter(contam_ASVs(), contaminant == "Yes"), 
                 shape = 21, 
                 size = 2, 
                 alpha = 0.2) +
      # Dark fill non-contaminants
      geom_point_interactive(data = filter(contam_ASVs(), contaminant == "No"), 
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
    girafe(code = print(abun_plot1))
  })
  
  # Prevalence Plot
  output$prev_plot <- renderGirafe({
    set.seed(1234) # for consistent jitter
    prev_plot1 <- contam_ASVs() %>% 
      ggplot(aes(x = `perc_prevalence_Neg-Ctrl`,
                 y = perc_prevalence_Good,
                 fill = keep,
                 tooltip = Genus)) +
      # Lightly fill contaminants
      geom_jitter_interactive(data = filter(contam_ASVs(), contaminant == "Yes"), 
                  width = 2, height = 0.5,
                  shape = 21,
                  size = 2,
                  alpha = 0.2) +
      # Dark fill non-contaminants
      geom_jitter_interactive(data = filter(contam_ASVs(), contaminant == "No"), 
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
    girafe(code = print(prev_plot1))
  })

  
  
  ## Rel-Abundance Stacked Barplots (Genus) ----------------------------------
  output$barplot <- renderPlot({
    
    # Get fill colors
    getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
    
    # Filter out bad otus from ps2
    ps2_good_relabun <- psmelt(ps2) %>% 
      # Remove contaminant ASVs
      filter(!OTU %in% filter(contam_ASVs(), keep == "No" & contaminant == "No")$OTU) %>% 
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
      head(input$N_genera)
    
    # Top N genera plot
    genus_abun_plot <- ps2_good_relabun %>% 
      filter(Genus %in% top_N_genera$Genus) %>% 
      group_by(Sample, Genus, Group) %>% 
      summarise(genus_abun = sum(Genus_rel_abun, na.rm = T)) %>%  
      ungroup() %>% 
      ggplot(aes(x = Sample, y = genus_abun, fill = Genus)) +
      geom_bar(stat = "identity", color = "black") +
      facet_wrap(~Group, scales = "free") +
      labs(title = paste0("Top ", 25, " Genera")) +
      scale_fill_manual(values = getPalette(10)) +
      labs(x = "Sample",
           y = "Relative Abundance",
           title = paste0("Genus (top ", input$N_genera, ") Relative Abundance")) +
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
  })
}


# Run App -----------------------------------------------------------------
#shinyApp(ui = ui, server = server)