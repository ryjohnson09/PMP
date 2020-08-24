# ASV Filtering


# Libraries ---------------------------------------------------------------
library(shiny)
library(tidyverse)
library(RColorBrewer)
library(ggiraph)


# UI ----------------------------------------------------------------------
ui <- fluidPage(
    
    fluidRow(
            column(3,
               # Prevalence
               sliderInput("n_prev",
                           "Prevalence in Neg-Ctrl Samples",
                           min = 0,
                           max = 100,
                           value = c(0, 100)),
               sliderInput("g_prev",
                           "Prevalence in Good Samples",
                           min = 0,
                           max = 100,
                           value = c(0, 100)),
               br(),
               #Abundance
               sliderInput("n_abun",
                           "Average Abundance (reads/sample) in Neg-Ctrl Samples",
                           min = 0,
                           max = 650,
                           value = c(0, 650)),
               sliderInput("g_abun",
                           "Average Abundance (reads/sample) in Good Samples",
                           min = 0,
                           max = 650,
                           value = c(5, 650)),
               br(),
               br(),
               # How many Genera
               sliderInput("N_genera",
                           "Show top N Genera",
                           min = 5,
                           max = 20,
                           value = 10)
            ),
            
            column(9,
                   mainPanel(
                       girafeOutput("prev_plot", width = "800px", height = "700px"),
                       girafeOutput("abun_plot", width = "800px", height = "700px"),
                       plotOutput("barplot", width = "1000px", height = "800px")
            )
        )
    ))



# Server ------------------------------------------------------------------
server <- function(input, output) {
    
    # Read in data
    ps2_tibble <- read_csv("ps2_tibble.csv")
    
    # Select contaminant ASVs in Neg-Ctrl Samples
    bad_asvs <- ps2_tibble %>% 
        filter(Gender == "Neg-Ctrl") %>% 
        filter(Abundance > 0)
    
    # Get mean abundance and percent prevalence of contaminant ASVs
    contam_ASVs <- reactive({
        ps2_tibble %>% 
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
            mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
            mutate(keep = ifelse(
                between(mean_abundance_Good, input$g_abun[1], input$g_abun[2]) &
                    between(`mean_abundance_Neg-Ctrl`, input$n_abun[1], input$n_abun[2]) &
                    between(perc_prevalence_Good, input$g_prev[1], input$g_prev[2]) &
                    between(`perc_prevalence_Neg-Ctrl`, input$n_prev[1], input$n_prev[2]), "Yes", "No"))          
        })
    
    # Prevalence Plot
    output$prev_plot <- renderGirafe({
        prev_plot1 <- contam_ASVs() %>% 
            ggplot(aes(x = `perc_prevalence_Neg-Ctrl`, y = perc_prevalence_Good, color = keep, tooltip = Genus)) +
            #geom_point(alpha = 0.7, size = 2.5) +
            geom_jitter_interactive(height = 1.5, width = 1.5, size = 2, alpha = 0.7) +
            theme_minimal() +
            labs(x = "Negative Controls",
                 y = "Real Samples",
                 title = "Percent Prevalence") +
            scale_color_manual(values=c("#d73027", "#4575b4"), 
                               name="Include or Remove\nContaminant ASV",
                               breaks=c("No", "Yes"),
                               labels=c("Remove", "Include")) +
            theme(
                axis.title = element_text(size = 15),
                plot.title = element_text(size = 18),
                axis.text = element_text(size = 12)
            )
        girafe(code = print(prev_plot1))
    })
    
    
    # Abundance Plot
    output$abun_plot <- renderGirafe({
        abund_plot1 <- contam_ASVs() %>% 
            ggplot(aes(x = `mean_abundance_Neg-Ctrl`, y = mean_abundance_Good, color = keep, tooltip = Genus)) +
            geom_point_interactive(size = 2, alpha = 0.7) +
            theme_minimal() +
            labs(x = "Negative Controls",
                 y = "Real Samples",
                 title = "Average ASV Abundance / Sample") +
            scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                               breaks = seq(0, 800, 200), limits = c(0, 800)) +
            scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), 
                               breaks = seq(0, 800, 200), limits = c(0,800)) +
            scale_color_manual(values=c("#d73027", "#4575b4"), 
                              name="Include or Remove\nContaminant ASV",
                              breaks=c("No", "Yes"),
                              labels=c("Remove", "Include")) +
            theme(
                axis.title = element_text(size = 15),
                plot.title = element_text(size = 18),
                axis.text = element_text(size = 12),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
            )
        girafe(code = print(abund_plot1))
    })
    
    
    # Rel-Abundance Plots (Genus)
    output$barplot <- renderPlot({
        
        # Get fill colors
        getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
        
        # Filter out bad otus from ps2
        ps2_good_relabun <- ps2_tibble %>% 
            # Remove contaminant ASVs
            filter(!OTU %in% filter(contam_ASVs(), keep == "No")$OTU) %>% 
            # Convert to Relative Abundance
            select(-OTU) %>% 
            group_by(Sample) %>% 
            mutate(ASV_rel_abun = (Abundance / sum(Abundance)) * 100) %>% 
            ungroup() %>% 
            # If no Genus taxon, convert to "Other"
            mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>% 
            group_by(Sample, Genus, Gender) %>% 
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
            group_by(Sample, Genus, Gender) %>% 
            summarise(genus_abun = sum(Genus_rel_abun, na.rm = T)) %>%  
            ungroup() %>% 
            group_by(Sample) %>% 
            ggplot(aes(x = Sample, y = genus_abun, fill = Genus)) +
            geom_bar(stat = "identity", color = "black") +
            facet_wrap(~Gender, scales = "free") +
            labs(title = paste0("Top ", input$N_genera, " Genera")) +
            scale_fill_manual(values = getPalette(input$N_genera)) +
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
shinyApp(ui = ui, server = server)


