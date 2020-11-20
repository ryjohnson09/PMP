# ASV Filtering


# Libraries ---------------------------------------------------------------
library(shiny)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ggiraph)
library(plotly)
library(decontam)
library(testthat)


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