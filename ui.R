library(shiny)
library(shinycssloaders)
library(plotly)
library(readxl)
library(bslib)

# make plot square

filepath = "./data/DE_protein.xlsx"
cohorts = c("BRCA","CCRCC","COAD","GBM","HCC","HNSCC","LSCC","LUAD","OV","PDA","UCEC")
cohortSheets <- excel_sheets(path = filepath)


shinyUI(navbarPage("ubiDash", theme = bs_theme(version = 5, bootswatch = "pulse"),
  includeCSS("./www/style.css"),
  
  # oncoplot tab ----
    tabPanel("Oncoplot", theme = bs_theme(version = 5, bootswatch = "flatly"),
           sidebarLayout(
             sidebarPanel(
               textInput("inputOncoGenesTXT", "Input text (separate genes with comma)", value = "", width = NULL, placeholder = NULL),
               sliderInput("maxOncoGenes", "Max Number of Genes to Display:",
                           min = 1, max = 50,
                           value = 10),
               h3("Cohorts"),
               checkboxGroupInput("oncoCohortSelect", "Select Cohorts:",
                                  cohorts, selected = cohorts),
               
               actionButton("generateOnco", "Generate Plot"),
               
               p("Download this data"),
               downloadButton(outputId = "downloadOncoMutations", label = "Download Mutations Data"),               
               downloadButton(outputId = "downloadOncoCohorts", label = "Download Cohort Data"),               
             ),
             mainPanel(
               textOutput("oncoError"),
               shinycssloaders::withSpinner(plotOutput("cohortOncoplot"), color = "#6A2798"),
               shinycssloaders::withSpinner(plotOutput("mutationOncoplot"), color = "#6A2798")
             )
          )
  ),
  # heatmap tab ----
  tabPanel("Heatmap",

    titlePanel("Heatmap"),

    sidebarLayout(
        sidebarPanel(
          h3("Genes"),
          # make sure these do not interfere
          textInput("inputGenesTXT", "Input text (separate genes with comma)", value = "", width = NULL, placeholder = NULL),
          p("OR"),
          fileInput("inputGenesCSV", "Upload CSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
          
          h3("Cohorts"),
          
          checkboxGroupInput("cohortSelect", "Select Cohorts:",
                             cohorts, selected = cohorts),
          checkboxInput("cluster", "Cluster columns",
                             value = TRUE),
          actionButton("generateHeatmap", "Generate Heatmap"),
          
          p("Download this data"),
          downloadButton(outputId = "downloadData", label = "Download"),
          
          # p("Download this plot"),
          # downloadButton(outputId = "downloadHeatmap", label = "Download")
          # 
          
        ),

        # Show a plot of the generated distribution
        mainPanel(
          textOutput("rowscut"),
          textOutput("err"),
          shinycssloaders::withSpinner(plotOutput("heatmap"), color = "#6A2798"),
            )
        
    )
  ),
  
  # volcano plot tab ----
  tabPanel("Volcano Plots", theme = bs_theme(version = 5, bootswatch = "pulse"),
    sidebarLayout(
      sidebarPanel(
        
        selectInput("cohort", "Cohort", cohortSheets),
        
        
          selectInput("PValCutoff", "Adjusted P Value Cutoff:", c(0.001, 0.01, 0.05)),
           sliderInput("logFCCutoff", "Log2(Tumor/Normal) Cutoff:",
                    min = .25, max = 2,
              value = 1, step = 0.05),
        
        actionButton("generateVolcano", "Generate Plot"),
        
        p("Download this data"),
        downloadButton(outputId = "downloadVolcanoData", label = "Download")
      ),
      mainPanel(
    #select which sheet from xlsx to read & generate graph
    shinycssloaders::withSpinner(plotlyOutput("distPlot"), color = "#6A2798"))
    )
  ),
  # correlation tab ----
  tabPanel("Correlations", theme = bs_theme(version = 5, bootswatch = "pulse"),
             sidebarLayout(
               sidebarPanel(
                 p("Download this data"),
                 downloadButton(outputId = "downloadCorrelationData", label = "Download")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Tumor Proteome - Normal Proteome",
                     shinycssloaders::withSpinner(plotlyOutput("pt_pn"), color = "#6A2798")
                   ),
                   tabPanel("Tumor Proteome - Tumor Transcriptome",
                     shinycssloaders::withSpinner(plotlyOutput("pt_rt"), color = "#6A2798")
                   )
                   # ,tabPanel("Normal Proteome - Tumor Transcriptome",
                   #   shinycssloaders::withSpinner(plotlyOutput("pn_rt"), color = "#6A2798")
                   # )
                   )
                )
            )
        )
    )
)
