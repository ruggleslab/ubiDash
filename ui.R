library(shiny)
library(shinycssloaders)
library(plotly)
library(readxl)
library(bslib)
library(DT)
library(corrplot)

# make plot square

filepath = "./data/DE_protein.xlsx"
cohorts = c("BRCA","CCRCC","COAD","GBM","HCC","HNSCC","LSCC","LUAD","OV","PDA","UCEC")
cohortSheets <- excel_sheets(path = filepath)
abstract =  paste(read.delim("./www/assets/abstract.txt"))

shinyUI(navbarPage("ubiDash", theme = bs_theme(version = 5, bootswatch = "pulse"),
  includeCSS("./www/style.css"),
  # home page tab ----
    tabPanel("Home", theme = bs_theme(version = 5, bootswatch = "flatly"),
             tags$div(
             h1("Integrated multi-omic analysis of Cullin-RING E3 ligases in cancer"),
              p(abstract),
            img(src='assets/OverviewSchematic.png', align = "right")%>% 
              tagAppendAttributes(class = 'home-img')
             )  %>%  
              tagAppendAttributes(class = 'title-page')
             
             ),
  # oncoplot tab ----
    tabPanel("Oncoplot", theme = bs_theme(version = 5, bootswatch = "flatly"),
           sidebarLayout(
             sidebarPanel(
               
               # gene inputs
               
               selectizeInput(inputId = "inputOncoGenesTXT",
                              choices = NULL,
                              label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
               p("OR"),
               fileInput("inputOncoGenesCSV", "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
                selectInput("oncoDelims", "Delimiter", c("tab", "comma")),               
               
               # options
               sliderInput("maxOncoGenes", "Max Number of Genes to Display:",
                           min = 1, max = 50,
                           value = 10),
               h3("Cohorts"),
               selectizeInput(inputId = "oncoCohortSelect",
                              choices = cohorts,selected = cohorts,
                              label = "Select cohorts",options = list(`live-search` = TRUE),  multiple = TRUE),
               
               actionButton("generateOnco", "Generate Plot"),
               
               # download options
               p("Download this data"),
               downloadButton(outputId = "downloadOncoMutations", label = "Download Mutations Data"),               
               downloadButton(outputId = "downloadOncoCohorts", label = "Download Cohort Data"),               
             ),
             mainPanel(
               textOutput("oncoError")%>% 
                 tagAppendAttributes(class = 'error'),
               # loading icons
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
          selectizeInput(inputId = "inputGenesTXT",
                         choices = NULL,
                         label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
          p("OR"),
          fileInput("inputGenesCSV", "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
          actionButton("resetHeatmapGenesCSV","clear"),
          selectInput("heatmapDelims", "Delimiter", c("tab", "comma")),               
          
          # cohort seleciton checkboxes
          h3("Cohorts"),
          selectizeInput(inputId = "cohortSelect",
                         choices = cohorts,selected = cohorts,
                         label = "Select cohorts",options = list(`live-search` = TRUE),  multiple = TRUE),
          checkboxInput("cluster", "Cluster columns",
                             value = TRUE),
          
          # toggle clustering and separation by tumor and normal
          checkboxInput("separateNandT", "Separate normal and tumor",
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
          textOutput("rowscut")%>% 
            tagAppendAttributes(class = 'error'),
          textOutput("err")%>% 
            tagAppendAttributes(class = 'error'),
          shinycssloaders::withSpinner(plotOutput("heatmap"), color = "#6A2798"),
            )
        
    )
  ),
  
  # volcano plot tab ----
  tabPanel("Volcano Plots", theme = bs_theme(version = 5, bootswatch = "pulse"),
    sidebarLayout(
      sidebarPanel(
        
        selectInput("cohort", "Cohort", cohortSheets),
       
        # input options
        selectizeInput(inputId = "inputVolcanoGenesTXT",
                       choices = NULL,
                       label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
        p("OR"),
        fileInput("inputVolcanoGenesCSV", "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c( ".tsv", "text/tsv",'.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
        selectInput("volcanoDelims", "Delimiter", c("tab", "comma")),               
        
        # cutoff options
        selectInput("PValCutoff", "Adjusted P Value Cutoff:", c(0.001, 0.01, 0.05)),
        sliderInput("logFCCutoff", "Log2(Tumor/Normal) Cutoff:",
              min = .25, max = 2, value = 1, step = 0.05),
        
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
                 
                 # gene input options
                 p("The underlying data for this tab is a pre-filtered table 
                   only containing genes correlated to CRLs with cor > | 0.5 |. 
                   So, input text should specify a CRL of interest to produce a 
                   scatterplot of all other proteins correlated to the specified 
                   CRL."),
                 textInput("inputCorrGenesTXT", "Input text (separate genes with comma)", value = "", width = NULL, placeholder = NULL),
                 
                 sliderInput("CorrelationCutoff", "Correlation Cutoff:",
                             min = .5, max = 1,
                             value = .5, step = 0.05),
                 actionButton("generateCorrelation", "Generate Plot"),
                 
                 p("Download this data"),
                 downloadButton(outputId = "downloadCorrelationData", label = "Download")
              
                  ),
               mainPanel(
                 tabsetPanel(
                   # three tabs to compare proteome and genome of normal and tumor
                   tabPanel("Tumor Proteome - Normal Proteome",
                   textOutput("corrErr1") %>% 
                     tagAppendAttributes(class = 'error'),
                     shinycssloaders::withSpinner(plotlyOutput("pt_pn"), color = "#6A2798")
                   ),
                   tabPanel("Tumor Proteome - Tumor Transcriptome",
                   textOutput("corrErr2") %>% 
                     tagAppendAttributes(class = 'error'),
                     shinycssloaders::withSpinner(plotlyOutput("pt_rt"), color = "#6A2798")
                   ),
                   tabPanel("Tumor Transcriptome - Normal Transcriptome",
                   textOutput("corrErr3") %>% 
                     tagAppendAttributes(class = 'error'),
                     shinycssloaders::withSpinner(plotlyOutput("rt_rn"), color = "#6A2798")
                   )
                  
                   )
                )
            )
        ),
  # violin plot tap ----
  # need to let users input genes of interest
  tabPanel("Violin Plot", theme = bs_theme(version = 5, bootswatch = "pulse"),
           
                p("Violin Plot"),
                shinycssloaders::withSpinner(plotlyOutput("violin"), color = "#6A2798")
           
        ),
  navbarMenu("Correlation Analyses",
    tabPanel("Gene vs Gene",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvGpanel",selected = "mRNA",
                             tabPanel("mRNA",
                                      selectizeInput(inputId = "inputGvgCorrRnaNameA",
                                                     choices = NULL,
                                                     label = "Choose protein A",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = "inputGvgCorrRnaNameB",
                                                     choices = NULL,
                                                     label = "Choose protein B",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      
                                      selectizeInput(inputId = "inputGvgCorrRnaCohort",
                                                     choices = cohorts,selected = "BRCA",
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvgCorrRnaStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvgCorrRna", label = "Start Analysis")
                                      ),
                             tabPanel("Protein",
                                      selectizeInput(inputId = "inputGvgCorrProteinNameA",
                                                     choices = NULL,
                                                     label = "Choose protein A",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = "inputGvgCorrProteinNameB",
                                                     choices = NULL,
                                                     label = "Choose protein B",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      
                                      selectizeInput(inputId = "inputGvgCorrProteinCohort",
                                                     choices = cohorts,selected = "BRCA",
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvgCorrProteinStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvgCorrProtein", label = "Start Analysis")
                                      ),
                             tabPanel("mRNA/Protein",
                                      selectizeInput(inputId = "inputGvgCorrBothNameA",
                                                     choices = NULL,
                                                     label = "Choose protein",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = "inputGvgCorrBothNameB",
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = "inputGvgCorrBothCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",selected = "BRCA",
                                                     options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvgCorrBothStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvgCorrBoth", label = "Start Analysis")
                                      ))
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvGpanel == 'mRNA'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvgCorrRna"),color = "#6A2798"),
                                  tableOutput(outputId = "outputTabGvgCorrRna")
                                  ),
                 conditionalPanel(condition = "input.GvGpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvgCorrProtein"),color = "#6A2798"),
                                  tableOutput(outputId = "outputTabGvgCorrProtein")
                                  ),
                 conditionalPanel(condition = "input.GvGpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvgCorrBoth"),color = "#6A2798"),
                                  tableOutput(outputId = "outputTabGvgCorrBoth")
                 )
                 
               )
             )),
    tabPanel("Gene vs All",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvApanel", selected = "mRNA",
                             tabPanel(title = "mRNA",
                                      selectizeInput(inputId = "inputGvACorrRnaGeneName",
                                                     choices = NULL,
                                                     label = "Choose a gene",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = "inputGvACorrRnaCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvACorrRnaType", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      sliderInput(inputId = "inputGvACorrRnaR",
                                                  label = "Set correlation coefficient cut-off",
                                                  value = c(0.5,1), min = -1, max = 1, step = 0.05),
                                      
                                      actionButton("generateGvACorrRna", label = "Start Analysis")
                                      ),
                             tabPanel(title = "Protein",
                                      selectizeInput(inputId = "inputGvACorrProteinGeneName",
                                                     choices = NULL,
                                                     label = "Choose a protein",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = "inputGvACorrProteinCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvACorrProteinType", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      sliderInput(inputId = "inputGvACorrProteinR",
                                                  label = "Set correlation coefficient cut-off",
                                                  value = c(0.5,1), min = -1, max = 1, step = 0.05),
                                      
                                      actionButton("generateGvACorrProtein", label = "Start Analysis")
                                      )
                             )
                 
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvApanel == 'Protein'",
                                  fluidRow(column(width = 4, 
                                                  withSpinner(dataTableOutput(outputId = "outputTabGvACorrProtein"),color = "#6A2798")),
                                           column(width = 8, 
                                                  withSpinner(plotOutput(outputId = "outputPlotGvACorrProtein"),color = "#6A2798")))
                                  ),
                 conditionalPanel(condition = "input.GvApanel == 'mRNA'",
                                  fluidRow(column(width = 4, 
                                                  withSpinner(dataTableOutput(outputId = "outputTabGvACorrRna"),color = "#6A2798")),
                                           column(width = 8, 
                                                  withSpinner(plotOutput(outputId = "outputPlotGvACorrRna"),color = "#6A2798")))
                                  )
               )
             )
    ),
    tabPanel("Gene Panel",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvPpanel",selected = "mRNA",
                             tabPanel("mRNA",
                                      selectizeInput(inputId = "inputGvpCorrRnaGeneName",
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = "inputGvpCorrRnaCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvpCorrRnaStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvpCorrRna", label = "Start Analysis")
                             ),
                             tabPanel("Protein",
                                      selectizeInput(inputId = "inputGvpCorrProteinGeneName",
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = "inputGvpCorrProteinCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvpCorrProteinStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvpCorrProtein", label = "Start Analysis")
                             ),
                             tabPanel("mRNA/Protein",
                                      selectizeInput(inputId = "inputGvpCorrBothGeneNameA",
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      selectizeInput(inputId = "inputGvpCorrBothGeneNameB",
                                                     choices = NULL,
                                                     label = "Choose protein",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = "inputGvpCorrBothCohort",
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = "inputGvpCorrBothStype", label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton("generateGvpCorrBoth", label = "Start Analysis")
                             ))
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvPpanel == 'mRNA'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrRna"),color = "#6A2798")
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrProtein"),color = "#6A2798")
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrBoth"),color = "#6A2798")
                 )
                 
               )
             ))
    
  )

  
    )
)
