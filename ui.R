library(shiny)
library(shinycssloaders)
library(plotly)
library(readxl)
library(bslib)
library(DT)
library(corrplot)
library(openxlsx)

source('./modules/oncoplot.R')
source('./modules/heatmap.R')
source('./modules/volcano.R')
source('./modules/correlations.R')
source('./modules/violin.R')
source('./modules/lollipop.R')


# make plot square

filepath = "./data/DE_protein.xlsx"
cohorts = c("BRCA","CCRCC","COAD","GBM","HCC","HNSCC","LSCC","LUAD","OV","PDA","UCEC")
cohortSheets <- excel_sheets(path = filepath)
abstract =  paste(read.delim("./www/assets/abstract.txt"))


# source: https://stackoverflow.com/questions/38558827/making-an-image-hyperlink-in-r-shiny-header
shinyUI(
  
  navbarPage(img(src="assets/ubiDash_logo.png", width="120"), theme = bs_theme(version = 5, bootswatch = "pulse"),
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
                   navbarMenu("Mutations",
                              
                              # oncoplot tab ----
                              tabPanel("Oncoplot", theme = bs_theme(version = 5, bootswatch = "flatly"),
                                       oncoplotUI("oncoTab", cohorts)
                              ),
                              # lollipop tab ----
                              tabPanel("Lollipop", theme = bs_theme(version = 5, bootswatch = "flatly"),
                                       lollipopUI("lollipopTab", cohorts)
                              )
                   ),
                   # heatmap tab ----
                   tabPanel("Heatmap",
                            heatmapUI("heatmapTab", cohorts)
                   ),
                   
                   # volcano plot tab ----
                   tabPanel("Volcano Plots", theme = bs_theme(version = 5, bootswatch = "pulse"),
                            volcanoUI("volcanoTab", cohortSheets)
                   ),
             
                   # correlation tab ----
                   tabPanel("Correlations", theme = bs_theme(version = 5, bootswatch = "pulse"),
                            correlationsUI("corrTab")
                            
                   ),
                   # violin plot tap ----
                   # need to let users input genes of interest
                   tabPanel("Violin Plot", theme = bs_theme(version = 5, bootswatch = "pulse"),
                            
                            violinUI("violinTab")
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
                                                            tableOutput(outputId = "outputTabGvgCorrRna"),
                                                            uiOutput(outputId = "downLoadButtonGvgCorrRna")
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
                                                            withSpinner(plotOutput(outputId = "outputPlotGvpCorrRna"),color = "#6A2798"),
                                                            uiOutput("dataDownLoadGvpCorrRna")
                                           ),
                                           conditionalPanel(condition = "input.GvPpanel == 'Protein'",
                                                            withSpinner(plotOutput(outputId = "outputPlotGvpCorrProtein"),color = "#6A2798"),
                                                            uiOutput("dataDownLoadGvpCorrProtein")
                                           ),
                                           conditionalPanel(condition = "input.GvPpanel == 'mRNA/Protein'",
                                                            withSpinner(plotOutput(outputId = "outputPlotGvpCorrBoth"),color = "#6A2798"),
                                                            uiOutput("dataDownLoadGvpCorrBoth")
                                           )
                                           
                                         )
                                       ))
                              
                   )
<<<<<<< HEAD
                  
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
                                      br(),
                                      p("By this feature one can explore the correlation of two selected genes in the selected tissue type and patient cohorts."),
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
                                      br(),
                                      p("By this feature one can explore the correlation of two selected proteins in the selected tissue type and patient cohorts."),
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
                                      br(),
                                      p("By this feature one can explore the correlation of a selected gene and protein in the selected tissue type and patient cohorts."),
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
                                  tableOutput(outputId = "outputTabGvgCorrRna"),
                                  uiOutput(outputId = "downLoadButtonGvgCorrRnaPlotUI")
                                  ),
                 conditionalPanel(condition = "input.GvGpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvgCorrProtein"),color = "#6A2798"),
                                  tableOutput(outputId = "outputTabGvgCorrProtein"),
                                  uiOutput(outputId = "downLoadButtonGvgCorrProteinPlotUI")
                                  ),
                 conditionalPanel(condition = "input.GvGpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvgCorrBoth"),color = "#6A2798"),
                                  tableOutput(outputId = "outputTabGvgCorrBoth"),
                                  uiOutput(outputId = "downLoadButtonGvgCorrBothPlotUI")
                 )
                 
               )
             )),
    tabPanel("Gene vs All",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvApanel", selected = "mRNA",
                             tabPanel(title = "mRNA",
                                      br(),
                                      p("By this feature users can explore the correlation profile of a selected gene based on a user defined correlation cutoff intervall. The algorithm perfomrs a correlation analysis with all available genes and filters based on the predefined parameters.",
                                        "As output users get a table with proteins meeting all the criteria, and a plot where the y axis is the -log 10 of the P values and the x axis is the correlation coefficient."),
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
                                      br(),
                                      p("By this feature users can explore the correlation profile of a selected protein based on a user defined correlation cutoff intervall. The algorithm perfomrs a correlation analysis with all available proteins and filters based on the predefined parameters.",
                                        "As output users get a table with proteins meeting all the criteria, and a plot where the y axis is the -log 10 of the P values and the x axis is the correlation coefficient."),
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
                                                  withSpinner(plotOutput(outputId = "outputPlotGvACorrProtein"),color = "#6A2798"),
                                                  uiOutput(outputId = "downLoadButtonGvACorrProteinPlotUI")
                                                  ))
                                  ),
                 conditionalPanel(condition = "input.GvApanel == 'mRNA'",
                                  fluidRow(column(width = 4, 
                                                  withSpinner(dataTableOutput(outputId = "outputTabGvACorrRna"),color = "#6A2798")),
                                           column(width = 8, 
                                                  withSpinner(plotOutput(outputId = "outputPlotGvACorrRna"),color = "#6A2798"),
                                                  uiOutput(outputId = "downLoadButtonGvACorrRnaPlotUI")))
                                  )
               )
             )
    ),
    tabPanel("Gene Panel",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvPpanel",selected = "mRNA",
                             tabPanel("mRNA",
                                      br(),
                                      p("By using the gene panel correlation feature, useres can generate a correlation matrix using the genes of interests in the selected patient cohorts and tissue types.",
                                        "The resulting plot is a bubble plot where the size of the dots are proportional to the negative 10 logarithm of the adjusted P value, the color of the dots represents the correletion coefficient.",
                                        "If the q value is smaller then 0.01 the dot is highlighted with a black outline. ",
                                        "Please note that the minimum number of input genes is five."),
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
                                      br(),
                                      p("By using the protein panel correlation feature, useres can generate a correlation matrix using the proteins of interests in the selected patient cohorts and tissue types.",
                                        "The resulting plot is a bubble plot where the size of the dots are proportional to the negative 10th logarithm of the adjusted P value, the color of the dots represents the correletion coefficient.",
                                        "If the q value is smaller then 0.01 the dot is highlighted with a black outline.",
                                        "Please note that the minimum number of input proteins is five."),
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
                                      br(),
                                      p("By using the mRNA/Protein panel correlation feature, useres can generate a correlation matrix using the proteins and genes of interests in the selected patient cohorts and tissue types.",
                                        "The resulting plot is a bubble plot where the size of the dots are proportional to the negative 10th logarithm of the adjusted P value, the color of the dots represents the correletion coefficient.",
                                        "If the q value is smaller then 0.01 the dot is highlighted with a black outline.",
                                        "Please note that the minimum number of input genes and proteins is five."),
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
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrRna"),color = "#6A2798"),
                                  uiOutput("dataDownLoadGvpCorrRna"),
                                  br(),
                                  uiOutput("GvpCorrRnaDownloadPlotUi")
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrProtein"),color = "#6A2798"),
                                  uiOutput("dataDownLoadGvpCorrProtein"),
                                  br(),
                                  uiOutput("GvpCorrProteinDownloadPlotUi")
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = "outputPlotGvpCorrBoth"),color = "#6A2798"),
                                  uiOutput("dataDownLoadGvpCorrBoth"),
                                  br(),
                                  uiOutput("GvpCorrBothDownloadPlotUi")
                 )
                 
               )
             ))
    
  )

  
    )
=======
                   
                   
)
>>>>>>> ae21f82513a02d11fb0c3f0f2b5f24f22eadf62d
)
