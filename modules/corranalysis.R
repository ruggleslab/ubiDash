corrAnalysisUI <- function(id, cohorts){
  ns <- NS(id)
  
  navbarMenu("Correlation Analysis",
    tabPanel("Gene vs Gene",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = ns("GvGpanel"),selected = "mRNA",
                             tabPanel("mRNA",
                                      selectizeInput(inputId =ns("inputGvgCorrRnaNameA"),
                                                                  choices = NULL,
                                                     label = "Choose protein A",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = ns("inputGvgCorrRnaNameB"),
                                                     choices = NULL,
                                                     label = "Choose protein B",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      
                                      selectizeInput(inputId = ns("inputGvgCorrRnaCohort"),
                                                     choices = cohorts,selected = "BRCA",
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvgCorrRnaStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvgCorrRna"), label = "Start Analysis")
                             ),
                             tabPanel("Protein",
                                      selectizeInput(inputId = ns("inputGvgCorrProteinNameA"),
                                                     choices = NULL,
                                                     label = "Choose protein A",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = ns("inputGvgCorrProteinNameB"),
                                                     choices = NULL,
                                                     label = "Choose protein B",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      
                                      selectizeInput(inputId = ns("inputGvgCorrProteinCohort"),
                                                     choices = cohorts,selected = "BRCA",
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvgCorrProteinStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvgCorrProtein"), label = "Start Analysis")
                             ),
                             tabPanel("mRNA/Protein",
                                      selectizeInput(inputId = ns("inputGvgCorrBothNameA"),
                                                     choices = NULL,
                                                     label = "Choose protein",options = list(`live-search` = TRUE),  multiple = F),
                                      selectizeInput(inputId = ns("inputGvgCorrBothNameB"),
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = ns("inputGvgCorrBothCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",selected = "BRCA",
                                                     options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvgCorrBothStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvgCorrBoth"), label = "Start Analysis")
                             ))
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvGpanel == 'mRNA'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvgCorrRna")),color = "#6A2798"),
                                  tableOutput(outputId = ns("outputTabGvgCorrRna")),
                                  uiOutput(outputId = ns("downLoadButtonGvgCorrRna"))
                 ),
                 conditionalPanel(condition = "input.GvGpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvgCorrProtein")),color = "#6A2798"),
                                  tableOutput(outputId = ns("outputTabGvgCorrProtein"))
                 ),
                 conditionalPanel(condition = "input.GvGpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvgCorrBoth")),color = "#6A2798"),
                                  tableOutput(outputId = ns("outputTabGvgCorrBoth"))
                 )
                 
               )
             )),
    tabPanel("Gene vs All",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvApanel", selected = "mRNA",
                             tabPanel(title = "mRNA",
                                      selectizeInput(inputId = ns("inputGvACorrRnaGeneName"),
                                                     choices = NULL,
                                                     label = "Choose a gene",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = ns("inputGvACorrRnaCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvACorrRnaType"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      sliderInput(inputId = ns("inputGvACorrRnaR"),
                                                  label = "Set correlation coefficient cut-off",
                                                  value = c(0.5,1), min = -1, max = 1, step = 0.05),
                                      
                                      actionButton(ns("generateGvACorrRna"), label = "Start Analysis")
                             ),
                             tabPanel(title = "Protein",
                                      selectizeInput(inputId = ns("inputGvACorrProteinGeneName"),
                                                     choices = NULL,
                                                     label = "Choose a protein",options = list(`live-search` = TRUE),  multiple = F),
                                      
                                      selectizeInput(inputId = ns("inputGvACorrProteinCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvACorrProteinType"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      sliderInput(inputId = ns("inputGvACorrProteinR"),
                                                  label = "Set correlation coefficient cut-off",
                                                  value = c(0.5,1), min = -1, max = 1, step = 0.05),
                                      
                                      actionButton(ns("generateGvACorrProtein"), label = "Start Analysis")
                             )
                 )
                 
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvApanel == 'Protein'",
                                  fluidRow(column(width = 4, 
                                                  withSpinner(dataTableOutput(outputId = ns("outputTabGvACorrProtein")),color = "#6A2798")),
                                           column(width = 8, 
                                                  withSpinner(plotOutput(outputId = ns("outputPlotGvACorrProtein")),color = "#6A2798")))
                 ),
                 conditionalPanel(condition = "input.GvApanel == 'mRNA'",
                                  fluidRow(column(width = 4, 
                                                  withSpinner(dataTableOutput(outputId = ns("outputTabGvACorrRna")),color = "#6A2798")),
                                           column(width = 8, 
                                                  withSpinner(plotOutput(outputId = ns("outputPlotGvACorrRna")),color = "#6A2798")))
                 )
               )
             )
    ),
    tabPanel("Gene Panel",
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(id = "GvPpanel",selected = "mRNA",
                             tabPanel("mRNA",
                                      selectizeInput(inputId = ns("inputGvpCorrRnaGeneName"),
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = ns("inputGvpCorrRnaCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvpCorrRnaStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvpCorrRna"), label = "Start Analysis")
                             ),
                             tabPanel("Protein",
                                      selectizeInput(inputId = ns("inputGvpCorrProteinGeneName"),
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = ns("inputGvpCorrProteinCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvpCorrProteinStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvpCorrProtein"), label = "Start Analysis")
                             ),
                             tabPanel("mRNA/Protein",
                                      selectizeInput(inputId = ns("inputGvpCorrBothGeneNameA"),
                                                     choices = NULL,
                                                     label = "Choose gene",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      selectizeInput(inputId = ns("inputGvpCorrBothGeneNameB"),
                                                     choices = NULL,
                                                     label = "Choose protein",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      selectizeInput(inputId = ns("inputGvpCorrBothCohort"),
                                                     choices = cohorts,
                                                     label = "Select cohort",options = list(`live-search` = TRUE),  multiple = TRUE),
                                      
                                      radioButtons(inputId = ns("inputGvpCorrBothStype"), label = "Choose sample type", choices = c("Normal", "Tumor"), selected = "Tumor"),
                                      
                                      actionButton(ns("generateGvpCorrBoth"), label = "Start Analysis")
                             ))
               ),
               mainPanel(
                 conditionalPanel(condition = "input.GvPpanel == 'mRNA'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvpCorrRna")),color = "#6A2798"),
                                  uiOutput(ns("dataDownLoadGvpCorrRna"))
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'Protein'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvpCorrProtein")),color = "#6A2798"),
                                  uiOutput(ns("dataDownLoadGvpCorrProtein"))
                 ),
                 conditionalPanel(condition = "input.GvPpanel == 'mRNA/Protein'",
                                  withSpinner(plotOutput(outputId = ns("outputPlotGvpCorrBoth")),color = "#6A2798"),
                                  uiOutput(ns("dataDownLoadGvpCorrBoth"))
                 )
                 
               )
             )))
}


corrAnalysisServer <- function(id, cohorts, cptac_protein,cptac_rna,cptac_samples) {
  moduleServer(id,function(input,output,session){
    # protein specific correlation ----
    source(file = "./RScripts/corrAnalysis.R")
    
    #Gene vs Gene correlation RNA
    updateSelectizeInput(session, 'inputGvgCorrRnaNameA', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
    updateSelectizeInput(session, 'inputGvgCorrRnaNameB', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
    
    resGvgCorrRna <- eventReactive(input$generateGvgCorrRna,{singleCorr(geneNameA = input$inputGvgCorrRnaNameA,
                                                                        geneNameB = input$inputGvgCorrRnaNameB,
                                                                        cohort = input$inputGvgCorrRnaCohort,
                                                                        tissueType = input$inputGvgCorrRnaStype, 
                                                                        cptac_protein = NULL,
                                                                        cptac_rna = cptac_rna,
                                                                        cptac_samples = cptac_samples)})
    
    output$outputPlotGvgCorrRna <- renderPlot({resGvgCorrRna()$singCorrPlot})
    output$outputTabGvgCorrRna <- renderTable({resGvgCorrRna()$singCorrTab},rownames = T, colnames = F)
    
    #Protein vs Protein correlation RNA
    updateSelectizeInput(session, 'inputGvgCorrProteinNameA', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
    updateSelectizeInput(session, 'inputGvgCorrProteinNameB', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
    
    resGvgCorrProtein <- eventReactive(input$generateGvgCorrProtein,{singleCorr(geneNameA = input$inputGvgCorrProteinNameA,
                                                                                geneNameB = input$inputGvgCorrProteinNameB,
                                                                                cohort = input$inputGvgCorrProteinCohort,
                                                                                tissueType = input$inputGvgCorrProteinStype, 
                                                                                cptac_protein = cptac_protein,
                                                                                cptac_rna = NULL, 
                                                                                cptac_samples = cptac_samples)})
    
    output$outputPlotGvgCorrProtein <- renderPlot({resGvgCorrProtein()$singCorrPlot})
    output$outputTabGvgCorrProtein <- renderTable({resGvgCorrProtein()$singCorrTab},rownames = T, colnames = F)
    
    #Protein vs Gene correlation
    updateSelectizeInput(session, 'inputGvgCorrBothNameA', choices = as.list(intersect(rownames(cptac_protein),rownames(cptac_rna))), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type Protein name'))
    updateSelectizeInput(session, 'inputGvgCorrBothNameB', choices = as.list(intersect(rownames(cptac_protein),rownames(cptac_rna))), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
    
    resGvgCorrBoth <- eventReactive(input$generateGvgCorrBoth,{singleCorr(geneNameA = input$inputGvgCorrBothNameA,
                                                                          geneNameB = input$inputGvgCorrBothNameB,
                                                                          cohort = input$inputGvgCorrBothCohort,
                                                                          tissueType = input$inputGvgCorrBothStype, 
                                                                          cptac_protein = cptac_protein,
                                                                          cptac_rna = cptac_rna,
                                                                          cptac_samples = cptac_samples)})
    
    output$outputPlotGvgCorrBoth <- renderPlot({resGvgCorrBoth()$singCorrPlot})
    output$outputTabGvgCorrBoth <- renderTable({resGvgCorrBoth()$singCorrTab},rownames = T, colnames = F)
    
    #Gene vs All correlation for proteins
    updateSelectizeInput(session, 'inputGvACorrProteinGeneName', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
    resGvACorrProtein <- eventReactive(input$generateGvACorrProtein,{multiCorr(geneName = input$inputGvACorrProteinGeneName,
                                                                               cutOff =  input$inputGvACorrProteinR, 
                                                                               cohort = input$inputGvACorrProteinCohort,
                                                                               tissueType =input$inputGvACorrProteinType,
                                                                               cptac_protein = cptac_protein,
                                                                               cptac_rna = NULL,
                                                                               cptac_samples = cptac_samples
    )
    })
    
    output$outputTabGvACorrProtein <- DT::renderDataTable(server = FALSE,{resGvACorrProtein()$result},rownames = F,extensions = 'Buttons', options = list(
      dom = 'frtBip',#Bfrtip '<<t>B>' B - buttons, f - search panel r - processing display element t - the table itself i - number of items p - paging
      buttons = list(
        list(extend = "excel", text = "Download results", filename = paste0(input$inputGvACorrProteinGeneName, "_ProteinCorrelations"),
             exportOptions = list(
               modifier = list(page = "all")
             )))))
    
    output$outputPlotGvACorrProtein <- renderPlot({resGvACorrProtein()$corrPlot})
    
    #Gene vs All correlation for RNA
    updateSelectizeInput(session, 'inputGvACorrRnaGeneName', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
    resGvACorrRna <- eventReactive(input$generateGvACorrRna,{multiCorr(geneName = input$inputGvACorrRnaGeneName,
                                                                       cutOff =  input$inputGvACorrRnaR,
                                                                       tissueType =input$inputGvACorrRnaType, 
                                                                       cohort = input$inputGvACorrRnaCohort,
                                                                       cptac_protein = NULL,
                                                                       cptac_rna = cptac_rna,
                                                                       cptac_samples = cptac_samples
    )
    })
    
    output$outputTabGvACorrRna <- DT::renderDataTable(server = FALSE,{resGvACorrRna()$result},rownames = F,extensions = 'Buttons', options = list(
      dom = 'frtBip',#Bfrtip '<<t>B>' B - buttons, f - search panel r - processing display element t - the table itself i - number of items p - paging
      buttons = list(
        list(extend = "excel", text = "Download results", filename = paste0(input$inputGvACorrRnaGeneName, "_RNAcorrelations"),
             exportOptions = list(
               modifier = list(page = "all")
             )))))
    
    output$outputPlotGvACorrRna <- renderPlot({resGvACorrRna()$corrPlot})
    
    #Gene panel correlation for RNA
    updateSelectizeInput(session, 'inputGvpCorrRnaGeneName', choices = as.list(rownames(cptac_rna)), server = TRUE,  options = list(placeholder = 'Type gene name'))
    
    resGvpCorrRna <- eventReactive(input$generateGvpCorrRna,{panelCorr(geneNamesRNA = input$inputGvpCorrRnaGeneName,
                                                                       geneNamesProtein = NULL, 
                                                                       tissueType = input$inputGvpCorrRnaStype, 
                                                                       cohort = input$inputGvpCorrRnaCohort, 
                                                                       cptac_protein = NULL, 
                                                                       cptac_rna = cptac_rna, 
                                                                       cptac_samples = cptac_samples)})
    
    output$outputPlotGvpCorrRna <- renderPlot({resGvpCorrRna()$corrPlot})
    
    output$dataDownLoadGvpCorrRna <- renderUI({
      req(input$generateGvpCorrRna)
      downloadButton("downLoadGvpCorrRna", label = "Download correlation data") })
    
    output$downLoadGvpCorrRna <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                                 content = function(file) {
                                                   write.xlsx(resGvpCorrRna()$corTab, file,row.names = TRUE)})
    
    
    #Gene panel correlation for Protein
    updateSelectizeInput(session, 'inputGvpCorrProteinGeneName', choices = as.list(rownames(cptac_protein)), server = TRUE, options = list(placeholder = 'Type gene name'))
    
    resGvpCorrProtein <- eventReactive(input$generateGvpCorrProtein,{panelCorr(geneNamesProtein = input$inputGvpCorrProteinGeneName,
                                                                               geneNamesRNA = NULL, 
                                                                               tissueType = input$inputGvpCorrProteinStype, 
                                                                               cohort = input$inputGvpCorrProteinCohort, 
                                                                               cptac_protein = cptac_protein, 
                                                                               cptac_rna = NULL, 
                                                                               cptac_samples = cptac_samples)})
    
    output$outputPlotGvpCorrProtein <- renderPlot({resGvpCorrProtein()$corrPlot})
    
    output$dataDownLoadGvpCorrProtein <- renderUI({
      req(input$generateGvpCorrProtein)
      downloadButton("downLoadGvpCorrProtein", label = "Download correlation data") })
    
    output$downLoadGvpCorrProtein <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                                     content = function(file) {
                                                       write.xlsx(resGvpCorrProtein()$corTab, file,row.names = TRUE)})
    
    #Gene panel correlation for Protein/RNA
    updateSelectizeInput(session, 'inputGvpCorrBothGeneNameA', choices = as.list(rownames(cptac_rna)), server = TRUE, options = list(placeholder = 'Type gene name'))
    updateSelectizeInput(session, 'inputGvpCorrBothGeneNameB', choices = as.list(rownames(cptac_protein)), server = TRUE, options = list(placeholder = 'Type protein name'))
    
    resGvpCorrBoth <- eventReactive(input$generateGvpCorrBoth,{panelCorr(geneNamesRNA = input$inputGvpCorrBothGeneNameA,
                                                                         geneNamesProtein = input$inputGvpCorrBothGeneNameB, 
                                                                         tissueType = input$inputGvpCorrBothStype, 
                                                                         cohort = input$inputGvpCorrBothCohort, 
                                                                         cptac_protein = cptac_protein, 
                                                                         cptac_rna = cptac_rna, 
                                                                         cptac_samples = cptac_samples)})
    
    output$outputPlotGvpCorrBoth <- renderPlot({resGvpCorrBoth()$corrPlot})
    
    output$dataDownLoadGvpCorrBoth <- renderUI({
      req(input$generateGvpCorrBoth)
      downloadButton("downLoadGvpCorrBoth", label = "Download correlation data") })
    
    output$downLoadGvpCorrBoth <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                                  content = function(file) {
                                                    write.xlsx(resGvpCorrBoth()$corTab, file,row.names = TRUE)})
    
  })
  
  }