correlationGeneUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Gene Panel Correlation Analysis"),
      p("(description paragraph)"),
      selectizeInput(
        inputId = ns("inputGvpCorrProteinGeneName"),
        choices = NULL,
        label = "Choose gene", 
        options = list(`live-search` = TRUE), 
        multiple = TRUE
      ),
      selectizeInput(
        inputId = ns("inputGvpCorrProteinCohort"),
        choices = cohorts, selected = "BRCA",
        label = "Select cohort", 
        options = list(`live-search` = TRUE), 
        multiple = TRUE
      ),
      radioButtons(
        inputId = ns("inputGvpCorrProteinStype"),
        label = "Choose sample type", 
        choices = c("Normal", "Tumor"), 
        selected = "Tumor"
      ),
      actionButton(ns("generateGvpCorrProtein"), label = "Start Analysis")
    ),
    mainPanel(
      withSpinner(plotlyOutput(outputId = ns("outputPlotGvpCorrProtein")), color = "#6A2798"),
      downloadButton(ns("downLoadGvpCorrProtein"), label = "Download correlation data")
    )
  )
}
correlationGeneServer <- function(id, cohorts, colorsTbl) {
  moduleServer(id, function(input, output, session) {
    # Load required functions and data
    source(file = "./scripts/corrAnalysis.R")
    cptac_samples <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")
    
    ### Gene Panel: Protein ###
    updateSelectizeInput(
      session, "inputGvpCorrProteinGeneName",
      selected = c("VHL", "A1BG", "A2M", "A2ML1", "A4GALT"),
      choices = cptac_protein_genes, server = TRUE,
      options = list(placeholder = "Type gene name")
    )
    
    resGvpCorrProtein <- eventReactive(input$generateGvpCorrProtein, {
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
      validate(
        need(length(isolate(input$inputGvpCorrProteinCohort)) >= 1, 'Please select at least one cohort.'),
        need(length(isolate(input$inputGvpCorrProteinGeneName)) >= 1, 'Please select at least one gene')
      )
      panelCorr(
        geneNamesProtein = input$inputGvpCorrProteinGeneName,
        geneNamesRNA     = NULL,
        tissueType       = input$inputGvpCorrProteinStype,
        cohort           = input$inputGvpCorrProteinCohort,
        cptac_protein    = cptac_protein,
        cptac_rna        = NULL,
        cptac_samples    = cptac_samples
      )
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    output$outputPlotGvpCorrProtein <- renderPlotly({
      resGvpCorrProtein()$corrPlot
    })
    
    output$downLoadGvpCorrProtein <- downloadHandler(
      filename = function() {
        paste0("GenePanelProtein-custom-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(resGvpCorrProtein()$corTab, file, row.names = TRUE)
      }
    )
    
  })
}