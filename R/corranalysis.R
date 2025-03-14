correlationGeneGeneUI <- function(id, cohorts) {
  ns <- NS(id)
  

    sidebarLayout(
      sidebarPanel(
        h3("Gene vs. Gene Correlation Analysis"),
        p("(description paragraph)"),
        tabsetPanel(
          id = ns("GvGpanel"), selected = "mRNA",
          tabPanel(
            title = "mRNA",
            selectizeInput(
              inputId = ns("inputGvgCorrRnaNameA"),
              choices = NULL,
              label = "Choose protein A", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrRnaNameB"),
              choices = NULL,
              label = "Choose protein B", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrRnaCohort"),
              choices = cohorts, selected = "BRCA",
              label = "Select cohort", 
              options = list(`live-search` = TRUE), 
              multiple = TRUE
            ),
            radioButtons(
              inputId = ns("inputGvgCorrRnaStype"),
              label = "Choose sample type", 
              choices = c("Normal", "Tumor"), 
              selected = "Tumor"
            ),
            actionButton(ns("generateGvgCorrRna"), label = "Start Analysis")
          ),
          tabPanel(
            title = "Protein",
            selectizeInput(
              inputId = ns("inputGvgCorrProteinNameA"),
              choices = NULL,
              label = "Choose protein A", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrProteinNameB"),
              choices = NULL,
              label = "Choose protein B", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrProteinCohort"),
              choices = cohorts, selected = "BRCA",
              label = "Select cohort", 
              options = list(`live-search` = TRUE), 
              multiple = TRUE
            ),
            radioButtons(
              inputId = ns("inputGvgCorrProteinStype"),
              label = "Choose sample type", 
              choices = c("Normal", "Tumor"), 
              selected = "Tumor"
            ),
            actionButton(ns("generateGvgCorrProtein"), label = "Start Analysis")
          ),
          tabPanel(
            title = "mRNA/Protein",
            selectizeInput(
              inputId = ns("inputGvgCorrBothNameA"),
              choices = NULL,
              label = "Choose protein", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrBothNameB"),
              choices = NULL,
              label = "Choose gene", 
              options = list(`live-search` = TRUE), 
              multiple = FALSE
            ),
            selectizeInput(
              inputId = ns("inputGvgCorrBothCohort"),
              choices = cohorts, selected = "BRCA",
              label = "Select cohort", 
              options = list(`live-search` = TRUE), 
              multiple = TRUE
            ),
            radioButtons(
              inputId = ns("inputGvgCorrBothStype"),
              label = "Choose sample type", 
              choices = c("Normal", "Tumor"), 
              selected = "Tumor"
            ),
            actionButton(ns("generateGvgCorrBoth"), label = "Start Analysis")
          )
        )
      ),
      mainPanel(
        conditionalPanel(
          condition = paste0("input['", ns("GvGpanel"), "\'] == 'mRNA'"),
          withSpinner(plotlyOutput(outputId = ns("outputPlotGvgCorrRna")), color = "#6A2798"),
          tableOutput(outputId = ns("outputTabGvgCorrRna")),
          downloadButton(ns("downLoadGvgCorrRna"), label = "Download correlation data")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("GvGpanel"), "\'] == 'Protein'"),
          withSpinner(plotlyOutput(outputId = ns("outputPlotGvgCorrProtein")), color = "#6A2798"),
          tableOutput(outputId = ns("outputTabGvgCorrProtein")),
          downloadButton(ns("downLoadGvgCorrProtein"), label = "Download correlation data")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("GvGpanel"), "\'] == 'mRNA/Protein'"),
          withSpinner(plotlyOutput(outputId = ns("outputPlotGvgCorrBoth")), color = "#6A2798"),
          tableOutput(outputId = ns("outputTabGvgCorrBoth")),
          downloadButton(ns("downLoadGvgCorrBoth"), label = "Download correlation data")
        )
      )
  )
}

correlationGeneGeneServer <- function(id, cohorts, colorsTbl) {
  moduleServer(id, function(input, output, session) {
    # Load required functions and data
    source(file = "./scripts/corrAnalysis.R")
    cptac_samples <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")
    
    ### mRNA Correlation ###
    updateSelectizeInput(
      session, "inputGvgCorrRnaNameA",
      choices = cptac_rna_genes, server = TRUE, selected = "BTRC",
      options = list(placeholder = "Type gene name")
    )
    updateSelectizeInput(
      session, "inputGvgCorrRnaNameB",
      choices = cptac_rna_genes, server = TRUE, selected = "CTNNB1",
      options = list(placeholder = "Type gene name")
    )
    
    resGvgCorrRna <- eventReactive(input$generateGvgCorrRna, {
      validate(
        need(length(isolate(input$inputGvgCorrRnaCohort)) >= 1, 'Please select at least one cohort.')
      )
      cptac_rna <- mem_get_from_drive("cptac_transformed_rna_expression_v1-1.rds")
      singleCorr(
        geneNameA     = input$inputGvgCorrRnaNameA,
        geneNameB     = input$inputGvgCorrRnaNameB,
        cohort        = input$inputGvgCorrRnaCohort,
        tissueType    = input$inputGvgCorrRnaStype,
        cptac_protein = NULL,
        cptac_rna     = cptac_rna,
        cptac_samples = cptac_samples,
        colorsTbl     = colorsTbl
      )
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    output$outputPlotGvgCorrRna <- renderPlotly({
      resGvgCorrRna()$singCorrPlot
    })
    output$outputTabGvgCorrRna <- renderTable({
      resGvgCorrRna()$singCorrTab
    }, rownames = TRUE, colnames = FALSE)
    
    
    ### Protein Correlation ###
    updateSelectizeInput(
      session, "inputGvgCorrProteinNameA",
      choices = cptac_protein_genes, server = TRUE, selected = "VHL",
      options = list(placeholder = "Type gene name")
    )
    updateSelectizeInput(
      session, "inputGvgCorrProteinNameB",
      choices = cptac_protein_genes, server = TRUE, selected = "HIF1A",
      options = list(placeholder = "Type gene name")
    )
    
    resGvgCorrProtein <- eventReactive(input$generateGvgCorrProtein, {
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
      validate(
        need(length(isolate(input$inputGvgCorrProteinCohort)) >= 1, 'Please select at least one cohort.')
      )
      singleCorr(
        geneNameA     = input$inputGvgCorrProteinNameA,
        geneNameB     = input$inputGvgCorrProteinNameB,
        cohort        = input$inputGvgCorrProteinCohort,
        tissueType    = input$inputGvgCorrProteinStype,
        cptac_protein = cptac_protein,
        cptac_rna     = NULL,
        cptac_samples = cptac_samples,
        colorsTbl     = colorsTbl
      )
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    output$outputPlotGvgCorrProtein <- renderPlotly({
      resGvgCorrProtein()$singCorrPlot
    })
    output$outputTabGvgCorrProtein <- renderTable({
      resGvgCorrProtein()$singCorrTab
    }, rownames = TRUE, colnames = FALSE)
    
    
    ### mRNA/Protein Correlation ###
    updateSelectizeInput(
      session, "inputGvgCorrBothNameA",
      choices = as.list(intersect(cptac_protein_genes, cptac_rna_genes)),
      server = TRUE, selected = "BTRC",
      options = list(placeholder = "Type Protein name")
    )
    updateSelectizeInput(
      session, "inputGvgCorrBothNameB",
      choices = as.list(intersect(cptac_protein_genes, cptac_rna_genes)),
      server = TRUE, selected = "CTNNB1",
      options = list(placeholder = "Type gene name")
    )
    
    resGvgCorrBoth <- eventReactive(input$generateGvgCorrBoth, {
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
      cptac_rna <- mem_get_from_drive("cptac_transformed_rna_expression_v1-1.rds")
      validate(
        need(length(isolate(input$inputGvgCorrBothCohort)) >= 1, 'Please select at least one cohort.')
      )
      singleCorr(
        geneNameA     = input$inputGvgCorrBothNameA,
        geneNameB     = input$inputGvgCorrBothNameB,
        cohort        = input$inputGvgCorrBothCohort,
        tissueType    = input$inputGvgCorrBothStype,
        cptac_protein = cptac_protein,
        cptac_rna     = cptac_rna,
        cptac_samples = cptac_samples,
        colorsTbl     = colorsTbl
      )
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    output$outputPlotGvgCorrBoth <- renderPlotly({
      resGvgCorrBoth()$singCorrPlot
    })
    output$outputTabGvgCorrBoth <- renderTable({
      resGvgCorrBoth()$singCorrTab
    }, rownames = TRUE, colnames = FALSE)
    
    
    ### Download Handlers ###
    output$downLoadGvgCorrRna <- downloadHandler(
      filename = function() {
        paste0("GeneVsGeneRNA-custom-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(resGvgCorrRna()$corTab, file, row.names = TRUE)
      }
    )
    output$downLoadGvgCorrProtein <- downloadHandler(
      filename = function() {
        paste0("GeneVsGeneProtein-custom-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(resGvgCorrProtein()$corTab, file, row.names = TRUE)
      }
    )
    output$downLoadGvgCorrBoth <- downloadHandler(
      filename = function() {
        paste0("GeneVsGeneProteinRNA-custom-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(resGvgCorrBoth()$corTab, file, row.names = TRUE)
      }
    )
  })
}