correlationGeneUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Gene Panel Correlation Analysis"),
      p(tab_summaries$`Pancan Gene Panel`),
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
      tags$small(
        "Note: Normal samples are limited in some cohorts. If analysis fails, try adding more cohorts or switching to Tumor.",
        style = "color: #6c757d; display: block; margin-bottom: 10px;"
      ),
      actionButton(ns("generateGvpCorrProtein"), label = "Generate Plot")
    ),
    mainPanel(
      withSpinner(plotlyOutput(outputId = ns("outputPlotGvpCorrProtein")), color = "#6A2798"),
      downloadButton(ns("downLoadGvpCorrProtein"), label = "Download correlation data"),
      downloadButton(ns("downloadGvpCorrProteinPDF"), label = "Download PDF")
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

    resGvpCorrProtein <- eventReactive(input$generateGvpCorrProtein,
      {
        cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
        if (length(isolate(input$inputGvpCorrProteinCohort)) < 1) {
          stop("Please select at least one cohort.")
        }
        if (length(isolate(input$inputGvpCorrProteinGeneName)) < 1) {
          stop("Please select at least one gene")
        }

        # Check sample availability for Normal tissue type
        if (input$inputGvpCorrProteinStype == "Normal") {
          # Note: column is "cohort" not "cohor" - matching corrAnalysis.R which uses "cohor"
          cptac_samples_filt <- cptac_samples[cptac_samples$cohort %in% input$inputGvpCorrProteinCohort, ]
          normal_count <- sum(cptac_samples_filt$type == "Normal", na.rm = TRUE)
          if (normal_count < 10) {
            # Get cohorts that have Normal samples
            cohorts_with_normal <- unique(cptac_samples$cohort[cptac_samples$type == "Normal"])
            stop(paste0(
              "Not enough Normal samples (", normal_count, " found, minimum 10 required). ",
              "Cohorts with Normal samples: ", paste(cohorts_with_normal, collapse = ", "), ". ",
              "Try selecting these cohorts or switch to 'Tumor' sample type."
            ))
          }
        }

        panelCorr(
          geneNamesProtein = input$inputGvpCorrProteinGeneName,
          geneNamesRNA     = NULL,
          tissueType       = input$inputGvpCorrProteinStype,
          cohort           = input$inputGvpCorrProteinCohort,
          cptac_protein    = cptac_protein,
          cptac_rna        = NULL,
          cptac_samples    = cptac_samples
        )
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

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

    output$downloadGvpCorrProteinPDF <- downloadHandler(
      filename = function() {
        paste0("GenePanelProtein-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        # Recreate the static ggplot for PDF export
        cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
        res <- panelCorr(
          geneNamesProtein = isolate(input$inputGvpCorrProteinGeneName),
          geneNamesRNA     = NULL,
          tissueType       = isolate(input$inputGvpCorrProteinStype),
          cohort           = isolate(input$inputGvpCorrProteinCohort),
          cptac_protein    = cptac_protein,
          cptac_rna        = NULL,
          cptac_samples    = cptac_samples
        )
        # Extract the ggplot object before ggplotly conversion
        pdf(file, width = 10, height = 8)
        print(res$corrPlot)
        dev.off()
      }
    )
  })
}
