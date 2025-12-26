oncoplotUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Oncoplot"),

      # gene inputs
      p(tab_summaries$Oncoplot),
      tags$small(
        "Genes are ranked by mutation frequency. Use custom gene set to focus on specific UPS components.",
        style = "color: #6c757d; display: block; margin-bottom: 10px;"
      ),
      selectInput(
        inputId = ns("geneSet"),
        label = "Default Gene Set",
        choices = c("Top E3 Ligases" = "e3", "Top UPS Genes" = "ups"),
        selected = "e3"
      ),
      sliderInput(ns("maxOncoGenes"), "Max Number of Genes to Display:",
        min = 1, max = 15,
        value = 15
      ),
      selectizeInput(
        inputId = ns("genes"),
        choices = NULL,
        label = "Custom Gene Set (overrides default)",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      # p("OR"),
      # fileInput(ns("inputOncoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      # selectInput(ns("oncoDelims"), "Delimiter", c("tab", "comma")),
      #
      # options
      selectizeInput(
        inputId = ns("cohorts"),
        choices = cohorts,
        selected = cohorts,
        label = "Select cohorts",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      checkboxInput(ns("sortByCohort"), "Sort by Cohort", value = FALSE),
      actionButton(ns("refresh"), "Generate Plot", class = "btn-primary"),

      # download options

      # downloadButton(outputId = ns("downloadOncoCohorts"), label = "Download Cohort Data"),
    ),
    mainPanel(
      textOutput(ns("error")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("mutationOncoplot")), color = "#6A2798"),
      br(),
      downloadButton(outputId = ns("downloadOncoMutations"), label = "Download Data"),
      downloadButton(outputId = ns("downloadOncoPDF"), label = "Download PDF"),
    )
  )
}


oncoplotServer <- function(id, top_e3_genes, top_ups_genes, cohorts, colorsTbl) {
  moduleServer(id, function(input, output, session) {
    # Use pre-loaded MAF data from global.R
    all_genes <- unique(getGeneSummary(pancan_maf)$Hugo_Symbol)
    updateSelectizeInput(session, "genes",
      choices = all_genes,
      server = TRUE
    )

    inputCohorts <- reactiveVal()

    filteredMaf <- eventReactive(input$refresh,
      {
        if (length(isolate(input$genes)) == 1) {
          stop("Two genes are needed to generate plot.")
        }

        # Use pre-loaded MAF data from global.R
        inputCohorts(isolate(input$cohorts))

        filteredMaf <- pancan_maf
        if (length(inputCohorts()) != 0 && length(inputCohorts()) != length(cohorts)) {
          filteredMaf <- subsetMaf(maf, tsb = maf@clinical.data[cohort_code %in% inputCohorts()]$Tumor_Sample_Barcode)
        }

        return(filteredMaf)
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )

    output$mutationOncoplot <- renderPlot({
      # cptac_samples <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")

      inputGenes <- isolate(input$genes)

      if (!isTruthy(inputGenes)) {
        # Use selected gene set (E3 or UPS)
        selectedGeneSet <- if (isolate(input$geneSet) == "ups") top_ups_genes else top_e3_genes
        inputGenes <- selectedGeneSet[1:isolate(input$maxOncoGenes)]
      }

      colors <- colorsTbl$color
      names(colors) <- colorsTbl$id

      colors <- colors[names(colors) %in% inputCohorts()]

      # vc_colors = c("#EF476F", "#F78C6B", "#FFD166", "#83D483", "#06D6A0", "#0CB0A9", "#118AB2","#073B4C", "#000000")
      vc_colors <- c("#0AA0BF", "#F6114A", "#F36E98", "#78B177", "#9862A2", "#F05006", "#25998F", "#FCA00C", "#000000")
      names(vc_colors) <- c(
        "Missense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Nonsense_Mutation",
        "Splice_Site",
        "In_Frame_Ins",
        "In_Frame_Del",
        "Nonstop_Mutation",
        "Multi_Hit"
      )

      return(oncoplot(
        maf = filteredMaf(),
        # top = isolate(input$maxOncoGenes),
        genes = inputGenes,
        annotationColor = list(cohort_code = colors),
        sortByAnnotation = isolate(input$sortByCohort),
        gene_mar = 8,
        # barcode_mar = 20,
        fontSize = 1,
        annotationFontSize = 1.5,
        legendFontSize = 1.5,
        colors = vc_colors,
        removeNonMutated = T, # this is new
        # drawColBar = F,
        # cBioPortal = T,
        clinicalFeatures = "cohort_code"
      ))
    })
    output$downloadOncoMutations <- downloadHandler(
      filename = function() {
        if (length(isolate(input$cohorts)) == length(cohorts) && !isTruthy(isolate(input$genes))) {
          paste("mutations-oncoplot-", Sys.Date(), ".csv", sep = "")
        } else {
          paste("mutations-oncoplot-custom", Sys.Date(), ".csv", sep = "")
        }
      },
      content = function(file) {
        write.csv(filteredMaf()@data, file)
      }
    )

    output$downloadOncoPDF <- downloadHandler(
      filename = function() {
        paste("oncoplot-", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        inputGenes <- isolate(input$genes)

        if (!isTruthy(inputGenes)) {
          selectedGeneSet <- if (isolate(input$geneSet) == "ups") top_ups_genes else top_e3_genes
          inputGenes <- selectedGeneSet[1:isolate(input$maxOncoGenes)]
        }

        colors <- colorsTbl$color
        names(colors) <- colorsTbl$id
        colors <- colors[names(colors) %in% inputCohorts()]

        vc_colors <- c("#0AA0BF", "#F6114A", "#F36E98", "#78B177", "#9862A2", "#F05006", "#25998F", "#FCA00C", "#000000")
        names(vc_colors) <- c(
          "Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
          "Splice_Site", "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Multi_Hit"
        )

        pdf(file, width = 12, height = 10)
        oncoplot(
          maf = filteredMaf(),
          genes = inputGenes,
          annotationColor = list(cohort_code = colors),
          sortByAnnotation = isolate(input$sortByCohort),
          gene_mar = 8,
          fontSize = 1,
          annotationFontSize = 1.5,
          legendFontSize = 1.5,
          colors = vc_colors,
          removeNonMutated = TRUE,
          clinicalFeatures = "cohort_code"
        )
        dev.off()
      }
    )
  })
}
