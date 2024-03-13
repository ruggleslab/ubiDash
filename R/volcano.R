volcanoUI <- function(id, cohortSheets) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Volcano Plot"),
      selectInput(ns("cohort"), "Cohort", cohortSheets),

      # input options
      selectizeInput(
        inputId = ns("inputVolcanoGenesTXT"),
        choices = NULL,
        label = "Input genes", options = list(`live-search` = TRUE), multiple = TRUE
      ),
      # p("OR"),
      # fileInput(ns("inputVolcanoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c( ".tsv", "text/tsv",'.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      # selectInput(ns("volcanoDelims"), "Delimiter", c("tab", "comma")),

      # cutoff options
      selectInput(ns("PValCutoff"), "Adjusted P Value Cutoff:", c(0.001, 0.01, 0.05)),
      sliderInput(ns("logFCCutoff"), "Log2(Tumor/Normal) Cutoff:",
        min = .25, max = 2, value = 1, step = 0.05
      ),
      actionButton(ns("generateVolcano"), "Generate Plot"),
      p("Download this data"),
      downloadButton(outputId = ns("downloadVolcanoData"), label = "Download")
    ),
    mainPanel(
      # select which sheet from xlsx to read & generate graph
      shinycssloaders::withSpinner(plotlyOutput(ns("distPlot")), color = "#6A2798")
    )
  )
}


volcanoServer <- function(id, cptac_protein, crl_genes) {
  moduleServer(id, function(input, output, session) {
    filepath <- "./data/DE_protein.xlsx"

    updateSelectizeInput(session, "inputVolcanoGenesTXT",
      choices = as.list(rownames(cptac_protein)),
      server = TRUE
    )

    output$distPlot <- renderPlotly({
      # allow user to save plots

      input$generateVolcano

      inputGenes <- getInputGenes(isolate(input$inputVolcanoGenesCSV), isolate(input$inputVolcanoGenesTXT), delim = isolate(input$volcanoDelims))
      if (length(inputGenes) == 0) {
        inputGenes <- c("")
      }
      # cutoffs
      logFCCutoff <- isolate(input$logFCCutoff)
      PValCutoff <- -log10(as.double(isolate(input$PValCutoff)))

      # prepare dataframe
      data <- read_excel(path = filepath, sheet = isolate(input$cohort))
      data$POI <- "Other"
      data$POI[data$gene %in% crl_genes] <- "CRL"
      data$POI[toupper(data$gene) %in% inputGenes] <- "Potential Substrate"
      data$log.scal.pval <- -log10(data$sca.adj.pval)

      # get significant genes
      labeledPoints <- filter(data, (gene %in% crl_genes & (abs(logFC) > logFCCutoff & abs(log.scal.pval) > PValCutoff)) | toupper(gene) %in% inputGenes)

      p <- ggplot(data, aes(x = logFC, y = log.scal.pval, color = POI)) +
        ylab("-Log10(adjusted p-val)") +
        xlab("Log2(Tumor/Normal)") +
        labs(title = isolate(input$cohort)) +
        geom_point(shape = 16, alpha = .5, mapping = aes(text = paste(
          "Gene: ", gene, "\n",
          "POI: ", POI, "\n",
          "LogFC: ", logFC, "\n",
          "P Value: ", log.scal.pval, "\n",
          sep = ""
        ))) +
        scale_color_manual(breaks = c("CRL", "Potential Substrate", "Other"), values = c("#FF0000", "#0000FF", "#888888")) +
        geom_vline(xintercept = logFCCutoff, linetype = "solid", color = "red") +
        geom_vline(xintercept = -logFCCutoff, linetype = "solid", color = "red") +
        geom_hline(yintercept = PValCutoff, linetype = "solid", color = "red")

      if (nrow(labeledPoints) > 0) {
        p <- p + geom_text(
          hjust = 1, vjust = 1, data = labeledPoints,
          mapping = aes(label = gene, fontface = "bold"),
          color = "#000000"
        )
      }

      ggplotly(p, tooltip = "text") %>% layout(dragmode = FALSE)
    })

    output$downloadVolcanoData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(read_excel(path = filepath, sheet = input$cohort), file)
      }
    )
  })
}
