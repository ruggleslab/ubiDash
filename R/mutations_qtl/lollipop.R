lollipopUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Lollipop Plot"),
      p(tab_summaries$Lollipop),
      selectizeInput(
        inputId = ns("gene"),
        choices = NULL,
        multiple = FALSE,
        label = "Gene",
        options = list(`live-search` = TRUE),
      ),
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("lollipop")),
        color = "#6A2798"
      ),
      br(),
      downloadButton(ns("downloadData"), "Download Data (CSV)"),
      downloadButton(ns("downloadPlot"), "Download Plot (PDF)")
    )
  )
}

lollipopServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    updateSelectizeInput(session, "gene",
      choices = maf_genes,
      server = TRUE, selected = "VHL"
    )

    output$lollipop <- renderPlot({
      # Use pre-loaded MAF data from global.R
      gene <- input$gene

      return(lollipopPlot(
        maf = pancan_maf,
        gene = gene,
        AACol = "Protein_Change",
        showMutationRate = TRUE,
        roundedRect = FALSE,
        # labelPos = 882
      ))
    })

    # Download handler for the plot (PDF format)
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(input$gene, "_lollipop_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        # Use pre-loaded MAF data from global.R
        gene <- input$gene

        pdf(file, width = 10, height = 6)
        lollipopPlot(
          maf = pancan_maf,
          gene = gene,
          AACol = "Protein_Change",
          showMutationRate = TRUE,
          roundedRect = FALSE
        )
        dev.off()
      }
    )

    # Download handler for the data (CSV format)
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0(input$gene, "_lollipop_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Use pre-loaded MAF data from global.R
        gene <- input$gene

        # Extract mutations for the selected gene
        gene_data <- pancan_maf@data[Hugo_Symbol == gene]
        write.csv(gene_data, file, row.names = FALSE)
      }
    )
  })
}
