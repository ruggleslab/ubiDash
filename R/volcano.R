volcanoUI <- function(id, cohorts){
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Volcano Plot"),
      p("(description paragraph)"),
      
      selectInput(ns("cohort"), "Cohort", cohorts),
      
      # input options
      selectizeInput(inputId = ns("genes"),
                     choices = NULL,
                     label =  "Genes",options = list(`live-search` = TRUE),  multiple = TRUE),
   
      # cutoff options
      selectInput(ns("PValCutoff"), "Adjusted P Value Cutoff:", c(0.001, 0.01, 0.05)),
      sliderInput(ns("logFCCutoff"), "Log2(Tumor/Normal) Cutoff:",
                  min = .25, max = 2, value = 1, step = 0.05),
      
      actionButton(ns("refresh"), "Generate Plot"),
      
      p("Download this data"),
      downloadButton(outputId = ns("downloadVolcanoData"), label = "Download")
    ),
    mainPanel(
      shinycssloaders::withSpinner(plotlyOutput(ns("distPlot")), color = "#6A2798"))
  )
}


volcanoServer <- function(id, crl_genes, DE_protein) {
  moduleServer(id,function(input,output,session){


        
    updateSelectizeInput(session, 'genes', 
                         choices = cptac_protein_genes, 
                         server = TRUE)
    
    data = eventReactive(input$refresh, {
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
      
      genes <- isolate(input$genes)
      
      # prepare dataframe
      data = as.data.frame(DE_protein[[isolate(input$cohort)]])
      data$POI <- "Other"
      data$POI[data$gene %in% crl_genes] = "CRL"
      data$POI[toupper(data$gene) %in% genes] = "Potential Substrate"
      data$log.scal.pval <- -log10(data$sca.adj.pval)
      
     return(data)
    
    },ignoreNULL = FALSE, ignoreInit=FALSE)
    
    
    output$distPlot <- renderPlotly({      

      genes <- isolate(input$genes)

      logFCCutoff <- isolate(input$logFCCutoff)
      PValCutoff <- -log10(as.double(isolate(input$PValCutoff)))
      
      # get significant genes
      labeledPoints <- filter(data(), (gene %in% crl_genes & (abs(logFC) > logFCCutoff & abs(log.scal.pval) > PValCutoff)) | toupper(gene) %in% genes)
      

   
      p <- ggplot(data(), aes(x = logFC, y = log.scal.pval, color = POI)) +
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
        req(data())
        write.csv(data(), file)
      }
    )
  })
}
