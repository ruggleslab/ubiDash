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
                     label =  "Genes of Interest",options = list(`live-search` = TRUE),  multiple = TRUE),
   
      # cutoff options
      textInput(ns("PValCutoff"), "Adjusted P Value Cutoff:",
                placeholder = "5e-8", value = "0.001"),
      
      sliderInput(ns("logFCCutoff"), "Log2(Tumor/Normal) Cutoff:",
                  min = .25, max = 2, value = 1, step = 0.05),
      
      actionButton(ns("refresh"), "Generate Plot"),
      
    
    ),
    mainPanel(
      shinycssloaders::withSpinner(plotlyOutput(ns("distPlot")), color = "#6A2798"),
      p("Download this data"),
      downloadButton(outputId = ns("downloadVolcanoData"), label = "Download")
      ),
  )
}


volcanoServer <- function(id, crl_genes, DE_protein) {
  moduleServer(id,function(input,output,session){


        
    updateSelectizeInput(session, 'genes', 
                         choices = cptac_protein_genes, 
                         server = TRUE)
    
    data = eventReactive(input$refresh, {
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
      
      genes <- isolate(input$genes)
      
      # prepare dataframe
      data = as.data.frame(DE_protein[[isolate(input$cohort)]])
      data$POI <- "Other"
      data$POI[data$gene %in% crl_genes] = "CRL"
      data$POI[toupper(data$gene) %in% genes] = "Gene of Interest"
      data$log.scal.pval <- -log10(data$sca.adj.pval)
      
     return(data)
    
    },ignoreNULL = FALSE, ignoreInit=FALSE)
    
    
    output$distPlot <- renderPlotly({      

      genes <- isolate(input$genes)

      logFCCutoff <- isolate(input$logFCCutoff)
      
      PValCutoff = -log10(as.numeric(gsub('\\*10\\^', 'e', isolate(input$PValCutoff))))
      print(PValCutoff)
      if(!isTruthy(PValCutoff)){
        PValCutoff = 0.0001
      }
      # print(PValCutoff)
      # PValCutoff <- -log10(as.double(isolate(input$PValCutoff)))
      
      # get significant genes
      labeledPoints <- filter(data(), (gene %in% crl_genes & (abs(logFC) > logFCCutoff & abs(log.scal.pval) > PValCutoff)) | toupper(gene) %in% genes)
      

   
      p <- ggplot(data(), aes(x = logFC, y = log.scal.pval, color = POI)) +
        ylab("-Log10(adjusted p-val)") +
        xlab("Log2(Tumor/Normal)") +
        labs(title = isolate(input$cohort)) +
        geom_point(shape = 16, alpha = .5, stroke=NA, size = 2, mapping = aes(text = paste(
          "Gene: ", gene, "\n",
          "POI: ", POI, "\n",
          "LogFC: ", logFC, "\n",
          "P Value: ", log.scal.pval, "\n",
          sep = ""
        ))) +
        scale_color_manual(breaks = c("CRL", "Gene of Interest", "Other"), values = c("#FF0000", "#0000FF", "#00000066")) +
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
