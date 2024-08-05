correlationsUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(

      # gene input options
      p("The underlying data for this tab is a pre-filtered table
                   only containing genes correlated to CRLs with cor > | 0.5 |.
                   So, input text should specify a CRL of interest to produce a
                   scatterplot of all other proteins correlated to the specified
                   CRL."),
      textInput(ns("inputCorrGenesTXT"), "Input text (separate genes with comma)", value = "", width = NULL, placeholder = NULL),
      sliderInput(ns("CorrelationCutoff"), "Correlation Cutoff:",
        min = .5, max = 1,
        value = .5, step = 0.05
      ),
      actionButton(ns("generateCorrelation"), "Generate Plot"),
      p("Download this data"),
      downloadButton(outputId = ns("downloadCorrelationData"), label = "Download")
    ),
    mainPanel(
      tabsetPanel(
        # three tabs to compare proteome and genome of normal and tumor
        tabPanel(
          "Tumor Proteome - Normal Proteome",
          textOutput(ns("corrErr1")) %>%
            tagAppendAttributes(class = "error"),
          shinycssloaders::withSpinner(plotlyOutput(ns("pt_pn")), color = "#6A2798")
        ),
        tabPanel(
          "Tumor Proteome - Tumor Transcriptome",
          textOutput(ns("corrErr2")) %>%
            tagAppendAttributes(class = "error"),
          shinycssloaders::withSpinner(plotlyOutput(ns("pt_rt")), color = "#6A2798")
        ),
        tabPanel(
          "Tumor Transcriptome - Normal Transcriptome",
          textOutput(ns("corrErr3")) %>%
            tagAppendAttributes(class = "error"),
          shinycssloaders::withSpinner(plotlyOutput(ns("rt_rn")), color = "#6A2798")
        )
      )
    )
  )
}

correlationsServer <- function(id, lines) {
  moduleServer(id,function(input,output,session){
    
    warning <- reactiveVal("")
    # proteome tumor - proteome normal
    output$pt_pn <- renderPlotly({
    
      
      # accessed separately three times  
      plot_tab <- mem_get_from_drive("cptac_multi_omic_correlations.rds")
      
      
      input$generateCorrelation
      output$corrErr1 <- renderText({
        warning()
      })
      tryCatch(
        {
          p <- generateCorrPlot(
            isolate(input$CorrelationCutoff), plot_tab, "pr_tumor.cor", "pr_normal.cor",
            isolate(input$inputCorrGenesTXT), warning
          )
          ggplotly(p, tooltip = "text")
        },
        error = function(e) {
          print(e)
          output$corrErr1 <- renderText({
            "Error: No data available for input genes"
          })
          return()
        }
      )
    })
    # proteome tunor - transcriptome tunor
    output$pt_rt <- renderPlotly({
      
      plot_tab <- mem_get_from_drive("cptac_multi_omic_correlations.rds")
      
      
      input$generateCorrelation
      output$corrErr2 <- renderText({
        warning()
      })
      tryCatch(
        {
          p <- generateCorrPlot(
            isolate(input$CorrelationCutoff), plot_tab, "pr_tumor.cor", "rn_tumor.cor",
            isolate(input$inputCorrGenesTXT), warning
          )
          ggplotly(p, tooltip = "text")
        },
        error = function(e) {
          output$corrErr2 <- renderText({
            "Error: No data available for input genes"
          })
          return()
        }
      )
    })

    # transcriptome tumor - transcriptome normal
    output$rt_rn <- renderPlotly({
      
      plot_tab <- mem_get_from_drive("cptac_multi_omic_correlations.rds")
      
      input$generateCorrelation
      output$corrErr3 <- renderText({
        warning()
      })
      tryCatch(
        {
          p <- generateCorrPlot(
            isolate(input$CorrelationCutoff), plot_tab, "rn_tumor.cor", "rn_normal.cor",
            isolate(input$inputCorrGenesTXT), warning
          )
          ggplotly(p, tooltip = "text")
        },
        error = function(e) {
          output$corrErr3 <- renderText({
            "Error: No data available for input genes"
          })
          return()
        }
      )
    })

    output$downloadCorrelationData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(plot_tab, file)
      }
    )
  })
}
