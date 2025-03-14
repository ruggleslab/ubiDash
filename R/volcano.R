volcanoUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Volcano Plot"),
      p("(description paragraph)"),
      selectInput(ns("class"), "UPS Class", choices = NULL),
      selectizeInput(ns("poi"), "Protein(s) of Interest", choices = NULL, multiple = TRUE),
      selectInput(ns("cohort"), "Cohort", cohorts),
      textInput(ns("PValCutoff"), "Adjusted P Value Cutoff:", placeholder = "5e-8", value = "0.001"),
      sliderInput(ns("logFCCutoff"), "Log2(Tumor/Normal) Cutoff:", min = 0.25, max = 2, value = 1, step = 0.05),
      actionButton(ns("refresh"), "Generate Plot")
    ),
    mainPanel(
      shinycssloaders::withSpinner(plotlyOutput(ns("distPlot")), color = "#6A2798"),
      p("Download data and plot:"),
      downloadButton(outputId = ns("downloadVolcanoData"), label = "Download Data"),
      downloadButton(outputId = ns("downloadSVG"), label = "Download SVG")
    )
  )
}

volcanoServer <- function(id, DE_protein) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load UPS data for class and gene information
    ups_data <- openxlsx::read.xlsx(
      xlsxFile = './data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx',
      sheet = 'dense'
    )
    
    # Update class choices
    observe({
      ups_class_selection <- ups_data %>%
        dplyr::filter(Branch == 'Ubiquitin Proteasome System') %>%
        dplyr::pull(Class) %>%
        unique()
      updateSelectInput(session, "class", choices = ups_class_selection)
    })
    
    # Update protein of interest choices based on selected class
    observeEvent(input$class, {
      req(input$class)
      class_genes <- ups_data %>%
        dplyr::filter(Branch == 'Ubiquitin Proteasome System', Class == input$class) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()
      updateSelectizeInput(session, "poi", choices = class_genes, server = TRUE)
    })
    
    data <- eventReactive(c(input$refresh,input$cohort), {
      selected_class <- isolate(input$class)
      poi_genes <- isolate(input$poi)
      
      # Get genes in the selected class
      class_genes <- ups_data %>%
        dplyr::filter(Branch == 'Ubiquitin Proteasome System', Class == selected_class) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()
      
      # Prepare dataframe with clear gene category assignment using case_when
      df <- as.data.frame(DE_protein[[isolate(input$cohort)]])
      df$POI <- dplyr::case_when(
        df$gene %in% poi_genes ~ "Protein of Interest",
        df$gene %in% class_genes ~ selected_class,
        TRUE ~ "Other"
      )
      df$log.scal.pval <- -log10(df$sca.adj.pval)
      
      return(list(data = df, selected_class = selected_class))
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    colors_tab <- read.csv("./data/colors.csv")
    colors_list <- setNames(colors_tab$color, colors_tab$id)

    # Reactive expression to generate the ggplot object
    plot_obj <- reactive({
      data_list <- data()
      df <- data_list$data
      selected_class <- data_list$selected_class
      
      logFCCutoff <- input$logFCCutoff
      PValCutoff <- as.numeric(gsub('\\*10\\^', 'e', input$PValCutoff))
      if (!isTruthy(PValCutoff)) {
        PValCutoff <- 0.001
      }
      PValCutoff_log <- -log10(PValCutoff)
      
  
      labeledPoints <- dplyr::filter(df, 
                                     (POI == "Protein of Interest") |
                                       (POI == selected_class & abs(logFC) >= logFCCutoff & log.scal.pval >= PValCutoff_log))
      
      # Define an improved color palette
      colors <- c("Other" = "grey70")
      colors[selected_class] <- colors_list[selected_class]
      
      p <- ggplot(df, aes(x = logFC, y = log.scal.pval)) +
        # Plot the "Other" points first
        geom_point(data = subset(df, POI == "Other"),
                   aes(text = paste("Gene: ", gene, "\n",
                                    "Category: ", POI, "\n",
                                    "LogFC: ", logFC, "\n",
                                    "Adj P: ", sca.adj.pval)),
                   color = "grey70", size = 2, alpha = 0.6) +
        # Then plot the non-"Other" points on top
        geom_point(data = subset(df, POI != "Other"),
                   aes(color = POI,
                       text = paste("Gene: ", gene, "\n",
                                    "Category: ", POI, "\n",
                                    "LogFC: ", logFC, "\n",
                                    "Adj P: ", sca.adj.pval)),
                   size = 2, alpha = 0.6) +
        scale_color_manual(values = {
          colors <- c("Protein of Interest" = "firebrick")
          colors[selected_class] <- colors_list[selected_class]
          colors
        }) +
        geom_vline(xintercept = logFCCutoff, linetype = "dashed", color = "black") +
        geom_vline(xintercept = -logFCCutoff, linetype = "dashed", color = "black") +
        geom_hline(yintercept = PValCutoff_log, linetype = "dashed", color = "black") +
        labs(title = input$cohort,
             x = "Log2(Tumor/Normal)",
             y = "-Log10(Adjusted P-value)") +
        theme_minimal()
      
      # Add text labels for the proteins of interest and those meeting the threshold
      if (nrow(labeledPoints) > 0) {
        p <- p + geom_text(data = labeledPoints, aes(label = gene),
                           vjust = -1, size = 3, check_overlap = TRUE)
      }
      
      return(p)
    })
    
    output$distPlot <- renderPlotly({
      ggplotly(plot_obj(), tooltip = "text") %>% layout(dragmode = FALSE)
    })
    
    output$downloadVolcanoData <- downloadHandler(
      filename = function() {
        paste("volcano_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(data())
        write.csv(data()$data, file, row.names = FALSE)
      }
    )
    
    output$downloadSVG <- downloadHandler(
      filename = function() {
        paste("volcano_plot_", Sys.Date(), ".svg", sep = "")
      },
      content = function(file) {
        # Save the ggplot object as an SVG file using ggsave
        ggsave(file, plot = plot_obj(), device = "svg", width = 7, height = 5)
      }
    )
  })
}