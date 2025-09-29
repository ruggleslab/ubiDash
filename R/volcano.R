volcanoUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Visualization"),
      p(tab_summaries$`Volcano Plots`),
      radioButtons(ns("plot_type"), "Plot Type:",
        choices = c("Volcano Plot" = "volcano", "Heatmap" = "heatmap"),
        selected = "volcano", inline = TRUE
      ),

      # Volcano plot inputs
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
        selectInput(ns("class_volcano"), "UPS Class",
          choices = NULL,
          selected = "E3 ubiquitin and UBL ligases"
        ),
        selectizeInput(ns("poi"), "Protein(s) of Interest", choices = NULL, multiple = TRUE),
        selectInput(ns("cohort"), "Cohort", cohorts),
        textInput(ns("PValCutoff"), "Adjusted P Value Cutoff:",
          placeholder = "5e-8", value = "0.001"
        ),
        sliderInput(ns("logFCCutoff"), "Log2(Tumor/Normal) Cutoff:",
          min = 0.25, max = 2, value = 1, step = 0.05
        )
      ),

      # Heatmap inputs
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
        selectizeInput(ns("class_heatmap"), "UPS Classes", choices = NULL, multiple = TRUE),
        selectizeInput(ns("filter_proteins_heatmap"),
          "Protein(s) of Interest",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Type or select one or more protein symbols",
            server = TRUE
          )
        )
      )
    ),
    mainPanel(
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] == 'volcano'"),
        shinycssloaders::withSpinner(plotlyOutput(ns("distPlot"), height = "800px"),
          color = "#6A2798"
        ),
        p("Download data and plot:"),
        downloadButton(outputId = ns("downloadVolcanoData"), label = "Download Data"),
        downloadButton(outputId = ns("downloadSVG"), label = "Download SVG")
      ),
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
        shinycssloaders::withSpinner(plotOutput(ns("heatmapPlot"), height = "800px"),
          color = "#6A2798"
        ),
        p("Download data and plot:"),
        downloadButton(outputId = ns("downloadHeatmapData"), label = "Download Data"),
        downloadButton(outputId = ns("downloadHeatmapPDF"), label = "Download PDF")
      )
    )
  )
}

volcanoServer <- function(id, DE_protein) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load UPS data and colors.
    ups_data <- openxlsx::read.xlsx(
      xlsxFile = "./data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
      sheet = "dense"
    )
    colors_tab <- read.csv("./data/colors.csv")
    colors <- setNames(colors_tab$color, colors_tab$id)

    # Update UPS class choices for volcano plot.
    observe({
      ups_class_selection <- ups_data %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System") %>%
        dplyr::pull(Class) %>%
        unique()
      updateSelectInput(session, "class_volcano",
        choices = ups_class_selection,
        selected = "E3 ubiquitin and UBL ligases"
      )
    })

    # Update UPS class choices for heatmap (select all by default).
    observe({
      ups_class_selection <- ups_data %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System") %>%
        dplyr::pull(Class) %>%
        unique()
      updateSelectizeInput(session, "class_heatmap",
        choices = ups_class_selection,
        selected = ups_class_selection
      )
    })

    # Update protein-of-interest choices for volcano plot based on selected class.
    observeEvent(input$class_volcano, {
      req(input$class_volcano)
      class_genes <- ups_data %>%
        dplyr::filter(
          Branch == "Ubiquitin Proteasome System",
          Class == input$class_volcano
        ) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()
      updateSelectizeInput(session, "poi", choices = class_genes, server = TRUE)
    })

    # Update protein search input for the heatmap.
    observe({
      if (input$plot_type == "heatmap") {
        ups_genes <- ups_data %>%
          dplyr::filter(Branch == "Ubiquitin Proteasome System") %>%
          dplyr::pull(Gene.Symbol) %>%
          unique() %>%
          sort()
        updateSelectizeInput(session, "filter_proteins_heatmap",
          choices = ups_genes,
          selected = NULL
        )
      }
    })

    # Volcano plot data processing.
    volcano_data <- reactive({
      if (input$plot_type != "volcano") {
        return(NULL)
      }

      selected_class <- input$class_volcano
      poi_genes <- input$poi

      class_genes <- ups_data %>%
        dplyr::filter(
          Branch == "Ubiquitin Proteasome System",
          Class == selected_class
        ) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()

      df <- as.data.frame(DE_protein[[input$cohort]])
      df$POI <- dplyr::case_when(
        df$gene %in% poi_genes ~ "Protein of Interest",
        df$gene %in% class_genes ~ selected_class,
        TRUE ~ "Other"
      )

      # Compute the traditional volcano plot metric as -log10(adjusted p-value).
      df$log.scal.pval <- -log10(df$sca.adj.pval)

      list(data = df, selected_class = selected_class)
    })

    # Heatmap data processing.
    heatmap_data <- reactive({
      if (input$plot_type != "heatmap") {
        return(NULL)
      }

      # Read the matrix file.
      out <- read.csv(file = "./data/EXP193_DEqMS_CPTAC_imputed_proteome_signed_p.adj_matrix.csv")

      selected_classes <- input$class_heatmap
      if (is.null(selected_classes)) {
        return(NULL)
      }

      ups_genes <- ups_data %>%
        dplyr::filter(
          Branch == "Ubiquitin Proteasome System",
          Class %in% selected_classes
        ) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()

      input_matrix <- out %>%
        dplyr::filter(gene %in% ups_genes) %>%
        column_to_rownames(var = "gene") %>%
        .[which(rowMeans(!is.na(.)) > 0.5), ] %>%
        as.matrix()

      # Extract UPS annotation.
      ups_annot <- ups_data %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System") %>%
        dplyr::select(Gene.Symbol, Class) %>%
        dplyr::filter(Gene.Symbol %in% rownames(input_matrix)) %>%
        dplyr::distinct(Gene.Symbol, .keep_all = TRUE) %>%
        column_to_rownames(var = "Gene.Symbol") %>%
        dplyr::select(Class) %>%
        .[rownames(input_matrix), , drop = FALSE]

      input_matrix <- input_matrix[rownames(input_matrix) %in% rownames(ups_annot), , drop = FALSE]

      # Create the display matrix.
      display_matrix <- input_matrix
      if (!is.null(input$filter_proteins_heatmap) && length(input$filter_proteins_heatmap) > 0) {
        rn <- rownames(display_matrix)
        rn[!(rn %in% input$filter_proteins_heatmap)] <- ""
        rownames(display_matrix) <- rn
        show_names <- TRUE
      } else {
        show_names <- FALSE
      }

      list(
        matrix = display_matrix,
        annotation = ups_annot,
        colors = colors,
        show_names = show_names
      )
    })

    # Volcano plot rendering.
    output$distPlot <- renderPlotly({
      req(volcano_data())
      df <- volcano_data()$data
      selected_class <- volcano_data()$selected_class

      logFCCutoff <- input$logFCCutoff
      PValCutoff <- as.numeric(input$PValCutoff)
      if (!isTruthy(PValCutoff)) PValCutoff <- 0.001
      PValCutoff_log <- -log10(PValCutoff)

      labeledPoints <- dplyr::filter(
        df,
        (POI == "Protein of Interest") |
          (POI == selected_class & abs(logFC) >= logFCCutoff & log.scal.pval >= PValCutoff_log)
      )

      p <- ggplot(df, aes(x = logFC, y = log.scal.pval)) +
        geom_point(
          data = subset(df, POI == "Other"),
          aes(text = paste(
            "Gene: ", gene, "\n",
            "Category: ", POI, "\n",
            "LogFC: ", logFC, "\n",
            "Adj P: ", round(log.scal.pval, 3)
          )),
          color = "grey70", size = 2, alpha = 0.6
        ) +
        geom_point(
          data = subset(df, POI != "Other"),
          aes(
            color = POI,
            text = paste(
              "Gene: ", gene, "\n",
              "Category: ", POI, "\n",
              "LogFC: ", logFC, "\n",
              "Adj P: ", round(log.scal.pval, 3)
            )
          ),
          size = 2, alpha = 0.6
        ) +
        scale_color_manual(values = {
          tmp <- c("Protein of Interest" = "firebrick")
          tmp[selected_class] <- colors[selected_class]
          tmp
        }) +
        geom_vline(xintercept = logFCCutoff, linetype = "dashed", color = "black") +
        geom_vline(xintercept = -logFCCutoff, linetype = "dashed", color = "black") +
        geom_hline(yintercept = c(PValCutoff_log), linetype = "dashed", color = "black") +
        labs(
          title = input$cohort,
          x = "Log2(Tumor/Normal)",
          y = "-Log10(Adjusted P-value)"
        ) +
        theme_minimal()

      if (nrow(labeledPoints) > 0) {
        p <- p + geom_text(
          data = labeledPoints, aes(label = gene),
          vjust = -1, size = 3, check_overlap = TRUE
        )
      }

      ggplotly(p, tooltip = "text") %>% layout(dragmode = FALSE, height = 700)
    })

    # Heatmap rendering.
    output$heatmapPlot <- renderPlot(
      {
        req(heatmap_data())
        data <- heatmap_data()
        row_annot <- rowAnnotation(
          Class = data$annotation$Class,
          col = list(Class = data$colors[levels(as.factor(data$annotation$Class))]),
          show_annotation_name = FALSE
        )

        Heatmap(data$matrix,
          column_title = NULL,
          row_title = NULL,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          show_column_names = TRUE,
          show_row_names = data$show_names,
          name = "-Log10(p.adj)",
          col = circlize::colorRamp2(
            c(-10, -2.5, 2.5, 10),
            c("blue", "white", "white", "red")
          ),
          left_annotation = row_annot,
          row_split = data$annotation$Class,
          show_row_dend = TRUE,
          show_heatmap_legend = TRUE
        )
      },
      height = 800
    )

    ### DOWNLOAD HANDLERS ###

    output$downloadVolcanoData <- downloadHandler(
      filename = function() {
        paste0("volcano_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(volcano_data())
        write.csv(volcano_data()$data, file, row.names = FALSE)
      }
    )

    output$downloadSVG <- downloadHandler(
      filename = function() {
        paste0("volcano_plot_", Sys.Date(), ".svg")
      },
      content = function(file) {
        req(volcano_data())
        df <- volcano_data()$data
        selected_class <- volcano_data()$selected_class

        logFCCutoff <- input$logFCCutoff
        PValCutoff <- as.numeric(input$PValCutoff)
        if (!isTruthy(PValCutoff)) PValCutoff <- 0.001
        PValCutoff_log <- -log10(PValCutoff)

        labeledPoints <- dplyr::filter(
          df,
          (POI == "Protein of Interest") |
            (POI == selected_class & abs(logFC) >= logFCCutoff & log.scal.pval >= PValCutoff_log)
        )

        p <- ggplot(df, aes(x = logFC, y = log.scal.pval)) +
          geom_point(
            data = subset(df, POI == "Other"),
            aes(text = paste(
              "Gene: ", gene, "\n",
              "Category: ", POI, "\n",
              "LogFC: ", logFC, "\n",
              "Adj P: ", round(log.scal.pval, 3)
            )),
            color = "grey70", size = 2, alpha = 0.6
          ) +
          geom_point(
            data = subset(df, POI != "Other"),
            aes(
              color = POI,
              text = paste(
                "Gene: ", gene, "\n",
                "Category: ", POI, "\n",
                "LogFC: ", logFC, "\n",
                "Adj P: ", round(log.scal.pval, 3)
              )
            ),
            size = 2, alpha = 0.6
          ) +
          scale_color_manual(values = {
            tmp <- c("Protein of Interest" = "firebrick")
            tmp[selected_class] <- colors[selected_class]
            tmp
          }) +
          geom_vline(xintercept = logFCCutoff, linetype = "dashed", color = "black") +
          geom_vline(xintercept = -logFCCutoff, linetype = "dashed", color = "black") +
          geom_hline(yintercept = c(PValCutoff_log), linetype = "dashed", color = "black") +
          labs(
            title = input$cohort,
            x = "Log2(Tumor/Normal)",
            y = "-Log10(Adjusted P-value)"
          ) +
          theme_minimal()

        if (nrow(labeledPoints) > 0) {
          p <- p + geom_text(
            data = labeledPoints, aes(label = gene),
            vjust = -1, size = 3, check_overlap = TRUE
          )
        }

        ggsave(file, plot = p, device = "svg", width = 8, height = 6)
      }
    )

    output$downloadHeatmapData <- downloadHandler(
      filename = function() {
        paste0("heatmap_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(heatmap_data())
        write.csv(heatmap_data()$matrix, file, row.names = TRUE)
      }
    )

    output$downloadHeatmapPDF <- downloadHandler(
      filename = function() {
        paste0("heatmap_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        req(heatmap_data())
        pdf(file, width = 10, height = 8)
        data <- heatmap_data()
        row_annot <- rowAnnotation(
          Class = data$annotation$Class,
          col = list(Class = data$colors[levels(as.factor(data$annotation$Class))]),
          show_annotation_name = FALSE
        )

        ht <- Heatmap(data$matrix,
          column_title = NULL,
          row_title = NULL,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          name = "-Log10(p.adj)",
          col = circlize::colorRamp2(
            c(-10, -2.5, 2.5, 10),
            c("blue", "white", "white", "red")
          ),
          left_annotation = row_annot,
          row_split = data$annotation$Class,
          show_row_dend = TRUE,
          show_heatmap_legend = TRUE
        )
        draw(ht)
        dev.off()
      }
    )
  })
}
