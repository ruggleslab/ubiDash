correlationGeneAllUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p(tab_summaries$`Pancan Gene vs all`),
        selectizeInput(
          inputId = ns("selectedProtein"),
          label = "Choose Protein(s) of Interest",
          choices = NULL,
          selected = NULL,
          multiple = TRUE # allow multiple gene selection
        ),
        sliderInput(
          inputId = ns("corThreshold"),
          label = "Set minimum correlation threshold (absolute value)",
          min = 0,
          max = 1,
          value = 0,
          step = 0.05
        ),
        selectInput(
          inputId = ns("colorPalette"),
          label = "Choose Color Palette",
          choices = c("Default" = "purple_yellow", "Yellow to White" = "default"),
          selected = "Default"
        ),
        actionButton(
          inputId = ns("plotButton"),
          label = "Generate Plot"
        )
      ),
      mainPanel(
        # This UI output will be populated with plots and download buttons
        uiOutput(ns("plots_ui"))
      )
    )
  )
}
# Server module: updated to include download handlers.
correlationGeneAllServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Read in your datasets (update the paths as needed)
    pancan_covariance <- readRDS("./data/Pancan_protein~protein_spearman_correlations_ups_variance.rds")
    biogrid <- readRDS("./data/BIOGRID-ALL-4.4.241.tab3.rds")

    # Update the protein-selection choices based on unique proteins.
    protein_choices <- sort(unique(pancan_covariance$ups))
    updateSelectizeInput(
      session,
      "selectedProtein",
      choices = protein_choices,
      selected = protein_choices[1]
    )

    # Reactive expression for the selected color palette.
    palette_used <- reactive({
      if (input$colorPalette == "default") {
        c("#ffe05f", "#b0b0b0", "white")
      } else if (input$colorPalette == "purple_yellow") {
        c("#ffe05f", "#6A2798")
      } else {
        c("#ffe05f", "#b0b0b0", "white")
      }
    })

    # Helper function to create the plot for a single gene.
    createPlotForGene <- function(selected_protein, threshold, palette_used) {
      # Filter the BioGRID dataset for interactions involving the current gene.
      biogrid_filtered <- biogrid %>%
        select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
        rename(
          Interactor_A = `Official Symbol Interactor A`,
          Interactor_B = `Official Symbol Interactor B`
        ) %>%
        filter(Interactor_A == selected_protein | Interactor_B == selected_protein) %>%
        mutate(Interactor = ifelse(Interactor_A == selected_protein, Interactor_B, Interactor_A)) %>%
        select(Interactor) %>%
        unique()

      # Filter and prepare the covariance data.
      df <- pancan_covariance %>%
        filter(!is.na(cor) & !is.na(std_dev)) %>%
        filter(ups == selected_protein | cs == selected_protein) %>%
        filter(abs(cor) >= threshold) %>% # apply threshold
        mutate(
          signed_fdr = -log10(fdr + 1e-10) * sign(cor),
          biogrid = ifelse(cs %in% biogrid_filtered$Interactor, "Yes", "No")
        ) %>%
        arrange(desc(signed_fdr))

      # For plotting, focus on rows where the selected protein is the query.
      tmp <- df %>%
        filter(ups == selected_protein) %>%
        mutate(biogrid = factor(biogrid, levels = c("No", "Yes"))) %>%
        arrange(biogrid == "Yes")

      # Identify top 5 and bottom 5 points (by correlation) to label.
      labeled_points <- tmp %>%
        slice_max(order_by = cor, n = 5, with_ties = FALSE) %>%
        bind_rows(
          tmp %>% slice_min(order_by = cor, n = 5, with_ties = FALSE)
        ) %>%
        mutate(text_color = ifelse(
          biogrid == "Yes", "#c70000",
          scales::col_numeric(
            palette = palette_used,
            domain = range(tmp$std_dev)
          )(std_dev)
        ))

      # Create the ggplot for this gene.
      p <- ggplot(tmp, aes(x = signed_fdr, y = cor)) +
        geom_point(
          aes(fill = std_dev, shape = biogrid),
          size = 3,
          color = ifelse(tmp$biogrid == "Yes", "#c70000", "white"),
          stroke = ifelse(tmp$biogrid == "Yes", 0.3, 0.01),
          alpha = 0.9
        ) +
        geom_text_repel(
          data = labeled_points,
          aes(label = cs, color = text_color),
          max.overlaps = 50,
          size = 3
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(
          x = expression("signed" - log[10]("FDR")),
          y = "Spearman, ρ",
          fill = "Standard Deviation",
          shape = "BioGRID"
        ) +
        scale_fill_gradientn(
          colors = palette_used,
          guide = guide_colorbar(barwidth = 1, barheight = 10)
        ) +
        scale_shape_manual(
          values = c("No" = 21, "Yes" = 21),
          labels = c("No" = "No", "Yes" = "Yes"),
          guide = guide_legend(
            override.aes = list(
              fill = c(palette_used[1], palette_used[1]),
              color = c("white", "#c70000"),
              stroke = c(0.01, 0.3)
            )
          )
        ) +
        scale_color_identity() +
        theme_minimal() +
        ggtitle(selected_protein)

      return(p)
    }

    # Reactive expression for generating a list of plots (one per selected gene).
    plotList <- eventReactive(input$plotButton, {
      req(input$selectedProtein, input$corThreshold, input$colorPalette)
      threshold <- input$corThreshold
      # Create a list of plots.
      plots <- lapply(input$selectedProtein, function(selected_protein) {
        createPlotForGene(selected_protein, threshold, palette_used())
      })
      return(plots)
    })

    # Reactive expression for the filtered data (for CSV download).
    filteredData <- eventReactive(input$plotButton, {
      req(input$selectedProtein, input$corThreshold)
      threshold <- input$corThreshold

      data_list <- lapply(input$selectedProtein, function(selected_protein) {
        biogrid_filtered <- biogrid %>%
          select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
          rename(
            Interactor_A = `Official Symbol Interactor A`,
            Interactor_B = `Official Symbol Interactor B`
          ) %>%
          filter(Interactor_A == selected_protein | Interactor_B == selected_protein) %>%
          mutate(Interactor = ifelse(Interactor_A == selected_protein, Interactor_B, Interactor_A)) %>%
          select(Interactor) %>%
          unique()

        df <- pancan_covariance %>%
          filter(!is.na(cor) & !is.na(std_dev)) %>%
          filter(ups == selected_protein | cs == selected_protein) %>%
          filter(abs(cor) >= threshold) %>% # apply threshold
          mutate(
            signed_fdr = -log10(fdr + 1e-10) * sign(cor),
            biogrid = ifelse(cs %in% biogrid_filtered$Interactor, "Yes", "No")
          ) %>%
          arrange(desc(signed_fdr)) %>%
          filter(ups == selected_protein)
        df$gene <- selected_protein
        return(df)
      })

      combined_df <- do.call(rbind, data_list)
      return(combined_df)
    })

    # Dynamically create UI output for each plot and its download buttons.
    # In the server module, replace the previous output$plots_ui with the following:
    output$plots_ui <- renderUI({
      req(plotList())
      ns <- session$ns
      n_genes <- length(plotList())

      # Decide on the number of columns per row
      cols_per_row <- if (n_genes <= 3) n_genes else 3
      col_width <- floor(12 / cols_per_row)

      # Split the indices into chunks (rows)
      idx <- seq_along(plotList())
      row_chunks <- split(idx, ceiling(seq_along(idx) / cols_per_row))

      # Build a fluidRow for each chunk
      fluidRows <- lapply(row_chunks, function(chunk) {
        fluidRow(
          lapply(chunk, function(i) {
            gene_name <- input$selectedProtein[i]
            column(
              width = col_width,
              h3(gene_name),
              withSpinner(
                plotOutput(ns(paste0("corPlot", i)), width = "100%", height = "600px"),
                color = "#6A2798"
              ),
              br(),
              downloadButton(ns(paste0("downloadPDF", i)),
                label = paste("Download", gene_name, "Plot as PDF")
              ),
              br(), br(),
              downloadButton(ns(paste0("downloadData", i)),
                label = paste("Download", gene_name, "Data as CSV")
              ),
              hr()
            )
          })
        )
      })

      # Return the list of fluidRows
      do.call(tagList, fluidRows)
    })

    # For each plot in plotList, create corresponding outputs and download handlers.
    observeEvent(plotList(), {
      lapply(seq_along(plotList()), function(i) {
        local({
          my_i <- i
          gene_name <- input$selectedProtein[my_i]

          # Render the plot.
          output[[paste0("corPlot", my_i)]] <- renderPlot(
            {
              plotList()[[my_i]]
            },
            height = 600
          )

          # Download handler for the PDF plot.
          output[[paste0("downloadPDF", my_i)]] <- downloadHandler(
            filename = function() {
              paste0("corPlot_", gene_name, "_", Sys.Date(), ".pdf")
            },
            content = function(file) {
              ggsave(file, plot = plotList()[[my_i]], device = "pdf", width = 9, height = 6, units = "in")
            }
          )

          # Download handler for the CSV data.
          output[[paste0("downloadData", my_i)]] <- downloadHandler(
            filename = function() {
              paste0("corPlot_data_", gene_name, "_", Sys.Date(), ".csv")
            },
            content = function(file) {
              df <- filteredData()
              # Subset the combined data for the current gene.
              gene_df <- df[df$gene == gene_name, ]
              write.csv(gene_df, file, row.names = FALSE)
            }
          )
        })
      })
    })
  })
}
