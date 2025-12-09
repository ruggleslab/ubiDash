# densityUI.R
densityUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    # Reduce the input panel's width to give more space to the plot
    sidebarPanel(
      width = 3,
      h3("Density Plot"),
      p(tab_summaries$`Density Plots`),
      selectizeInput(
        inputId = ns("gene"),
        label = "Gene",
        choices = NULL,
        multiple = FALSE,
        options = list(`live-search` = TRUE)
      )
    ),
    mainPanel(
      textOutput(ns("err")) %>% tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("plot")), color = "#6A2798"),
      br(),
      downloadButton(outputId = ns("download_data"), label = "Download CSV"),
      downloadButton(outputId = ns("download_plot"), label = "Download PDF")
    )
  )
}
# densityServer.R
densityServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Load protein genes if not already available
    if (!exists("cptac_protein_genes")) {
      load("./data/searchlists.RData")
    }
    dataset <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
    protein_genes <- rownames(dataset)

    # Update the gene selectize input with protein genes
    updateSelectizeInput(
      session,
      inputId = "gene",
      choices = protein_genes,
      server = TRUE,
      selected = if ("CRBN" %in% protein_genes) "CRBN" else protein_genes[1]
    )
    # Load and process the data when the user clicks refresh
    sample <- eventReactive(input$gene,
      {
        if (!isTruthy(input$gene)) {
          stop("Please select a gene.")
        }

        # Get the gene and datatype selections
        ups_gn <- isolate(input$gene)
        xlab <- "log2 Normalized Protein Expression"
        # Load the appropriate dataset using mem_get_from_drive()
        dataset <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
        # Load sample metadata from the drive
        sample <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")

        # Subset to the desired datatype and only include samples that passed QC
        sample <- subset(
          sample,
          datatype == "protein" &
            harmonized_data_id %in% colnames(dataset) &
            include == "TRUE"
        ) %>%
          dplyr::select(c(harmonized_data_id, cohort, type, cancer_subtype)) %>%
          unique()

        # Check if the gene exists in the dataset
        if (!(ups_gn %in% rownames(dataset))) {
          stop(paste("Gene", ups_gn, "not found in protein dataset."))
        }

        # Add the gene expression values to the sample data frame
        sample[[ups_gn]] <- as.numeric(dataset[ups_gn, sample$harmonized_data_id])

        # Filter out samples with missing expression values
        sample <- sample %>%
          dplyr::filter(!is.na(!!rlang::sym(ups_gn)))

        # Check if we have enough data
        if (nrow(sample) < 10) {
          stop(paste("Not enough data for gene", ups_gn, ". Need at least 10 samples with expression data."))
        }

        # Reorder cohorts based on the median expression of the gene of interest
        # Only reorder if we have data for all cohorts
        tryCatch(
          {
            sample <- sample %>%
              dplyr::mutate(cohort = forcats::fct_reorder(cohort, !!rlang::sym(ups_gn), .fun = median, .na_rm = TRUE))
          },
          error = function(e) {
            # If reordering fails (due to missing data), keep original cohort order
            warning("Could not reorder cohorts due to missing data for some cohorts. Using original order.")
          }
        )

        sample
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )


    # Render the density ridge plot
    output$plot <- renderPlot(
      {
        req(sample())
        ups_gn <- isolate(input$gene)
        xlab <- "log2 Normalized Protein Expression"

        # Calculate annotations (sample sizes and median positions)
        annotations <- sample() %>%
          dplyr::group_by(cohort, type) %>%
          dplyr::summarise(
            n = dplyr::n(),
            label_x = median(!!rlang::sym(ups_gn), na.rm = TRUE)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(label = paste0("N = ", n))

        # Define the x-axis range and ensure label positions fall within it
        x_range <- range(sample()[[ups_gn]], na.rm = TRUE)
        annotations <- annotations %>%
          dplyr::mutate(label_x = pmax(pmin(label_x, x_range[2]), x_range[1]))

        tryCatch(
          {
            ggplot(sample()) +
              ggridges::geom_density_ridges(
                aes(x = !!rlang::sym(ups_gn), y = cohort, fill = type),
                alpha = 0.5,
                bandwidth = 0.3,
                scale = 1,
                quantiles = 2,
                quantile_lines = TRUE
              ) +
              ggrepel::geom_text_repel(
                data = annotations,
                aes(x = label_x, y = cohort, label = label, color = type),
                nudge_y = 0.7,
                segment.color = NA
              ) +
              labs(
                x = glue::glue("{ups_gn}\n{xlab}"),
                y = NULL,
                title = glue::glue("{ups_gn} Protein Expression Distribution")
              ) +
              ggprism::theme_prism(border = TRUE) +
              scale_fill_manual(values = c("Normal" = "orange", "Tumor" = "purple")) +
              scale_color_manual(values = c("Normal" = "orange", "Tumor" = "purple"), guide = FALSE)
          },
          error = function(e) {
            # If plotting fails, show an error message
            ggplot() +
              annotate("text",
                x = 0.5, y = 0.5,
                label = paste("Error creating plot:", e$message),
                size = 6, hjust = 0.5, vjust = 0.5
              ) +
              theme_void()
          }
        )
      },
      width = 1000,
      height = 800
    )


    # Download handler for the processed sample data as CSV
    output$download_data <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(sample())
        write.csv(sample(), file, row.names = FALSE)
      }
    )


    # Download handler for saving the plot as a PDF file
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(isolate(input$gene), "_distribution.pdf")
      },
      content = function(file) {
        req(sample())
        ups_gn <- isolate(input$gene)
        xlab <- "log2 Normalized Protein Expression"

        annotations <- sample() %>%
          dplyr::group_by(cohort, type) %>%
          dplyr::summarise(
            n = dplyr::n(),
            label_x = median(!!rlang::sym(ups_gn), na.rm = TRUE)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(label = paste0("N = ", n))

        x_range <- range(sample()[[ups_gn]], na.rm = TRUE)
        annotations <- annotations %>%
          dplyr::mutate(label_x = pmax(pmin(label_x, x_range[2]), x_range[1]))

        tryCatch(
          {
            p <- ggplot(sample()) +
              ggridges::geom_density_ridges(
                aes(x = !!rlang::sym(ups_gn), y = cohort, fill = type),
                alpha = 0.5,
                bandwidth = 0.3,
                scale = 1,
                quantiles = 2,
                quantile_lines = TRUE
              ) +
              ggrepel::geom_text_repel(
                data = annotations,
                aes(x = label_x, y = cohort, label = label, color = type),
                nudge_y = 0.7,
                segment.color = NA
              ) +
              labs(
                x = glue::glue("{ups_gn}\n{xlab}"),
                y = NULL,
                title = glue::glue("{ups_gn} Protein Expression Distribution")
              ) +
              ggprism::theme_prism(border = TRUE) +
              scale_fill_manual(values = c("Normal" = "orange", "Tumor" = "purple")) +
              scale_color_manual(values = c("Normal" = "orange", "Tumor" = "purple"), guide = FALSE)

            ggsave(filename = file, plot = p, width = 10, height = 7, device = "pdf")
          },
          error = function(e) {
            # If plotting fails, create a simple error plot
            error_plot <- ggplot() +
              annotate("text",
                x = 0.5, y = 0.5,
                label = paste("Error creating plot:", e$message),
                size = 6, hjust = 0.5, vjust = 0.5
              ) +
              theme_void()
            ggsave(filename = file, plot = error_plot, width = 10, height = 7, device = "pdf")
          }
        )
      }
    )
  })
}
