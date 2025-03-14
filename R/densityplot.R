# densityUI.R
densityUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    # Reduce the input panel's width to give more space to the plot
    sidebarPanel(
      width = 3,
      h3("Density Plot"),
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
      p("Download the sample data:"),
      downloadButton(outputId = ns("download_data"), label = "Download CSV"),
      p("Download the plot as SVG:"),
      downloadButton(outputId = ns("download_plot"), label = "Download SVG")
    )
  )
}
# densityServer.R
densityServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Update the gene selectize input (assuming 'maf_genes' is defined globally)
    updateSelectizeInput(
      session,
      inputId = "gene",
      choices = maf_genes,
      server = TRUE,
      selected = "KBTBD11"
    )
    # Load and process the data when the user clicks refresh
    sample <- eventReactive(input$gene, {
      validate(
        need(input$gene, "Please select a gene."))
      
      # Get the gene and datatype selections
      ups_gn   <- isolate(input$gene)
      xlab     <- "Log2(TPM+1)"
      # Load the appropriate dataset using mem_get_from_drive()
      dataset <- mem_get_from_drive("cptac_transformed_rna_expression_v1-1.rds")
      # Load sample metadata from the drive
      sample <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")
      
      # Subset to the desired datatype and only include samples that passed QC
      sample <- subset(sample,
                       datatype == 'rna' &
                         harmonized_data_id %in% colnames(dataset) &
                         include == "TRUE") %>%
        dplyr::select(c(harmonized_data_id, cohort, type, cancer_subtype)) %>%
        unique()
      
      # Add the gene expression values to the sample data frame
      sample[[ups_gn]] <- as.numeric(dataset[ups_gn, sample$harmonized_data_id])
      
      # Reorder cohorts based on the median expression of the gene of interest
      sample <- sample %>%
        dplyr::mutate(cohort = forcats::fct_reorder(cohort, !!rlang::sym(ups_gn), .fun = median, .na_rm = TRUE))
      
      sample
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    
    # Render the density ridge plot
    output$plot <- renderPlot({
      req(sample())
      ups_gn   <- isolate(input$gene)
      xlab     <- "Log2(TPM+1)"
      
      # Calculate annotations (sample sizes and median positions)
      annotations <- sample() %>%
        dplyr::group_by(cohort, type) %>%
        dplyr::summarise(n = dplyr::n(),
                         label_x = median(!!rlang::sym(ups_gn), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(label = paste0("N = ", n))
      
      # Define the x-axis range and ensure label positions fall within it
      x_range <- range(sample()[[ups_gn]], na.rm = TRUE)
      annotations <- annotations %>%
        dplyr::mutate(label_x = pmax(pmin(label_x, x_range[2]), x_range[1]))
      
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
          title = glue::glue("{ups_gn} Global protein expression")
        ) +
        ggprism::theme_prism(border = TRUE) +
        scale_fill_manual(values = c("Normal" = "orange", "Tumor" = "purple")) +
        scale_color_manual(values = c("Normal" = "orange", "Tumor" = "purple"), guide = FALSE)
    }, width = 1000, height = 800)
    
    
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
    
    
    # Download handler for saving the plot as an SVG file
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(isolate(input$gene), "_distribution.svg")
      },
      content = function(file) {
        req(sample())
        ups_gn   <- isolate(input$gene)
        xlab <- "Log2(TPM+1)"
        
        annotations <- sample() %>%
          dplyr::group_by(cohort, type) %>%
          dplyr::summarise(n = dplyr::n(),
                           label_x = median(!!rlang::sym(ups_gn), na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(label = paste0("N = ", n))
        
        x_range <- range(sample()[[ups_gn]], na.rm = TRUE)
        annotations <- annotations %>%
          dplyr::mutate(label_x = pmax(pmin(label_x, x_range[2]), x_range[1]))
        
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
            title = glue::glue("{ups_gn} Global protein expression")
          ) +
          ggprism::theme_prism(border = TRUE) +
          scale_fill_manual(values = c("Normal" = "orange", "Tumor" = "purple")) +
          scale_color_manual(values = c("Normal" = "orange", "Tumor" = "purple"), guide = FALSE)
        
        ggsave(filename = file, plot = p, width = 5, height = 5, device = "svg")
      }
    )
    
  })
}