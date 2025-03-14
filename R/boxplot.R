boxplotUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        # Select between Proteome or TSI data
        # Add static significance legend
        p(HTML("<b>Significance Levels:</b><br>",
               "ns: p > 0.05 (not significant)<br>",
               "*: p ≤ 0.05<br>",
               "**: p ≤ 0.01<br>",
               "***: p ≤ 0.001")),
        selectInput(ns("data_type"), "Data Type",
                    choices = c("Proteome", "TSI"), selected = "Proteome"),
        # Protein of Interest (auto-populated)
        selectInput(ns("protein"), "Protein of Interest", choices = NULL, selected = "CRBN"),
        # Sample type selection: allow Tumor, Normal, or both
        selectInput(ns("sample_type"), "Sample Type",
                    choices = c("Tumor", "Normal"), 
                    selected = c("Tumor", "Normal"), multiple = TRUE),
        # Cohort filter (auto-populated)
        selectInput(ns("cohort_filter"), "Select Cohort(s)",
                    choices = NULL, multiple = TRUE),
        # Split by: Cohort or Tissue
        selectInput(ns("split_by"), "Split Plot By",
                    choices = c("Cohort", "Tissue"), selected = "Cohort")
      ),
      mainPanel(
        plotOutput(ns("boxplot")),
        br(),
        downloadButton(ns("downloadData"), "Download Data (CSV)"),
        downloadButton(ns("downloadPlot"), "Download Plot (SVG)")
      )
    )
  )
}

# boxplotServer.R
boxplotServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load required datasets and color mapping.
    sampl.df <- readRDS("./data/sampl.df.rds")
    cptac_sampl <- readRDS("./data/cptac_samples_info_across_datatypes_v1.1.rds") %>%
      filter(include == "TRUE" & datatype == "protein") %>%
      select(harmonized_data_id, tissue, type, cohort) %>%
      unique()
    mat <- readRDS("./data/mat.rds")
    mat_prot <- readRDS("./data/mat_prot.rds")
    
    colors_tab <- read.csv("./data/colors.csv")
    colors <- setNames(colors_tab$color, colors_tab$id)

    
    tissue_colors <- colors_tab %>%
      filter(description == "tissue") %>%
      select(id, color) %>%
      deframe()
    # Update cohort choices based on available cohorts.
    observe({
      cohorts <- unique(cptac_sampl$cohort)
      updateSelectInput(session, "cohort_filter", choices = cohorts, selected = cohorts)
    })
    
    # Auto-populate protein choices based on available proteins in the TSI matrix.
    observe({
      protein_choices <- sort(unique(rownames(mat)))
      updateSelectInput(session, "protein", choices = protein_choices, selected = "CRBN")
    })
    
    # Reactive expression for TSI data.
    tsi_data <- reactive({
      req(input$protein)
      gene_data <- mat %>%
        rankGenes() %>%
        apply(., 2, function(x) x / max(x, na.rm = TRUE)) %>%
        as.data.frame() %>%
        rownames_to_column(var = "gene") %>%
        filter(gene == input$protein) %>%
        pivot_longer(cols = -gene,
                     names_to = "harmonized_data_id",
                     values_to = "value") %>%
        inner_join(sampl.df, by = "harmonized_data_id")
      gene_data
    })
    
    # Reactive expression for Proteome data.
    proteome_data <- reactive({
      req(input$protein)
      gene_data_prot <- as.data.frame(mat_prot) %>%
        rownames_to_column(var = "gene") %>%
        filter(gene == input$protein) %>%
        pivot_longer(cols = -gene,
                     names_to = "harmonized_data_id",
                     values_to = "value") %>%
        inner_join(cptac_sampl, by = "harmonized_data_id")
      gene_data_prot
    })
    
    # Choose the appropriate data based on the data type selection.
    selected_data <- reactive({
      if (input$data_type == "TSI") {
        tsi_data()
      } else {
        proteome_data()
      }
    })
    
    # Filter the selected data by sample type and cohort selection.
    filtered_data <- reactive({
      data <- selected_data()
      if (!is.null(input$sample_type) && length(input$sample_type) > 0) {
        data <- data %>% filter(type %in% input$sample_type)
      }
      if (!is.null(input$cohort_filter) && length(input$cohort_filter) > 0) {
        data <- data %>% filter(cohort %in% input$cohort_filter)
      }
      data
    })
    
    # Reactive expression for the ggplot object.
    # In your plot_obj reactive expression:
    plot_obj <- reactive({
      data <- filtered_data()
      req(nrow(data) > 0)
      
      y_label <- ifelse(input$data_type == "TSI",
                        "Tissue Specificity Index",
                        "Normalized Proteome Expression")
      plot_title <- paste(input$protein, paste(input$sample_type, collapse = ", "),
                          "(", input$data_type, ")")
      
      y_min <- min(data$value, na.rm = TRUE)
      y_max <- max(data$value, na.rm = TRUE)
      
      if (input$split_by == "Cohort") {
        # Group by cohort and type for sample sizes
        sample_sizes <- data %>%
          group_by(cohort, type) %>%
          summarise(n = n(), min_value = min(value, na.rm = TRUE), .groups = "drop")
        
        # Create the base boxplot
        p <- ggplot(data, aes(x = cohort, y = value, fill = type)) +
          geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(aes(color = type),
                      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                      alpha = 0.6, size = 1.2) +
          geom_text(data = sample_sizes,
                    aes(x = cohort, y = min_value, label = glue("n={n}"), group = type),
                    position = position_dodge(width = 0.75), size = 3, vjust = 0) +
          labs(x = "Cohort", y = y_label, title = plot_title) +
          scale_fill_manual(values = colors[unique(data$type)]) +
          scale_color_manual(values = colors[unique(data$type)])
        
        # Identify cohorts with both Tumor and Normal samples for statistical tests
        valid_cohorts <- data %>%
          group_by(cohort) %>%
          filter(all(c("Tumor", "Normal") %in% type)) %>%
          pull(cohort) %>%
          unique()
        
        # Add statistical annotations for valid cohorts
        if (length(valid_cohorts) > 0 && length(unique(data$type)) > 1) {
          valid_data <- data %>% filter(cohort %in% valid_cohorts)
          
          # Kruskal-Wallis test with significance stars
          p <- p +
            stat_compare_means(data = valid_data,
                               method = "kruskal.test",
                               label = "p.signif",  # Use "p.format" for p-value instead
                               label.y = y_max)  # Position at the top of the plot
          
          # Pairwise Wilcoxon tests (already using stars)
          comparisons <- list(c("Tumor", "Normal"))
          if (all(comparisons[[1]] %in% unique(valid_data$type))) {
            p <- p + stat_compare_means(data = valid_data,
                                        comparisons = comparisons,
                                        aes(group = type),
                                        label = "p.signif",
                                        method = "wilcox.test",
                                        label.y = y_max * 0.9)  # Slightly below Kruskal-Wallis
          }
        }
      }
      
      # Finalize the plot
      p + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_continuous(limits = c(y_min, y_max * 1.05))
    })
    
    # Render the plot.
    output$boxplot <- renderPlot({
      print(plot_obj())
    })
    
    # Download handler for the filtered data.
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0(input$protein, "_", input$data_type, "_", 
               paste(input$sample_type, collapse = "_"), "_", 
               paste(input$cohort_filter, collapse = "_"), ".csv")
      },
      content = function(file) {
        data <- filtered_data()
        write.csv(data, file, row.names = FALSE)
      }
    )
    
    # Download handler for the plot (SVG format).
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(input$protein, "_", input$data_type, "_", 
               paste(input$sample_type, collapse = "_"), "_", 
               paste(input$cohort_filter, collapse = "_"), ".svg")
      },
      content = function(file) {
        ggsave(file, plot = plot_obj(), device = "svg", width = 8, height = 6)
      }
    )
    
  })
}