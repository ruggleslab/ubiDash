boxplotUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p(tab_summaries$Boxplot),
        # Add static significance legend
        p(HTML(
          "<b>Significance Levels:</b><br>",
          "ns: p > 0.05 (not significant)<br>",
          "*: p ≤ 0.05<br>",
          "**: p ≤ 0.01<br>",
          "***: p ≤ 0.001"
        )),
        # Protein of Interest (auto-populated)
        selectInput(ns("protein"), "Protein of Interest", choices = NULL, selected = "CRBN"),
        # Sample type selection: allow Tumor, Normal, or both
        selectInput(ns("sample_type"), "Sample Type",
          choices = c("Tumor", "Normal"),
          selected = c("Tumor", "Normal"), multiple = TRUE
        ),
        # Cohort filter (auto-populated)
        selectInput(ns("cohort_filter"), "Select Cohort(s)",
          choices = NULL, multiple = TRUE
        ),
        # Split by: Cohort or Tissue
        selectInput(ns("split_by"), "Split Plot By",
          choices = c("Cohort", "Tissue"), selected = "Cohort"
        )
      ),
      mainPanel(
        plotOutput(ns("boxplot")),
        br(),
        downloadButton(ns("downloadData"), "Download Data (CSV)"),
        downloadButton(ns("downloadPlot"), "Download Plot (PDF)")
      )
    )
  )
}

# boxplotServer.R
boxplotServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load required datasets and color mapping.
    cptac_sampl <- readRDS("./data/cptac_samples_info_across_datatypes_v1.1.rds") %>%
      filter(include == "TRUE" & datatype == "protein") %>%
      select(harmonized_data_id, tissue, type, cohort) %>%
      unique()
    mat_prot <- readRDS("./data/mat_prot.rds")

    colors_tab <- read.csv("./data/colors.csv")
    colors <- setNames(colors_tab$color, colors_tab$id)


    tissue_colors <- colors_tab %>%
      filter(description == "tissue") %>%
      select(id, color) %>%
      deframe()
    # Update cohort choices based on available cohorts that have data in mat_prot.
    observe({
      # Find samples that exist in both cptac_sampl and mat_prot
      available_samples <- intersect(cptac_sampl$harmonized_data_id, colnames(mat_prot))
      # Get cohorts that have at least one available sample
      available_cohorts <- cptac_sampl %>%
        filter(harmonized_data_id %in% available_samples) %>%
        pull(cohort) %>%
        unique() %>%
        sort()
      updateSelectInput(session, "cohort_filter", choices = available_cohorts, selected = available_cohorts)
    })

    # Auto-populate protein choices based on available proteins in the protein matrix.
    observe({
      protein_choices <- sort(unique(rownames(mat_prot)))
      updateSelectInput(session, "protein", choices = protein_choices, selected = "CRBN")
    })

    # Reactive expression for Protein data.
    selected_data <- reactive({
      req(input$protein)
      gene_data_prot <- as.data.frame(mat_prot) %>%
        rownames_to_column(var = "gene") %>%
        filter(gene == input$protein) %>%
        pivot_longer(
          cols = -gene,
          names_to = "harmonized_data_id",
          values_to = "value"
        ) %>%
        inner_join(cptac_sampl, by = "harmonized_data_id")
      gene_data_prot
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

      y_label <- "log2 Normalized Protein Expression"
      plot_title <- paste(
        input$protein, paste(input$sample_type, collapse = ", "),
        "(Protein)"
      )

      y_min <- min(data$value, na.rm = TRUE)
      y_max <- max(data$value, na.rm = TRUE)

      if (input$split_by == "Cohort") {
        # Group by cohort and type for sample sizes
        sample_sizes <- data %>%
          group_by(cohort, type) %>%
          summarise(
            n = sum(!is.na(value)), # Count only non-missing values
            min_value = ifelse(sum(!is.na(value)) > 0,
              min(value, na.rm = TRUE),
              y_min
            ), # Use y_min if all values are missing
            .groups = "drop"
          )

        # Create the base boxplot
        p <- ggplot(data, aes(x = cohort, y = value, fill = type)) +
          geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(aes(color = type),
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
            alpha = 0.6, size = 1.2
          ) +
          geom_text(
            data = sample_sizes,
            aes(x = cohort, y = min_value, label = glue("n={n}"), group = type),
            position = position_dodge(width = 0.75), size = 3, vjust = 0
          ) +
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
          # Calculate p-values manually for each cohort and add as text annotations
          stat_results <- data %>%
            filter(cohort %in% valid_cohorts) %>%
            group_by(cohort) %>%
            do({
              # Filter out missing values first
              tumor_vals <- .$value[.$type == "Tumor" & !is.na(.$value)]
              normal_vals <- .$value[.$type == "Normal" & !is.na(.$value)]

              # Only perform test if both groups have at least 3 non-missing samples
              if (length(tumor_vals) >= 3 & length(normal_vals) >= 3) {
                tryCatch(
                  {
                    # Wilcoxon rank-sum test (Mann-Whitney U test)
                    test_result <- wilcox.test(tumor_vals, normal_vals, exact = FALSE)
                    p_val <- test_result$p.value

                    # Convert p-value to significance stars
                    sig_label <- case_when(
                      is.na(p_val) ~ "",
                      p_val > 0.05 ~ "ns",
                      p_val <= 0.001 ~ "***",
                      p_val <= 0.01 ~ "**",
                      p_val <= 0.05 ~ "*",
                      TRUE ~ ""
                    )
                    data.frame(p_value = p_val, label = sig_label, y_pos = y_max * 1.02)
                  },
                  error = function(e) {
                    # If test fails for any reason, return empty result
                    data.frame(p_value = NA, label = "", y_pos = y_max * 1.02)
                  }
                )
              } else {
                # Not enough samples for reliable test
                data.frame(p_value = NA, label = "", y_pos = y_max * 1.02)
              }
            })

          # Add text annotations for significance
          if (nrow(stat_results) > 0 && any(stat_results$label != "")) {
            p <- p + geom_text(
              data = stat_results %>% filter(label != ""),
              aes(x = cohort, y = y_pos, label = label),
              size = 4, vjust = 0, inherit.aes = FALSE
            )
          }
        }
      } else if (input$split_by == "Tissue") {
        # Group by tissue and type for sample sizes
        sample_sizes <- data %>%
          group_by(tissue, type) %>%
          summarise(
            n = sum(!is.na(value)), # Count only non-missing values
            min_value = ifelse(sum(!is.na(value)) > 0,
              min(value, na.rm = TRUE),
              y_min
            ), # Use y_min if all values are missing
            .groups = "drop"
          )

        # Create the base boxplot
        p <- ggplot(data, aes(x = tissue, y = value, fill = type)) +
          geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(aes(color = type),
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
            alpha = 0.6, size = 1.2
          ) +
          geom_text(
            data = sample_sizes,
            aes(x = tissue, y = min_value, label = glue("n={n}"), group = type),
            position = position_dodge(width = 0.75), size = 3, vjust = 0
          ) +
          labs(x = "Tissue", y = y_label, title = plot_title) +
          scale_fill_manual(values = colors[unique(data$type)]) +
          scale_color_manual(values = colors[unique(data$type)])

        # Identify tissues with both Tumor and Normal samples for statistical tests
        valid_tissues <- data %>%
          group_by(tissue) %>%
          filter(all(c("Tumor", "Normal") %in% type)) %>%
          pull(tissue) %>%
          unique()

        # Add statistical annotations for valid tissues
        if (length(valid_tissues) > 0 && length(unique(data$type)) > 1) {
          # Calculate p-values manually for each tissue and add as text annotations
          stat_results <- data %>%
            filter(tissue %in% valid_tissues) %>%
            group_by(tissue) %>%
            do({
              # Filter out missing values first
              tumor_vals <- .$value[.$type == "Tumor" & !is.na(.$value)]
              normal_vals <- .$value[.$type == "Normal" & !is.na(.$value)]

              # Only perform test if both groups have at least 3 non-missing samples
              if (length(tumor_vals) >= 3 & length(normal_vals) >= 3) {
                tryCatch(
                  {
                    # Wilcoxon rank-sum test (Mann-Whitney U test)
                    test_result <- wilcox.test(tumor_vals, normal_vals, exact = FALSE)
                    p_val <- test_result$p.value

                    # Convert p-value to significance stars
                    sig_label <- case_when(
                      is.na(p_val) ~ "",
                      p_val > 0.05 ~ "ns",
                      p_val <= 0.001 ~ "***",
                      p_val <= 0.01 ~ "**",
                      p_val <= 0.05 ~ "*",
                      TRUE ~ ""
                    )
                    data.frame(p_value = p_val, label = sig_label, y_pos = y_max * 1.02)
                  },
                  error = function(e) {
                    # If test fails for any reason, return empty result
                    data.frame(p_value = NA, label = "", y_pos = y_max * 1.02)
                  }
                )
              } else {
                # Not enough samples for reliable test
                data.frame(p_value = NA, label = "", y_pos = y_max * 1.02)
              }
            })

          # Add text annotations for significance
          if (nrow(stat_results) > 0 && any(stat_results$label != "")) {
            p <- p + geom_text(
              data = stat_results %>% filter(label != ""),
              aes(x = tissue, y = y_pos, label = label),
              size = 4, vjust = 0, inherit.aes = FALSE
            )
          }
        }
      }

      # Finalize the plot
      p + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_continuous(limits = c(y_min, y_max * 1.1)) # More space for stars
    })

    # Render the plot.
    output$boxplot <- renderPlot({
      print(plot_obj())
    })

    # Download handler for the filtered data.
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0(
          input$protein, "_Protein_",
          paste(input$sample_type, collapse = "_"), "_",
          paste(input$cohort_filter, collapse = "_"), ".csv"
        )
      },
      content = function(file) {
        data <- filtered_data()
        write.csv(data, file, row.names = FALSE)
      }
    )

    # Download handler for the plot (PDF format).
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(
          input$protein, "_Protein_",
          paste(input$sample_type, collapse = "_"), "_",
          paste(input$cohort_filter, collapse = "_"), ".pdf"
        )
      },
      content = function(file) {
        ggsave(file, plot = plot_obj(), device = "pdf", width = 8, height = 6)
      }
    )
  })
}
