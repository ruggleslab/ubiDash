scatterPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          # Left column for inputs:
          # Display descriptive text first
          column(
            12,
            p(tab_summaries$`Pancan Gene vs Gene`),
            tags$small(
              tags$strong("Note:"),
              " Each substrate generates a separate plot showing Spearman correlation with statistical significance.",
              style = "color: #6c757d; display: block; margin-bottom: 10px;"
            )
          )
        ),
        fluidRow(
          tags$img(
            src = "assets/E3~substrate_cartoon.png",
            style = "width:70%; height:auto; display:block; margin:1px 0px 0px 0px;"
          )
        ),
        fluidRow(
          column(
            12,
            selectizeInput(
              ns("ups_class"),
              "Class of UPS",
              choices = NULL
            ),
            # Static label for the second input.
            selectizeInput(
              ns("e3"),
              "Protein of Interest",
              choices = NULL # Updated on the server based on UPS class.
            ),
            selectizeInput(
              ns("cs"),
              "Candidate Substrate Genes",
              choices = NULL, # Updated on the server from the full UPS list.
              multiple = TRUE # Allow multiple selections.
            ),
            # New cohort filter input; default: all cohorts selected.
            selectizeInput(
              ns("cohorts"),
              "Cohort(s)",
              choices = NULL,
              multiple = TRUE
            )
          )
        )
      ),
      mainPanel(
        conditionalPanel(
          condition = sprintf("output['%s']", ns("data_loaded")),
          uiOutput(ns("scatterPlotsUI"))
        ),
        conditionalPanel(
          condition = sprintf("!output['%s']", ns("data_loaded")),
          div(
            style = "text-align: center; padding: 50px;",
            shinycssloaders::withSpinner(
              div(style = "height: 200px;"),
              color = "#6A2798",
              type = 4,
              size = 2
            ),
            h4("Loading data and initializing inputs...", style = "color: #6A2798;")
          )
        )
      )
    )
  )
}
scatterPlotServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Lazy load data containers
    sampl_data <- reactiveVal(NULL)
    prot_data <- reactiveVal(NULL)
    ups_data_val <- reactiveVal(NULL)
    processed_prot <- reactiveVal(NULL)
    processed_samples <- reactiveVal(NULL)
    processed_colors <- reactiveVal(NULL)

    # Load data only when first needed
    load_scatter_data <- function() {
      if (is.null(sampl_data())) {
        sampl_ <- readRDS("./data/cptac_plus_samples_info_across_datatypes_v1.1.rds")
        sampl_data(sampl_)
      }
      if (is.null(prot_data())) {
        prot_ <- readRDS("./data/CPTAC_proteome.rds")
        prot_data(prot_)
      }
      if (is.null(ups_data_val())) {
        ups_data <- readRDS("./data/Human_Proteostasis_Network.rds")
        ups_data_val(ups_data)
      }
    }

    # Process data only when needed
    get_processed_data <- function() {
      if (is.null(processed_prot())) {
        # Load pre-scaled proteome data (now correctly formatted)
        prot <- readRDS("./data/CPTAC_proteome_scaled.rds")
        processed_prot(prot)
      }

      if (is.null(processed_samples())) {
        # Load sample data first
        if (is.null(sampl_data())) {
          sampl_ <- readRDS("./data/cptac_plus_samples_info_across_datatypes_v1.1.rds")
          sampl_data(sampl_)
        }

        samples <- sampl_data() %>%
          dplyr::filter(
            datatype == "protein",
            harmonized_data_id %in% rownames(processed_prot())
          ) %>%
          dplyr::select(harmonized_data_id, cohort) %>%
          unique()

        processed_samples(samples)
      }

      if (is.null(processed_colors())) {
        colors <- read.csv("./data/colors.csv") %>%
          dplyr::filter(id %in% unique(processed_samples()$cohort)) %>%
          {
            setNames(.$color, .$id)
          }

        processed_colors(colors)
      }

      list(
        prot = processed_prot(),
        samples = processed_samples(),
        colors = processed_colors()
      )
    }
    ## ----------------------------
    ## Update UPS Class, Protein of Interest, and Cohort Filter Inputs
    ## ----------------------------
    # Track loading state
    data_loaded <- reactiveVal(FALSE)

    output$data_loaded <- reactive({
      data_loaded()
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

    # Load UPS data and update inputs when module is first accessed
    observe({
      tryCatch(
        {
          load_scatter_data()
          ups_data <- ups_data_val()

          ups_class_selection <- ups_data %>%
            dplyr::filter(Branch == "Ubiquitin Proteasome System") %>%
            dplyr::pull(Class) %>%
            unique()

          updateSelectizeInput(
            session, "ups_class",
            choices = ups_class_selection,
            selected = "E3 ubiquitin and UBL ligases",
            server = TRUE
          )

          # Also initialize cohort choices here to avoid separate observe
          data <- get_processed_data()
          cohort_choices <- sort(unique(data$samples$cohort))
          updateSelectizeInput(
            session, "cohorts",
            choices = cohort_choices,
            selected = cohort_choices,
            server = TRUE
          )

          # Mark data as loaded
          data_loaded(TRUE)
        },
        error = function(e) {
          showNotification(paste("Error loading data:", e$message), type = "error")
          data_loaded(FALSE)
        }
      )
    })

    # Update the "Protein of Interest" input based on the selected UPS class.
    observeEvent(input$ups_class, {
      req(input$ups_class)
      load_scatter_data()
      ups_data <- ups_data_val()

      current_genes <- ups_data %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System", Class == input$ups_class) %>%
        dplyr::pull(Gene.Symbol) %>%
        unique()

      updateSelectizeInput(
        session, "e3",
        choices = current_genes,
        selected = if ("SKP2" %in% current_genes) "SKP2" else if (length(current_genes) > 0) current_genes[1] else NULL,
        server = TRUE
      )

      # Update candidate substrate genes whenever the UPS class changes.
      full_candidates <- sort(unique(ups_data %>% dplyr::pull(Gene.Symbol)))
      default_cs <- if ("CDKN1B" %in% full_candidates) {
        "CDKN1B"
      } else if (length(full_candidates) > 0) {
        full_candidates[1]
      } else {
        character(0)
      }
      updateSelectizeInput(
        session, "cs",
        choices = full_candidates,
        selected = default_cs
      )
    })

    # Remove the protein of interest gene from the candidate list when selected.
    observeEvent(input$e3, {
      req(input$e3)
      load_scatter_data()
      ups_data <- ups_data_val()

      full_candidates <- sort(unique(ups_data %>% dplyr::pull(Gene.Symbol)))
      candidate_choices <- setdiff(full_candidates, input$e3)
      new_selected <- setdiff(input$cs, input$e3)

      updateSelectizeInput(
        session, "cs",
        choices = candidate_choices,
        selected = new_selected,
        server = TRUE
      )
    })

    ## ----------------------------
    ## REACTIVE PLOTS GENERATION
    ## ----------------------------
    plots_list <- reactiveVal(list())

    observe({
      req(input$e3, input$cs)

      # Get processed data (isolate to prevent reactive cascade)
      data <- isolate(get_processed_data())
      prot <- data$prot
      samples <- data$samples
      colors <- data$colors

      e3_gene <- trimws(input$e3)
      cs_genes <- input$cs
      target_genes <- c(e3_gene, cs_genes)

      df <- prot %>%
        tibble::rownames_to_column("harmonized_data_id") %>%
        dplyr::right_join(samples, by = "harmonized_data_id") %>%
        tibble::column_to_rownames("harmonized_data_id") %>%
        dplyr::select(cohort, dplyr::matches(paste0("^(", paste(target_genes, collapse = "|"), ")")))

      # Use only the cohorts selected by the user.
      selected_cohorts <- if (is.null(input$cohorts) || length(input$cohorts) == 0) {
        sort(unique(samples$cohort))
      } else {
        input$cohorts
      }

      plot_list <- list()

      for (gene in cs_genes) {
        # Check if gene exists in dataframe
        if (!gene %in% colnames(df)) {
          plot_list[[gene]] <- ggplot2::ggplot() +
            ggplot2::annotate("text",
              x = 0.5, y = 0.5,
              label = paste("Gene", gene, "not found in dataset"),
              size = 5, hjust = 0.5, vjust = 0.5
            ) +
            ggplot2::theme_void()
          next
        }

        plots_cohort <- lapply(selected_cohorts, function(current_cohort) {
          df_cohort <- df %>% dplyr::filter(cohort == current_cohort)

          if (nrow(df_cohort) == 0 ||
            !e3_gene %in% colnames(df_cohort) ||
            !gene %in% colnames(df_cohort) ||
            sum(is.finite(df_cohort[[e3_gene]]) & is.finite(df_cohort[[gene]])) < 3) {
            ggplot2::ggplot() +
              ggplot2::annotate("text",
                x = 0.5, y = 0.5,
                label = paste("Cohort", current_cohort, ":\nNot enough data available"),
                size = 4, hjust = 0.5, vjust = 0.5
              ) +
              ggplot2::theme_void() +
              ggplot2::ggtitle(current_cohort) +
              ggplot2::theme(aspect.ratio = 1)
          } else {
            cor_value <- tryCatch(
              {
                cor.test(df_cohort[[e3_gene]], df_cohort[[gene]], method = "spearman")$estimate
              },
              error = function(e) NA
            )

            tryCatch(
              {
                ggplot2::ggplot(df_cohort, ggplot2::aes_string(x = e3_gene, y = gene)) +
                  ggplot2::geom_point(color = colors[current_cohort], size = 3, alpha = 0.7) +
                  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black") +
                  ggpubr::stat_cor(
                    method = "spearman",
                    cor.coef.name = "rho",
                    label.y = max(df_cohort[[gene]], na.rm = TRUE) * 0.9
                  ) +
                  ggplot2::ggtitle(glue::glue("{current_cohort} (rho: {ifelse(is.na(cor_value), 'NA', round(cor_value,2))})")) +
                  ggplot2::theme(aspect.ratio = 1)
              },
              error = function(e) {
                ggplot2::ggplot() +
                  ggplot2::annotate("text",
                    x = 0.5, y = 0.5,
                    label = paste("Error plotting", current_cohort),
                    size = 4, hjust = 0.5, vjust = 0.5
                  ) +
                  ggplot2::theme_void() +
                  ggplot2::ggtitle(current_cohort)
              }
            )
          }
        })

        # Filter out NULL plots and validate before wrapping
        plots_cohort <- Filter(function(x) !is.null(x) && inherits(x, "gg"), plots_cohort)

        if (length(plots_cohort) > 0) {
          p_final <- tryCatch(
            {
              p <- patchwork::wrap_plots(plots_cohort, ncol = ceiling(length(plots_cohort) / 2))
              p + patchwork::plot_annotation(title = glue::glue("{e3_gene} vs {gene}"))
            },
            error = function(e) {
              ggplot2::ggplot() +
                ggplot2::annotate("text",
                  x = 0.5, y = 0.5,
                  label = paste("Error combining plots for", gene),
                  size = 5, hjust = 0.5
                ) +
                ggplot2::theme_void()
            }
          )
          plot_list[[gene]] <- p_final
        } else {
          plot_list[[gene]] <- ggplot2::ggplot() +
            ggplot2::annotate("text",
              x = 0.5, y = 0.5,
              label = paste("No valid plots for", gene),
              size = 5
            ) +
            ggplot2::theme_void()
        }
      }

      plots_list(plot_list)

      output$scatterPlotsUI <- renderUI({
        plot_output_list <- lapply(names(plots_list()), function(gene) {
          plotname <- paste0("plot_", gene)
          # For each candidate gene we output the plot plus two download buttons.
          tagList(
            shinycssloaders::withSpinner(
              plotOutput(ns(plotname)),
              color = "#6A2798"
            ),
            downloadButton(ns(paste0("downloadData_", gene)), label = "Download Data"),
            downloadButton(ns(paste0("downloadPDF_", gene)), label = "Download PDF Plot"),
            tags$hr()
          )
        })
        do.call(tagList, plot_output_list)
      })

      for (gene in names(plots_list())) {
        local({
          my_gene <- gene
          # Render the plot
          output[[paste0("plot_", my_gene)]] <- renderPlot({
            plots_list()[[my_gene]]
          })
          # Download handler for data (CSV)
          output[[paste0("downloadData_", my_gene)]] <- downloadHandler(
            filename = function() {
              paste0("data_", my_gene, "_", Sys.Date(), ".csv")
            },
            content = function(file) {
              # Extract the columns for the protein of interest and the candidate gene.
              # Note: We assume that the column names in df match the gene symbols.
              data_to_save <- df %>% dplyr::select(cohort, !!e3_gene, !!my_gene)
              write.csv(data_to_save, file, row.names = FALSE)
            }
          )
          # Download handler for SVG plot
          output[[paste0("downloadPDF_", my_gene)]] <- downloadHandler(
            filename = function() {
              paste0("plot_", my_gene, "_", Sys.Date(), ".pdf")
            },
            content = function(file) {
              ggplot2::ggsave(file, plot = plots_list()[[my_gene]], device = "pdf")
            }
          )
        })
      }
    })
  })
}
