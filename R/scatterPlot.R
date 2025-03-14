scatterPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          # Left column for inputs:
          # Display the full PNG image (adjust style as needed)
          tags$img(src = "assets/E3~substrate_cartoon.png",
                   style = "width:100%; height:auto; display:block; margin-top:15px;"),
          column(8,
                 helpText("Select a UPS class, a protein of interest, candidate substrate genes, and cohort(s) from the available list."),
                 selectizeInput(
                   ns("ups_class"), 
                   "Class of UPS",
                   choices = NULL
                 ),
                 # Static label for the second input.
                 selectizeInput(
                   ns("e3"), 
                   "Protein of Interest", 
                   choices = NULL  # Updated on the server based on UPS class.
                 ),
                 selectizeInput(
                   ns("cs"), 
                   "Candidate Substrate Genes", 
                   choices = NULL,      # Updated on the server from the full UPS list.
                   multiple = TRUE      # Allow multiple selections.
                 ),
                 # New cohort filter input; default: all cohorts selected.
                 selectizeInput(
                   ns("cohorts"),
                   "Cohort(s)",
                   choices = NULL,
                   multiple = TRUE
                 ),
                 actionButton(ns("generatePlots"), "Generate Scatter Plots")
          )
        )
      ),
      mainPanel(
        uiOutput(ns("scatterPlotsUI"))
      )
    )
  )
}
scatterPlotServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    ## ----------------------------
    ## Data Loading & Pre-processing
    ## ----------------------------
    sampl_ <- readRDS("./data/cptac_plus_samples_info_across_datatypes_v1.1.rds")
    
    prot_ <- read.csv("./data/Combine_PanCan_CPTAC_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed.tsv",
                      sep = "\t", check.names = FALSE)
    
    ups_data <- openxlsx::read.xlsx(
      xlsxFile = './data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx',
      sheet = 'dense'
    )
    
    prot <- prot_ %>% 
      dplyr::select(-c(Gene, GENECODE34_Symbol)) %>%
      tibble::column_to_rownames(var = "HGNC_Approved_Symbol") %>% 
      sweep(., 2, apply(., 2, median, na.rm = TRUE)) %>%  
      t() %>% 
      scale() %>%   
      as.data.frame()
    
    samples <- sampl_ %>% 
      dplyr::filter(datatype == "protein",
                    harmonized_data_id %in% rownames(prot)) %>% 
      dplyr::select(harmonized_data_id, cohort) %>% 
      unique()
    
    colors <- read.csv("./data/colors.csv") %>% 
      dplyr::filter(id %in% unique(samples$cohort)) %>% 
      { setNames(.$color, .$id) }
    print("colors scatter")
    print(colors)
    ## ----------------------------
    ## Update UPS Class, Protein of Interest, and Cohort Filter Inputs
    ## ----------------------------
    ups_class_selection <- ups_data %>%
      dplyr::filter(Branch == 'Ubiquitin Proteasome System') %>%
      dplyr::pull(Class) %>%
      unique()
    
    updateSelectizeInput(
      session, "ups_class",
      choices = ups_class_selection,
      selected = "E3 ubiquitin and UBL ligases",  
      server = TRUE
    )
    
    # Update the "Protein of Interest" input based on the selected UPS class.
    observeEvent(input$ups_class, {
      req(input$ups_class)
      
      current_genes <- ups_data %>%
        dplyr::filter(Branch == 'Ubiquitin Proteasome System', Class == input$ups_class) %>%
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
      updateSelectizeInput(
        session, "cs",
        choices = full_candidates,
        selected = if ("CDKN1B" %in% full_candidates) "CDKN1B" else character(0)
      )
    })
    
    # Remove the protein of interest gene from the candidate list when selected.
    observeEvent(input$e3, {
      req(input$e3)
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
    
    # Update cohort filter choices based on the samples data.
    observe({
      cohort_choices <- sort(unique(samples$cohort))
      updateSelectizeInput(
        session, "cohorts",
        choices = cohort_choices,
        selected = cohort_choices,  # Default selection: all cohorts
        server = TRUE
      )
    })
    
    ## ----------------------------
    ## REACTIVE PLOTS GENERATION
    ## ----------------------------
    plots_list <- reactiveVal(list())
    
    observeEvent(input$generatePlots, {
      req(input$e3, input$cs)
      
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
        plots_cohort <- lapply(selected_cohorts, function(current_cohort) {
          df_cohort <- df %>% dplyr::filter(cohort == current_cohort)
          
          if (nrow(df_cohort) == 0 || 
              sum(is.finite(df_cohort[[e3_gene]]) & is.finite(df_cohort[[gene]])) < 3) {
            ggplot2::ggplot() +
              ggplot2::annotate("text", x = 0.5, y = 0.5, 
                                label = paste("Cohort", current_cohort, ":\nNot enough data available"),
                                size = 4, hjust = 0.5, vjust = 0.5) +
              ggplot2::theme_void() +
              ggplot2::ggtitle(current_cohort) +
              ggplot2::theme(aspect.ratio = 1)
          } else {
            cor_value <- tryCatch({
              cor.test(df_cohort[[e3_gene]], df_cohort[[gene]], method = "spearman")$estimate
            }, error = function(e) NA)
            
            ggplot2::ggplot(df_cohort, ggplot2::aes_string(x = e3_gene, y = gene)) +
              ggplot2::geom_point(color = colors[current_cohort], size = 3, alpha = 0.7) +
              ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black") +
              ggpubr::stat_cor(method = "spearman",
                               cor.coef.name = "rho",
                               label.y = max(df_cohort[[gene]], na.rm = TRUE) * 0.9) +
              ggplot2::ggtitle(glue::glue("{current_cohort} (rho: {ifelse(is.na(cor_value), 'NA', round(cor_value,2))})")) +
              ggplot2::theme(aspect.ratio = 1)
          }
        })
        
        p_final <- patchwork::wrap_plots(plots_cohort, ncol = ceiling(length(plots_cohort)/2))
        p_final <- p_final + patchwork::plot_annotation(
          title = glue::glue("{e3_gene} vs {gene}")
        )
        plot_list[[gene]] <- p_final
      }
      
      plots_list(plot_list)
      
      output$scatterPlotsUI <- renderUI({
        plot_output_list <- lapply(names(plots_list()), function(gene) {
          plotname <- paste0("plot_", gene)
          # For each candidate gene we output the plot plus two download buttons.
          tagList(
            plotOutput(ns(plotname)),
            downloadButton(ns(paste0("downloadData_", gene)), label = "Download Data"),
            downloadButton(ns(paste0("downloadSVG_", gene)), label = "Download SVG Plot"),
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
          output[[paste0("downloadSVG_", my_gene)]] <- downloadHandler(
            filename = function() {
              paste0("plot_", my_gene, "_", Sys.Date(), ".svg")
            },
            content = function(file) {
              ggplot2::ggsave(file, plot = plots_list()[[my_gene]], device = "svg")
            }
          )
        })
      }
    })
  })
}