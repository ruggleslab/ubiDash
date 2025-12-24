## ------------------------------------------------------------------------
## UI
## ------------------------------------------------------------------------
drugResponseUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Response Analysis"),
      p(tab_summaries$`CCLE Drug Sensitivity`),
      # Display the schematic explanation figure
      tags$img(
        src = "assets/Gonzalez-Robles-Figures-5-02.png",
        style = "width:100%; height:auto; display:block; margin:15px 0px;"
      ),
      radioButtons(
        inputId = ns("mode"),
        label = "Analysis mode:",
        choices = c(
          "Drug-centric"             = "drug",
          "UPS-centric"              = "ups",
          "Expression (box+density)" = "dist"
        )
      ),

      # Stratification option for Drug-centric analysis
      conditionalPanel(
        condition = sprintf("input['%s'] == 'drug'", ns("mode")),
        radioButtons(
          inputId = ns("stratify_by"),
          label = "Stratify by:",
          choices = c("Tissue" = "tissue", "Cancer Type" = "Cancer_type"),
          selected = "tissue"
        )
      ),

      ## drug-centric inputs
      conditionalPanel(
        condition = sprintf("input['%s'] == 'drug'", ns("mode")),
        selectizeInput(ns("drug"), "Drug", choices = NULL, options = list(`live-search` = TRUE))
      ),

      ## ups-centric inputs
      conditionalPanel(
        condition = sprintf("input['%s'] == 'ups'", ns("mode")),
        selectizeInput(ns("e3"), "E3 (UPS) gene", choices = NULL, options = list(`live-search` = TRUE)),
        selectizeInput(ns("tissue_ups"), "Tissue (optional)", choices = NULL, options = list(
          `live-search`     = TRUE,
          placeholder       = "All tissues (leave blank)",
          allowEmptyOption  = TRUE
        ), selected = NULL)
      ),

      ## distribution inputs
      conditionalPanel(
        condition = sprintf("input['%s'] == 'dist'", ns("mode")),
        selectizeInput(
          ns("e3_dist"),
          "E3 gene",
          choices  = NULL,
          options  = list(`live-search` = TRUE),
          selected = NULL
        ),
        selectizeInput(
          ns("tissue_dist"),
          "Tissue",
          choices  = NULL,
          options  = list(`live-search` = TRUE),
          selected = NULL
        ),
        selectizeInput(
          ns("cs"),
          "Co-factor (optional)",
          choices = NULL,
          options = list(
            `live-search`     = TRUE,
            placeholder       = "None (leave blank)",
            allowEmptyOption  = TRUE
          ),
          selected = NULL
        ),
        selectizeInput(
          ns("mut"),
          "Mutation (optional)",
          choices = NULL,
          options = list(
            `live-search`     = TRUE,
            placeholder       = "None (leave blank)",
            allowEmptyOption  = TRUE
          ),
          selected = NULL
        ),
        selectizeInput(
          ns("drug_dist"),
          "Drug",
          choices  = NULL,
          options  = list(`live-search` = TRUE),
          selected = NULL
        )
      )

      # -- note: no actionButton any more --
    ),
    mainPanel(
      shinycssloaders::withSpinner(
        plotOutput(ns("result_plot"), height = "800px"),
        color = "#6A2798"
      ),
      br(),
      downloadButton(ns("downloadData"), "Download Data (CSV)"),
      downloadButton(ns("downloadPlot"), "Download Plot (PDF)")
    )
  )
}


## ------------------------------------------------------------------------
## Server
## ------------------------------------------------------------------------
drugResponseServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # 1) Load data once
    print("Loading data...")
    cell_tx_LFC_collapsed <- readRDS("./data/Repurposing_Public_24Q2_LFC_COLLAPSED_tgr.rds") %>%
      select(depmap_id, IDs, LFC, Drug.Name, dose, tissue, Cancer_type) %>%
      distinct()
    sample_lfq <- readRDS("./data/ccle_lfq_sample_info.rds") %>%
      select(BROAD_ID, CCLE_ID) %>%
      rename(depmap_id = BROAD_ID)
    protein <- readRDS("./data/ccle_lfq_protein_expression.rds") %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("CCLE_ID") %>%
      left_join(sample_lfq, by = "CCLE_ID") %>%
      select(depmap_id, everything())
    protein_genes <- setdiff(names(protein), "depmap_id")
    ups <- readRDS("./data/Human_Proteostasis_Network.rds") %>%
      filter(Branch == "Ubiquitin Proteasome System") %>%
      select(Gene.Symbol, Class) %>%
      group_by(Gene.Symbol) %>%
      slice(1) %>%
      ungroup()
    mutations_df <- readRDS("./data/OmicsSomaticMutationsMatrixHotspot.rds") %>%
      as.data.frame() %>%
      rownames_to_column("depmap_id")
    print("Data loaded.")
    # 2) Populate selectors
    updateSelectizeInput(
      session, "drug",
      choices = sort(unique(cell_tx_LFC_collapsed$Drug.Name)),
      selected = "ONCRASIN-1",
      server = TRUE
    )
    updateSelectizeInput(
      session, "e3",
      choices = sort(protein_genes),
      server = TRUE
    )
    updateSelectizeInput(
      session, "tissue_ups",
      choices = sort(unique(cell_tx_LFC_collapsed$tissue)),
      selected = character(0),
      options = list(
        `live-search`     = TRUE,
        placeholder       = "All tissues (leave blank)",
        allowEmptyOption  = TRUE
      ),
      server = TRUE
    )
    updateSelectizeInput(
      session, "e3_dist",
      choices = sort(protein_genes),
      selected = "CLPX",
      server = TRUE
    )
    updateSelectizeInput(
      session, "cs",
      choices = sort(protein_genes),
      selected = character(0),
      options = list(
        `live-search`     = TRUE,
        placeholder       = "None (leave blank)",
        allowEmptyOption  = TRUE
      ),
      server = TRUE
    )
    updateSelectizeInput(
      session, "mut",
      choices = setdiff(names(mutations_df), "depmap_id"),
      selected = "PIK3CD (5293)",
      options = list(
        `live-search`     = TRUE,
        placeholder       = "None (leave blank)",
        allowEmptyOption  = TRUE
      ),
      server = TRUE
    )
    updateSelectizeInput(
      session, "drug_dist",
      choices = sort(unique(cell_tx_LFC_collapsed$Drug.Name)),
      selected = "ONCRASIN-1",
      server = TRUE
    )
    updateSelectizeInput(
      session, "tissue_dist",
      choices = sort(unique(cell_tx_LFC_collapsed$tissue)),
      selected = "AUTONOMIC_GANGLIA",
      server = TRUE
    )

    # 3) Reactive data object for downloads
    data_obj <- reactive({
      req(input$mode)

      if (input$mode == "drug") {
        req(input$drug)
        strat_col <- if (!is.null(input$stratify_by)) input$stratify_by else "tissue"
        protein %>%
          select(depmap_id, any_of(ups$Gene.Symbol)) %>%
          pivot_longer(cols = -depmap_id, names_to = "hgnc_symbol", values_to = "prot_expression") %>%
          inner_join(filter(cell_tx_LFC_collapsed, Drug.Name == input$drug), by = "depmap_id") %>%
          filter(is.finite(prot_expression), is.finite(LFC)) %>%
          group_by(IDs, hgnc_symbol, .data[[strat_col]]) %>%
          summarise(
            cor = if (n() >= 3) cor(prot_expression, LFC, method = "spearman") else NA_real_,
            p = if (n() >= 3) cor.test(prot_expression, LFC, method = "spearman")$p.value else NA_real_,
            .groups = "drop"
          ) %>%
          filter(!is.na(cor)) %>%
          mutate(fdr = p.adjust(p, "BH"))
      } else if (input$mode == "ups") {
        req(input$e3)
        tissue_filter <- if (nzchar(input$tissue_ups)) input$tissue_ups else NULL

        df_ups <- protein %>%
          select(depmap_id, !!sym(input$e3)) %>%
          inner_join(cell_tx_LFC_collapsed, by = "depmap_id")

        # Filter by tissue if selected
        if (!is.null(tissue_filter)) {
          df_ups <- df_ups %>% filter(tissue == tissue_filter)
        }

        df_ups %>%
          group_by(IDs, Drug.Name) %>%
          filter(n() >= 3) %>%
          summarise(
            ups = input$e3,
            cor = if (sum(is.finite(!!sym(input$e3)) & is.finite(LFC)) >= 3) {
              cor.test(!!sym(input$e3), LFC, method = "spearman", use = "complete.obs")$estimate
            } else {
              NA_real_
            },
            p = if (sum(is.finite(!!sym(input$e3)) & is.finite(LFC)) >= 3) {
              cor.test(!!sym(input$e3), LFC, method = "spearman", use = "complete.obs")$p.value
            } else {
              NA_real_
            },
            n = sum(is.finite(!!sym(input$e3)) & is.finite(LFC)),
            .groups = "drop"
          ) %>%
          filter(!is.na(cor)) %>%
          arrange(desc(abs(cor))) %>%
          mutate(
            fdr = p.adjust(p, method = "BH"),
            fdr_sig = if_else(fdr < 0.1, "*", "")
          )
      } else if (input$mode == "dist") {
        req(input$e3_dist, input$drug_dist, input$tissue_dist)
        drg <- input$drug_dist
        tis <- input$tissue_dist
        cof <- if (nzchar(input$cs)) input$cs else NULL
        mut <- if (nzchar(input$mut)) input$mut else NULL

        # 1) pick columns & filter to single drug + single tissue
        trgt <- c(input$e3_dist, cof) %>% intersect(names(protein))
        df0 <- cell_tx_LFC_collapsed %>%
          filter(Drug.Name == drg, tissue == tis) %>%
          inner_join(select(protein, depmap_id, all_of(trgt)), by = "depmap_id") %>%
          mutate(across(all_of(trgt),
            ~ if_else(.x > median(.x, na.rm = TRUE), "High", "Low"),
            .names = "{.col}_status"
          ))

        # 2) optionally merge mutation and rename to "Mutation"
        if (!is.null(mut) && mut %in% names(mutations_df)) {
          df0 <- df0 %>%
            left_join(select(mutations_df, depmap_id, !!sym(mut)), by = "depmap_id") %>%
            rename(Mutation = !!sym(mut))
          group_mut <- "Mutation"
        } else {
          group_mut <- NULL
        }

        # 3) build contrast name & grouping vars
        e3_col <- paste0(input$e3_dist, "_status")
        if (!is.null(cof)) {
          cs_col <- paste0(cof, "_status")
          combined_col <- paste0(input$e3_dist, "_", cof, "_status")
          df0 <- df0 %>% mutate(!!combined_col := paste0(.data[[e3_col]], "_", .data[[cs_col]]))
          contrast <- combined_col
          group_vars <- c(contrast, group_mut)
          df0 <- df0 %>%
            filter(
              if_all(all_of(group_vars), ~ !is.na(.)), # drop any NA in status or Mutation
              !is.na(LFC)
            )
        } else {
          contrast <- e3_col
          group_vars <- c(contrast, group_mut)
        }

        df0
      }
    })

    # 4) Reactive plot object
    plot_obj <- reactive({
      req(input$mode)

      if (input$mode == "drug") {
        req(input$drug)
        strat_col <- if (!is.null(input$stratify_by)) input$stratify_by else "tissue"
        df <- protein %>%
          select(depmap_id, any_of(ups$Gene.Symbol)) %>%
          pivot_longer(cols = -depmap_id, names_to = "hgnc_symbol", values_to = "prot_expression") %>%
          inner_join(filter(cell_tx_LFC_collapsed, Drug.Name == input$drug), by = "depmap_id") %>%
          filter(is.finite(prot_expression), is.finite(LFC)) %>%
          group_by(IDs, hgnc_symbol, .data[[strat_col]]) %>%
          summarise(
            cor = if (n() >= 3) cor(prot_expression, LFC, method = "spearman") else NA_real_,
            p = if (n() >= 3) cor.test(prot_expression, LFC, method = "spearman")$p.value else NA_real_,
            .groups = "drop"
          ) %>%
          filter(!is.na(cor)) %>%
          mutate(fdr = p.adjust(p, "BH"))
        print(head(df))
        top_pts <- filter(df, p < 0.05)
        # Use colors for tissue or cancer type
        color_desc <- if (strat_col == "tissue") "tissue" else "Cancer_type"
        colors <- read.csv("./data/colors.csv") %>%
          filter(description == color_desc) %>%
          pull(color) %>%
          setNames(unique(df[[strat_col]]))

        ggplot(df, aes(x = cor, y = -log10(p), color = .data[[strat_col]])) +
          geom_point(alpha = abs(df$cor)) +
          geom_text_repel(data = top_pts, aes(label = hgnc_symbol), max.overlaps = 8) +
          geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
          geom_vline(xintercept = 0) +
          scale_color_manual(values = colors) +
          facet_wrap(~IDs) +
          labs(
            title = paste("Drug-centric:", input$drug),
            x = "Spearman ρ", y = "-log10(P-value)",
            color = if (strat_col == "tissue") "Tissue" else "Cancer Type"
          ) +
          theme_minimal() +
          annotate("text", x = -0.35, y = 0, label = "bold(Sensitive)", parse = TRUE, size = 5) +
          annotate("text", x = 0.35, y = 0, label = "bold(Resistant)", parse = TRUE, size = 5)
      } else if (input$mode == "ups") {
        req(input$e3)
        tissue_filter <- if (nzchar(input$tissue_ups)) input$tissue_ups else NULL

        df_ups <- protein %>%
          select(depmap_id, !!sym(input$e3)) %>%
          inner_join(cell_tx_LFC_collapsed, by = "depmap_id")

        # Filter by tissue if selected
        if (!is.null(tissue_filter)) {
          df_ups <- df_ups %>% filter(tissue == tissue_filter)
        }

        df <- df_ups %>%
          group_by(IDs, Drug.Name) %>%
          filter(n() >= 3) %>%
          summarise(
            ups = input$e3,
            cor = if (sum(is.finite(!!sym(input$e3)) & is.finite(LFC)) >= 3) {
              cor.test(!!sym(input$e3), LFC, method = "spearman", use = "complete.obs")$estimate
            } else {
              NA_real_
            },
            p = if (sum(is.finite(!!sym(input$e3)) & is.finite(LFC)) >= 3) {
              cor.test(!!sym(input$e3), LFC, method = "spearman", use = "complete.obs")$p.value
            } else {
              NA_real_
            },
            n = sum(is.finite(!!sym(input$e3)) & is.finite(LFC)),
            .groups = "drop"
          ) %>%
          filter(!is.na(cor)) %>%
          arrange(desc(abs(cor))) %>%
          mutate(
            fdr = p.adjust(p, method = "BH"),
            fdr_sig = if_else(fdr < 0.1, "*", "")
          )

        top_points <- filter(df, p < 0.05)

        ggplot(df, aes(x = cor, y = -log10(p))) +
          annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
          annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
          geom_point(alpha = abs(df$cor)) +
          geom_text_repel(data = top_points, aes(label = Drug.Name), max.overlaps = 10) +
          geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
          geom_vline(xintercept = 0) +
          labs(
            title = if (!is.null(tissue_filter)) paste(input$e3, "-", tissue_filter) else input$e3,
            x = "Spearman, ρ",
            y = expression(-log[10]("P-value"))
          ) +
          theme_minimal() +
          annotate("text", x = -0.25, y = max(-log10(df$p), na.rm = TRUE) * 0.95, label = "bold(Sensitive)", parse = TRUE, size = 5, color = "blue") +
          annotate("text", x = 0.25, y = max(-log10(df$p), na.rm = TRUE) * 0.95, label = "bold(Resistant)", parse = TRUE, size = 5, color = "red")
      } else if (input$mode == "dist") {
        req(input$e3_dist, input$drug_dist, input$tissue_dist)
        drg <- input$drug_dist
        tis <- input$tissue_dist
        cof <- if (nzchar(input$cs)) input$cs else NULL
        mut <- if (nzchar(input$mut)) input$mut else NULL

        # 1) pick columns & filter to single drug + single tissue
        trgt <- c(input$e3_dist, cof) %>% intersect(names(protein))
        df0 <- cell_tx_LFC_collapsed %>%
          filter(Drug.Name == drg, tissue == tis) %>%
          inner_join(select(protein, depmap_id, all_of(trgt)), by = "depmap_id") %>%
          mutate(across(all_of(trgt),
            ~ if_else(.x > median(.x, na.rm = TRUE), "High", "Low"),
            .names = "{.col}_status"
          ))

        # 2) optionally merge mutation and rename to “Mutation”
        if (!is.null(mut) && mut %in% names(mutations_df)) {
          df0 <- df0 %>%
            left_join(select(mutations_df, depmap_id, !!sym(mut)), by = "depmap_id") %>%
            rename(Mutation = !!sym(mut))
          group_mut <- "Mutation"
        } else {
          group_mut <- NULL
        }

        # 3) build contrast name & grouping vars
        e3_col <- paste0(input$e3_dist, "_status")
        if (!is.null(cof)) {
          cs_col <- paste0(cof, "_status")
          combined_col <- paste0(input$e3_dist, "_", cof, "_status")
          df0 <- df0 %>% mutate(!!combined_col := paste0(.data[[e3_col]], "_", .data[[cs_col]]))
          contrast <- combined_col
          group_vars <- c(contrast, group_mut)
          df0 <- df0 %>%
            filter(
              if_all(all_of(group_vars), ~ !is.na(.)), # drop any NA in status or Mutation
              !is.na(LFC)
            )
        } else {
          contrast <- e3_col
          group_vars <- c(contrast, group_mut)
        }

        # 4) shared settings
        lfc_lim <- range(df0$LFC, na.rm = TRUE)
        colors <- c(
          "High" = "blue", "Low" = "orange",
          "High_High" = "#4E79A7", "Low_Low" = "#F28E2B",
          "High_Low" = "#76B7B2", "Low_High" = "#E15759"
        )

        # 5) compute medians for dashed lines
        df_meds <- df0 %>%
          group_by_at(group_vars) %>%
          summarise(mean_val = median(LFC, na.rm = TRUE), .groups = "drop")

        # 6) Density plot
        p1 <- ggplot(df0, aes_string(x = "LFC", color = contrast)) +
          geom_density(alpha = 0.3) +
          geom_vline(
            data = df_meds,
            aes_string(xintercept = "mean_val", color = contrast),
            linetype = "dashed", size = 1
          ) +
          scale_color_manual(values = colors) +
          labs(
            title = if (!is.null(mut)) glue("{tis} | {mut}") else tis,
            x = NULL, y = "Density"
          ) +
          xlim(lfc_lim) +
          facet_wrap(if (!is.null(mut)) vars(Drug.Name, Mutation) else vars(Drug.Name), scales = "fixed") +
          theme_minimal() +
          theme(legend.position = "right")

        # 7) Boxplot + jitter + sample-size annotation
        df_counts <- df0 %>% count(across(all_of(group_vars)), name = "n")
        contrast_lvls <- unique(df0[[contrast]])

        p2 <- ggplot(df0, aes_string(x = contrast, y = "LFC", color = contrast)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2) +
          geom_text(
            data = df_counts,
            aes_string(
              x     = contrast,
              y     = lfc_lim[1] + diff(lfc_lim) * 0.1,
              label = "paste0('n=', n)"
            ),
            inherit.aes = FALSE, size = 3.5
          ) +
          scale_color_manual(values = colors) +
          labs(x = NULL, y = expression(log[2] * "(cell viability)")) +
          ylim(lfc_lim) +
          facet_wrap(if (!is.null(mut)) vars(Drug.Name, Mutation) else vars(Drug.Name), scales = "fixed") +
          coord_flip() +
          theme_minimal() +
          theme(legend.position = "right")

        # 8) add Wilcoxon pairwise comparisons
        if (length(contrast_lvls) > 1) {
          p2 <- p2 +
            stat_compare_means(
              comparisons = combn(contrast_lvls, 2, simplify = FALSE),
              label = "p.format",
              method = "wilcox.test",
              label.y = seq(
                from = max(df0$LFC, na.rm = TRUE) - 0.5,
                to   = max(df0$LFC, na.rm = TRUE) - 2,
                by   = -0.2
              )
            )
        }

        # 9) combine density + boxplot
        ggpubr::ggarrange(
          p1, p2,
          ncol = 1,
          common.legend = TRUE,
          legend = "bottom",
          align = "hv"
        )
      }
    })

    # 5) Render the reactive plot
    output$result_plot <- renderPlot(
      {
        plot_obj()
      },
      height = 600
    )

    # 6) Download handlers
    output$downloadData <- downloadHandler(
      filename = function() {
        if (input$mode == "drug") {
          paste0("drug_response_", input$drug, "_", Sys.Date(), ".csv")
        } else if (input$mode == "ups") {
          paste0("ups_response_", input$e3, "_", Sys.Date(), ".csv")
        } else if (input$mode == "dist") {
          paste0("expression_distribution_", input$e3_dist, "_", input$drug_dist, "_", input$tissue_dist, "_", Sys.Date(), ".csv")
        }
      },
      content = function(file) {
        data <- data_obj()
        write.csv(data, file, row.names = FALSE)
      }
    )

    output$downloadPlot <- downloadHandler(
      filename = function() {
        if (input$mode == "drug") {
          paste0("drug_response_plot_", input$drug, "_", Sys.Date(), ".pdf")
        } else if (input$mode == "ups") {
          paste0("ups_response_plot_", input$e3, "_", Sys.Date(), ".pdf")
        } else if (input$mode == "dist") {
          paste0("expression_distribution_plot_", input$e3_dist, "_", input$drug_dist, "_", input$tissue_dist, "_", Sys.Date(), ".pdf")
        }
      },
      content = function(file) {
        pdf(file, width = 12, height = 8)
        print(plot_obj())
        dev.off()
      }
    )
  })
}
