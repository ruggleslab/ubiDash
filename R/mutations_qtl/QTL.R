QTLUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Quantitative Trait Locus (QTL) Analysis"),
      p(tab_summaries$`QTL Analysis`),
      p("Plots may take over a minute to load."),
      p("The data has been filtered by the hDC, hChR, and hTF genesets to fit within the application. For this reason, the volcano plot may be blank for some inputs."),
      selectizeInput(
        inputId = ns("mutgene"),
        choices = NULL,
        multiple = FALSE,
        label = "Mutant Gene",
        options = list(`live-search` = TRUE)
      ),
      selectizeInput(
        inputId = ns("gene"),
        choices = NULL,
        multiple = FALSE,
        label = "UPS Protein",
        options = list(`live-search` = TRUE)
      ),
      selectizeInput(
        inputId = ns("genelist"),
        choices = NULL,
        selected = "hCD",
        multiple = FALSE,
        label = "Gene List (for highlighting)",
        options = list(`live-search` = TRUE)
      ),
      fluidRow(
        column(6, sliderInput(ns("pvalcutoff"), "-Log10(P-value):",
          min = 1, max = 15, value = 7, step = 0.5
        )),
        column(6, sliderInput(ns("n_genes"), "Min Samples/Group:",
          min = 0, max = 1000, value = 1, step = 1
        ))
      ),
      actionButton(ns("refresh"), "Generate Plots", class = "btn-primary"),
      hr(),
      tags$div(
        style = "font-size: 11px; color: #666; margin-top: 10px;",
        tags$strong("Boxplot Significance Legend:"),
        tags$br(),
        "ns: p > 0.05", tags$br(),
        "*: p ≤ 0.05", tags$br(),
        "**: p ≤ 0.01", tags$br(),
        "***: p ≤ 0.001", tags$br(),
        "****: p ≤ 0.0001"
      )
    ),
    mainPanel(
      tags$img(
        src = "assets/Gonzalez-Robles-Figures-3-01.png",
        style = "max-width:60%;height:auto;margin-bottom:20px;",
        alt = "QTL Analysis Context"
      ),
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      h4("Volcano Plots"),
      fluidRow(
        column(
          6,
          shinycssloaders::withSpinner(plotlyOutput(ns("volcano_mutation")),
            color = "#6A2798"
          )
        ),
        column(
          6,
          shinycssloaders::withSpinner(plotlyOutput(ns("volcano_ups")),
            color = "#6A2798"
          )
        )
      ),
      h4("Boxplots by Cohort"),
      div(
        shinycssloaders::withSpinner(plotlyOutput(ns("boxplot")),
          color = "#6A2798"
        )
      ),
      br(),
      fluidRow(
        column(6, downloadButton(outputId = ns("downloadDat"), label = "Download CSV")),
        column(6, downloadButton(outputId = ns("downloadPDF"), label = "Download PDF"))
      )
    )
  )
}
QTLServer <- function(id, genelists, colorsTbl, cohorts) {
  moduleServer(id, function(input, output, session) {
    # Lazy load data once per session
    qtl_cache <- reactiveVal(NULL)
    sample_mapping_cache <- reactiveVal(NULL)
    protein_cache <- reactiveVal(NULL)
    mut_cache <- reactiveVal(NULL)

    get_qtl_data <- function() {
      if (is.null(qtl_cache())) {
        qtl_cache(readRDS("./data/CPTAC_Pancan_pQTLs_with_class.rds"))
      }
      qtl_cache()
    }

    get_sample_mapping <- function() {
      if (is.null(sample_mapping_cache())) {
        sample_mapping_cache(
          readRDS("./data/UbiDash_CPTAC_PanCan_sample_ID_mapping_across_datatypes_v1.1.rds") %>%
            filter(datatype == "protein", include == TRUE) %>%
            distinct(harmonized_data_id, .keep_all = TRUE)
        )
      }
      sample_mapping_cache()
    }

    get_protein_data <- function() {
      if (is.null(protein_cache())) {
        protein_cache(
          readRDS("./data/UbiDash_CPTAC_PanCan_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed_clean.rds") %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column("harmonized_data_id")
        )
      }
      protein_cache()
    }

    get_mut_data <- function() {
      if (is.null(mut_cache())) {
        mut_cache(mem_get_from_drive("cptac_v1.1_mut.rds"))
      }
      mut_cache()
    }

    # Get data references
    qtl_data <- get_qtl_data()
    sample_mapping <- get_sample_mapping()
    protein_data <- get_protein_data()
    mut <- get_mut_data()

    # Precompute static choices
    ups_choices <- sort(unique(qtl_data$ups_protein))
    mut_choices <- colnames(mut)
    genelist_choices <- names(genelists)

    # Update input choices
    updateSelectizeInput(session, "gene",
      choices = ups_choices,
      server = TRUE, selected = "TRIM29"
    )
    updateSelectizeInput(session, "genelist",
      choices = genelist_choices,
      server = TRUE, selected = "hCD"
    )
    updateSelectizeInput(session, "mutgene",
      choices = mut_choices,
      server = TRUE, selected = "TP53.mut"
    )

    # Boxplot output
    output$boxplot <- renderPlotly({
      req(input$refresh)
      gn <- isolate(input$gene)
      m <- isolate(input$mutgene)
      if (!isTruthy(gn)) gn <- "ABL1"
      if (!isTruthy(m)) m <- "VHL"
      colors <- colorsTbl
      colors <- setNames(object = colors$color, nm = colors$id)
      selected_cohorts <- cohorts
      df <- protein_data %>%
        left_join(sample_mapping, by = "harmonized_data_id") %>%
        filter(cohort %in% selected_cohorts, !is.na(harmonized_data_id)) %>%
        select(harmonized_data_id, cohort, all_of(gn)) %>%
        rename(protein_expr = all_of(gn)) %>%
        left_join(
          data.frame(
            harmonized_data_id = rownames(mut),
            mut_status = mut[, m, drop = TRUE]
          ),
          by = "harmonized_data_id"
        ) %>%
        filter(!is.na(protein_expr), !is.na(mut_status))
      df$mut_status <- factor(df$mut_status, levels = c("WT", "Mutant"))
      n_req <- isolate(input$n_genes)
      n_mutant <- sum(df$mut_status == "Mutant", na.rm = TRUE)
      n_wt <- sum(df$mut_status == "WT", na.rm = TRUE)
      if (n_mutant < n_req | n_wt < n_req) {
        plot <- ggplot() +
          annotate("text",
            x = 0.5, y = 0.5,
            label = paste(
              "Not enough samples in one of the groups to meet your specified minimum:",
              "\n- WT samples:", n_wt,
              "\n- Mutant samples:", n_mutant,
              "\n- Your specified minimum per group:", n_req
            ),
            size = 6, hjust = 0.5
          ) +
          theme_void()
        return(ggplotly(plot))
      }
      p <- ggplot(df, aes(x = mut_status, y = protein_expr, color = cohort)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(position = position_jitter(0.2), alpha = 0.6, size = 1.2) +
        scale_color_manual(values = colors) +
        facet_wrap(~cohort, scales = "free_y", labeller = label_wrap_gen(width = 10), nrow = 2) +
        theme_minimal() +
        labs(
          x = m, y = gn,
          title = "CPTAC Pancan Protein Expression by Mutation Status",
          color = "Cohort"
        ) +
        theme(
          legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 9)
        )
      cohort_levels <- unique(df$cohort)
      for (coh in cohort_levels) {
        subdf <- df[df$cohort == coh, ]
        if (length(unique(subdf$mut_status)) == 2 &&
          sum(subdf$mut_status == "WT") > 0 &&
          sum(subdf$mut_status == "Mutant") > 0) {
          p <- p + stat_compare_means(
            data = subdf,
            aes(x = mut_status, y = protein_expr),
            method = "wilcox.test",
            # label = "p.format",
            label = "p.signif",
            label.x.npc = "center",
            label.y.npc = "top",
            na.rm = TRUE,
          )
        }
      }
      return(ggplotly(p, height = 800))
    })

    # Reactive data frames for volcano plots
    df_mutation <- eventReactive(input$refresh,
      {
        n <- isolate(input$n_genes)
        gn <- gsub(".mut", "", isolate(input$mutgene))
        if (!isTruthy(gn)) gn <- "VHL"
        df <- qtl_data %>%
          filter(n_mut > n &
            is.numeric(estimate) &
            !is.nan(estimate) &
            !is.infinite(estimate) &
            mutation == gn) %>%
          mutate(Label = paste(ups_protein, mutation, sep = "-")) %>%
          mutate(Class = as.factor(Class))
        return(df)
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )
    df_ups <- eventReactive(input$refresh,
      {
        n <- isolate(input$n_genes)
        gn <- isolate(input$gene)
        if (!isTruthy(gn)) gn <- "ABL1"
        df <- qtl_data %>%
          filter(ups_protein == gn & n_mut > n) %>%
          mutate(Label = paste(ups_protein, mutation, sep = "-")) %>%
          mutate(Class = as.factor(mutation))
        return(df)
      },
      ignoreNULL = FALSE,
      ignoreInit = FALSE
    )

    # Mutation-centric volcano plot
    output$volcano_mutation <- renderPlotly({
      req(input$refresh)

      n <- isolate(input$n_genes)
      p_adj_log <- isolate(input$pvalcutoff) # Already in -log10 scale
      p_adj <- 10^(-p_adj_log) # Convert back for display

      df_data <- df_mutation()

      if (nrow(df_data) == 0) {
        plot <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No data available for this mutation", size = 5) +
          theme_void()
        return(ggplotly(plot))
      }
      class_levels <- levels(df_data$Class)
      color_vec <- colorsTbl[colorsTbl$id %in% class_levels, ]$color
      if (length(color_vec) == 0 || length(class_levels) == 0 || length(color_vec) != length(class_levels)) {
        colors <- setNames(rep("#555555", length(class_levels)), class_levels)
      } else {
        colors <- setNames(color_vec, class_levels)
      }

      gn <- gsub(".mut", "", isolate(input$mutgene))
      if (!isTruthy(gn)) gn <- "VHL"

      p <- ggplot(df_data, aes(x = estimate, y = -log10(fdr), label = ups_protein, color = Class)) +
        geom_point(alpha = 0.5, size = 2.5) +
        theme_prism(border = TRUE) +
        scale_color_manual(values = colors, name = "UPS Protein Class") +
        labs(
          x = glue("Estimate ({gn}-mut/wt)"), y = "-log10(p-adj)",
          title = glue("Mutation Centric: {gn}"),
          caption = glue("n_mut > {n} | p-adj cutoff = {format(p_adj, scientific = TRUE)}")
        ) +
        geom_hline(yintercept = p_adj_log, linetype = "dashed", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black")

      return(ggplotly(p))
    })

    # UPS-centric volcano plot
    output$volcano_ups <- renderPlotly({
      req(input$refresh)

      n <- isolate(input$n_genes)
      p_adj_log <- isolate(input$pvalcutoff) # Already in -log10 scale
      p_adj <- 10^(-p_adj_log) # Convert back for display

      df_data <- df_ups()

      if (nrow(df_data) == 0) {
        plot <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No data available for this UPS protein", size = 5) +
          theme_void()
        return(ggplotly(plot))
      }

      selected_genelist <- genelists[[isolate(input$genelist)]]
      if (!isTruthy(selected_genelist)) {
        selected_genelist <- genelists[["hCD"]]
      }

      gn <- isolate(input$gene)
      if (!isTruthy(gn)) gn <- "ABL1"

      p <- ggplot(df_data, aes(x = estimate, y = -log10(fdr), label = mutation, color = mutation %in% selected_genelist)) +
        geom_point(alpha = 0.5, size = 2.5) +
        scale_color_manual(name = glue("{isolate(input$genelist)} Gene"), values = c("TRUE" = "#44BBFF", "FALSE" = "#55555555")) +
        theme_minimal() +
        labs(
          x = glue("Estimate (mut/wt) {gn}"), y = "-log10(p-adj)",
          title = glue("UPS Centric: {gn}"),
          caption = glue("n_mut > {n} | p-adj cutoff = {format(p_adj, scientific = TRUE)}")
        ) +
        geom_hline(yintercept = p_adj_log, linetype = "dashed", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black")

      return(ggplotly(p))
    })

    # Download handler - combine both datasets
    output$downloadDat <- downloadHandler(
      filename = function() {
        paste("QTL-analysis-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        mutation_data <- df_mutation() %>% mutate(analysis_type = "mutation_centric")
        ups_data <- df_ups() %>% mutate(analysis_type = "ups_centric")
        combined <- bind_rows(mutation_data, ups_data)
        write.csv(combined, file)
      }
    )

    # PDF download handler
    output$downloadPDF <- downloadHandler(
      filename = function() {
        paste("QTL-boxplot-", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        gn <- isolate(input$gene)
        m <- isolate(input$mutgene)
        if (!isTruthy(gn)) gn <- "ABL1"
        if (!isTruthy(m)) m <- "VHL"
        colors <- colorsTbl
        colors <- setNames(object = colors$color, nm = colors$id)
        selected_cohorts <- cohorts
        df <- protein_data %>%
          left_join(sample_mapping, by = "harmonized_data_id") %>%
          filter(cohort %in% selected_cohorts, !is.na(harmonized_data_id)) %>%
          select(harmonized_data_id, cohort, all_of(gn)) %>%
          rename(protein_expr = all_of(gn)) %>%
          left_join(
            data.frame(
              harmonized_data_id = rownames(mut),
              mut_status = mut[, m, drop = TRUE]
            ),
            by = "harmonized_data_id"
          ) %>%
          filter(!is.na(protein_expr), !is.na(mut_status))
        df$mut_status <- factor(df$mut_status, levels = c("WT", "Mutant"))

        p <- ggplot(df, aes(x = mut_status, y = protein_expr, color = cohort)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position = position_jitter(0.2), alpha = 0.6, size = 1.2) +
          scale_color_manual(values = colors) +
          facet_wrap(~cohort, scales = "free_y", labeller = label_wrap_gen(width = 10), nrow = 2) +
          theme_minimal() +
          labs(
            x = m, y = gn,
            title = "CPTAC Pancan Protein Expression by Mutation Status",
            subtitle = "Significance: ns (p>0.05), * (p≤0.05), ** (p≤0.01), *** (p≤0.001), **** (p≤0.0001)",
            color = "Cohort"
          ) +
          theme(
            legend.position = "bottom",
            strip.text = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9, color = "gray40")
          )

        cohort_levels <- unique(df$cohort)
        for (coh in cohort_levels) {
          subdf <- df[df$cohort == coh, ]
          if (length(unique(subdf$mut_status)) == 2 &&
            sum(subdf$mut_status == "WT") > 0 &&
            sum(subdf$mut_status == "Mutant") > 0) {
            p <- p + stat_compare_means(
              data = subdf,
              aes(x = mut_status, y = protein_expr),
              method = "wilcox.test",
              label = "p.signif",
              label.x.npc = "center",
              label.y.npc = "top",
              na.rm = TRUE
            )
          }
        }

        ggsave(file, plot = p, width = 14, height = 12, device = "pdf")
      }
    )
  })
}
