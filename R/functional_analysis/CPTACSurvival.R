# UI Module: CPTACSurvivalUI
CPTACSurvivalUI <- function(id) {
    ns <- NS(id)

    sidebarLayout(
        sidebarPanel(
            h3("CPTAC Survival Explorer"),
            p("Visualize CPTAC pan-cancer survival curves based on a protein of interest."),
            # Illustration image
            tags$div(
                style = "text-align: center; margin-bottom: 15px;",
                tags$img(
                    src = "assets/Gonzalez-Robles-Figures-4-01.png",
                    alt = "Survival Analysis Illustration",
                    style = "max-width: 100%; height: auto; border-radius: 5px;"
                )
            ),
            selectizeInput(
                ns("gene_symbol"),
                label = "Protein / gene symbol:",
                choices = NULL,
                multiple = FALSE,
                options = list(
                    placeholder = "Type or select one protein symbol",
                    maxOptions = 20
                )
            ),
            radioButtons(
                ns("analysis_type"),
                label = "Analysis type:",
                choices = c("Pan-Cancer" = "pancan", "Per-Cancer" = "percan"),
                selected = "pancan",
                inline = TRUE
            ),
            conditionalPanel(
                condition = sprintf("input['%s'] == 'percan'", ns("analysis_type")),
                selectizeInput(
                    ns("cohort_selection"),
                    label = "Select cohort(s):",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE,
                    options = list(
                        placeholder = "Select one or more CPTAC cohorts",
                        maxItems = 5
                    )
                )
            ),
            sliderInput(
                ns("quantile_cutoff"),
                label = "Quantile threshold for Low / High groups (%):",
                min = 10,
                max = 45,
                step = 5,
                value = 25
            ),
            tags$div(
                class = "help-topic alert alert-secondary",
                style = "padding: 0.65rem 0.85rem; font-size: 0.85rem;",
                tags$strong("Help:"),
                "The threshold splits samples based on protein expression: samples in the bottom X% are 'Low', top X% are 'High', and the middle samples are excluded. For example, at 25%, the bottom 25% percentile = Low group, top 75% percentile = High group. Quantiles are calculated per cohort independently.",
                tags$br()
            ),
            tags$small(
                "Low and High groups are defined by the cohort-specific quantile cutoffs. Middle samples are excluded to emphasize expression extremes.",
                style = "color: #6c757d;"
            ),
            tags$small(
                tags$strong("Per-Cancer mode:"),
                " Analyze survival within specific cohorts to identify lineage-specific prognostic markers.",
                style = "color: #6c757d; display: block; margin-top: 10px;"
            )
        ),
        mainPanel(
            uiOutput(ns("cohort_summary")),
            plotOutput(ns("survival_plot"), height = "800px"),
            br(),
            downloadButton(ns("downloadData"), "Download Data (CSV)"),
            downloadButton(ns("downloadPlot"), "Download Plot (PDF)"),
            br(),
            br(),
            verbatimTextOutput(ns("pvalue_text")),
            h4("Grouped sample counts"),
            tableOutput(ns("level_counts"))
        )
    )
}

# Server Module: CPTACSurvivalServer
CPTACSurvivalServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        proteome_file <- "./data/UbiDash_CPTAC_PanCan_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed_clean.rds"
        sample_file <- "./data/UbiDash_CPTAC_PanCan_sample_ID_mapping_across_datatypes_v1.1.rds"
        overall_col <- "Overall.survival,.days"
        status_col <- "Survival.status.(1,.dead;.0,.alive)"

        protein_scaled <- readRDS(proteome_file) %>%
            sweep(., 1, rowMeans(., na.rm = TRUE), "-") %>%
            scale(
                .,
                center = apply(., 2, median, na.rm = TRUE),
                scale = apply(., 2, mad, na.rm = TRUE)
            ) %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "harmonized_data_id")

        sample_meta <- readRDS(sample_file) %>%
            filter(datatype == "protein" & include == "TRUE") %>%
            filter(harmonized_data_id %in% protein_scaled$harmonized_data_id) %>%
            select(
                harmonized_data_id,
                cohort,
                tissue,
                type,
                cancer_subtype,
                age,
                sex,
                histologic_grade,
                ESTIMATE_TumorPurity,
                all_of(overall_col),
                all_of(status_col)
            ) %>%
            distinct() %>%
            mutate(cohort = as.character(cohort))

        gene_choices <- sort(setdiff(colnames(protein_scaled), "harmonized_data_id"))
        cohort_values <- sample_meta$cohort
        cohort_values <- cohort_values[!is.na(cohort_values) & cohort_values != ""]
        cohort_choices <- unique(cohort_values)
        cohort_choices <- sort(cohort_choices)

        updateSelectizeInput(
            session,
            "gene_symbol",
            choices = gene_choices,
            selected = if (length(gene_choices) > 0) gene_choices[1] else NULL
        )

        updateSelectizeInput(
            session,
            "cohort_selection",
            choices = cohort_choices,
            selected = NULL
        )

        selected_cohorts <- reactive({
            if (input$analysis_type == "pancan") {
                return(NULL) # NULL means all cohorts (Pan-Cancer)
            } else {
                choices <- input$cohort_selection
                if (is.null(choices) || length(choices) == 0) {
                    return(NULL)
                }
                return(choices)
            }
        })

        survival_df <- reactive({
            req(input$gene_symbol)
            gene <- input$gene_symbol
            quantile_prob <- input$quantile_cutoff / 100

            df <- sample_meta %>%
                filter(type == "Tumor") %>%
                filter(!is.na(cohort) & cohort != "")

            cohort_filter <- selected_cohorts()
            if (!is.null(cohort_filter)) {
                df <- df %>%
                    filter(cohort %in% cohort_filter)
            }

            expr <- protein_scaled %>%
                select(harmonized_data_id, all_of(gene))

            df <- df %>%
                left_join(expr, by = "harmonized_data_id") %>%
                mutate(expression = .data[[gene]]) %>%
                filter(!is.na(expression)) %>%
                mutate(
                    survival_days = as.numeric(.data[[overall_col]]),
                    censoring_status = as.numeric(.data[[status_col]])
                ) %>%
                filter(!is.na(survival_days) & !is.na(censoring_status))
            df <- df %>%
                group_by(cohort) %>%
                mutate(
                    q_low = stats::quantile(expression, probs = quantile_prob, na.rm = TRUE),
                    q_high = stats::quantile(expression, probs = 1 - quantile_prob, na.rm = TRUE),
                    level = case_when(
                        expression <= q_low ~ "Low",
                        expression >= q_high ~ "High",
                        TRUE ~ NA_character_
                    )
                ) %>%
                ungroup() %>%
                mutate(level = factor(level, levels = c("Low", "High"))) %>%
                filter(!is.na(level))

            if (nrow(df) == 0) {
                return(NULL)
            }

            df
        })

        km_fit <- reactive({
            df <- survival_df()
            req(df)
            survfit(Surv(survival_days, censoring_status) ~ level, data = df)
        })

        cohort_label <- reactive({
            df <- survival_df()
            if (is.null(df)) {
                return("Pan-Cancer")
            }
            if (is.null(selected_cohorts())) {
                return("Pan-Cancer")
            }
            paste(unique(df$cohort), collapse = ", ")
        })

        output$cohort_summary <- renderUI({
            df <- survival_df()
            req(df)
            n_low <- sum(df$level == "Low", na.rm = TRUE)
            n_high <- sum(df$level == "High", na.rm = TRUE)
            total <- nrow(df)

            tagList(
                p(glue::glue("Gene: {input$gene_symbol}")),
                p(glue::glue("Cohorts: {cohort_label()}")),
                p(glue::glue(
                    "Samples: {total} (Low: {n_low}, High: {n_high})."
                )),
                p(glue::glue("Quantiles: {input$quantile_cutoff}% (Low) / {100 - input$quantile_cutoff}% (High)."))
            )
        })

        output$survival_plot <- renderPlot(
            {
                df <- survival_df()
                shiny::validate(
                    shiny::need(!is.null(df) && nrow(df) > 0, "No eligible samples after filtering. Try another gene or cohort."),
                    shiny::need(dplyr::n_distinct(df$level) == 2, "Select a gene with enough samples in both Low and High groups.")
                )

                fit <- km_fit()
                cohort_title <- cohort_label()

                surv_plot <- survminer::ggsurvplot(
                    fit,
                    data = df,
                    size = 0.8,
                    palette = c("#4D2D8C", "#E4004B"),
                    conf.int = FALSE,
                    pval = TRUE,
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    risk.table.height = 0.2,
                    risk.table.fontsize = 4.2,
                    tables.theme = ggplot2::theme(
                        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 12), size = 10)
                    ),
                    ggtheme = ggpubr::theme_pubr(),
                    title = glue::glue("{input$gene_symbol} ({cohort_title})")
                )

                surv_plot$plot <- surv_plot$plot +
                    ggplot2::theme(
                        plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
                    )

                print(surv_plot)
            },
            height = 800,
            res = 96
        )

        output$pvalue_text <- renderText({
            df <- survival_df()
            req(df)
            fit <- km_fit()
            pval <- tryCatch(
                survminer::surv_pvalue(fit, data = df)$pval,
                error = function(e) NA_real_
            )
            if (is.na(pval)) {
                "Log-rank p-value: not available"
            } else {
                glue::glue("Log-rank p-value: {formatC(pval, digits = 3, format = 'g')}")
            }
        })

        output$level_counts <- renderTable(
            {
                df <- survival_df()
                req(df)
                df %>%
                    count(cohort, level) %>%
                    tidyr::pivot_wider(
                        names_from = level,
                        values_from = n,
                        values_fill = 0
                    ) %>%
                    dplyr::arrange(cohort)
            },
            rownames = FALSE
        )

        # Download handler for the plot (PDF format)
        output$downloadPlot <- downloadHandler(
            filename = function() {
                paste0(input$gene_symbol, "_survival_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file, width = 10, height = 8)
                df <- survival_df()
                req(df)
                fit <- km_fit()
                cohort_title <- cohort_label()
                surv_plot <- survminer::ggsurvplot(
                    fit,
                    data = df,
                    size = 0.8,
                    palette = c("#4D2D8C", "#E4004B"),
                    conf.int = FALSE,
                    pval = TRUE,
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    risk.table.height = 0.2,
                    risk.table.fontsize = 4.2,
                    tables.theme = ggplot2::theme(
                        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 12), size = 10)
                    ),
                    ggtheme = ggpubr::theme_pubr(),
                    title = glue::glue("{input$gene_symbol} ({cohort_title})")
                )
                print(surv_plot)
                dev.off()
            }
        )

        # Download handler for the data (CSV format)
        output$downloadData <- downloadHandler(
            filename = function() {
                paste0(input$gene_symbol, "_survival_data_", Sys.Date(), ".csv")
            },
            content = function(file) {
                df <- survival_df()
                req(df)
                write.csv(df, file, row.names = FALSE)
            }
        )
    })
}
