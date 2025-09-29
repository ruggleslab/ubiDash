# UI for Ubiquityl Scatter Plot Module
LSCCUbiquitylUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Ubiquityl Scatter Plot"),
      p(tab_summaries$`LSCC Ubiquityl`),
      selectInput(ns("e3"), "Select E3 protein:", choices = NULL),
      selectInput(ns("substrate"), "Select Candidate Substrate:", choices = NULL),
      # Allow multiple selection for ubiquityl sites
      selectizeInput(ns("ub_site"), "Select Candidate Ubiquityl Site(s):", choices = NULL, multiple = TRUE),
      actionButton(ns("generate_plot"), "Generate Scatter Plot")
    ),
    mainPanel(
      conditionalPanel(
        condition = sprintf("!output['%s']", ns("no_sites_trigger")),
        shinycssloaders::withSpinner(plotOutput(ns("ubiScatterPlot")), color = "#6A2798")
      ),
      conditionalPanel(
        condition = sprintf("output['%s']", ns("no_sites_trigger")),
        uiOutput(ns("no_sites_msg"))
      )
    )
  )
}

# Server for Ubiquityl Scatter Plot Module
LSCCUbiquitylServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load data frames.
    # prot is loaded from final_data_2.csv (contains E3 proteins and candidate substrate values),
    # and ub.prot is loaded from final_data.csv (contains candidate ubiquityl site measurements in "CandidateSubstrate_site" format).
    prot <- read.csv("./data/LSCC_prot_former_Combine_PanCan_CPTAC_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed.csv", stringsAsFactors = FALSE)
    ub.prot <- read.csv("./data/LSCC_ub_prot_former_1-s2.0-S0092867421008576-mmc3.csv", stringsAsFactors = FALSE)

    ### Update E3 Input
    # Candidate E3 proteins: from prot, use columns (excluding sample_id) that do not have an underscore.
    observe({
      req(prot)
      candidate_e3 <- setdiff(colnames(prot), "sample_id")
      candidate_e3 <- candidate_e3[!grepl("_", candidate_e3)]
      current_e3 <- isolate(input$e3)
      selected_e3 <- if (!is.null(current_e3) && current_e3 %in% candidate_e3) current_e3 else candidate_e3[1]
      updateSelectInput(session, "e3", choices = candidate_e3, selected = selected_e3)
    })

    ### Update Candidate Substrate Input
    # Candidate substrates: from prot, columns without underscores and not equal to the selected E3.
    observe({
      req(prot)
      candidate_sub <- setdiff(colnames(prot), "sample_id")
      candidate_sub <- candidate_sub[!grepl("_", candidate_sub)]
      if (!is.null(input$e3)) {
        candidate_sub <- candidate_sub[!candidate_sub %in% input$e3]
      }
      current_sub <- isolate(input$substrate)
      selected_sub <- if (!is.null(current_sub) && current_sub %in% candidate_sub) current_sub else candidate_sub[1]
      updateSelectInput(session, "substrate", choices = candidate_sub, selected = selected_sub)
    })

    ### Update Candidate Ubiquityl Site Input
    # Candidate ubiquityl sites: from ub.prot, choose columns (excluding sample_id) that contain an underscore
    # and start with the candidate substrate.
    observe({
      req(ub.prot, input$substrate)
      candidate_ub <- setdiff(colnames(ub.prot), "sample_id")
      candidate_ub <- candidate_ub[grepl("_", candidate_ub)]
      candidate_ub <- candidate_ub[grep(paste0("^", input$substrate, "_"), candidate_ub)]
      current_site <- isolate(input$ub_site)
      # When multiple selections are allowed, preserve any that are still valid.
      if (!is.null(current_site) && length(current_site) > 0) {
        valid_sites <- current_site[current_site %in% candidate_ub]
        selected_site <- if (length(valid_sites) > 0) valid_sites else candidate_ub[1]
      } else {
        selected_site <- candidate_ub[1]
      }
      updateSelectizeInput(session, "ub_site", choices = candidate_ub, selected = selected_site)
    })
    # If no site available
    no_sites_trigger <- reactiveVal(FALSE)
    output$no_sites_trigger <- reactive({
      no_sites_trigger()
    })
    outputOptions(output, "no_sites_trigger", suspendWhenHidden = FALSE)
    output$no_sites_msg <- renderUI({
      if (no_sites_trigger()) {
        div(
          style = " display: flex;
        align-items: center;
        justify-content: center;
        height: 100px;
        width: 100%;
        color: #6A2798;
        font-weight: bold;
        font-size: 1.2em;",
          paste0("No ubiquityl sites available for ", input$substrate)
        )
      } else {
        NULL
      }
    })
    ### Create Joined Data
    # Join the two data frames on sample_id.
    joined_data <- reactive({
      req(prot, ub.prot)
      prot %>% inner_join(ub.prot, by = "sample_id")
    })

    ### Generate Scatter Plot
    # For each selected ubiquityl site (target column), generate a scatter plot:
    # - X-axis: the selected E3 protein measurement (from prot)
    # - Y-axis: the candidate ubiquityl site measurement (target column from ub.prot)
    scatter_plot <- eventReactive(input$generate_plot, {
      if (is.null(input$ub_site) || length(input$ub_site) == 0) {
        no_sites_trigger(TRUE)
        return(NULL)
      }
      no_sites_trigger(FALSE)
      # Build a list of plots for each selected candidate ubiquityl site.
      plots <- lapply(input$ub_site, function(site) {
        target_col <- site # since each selected site already is in "CandidateSubstrate_site" format
        df <- joined_data() %>%
          select(sample_id, matches(paste0("^(", input$e3, "|", target_col, ")")))

        if (!(target_col %in% colnames(df))) {
          showNotification(glue("Candidate ubiquityl site {target_col} not found in joined data"), type = "error")
          return(NULL)
        }
        if (!(input$e3 %in% colnames(df))) {
          showNotification(glue("E3 protein {input$e3} not found in joined data"), type = "error")
          return(NULL)
        }

        # Generate scatter plot with X-axis = E3 measurement and Y-axis = candidate ubiquityl site.
        cor_value <- tryCatch(
          {
            cor.test(df[[input$e3]], df[[target_col]], method = "spearman")$estimate
          },
          error = function(e) NA
        )
        ggplot(df, aes(x = !!rlang::sym(input$e3), y = !!rlang::sym(target_col))) +
          geom_point(color = "slategray3", size = 3, alpha = 0.7) +
          geom_smooth(
            method = "lm", se = TRUE,
            color = ifelse(is.na(cor_value) || cor_value < 0, "blue", "red")
          ) +
          stat_cor(
            method = "spearman",
            cor.coef.name = "rho",
            aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
            label.y = max(df[[target_col]], na.rm = TRUE) * 0.9
          ) +
          ggtitle(glue("{input$e3} vs {target_col}")) +
          theme_prism() +
          theme(legend.position = "none")
      })

      # Remove any NULL plots
      plots <- plots[!sapply(plots, is.null)]
      if (length(plots) == 0) {
        return(NULL)
      }
      wrap_plots(plots, ncol = 2)
    })

    output$ubiScatterPlot <- renderPlot({
      scatter_plot()
    })
  })
}
