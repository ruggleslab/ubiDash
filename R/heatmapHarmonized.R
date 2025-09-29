# UI Module: heatmapHarmonizedUI
heatmapHarmonizedUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      p(tab_summaries$Heatmap),

      # Filter inputs (method filter removed)
      selectizeInput(
        ns("filter_dataset"),
        label    = "Dataset:",
        choices  = NULL,
        multiple = TRUE
      ),
      selectizeInput(
        ns("filter_type"),
        label    = "Type:",
        choices  = NULL,
        multiple = TRUE
      ),
      selectizeInput(
        ns("filter_tissue"),
        label    = "Tissue:",
        choices  = NULL,
        multiple = TRUE
      ),

      # Gene search input with server-side search
      selectizeInput(
        ns("filter_genes"),
        label = "Gene(s) of interest:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Type or select one or more gene symbols",
          server = TRUE # Enable server-side search
        )
      ),

      # Clustering and splitting options
      checkboxInput(ns("cluster_columns"), "Cluster columns", value = TRUE),
      selectizeInput(
        ns("split_columns_by"),
        label = "Split columns by:",
        choices = c("None", "Dataset", "Type", "Tissue"),
        selected = "Dataset"
      ),

      # Buttons
      actionButton(ns("reset_filters"), "Reset Filters", class = "btn-primary")
    ),
    mainPanel(
      h3("UPS Protein Expression Profile"),
      plotOutput(ns("static_heatmap")),
      p("Download this data or plot:"),
      downloadButton(outputId = ns("download_data"), label = "Download Data"),
      downloadButton(outputId = ns("download_plot"), label = "Download Plot")
    )
  )
}
# Server Module: heatmapHarmonizedServer
heatmapHarmonizedServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ============
    ## Data Loading
    ## ============
    mat_ups_file <- "./data/mat_ups.rds"
    sampl_df_file <- "./data/sampl_df.rds"
    ups_df_file <- "./data/ups_df.rds"
    colors_file <- "./data/colors.rds"

    mat.ups <- readRDS(mat_ups_file)
    sampl.df <- readRDS(sampl_df_file)
    ups.df <- readRDS(ups_df_file)
    colors <- readRDS(colors_file)

    ## =============================
    ## Initialize Filter Input Choices
    ## =============================
    updateSelectizeInput(session, "filter_dataset",
      choices  = sort(unique(sampl.df$dataset)),
      selected = sort(unique(sampl.df$dataset))
    )
    updateSelectizeInput(session, "filter_type",
      choices  = sort(unique(sampl.df$type)),
      selected = sort(unique(sampl.df$type))
    )
    updateSelectizeInput(session, "filter_tissue",
      choices  = sort(unique(sampl.df$tissue)),
      selected = sort(unique(sampl.df$tissue))
    )
    updateSelectizeInput(session, "filter_genes",
      choices  = sort(rownames(mat.ups)),
      selected = NULL
    )

    ## ===============================
    ## Debounced Filter State
    ## ===============================
    filter_state <- reactive({
      list(
        dataset = input$filter_dataset,
        type = input$filter_type,
        tissue = input$filter_tissue,
        genes = input$filter_genes,
        cluster_columns = input$cluster_columns,
        split_columns_by = input$split_columns_by
      )
    })

    debounced_state <- shiny::debounce(filter_state, millis = 300)

    ## Update type when dataset changes (method filter removed)
    observeEvent(input$filter_dataset, {
      req(input$filter_dataset)
      filtered_sampl <- sampl.df[sampl.df$dataset %in% input$filter_dataset, ]
      updateSelectizeInput(session, "filter_type",
        choices  = sort(unique(filtered_sampl$type)),
        selected = sort(unique(filtered_sampl$type))
      )
    })

    ## ==========================
    ## Reactive: Filtered Data
    ## ==========================
    filteredData <- reactive({
      req(sampl.df, mat.ups, ups.df)

      state <- debounced_state()
      shiny::req(state)

      selected_dataset <- state$dataset
      if (is.null(selected_dataset)) selected_dataset <- character(0)
      selected_type <- state$type
      if (is.null(selected_type)) selected_type <- character(0)
      selected_tissue <- state$tissue
      if (is.null(selected_tissue)) selected_tissue <- character(0)

      # Filter columns (samples) based on selected dataset, type, and tissue.
      keep_cols <- sampl.df$dataset %in% selected_dataset &
        sampl.df$type %in% selected_type &
        sampl.df$tissue %in% selected_tissue

      if (!any(keep_cols)) {
        return(list(
          mat   = matrix(nrow = 0, ncol = 0),
          sampl = sampl.df[FALSE, ],
          ups   = ups.df
        ))
      }

      sampl_sub <- sampl.df[keep_cols, , drop = FALSE]
      mat_sub <- mat.ups[, keep_cols, drop = FALSE]

      # Filter rows (genes) based on selected gene symbols
      selected_genes <- state$genes
      if (!is.null(selected_genes) && length(selected_genes) > 0) {
        keep_rows <- rownames(mat_sub) %in% selected_genes
      } else {
        keep_rows <- rep(TRUE, nrow(mat_sub))
      }
      mat_sub <- mat_sub[keep_rows, , drop = FALSE]
      ups_sub <- ups.df[keep_rows, , drop = FALSE]

      ## --- Data Cleaning ---
      # Instead of removing rows with Inf, NA, or constant values,
      # replace any Inf or NA values with 0 so that all genes are retained.
      mat_sub[is.infinite(mat_sub)] <- 0
      mat_sub[is.na(mat_sub)] <- 0

      list(
        mat   = mat_sub,
        sampl = sampl_sub,
        ups   = ups_sub
      )
    })

    ## =====================
    ## Render Heatmap Plot
    ## =====================
    output$static_heatmap <- renderPlot(
      {
        state <- debounced_state()
        shiny::req(state)

        selected_genes <- state$genes
        # If exactly one gene is selected, show a friendly message
        if (!is.null(selected_genes) && length(selected_genes) == 1) {
          plot.new()
          text(
            x = 0.5, y = 0.5,
            "Please select at least 2 genes for a meaningful heatmap.",
            cex = 1.5, col = "#6A2798"
          )
          return(NULL)
        }

        data_list <- filteredData()
        mat_sub <- data_list$mat
        sampl_sub <- data_list$sampl
        ups_sub <- data_list$ups

        # Check for empty data
        if (nrow(mat_sub) == 0 || ncol(mat_sub) == 0) {
          plot.new()
          text(x = 0.5, y = 0.5, "No data to display. Adjust filters.", cex = 1.2, col = "#6A2798")
          return(NULL)
        }

        ## ----- Annotations for Heatmap -----
        # Subset sampl_sub to remove the "method" column from the annotations/legend
        annot_df <- sampl_sub[, c("dataset", "type", "tissue"), drop = FALSE]

        # Column annotation (using subset of sample data)
        annot <- ComplexHeatmap::columnAnnotation(
          df = annot_df,
          na_col = "grey",
          col = list(
            dataset = colors[levels(as.factor(annot_df$dataset))],
            type    = colors[levels(as.factor(annot_df$type))],
            tissue  = colors[levels(as.factor(annot_df$tissue))]
          ),
          show_annotation_name = FALSE,
          annotation_legend_param = list(direction = "horizontal", ncol = 1)
        )

        # Row annotation (using UPS classification)
        gn.annot <- ComplexHeatmap::rowAnnotation(
          df = ups_sub,
          na_col = "grey",
          col = list(
            UPS_Class = colors[levels(as.factor(ups_sub$UPS_Class))]
          ),
          show_annotation_name = FALSE,
          annotation_legend_param = list(direction = "horizontal")
        )

        ## ----- Create the Heatmap -----
        split_by <- state$split_columns_by
        if (is.null(split_by)) {
          split_by <- "None"
        }
        split_col_name <- tolower(split_by)
        column_split_var <- if (!identical(split_by, "None") && split_col_name %in% colnames(sampl_sub)) {
          sampl_sub[[split_col_name]]
        } else {
          NULL
        }
        heatmap_args <- list(
          mat_sub,
          name = "Relative Ranked Score",
          na_col = "#C0C0C2",
          col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#5B99C2", "#1F316F")),
          cluster_rows = if (nrow(mat_sub) > 1) TRUE else FALSE,
          cluster_columns = isTRUE(state$cluster_columns),
          show_row_names = if (nrow(mat_sub) < 30) TRUE else FALSE,
          show_column_names = FALSE,
          heatmap_legend_param = list(
            title = "Relative Ranked Score",
            title_position = "topleft",
            legend_direction = "vertical"
          ),
          column_split = column_split_var,
          top_annotation = annot,
          left_annotation = gn.annot,
          row_split = ups_sub$UPS_Class,
          row_title = NULL
        )
        if (!is.null(column_split_var) && identical(split_by, "Tissue")) {
          heatmap_args$column_title <- rep("", length(unique(column_split_var)))
        }
        ht <- do.call(ComplexHeatmap::Heatmap, heatmap_args)

        ## ----- Draw the Heatmap -----
        ComplexHeatmap::draw(
          ht,
          heatmap_legend_side = "right",
          annotation_legend_side = "right",
          padding = grid::unit(c(1, 1, 1, 1), "mm")
        )
      },
      height = 700
    )

    ## ======================
    ## Reset Filters Button
    ## ======================
    observeEvent(input$reset_filters, {
      updateSelectizeInput(session, "filter_dataset",
        selected = sort(unique(sampl.df$dataset))
      )
      updateSelectizeInput(session, "filter_type",
        selected = sort(unique(sampl.df$type))
      )
      updateSelectizeInput(session, "filter_tissue",
        selected = sort(unique(sampl.df$tissue))
      )
      updateSelectizeInput(session, "filter_genes",
        selected = NULL
      )
      updateCheckboxInput(session, "cluster_columns", value = TRUE)
      updateSelectizeInput(session, "split_columns_by", selected = "Dataset")
    })

    ## ========================
    ## Download Data and Plot
    ## ========================
    output$download_data <- downloadHandler(
      filename = function() {
        paste("summarized-profile-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        filtered <- filteredData()$mat
        write.csv(as.data.frame(filtered), file)
      }
    )
    # Make sure svglite is installed: install.packages("svglite")
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("heatmap-", Sys.Date(), ".svg", sep = "")
      },
      content = function(file) {
        svglite::svglite(file, width = 12, height = 8)

        data_list <- filteredData()
        mat_sub <- data_list$mat
        sampl_sub <- data_list$sampl
        ups_sub <- data_list$ups

        state <- debounced_state()
        shiny::req(state)

        annot_df <- sampl_sub[, c("dataset", "type", "tissue"), drop = FALSE]

        annot <- ComplexHeatmap::columnAnnotation(
          df = annot_df,
          na_col = "grey",
          col = list(
            dataset = colors[levels(as.factor(annot_df$dataset))],
            type    = colors[levels(as.factor(annot_df$type))],
            tissue  = colors[levels(as.factor(annot_df$tissue))]
          ),
          show_annotation_name = FALSE
        )

        gn.annot <- ComplexHeatmap::rowAnnotation(
          df = ups_sub,
          na_col = "grey",
          col = list(
            UPS_Class = colors[levels(as.factor(ups_sub$UPS_Class))]
          ),
          show_annotation_name = FALSE
        )

        split_by <- state$split_columns_by
        if (is.null(split_by)) {
          split_by <- "None"
        }
        split_col_name <- tolower(split_by)
        column_split_var <- if (!identical(split_by, "None") && split_col_name %in% colnames(sampl_sub)) {
          sampl_sub[[split_col_name]]
        } else {
          NULL
        }
        heatmap_args <- list(
          mat_sub,
          name = "Relative Ranked Score",
          na_col = "#C0C0C2",
          col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#5B99C2", "#1F316F")),
          cluster_rows = if (nrow(mat_sub) > 1) TRUE else FALSE,
          cluster_columns = isTRUE(state$cluster_columns),
          show_row_names = if (nrow(mat_sub) < 30) TRUE else FALSE,
          show_column_names = FALSE,
          heatmap_legend_param = list(
            title = "Relative Ranked Score",
            title_position = "topcenter",
            legend_direction = "horizontal",
            legend_position = "bottom"
          ),
          column_split = column_split_var,
          top_annotation = annot,
          left_annotation = gn.annot,
          row_split = ups_sub$UPS_Class,
          row_title = NULL
        )
        if (!is.null(column_split_var) && identical(split_by, "Tissue")) {
          heatmap_args$column_title <- rep("", length(unique(column_split_var)))
        }
        ht <- do.call(ComplexHeatmap::Heatmap, heatmap_args)

        ComplexHeatmap::draw(ht)
        dev.off()
      }
    )
  })
}
