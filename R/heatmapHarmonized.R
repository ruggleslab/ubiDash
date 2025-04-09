# UI Module: heatmapHarmonizedUI
heatmapHarmonizedUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      
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
        label    = "Gene(s) of interest:",
        choices  = NULL,
        multiple = TRUE,
        options  = list(
          placeholder = "Type or select one or more gene symbols",
          server = TRUE  # Enable server-side search
        )
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
    mat_ups_file  <- "./data/mat_ups.rds"
    sampl_df_file <- "./data/sampl_df.rds"
    ups_df_file   <- "./data/ups_df.rds"
    colors_file   <- "./data/colors.rds"
    
    mat.ups  <- readRDS(mat_ups_file)
    sampl.df <- readRDS(sampl_df_file)
    ups.df   <- readRDS(ups_df_file)
    colors   <- readRDS(colors_file)
    
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
      
      # Filter columns (samples) based on selected dataset, type, and tissue.
      keep_cols <- sampl.df$dataset %in% input$filter_dataset &
        sampl.df$type    %in% input$filter_type &
        sampl.df$tissue  %in% input$filter_tissue
      
      if (!any(keep_cols)) {
        return(list(
          mat   = matrix(nrow = 0, ncol = 0),
          sampl = sampl.df[FALSE, ],
          ups   = ups.df
        ))
      }
      
      sampl_sub <- sampl.df[keep_cols, , drop = FALSE]
      mat_sub   <- mat.ups[, keep_cols, drop = FALSE]
      
      # Filter rows (genes) based on selected gene symbols
      if (length(input$filter_genes) > 0) {
        keep_rows <- rownames(mat_sub) %in% input$filter_genes
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
    output$static_heatmap <- renderPlot({
      # If exactly one gene is selected, show a friendly message
      if (!is.null(input$filter_genes) && length(input$filter_genes) == 1) {
        plot.new()
        text(x = 0.5, y = 0.5,
             "Please select at least 2 genes for a meaningful heatmap.",
             cex = 1.5, col = "#6A2798")
        return(NULL)
      }
      
      data_list <- filteredData()
      mat_sub   <- data_list$mat
      sampl_sub <- data_list$sampl
      ups_sub   <- data_list$ups
      
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
      ht <- ComplexHeatmap::Heatmap(
        mat_sub,
        name = "Tissue Specificity Score",
        na_col = "#C0C0C2",
        col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#5B99C2", "#1F316F")),
        cluster_rows    = if (nrow(mat_sub) > 1) TRUE else FALSE,
        show_row_names  = if (nrow(mat_sub) < 30) TRUE else FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(
          title = "Tissue Specificity Index",
          title_position = "topleft",
          legend_direction = "vertical"
        ),
        column_split    = sampl_sub$dataset,
        top_annotation  = annot,
        left_annotation = gn.annot,
        row_split       = ups_sub$UPS_Class,
        row_title       = NULL
      )
      
      ## ----- Draw the Heatmap -----
      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side     = "right",
        annotation_legend_side  = "right",
        padding = grid::unit(c(1, 1, 1, 1), "mm")
      )
    }, height = 700)
    
    ## ======================
    ## Reset Filters Button
    ## ======================
    observeEvent(input$reset_filters, {
      updateSelectizeInput(session, "filter_dataset",
                           selected = sort(unique(sampl.df$dataset)))
      updateSelectizeInput(session, "filter_type",
                           selected = sort(unique(sampl.df$type)))
      updateSelectizeInput(session, "filter_tissue",
                           selected = sort(unique(sampl.df$tissue)))
      updateSelectizeInput(session, "filter_genes",
                           selected = NULL)
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
        mat_sub   <- data_list$mat
        sampl_sub <- data_list$sampl
        ups_sub   <- data_list$ups
        
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
        
        ht <- ComplexHeatmap::Heatmap(
          mat_sub,
          name = "Tissue Specificity Score",
          na_col = "#C0C0C2",
          col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#5B99C2", "#1F316F")),
          cluster_rows    = if (nrow(mat_sub) > 1) TRUE else FALSE,
          show_row_names  = if (nrow(mat_sub) < 30) TRUE else FALSE,
          show_column_names = FALSE,
          heatmap_legend_param = list(
            title = "Tissue Specificity Index",
            title_position = "topcenter",
            legend_direction = "horizontal",
            legend_position = "bottom"
          ),
          column_split    = sampl_sub$dataset,
          top_annotation  = annot,
          left_annotation = gn.annot,
          row_split       = ups_sub$UPS_Class,
          row_title       = NULL
        )
        
        ComplexHeatmap::draw(ht)
        dev.off()
      }
    )
    
  })
}