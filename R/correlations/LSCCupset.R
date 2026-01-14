LSCCUpsetUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Upset Plot"),
      p(tab_summaries$`LSCC Upset`),
      checkboxInput(ns("filter_all_true"),
        "Show only intersections where all sets are TRUE in table",
        value = FALSE
      )
    ),
    mainPanel(
      shinycssloaders::withSpinner(plotOutput(ns("distPlot")), color = "#6A2798"),
      br(),
      DT::dataTableOutput(ns("upsetTable")),
      br(),
      p("Download Upset data and plot:"),
      downloadButton(outputId = ns("downloadUpsetData"), label = "Download Upset Data"),
      downloadButton(outputId = ns("downloadPDF"), label = "Download Upset PDF")
    )
  )
}

LSCCUpsetServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    sets <- c("neg_prot_prot", "pos_prot_ub", "yes_biogrid")

    full_data <- reactive({
      req(file.exists("./data/EXP200_Percan_LSCC_upset_plot.xlsx"))
      df <- readxl::read_excel("./data/EXP200_Percan_LSCC_upset_plot.xlsx", skip = 1, col_names = TRUE)
      if ("%" %in% names(df)) {
        df <- df[, -which(names(df) == "%")]
      }
      for (set in sets) {
        df[[set]] <- as.logical(df[[set]])
        df[[set]][is.na(df[[set]])] <- FALSE
      }
      df
    })

    table_data <- reactive({
      df <- full_data()
      if (isTRUE(input$filter_all_true)) {
        df <- df[rowSums(df[, sets]) == length(sets), ]
      }
      df
    })

    output$upsetTable <- DT::renderDataTable({
      DT::datatable(table_data(), selection = "single", options = list(pageLength = 10))
    })

    selected_sets <- reactive({
      req(input$upsetTable_rows_selected)
      df <- table_data()
      selected_row <- df[input$upsetTable_rows_selected, ]
      sets[as.logical(selected_row[, sets])]
    })

    upset_plot <- reactive({
      df <- full_data()
      queries_arg <- NULL
      if (!is.null(input$upsetTable_rows_selected) && length(input$upsetTable_rows_selected) > 0) {
        true_sets <- selected_sets()
        if (length(true_sets) > 0) {
          queries_arg <- list(
            ComplexUpset::upset_query(
              intersect = true_sets,
              color = "red",
              fill = scales::alpha("lightcoral", 0.5),
              only_components = c("intersections_matrix", "Intersection size")
            )
          )
        }
      }
      ComplexUpset::upset(
        df,
        intersect = sets,
        name = "Gene Pair Intersections",
        sort_intersections_by = c("degree"),
        queries = queries_arg,
        base_annotations = list("Intersection size" = intersection_size(counts = TRUE)),
        set_sizes = upset_set_size(position = "right"),
        width_ratio = 0.15
      )
    })

    output$distPlot <- renderPlot({
      upset_plot()
    })

    output$downloadUpsetData <- downloadHandler(
      filename = function() {
        paste("upset_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(full_data(), file, row.names = FALSE)
      }
    )

    output$downloadPDF <- downloadHandler(
      filename = function() {
        paste("upset_plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = upset_plot(), device = "pdf", width = 7, height = 5)
      }
    )
  })
}
