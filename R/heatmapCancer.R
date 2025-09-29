heatmapCancerUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Customize UPS Protein Heatmap"),
      p(tab_summaries$`Heatmap Cancer Subtypes`),
      selectInput(ns("classSelect"), "UPS Class",
        choices = NULL, # will be updated dynamically
        multiple = TRUE
      ),
      selectInput(ns("cohortSelect"), "Cohort",
        choices = NULL, # will be updated dynamically
        multiple = TRUE
      ),
      selectInput(ns("cancerSubtypeSelect"), "Cancer Subtype",
        choices = NULL, # will be updated dynamically
        multiple = TRUE
      ),
      # NEW: Protein search input with server‐side search.
      selectizeInput(
        ns("filter_proteins"),
        label = "Protein(s) of interest:",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Type or select one or more protein symbols",
          server = TRUE
        )
      ),
      radioButtons(ns("displayMode"), "Display Mode",
        choices = c("Show as annotations" = "all", "Show only selected proteins" = "selected"),
        selected = "all"
      )
    ),
    mainPanel(
      shinycssloaders::withSpinner(
        plotOutput(ns("heatmapPlot"), height = "800px"),
        color = "#6A2798"
      ),
      br(),
      downloadButton(ns("download"), "Download PDF"),
      downloadButton(ns("downloadData"), "Download Data")
    )
  )
}
heatmapCancerServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load raw data once
    data_raw <- reactive({
      prot <- read.csv("./data/DEqMS.results_cancer_subtypes_with_tumotPurity.csv") # adjust file path
      pn <- openxlsx::read.xlsx(
        "./data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
        sheet = "dense"
      )[, -1] %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System")
      list(prot = prot, pn = pn)
    })

    # Base merged matrix for available proteins
    base_merged_mat <- reactive({
      dat <- data_raw()
      prot <- dat$prot
      pn <- dat$pn

      # Filter pn by selected UPS class
      if (!is.null(input$classSelect) && length(input$classSelect) > 0) {
        pn <- pn %>% dplyr::filter(Class %in% input$classSelect)
      }

      # Filter prot by selected cohorts
      if (!is.null(input$cohortSelect) && length(input$cohortSelect) > 0) {
        prot <- prot %>% dplyr::filter(cohort %in% input$cohortSelect)
      }

      # Merge network and differential expression data
      proteins <- pn %>%
        dplyr::group_by(Gene.Symbol) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(prot, by = c("Gene.Symbol" = "gene"))

      # Build matrix for each cohort
      mat <- list()
      for (i in unique(prot$cohort)) {
        tmp <- proteins %>%
          dplyr::filter(cohort == i) %>%
          dplyr::select(Gene.Symbol, logFC, contrast) %>%
          dplyr::group_by(Gene.Symbol) %>%
          tidyr::spread(key = contrast, value = logFC)
        colnames(tmp)[-1] <- paste0(i, "_", colnames(tmp)[-1])
        mat[[i]] <- tmp
      }

      merged_mat <- purrr::reduce(mat, dplyr::full_join, by = "Gene.Symbol") %>%
        tibble::column_to_rownames(var = "Gene.Symbol") %>%
        as.data.frame() %>%
        dplyr::filter(rowMeans(is.na(.)) < 0.2) %>%
        as.matrix()

      list(merged_mat = merged_mat, proteins = proteins)
    })

    # Update select input choices based on loaded data (only once at startup)
    observe({
      dat <- data_raw()

      updateSelectInput(session, "classSelect",
        choices = unique(dat$pn$Class),
        selected = unique(dat$pn$Class)
      )

      updateSelectInput(session, "cohortSelect",
        choices = unique(dat$prot$cohort),
        selected = unique(dat$prot$cohort)
      )

      cancer_subtypes <- unique(sapply(strsplit(dat$prot$contrast, "_vs_"), function(x) x[1]))
      updateSelectInput(session, "cancerSubtypeSelect",
        choices = cancer_subtypes,
        selected = cancer_subtypes
      )

      # Update the protein search input (using all available proteins from network data)
      updateSelectizeInput(session, "filter_proteins",
        choices = sort(unique(dat$pn$Gene.Symbol)),
        selected = NULL
      )
    })

    # Reactive for available proteins based on current filters
    available_proteins <- reactive({
      base_mat <- base_merged_mat()
      rownames(base_mat$merged_mat)
    })

    # Update protein choices when class/cohort filters change (but not when proteins are selected)
    observeEvent(c(input$classSelect, input$cohortSelect),
      {
        current_proteins <- available_proteins()
        current_selection <- input$filter_proteins

        # Only update if the available proteins have actually changed
        if (!identical(sort(current_proteins), sort(attr(current_selection, "choices")))) {
          # Preserve valid selections
          valid_selection <- current_selection[current_selection %in% current_proteins]

          updateSelectizeInput(session, "filter_proteins",
            choices = sort(current_proteins),
            selected = valid_selection
          )
        }
      },
      ignoreNULL = FALSE
    )

    # Build the heatmap object
    heatmap_obj <- reactive({
      base <- base_merged_mat()
      merged_mat <- base$merged_mat
      proteins <- base$proteins

      # Filter based on protein selection and display mode
      if (!is.null(input$filter_proteins) && length(input$filter_proteins) > 0) {
        if (input$displayMode == "selected") {
          # Show only selected proteins that are present in the data
          valid_selected <- input$filter_proteins[input$filter_proteins %in% rownames(merged_mat)]
          if (length(valid_selected) > 0) {
            merged_mat <- merged_mat[valid_selected, , drop = FALSE]
            display_mat <- merged_mat
            show_names <- TRUE
          } else {
            # No valid selected proteins, show all
            display_mat <- merged_mat
            show_names <- (nrow(display_mat) < 30)
          }
        } else {
          # Show all in UPS class with annotations
          gene_class_lookup <- proteins %>%
            dplyr::select(Gene.Symbol, Class) %>%
            dplyr::distinct()
          selected_classes <- unique(gene_class_lookup$Class[
            gene_class_lookup$Gene.Symbol %in% input$filter_proteins
          ])
          gene_info <- data.frame(gene = rownames(merged_mat), stringsAsFactors = FALSE) %>%
            dplyr::inner_join(gene_class_lookup, by = c("gene" = "Gene.Symbol"))
          keep_rows <- gene_info$Class %in% selected_classes
          merged_mat <- merged_mat[keep_rows, , drop = FALSE]
          display_mat <- merged_mat
          rn <- rownames(display_mat)
          rn[!(rn %in% input$filter_proteins)] <- ""
          rownames(display_mat) <- rn
          show_names <- TRUE
        }
      } else {
        display_mat <- merged_mat
        show_names <- (nrow(display_mat) < 30)
      }

      # Create column annotations
      annotations <- data.frame(
        Cohort = gsub("_.*", "", colnames(display_mat)),
        Cancer_subtype = gsub("^[^_]*_|_vs_Others$", "", colnames(display_mat)),
        row.names = colnames(display_mat)
      )

      # Subset columns based on cancer subtype filter
      if (!is.null(input$cancerSubtypeSelect) && length(input$cancerSubtypeSelect) > 0) {
        extracted_subtypes <- gsub("^[^_]*_|_vs_Others$", "", colnames(display_mat))
        keep_cols <- extracted_subtypes %in% input$cancerSubtypeSelect
        display_mat <- display_mat[, keep_cols, drop = FALSE]
        annotations <- annotations[keep_cols, , drop = FALSE]
      }

      # Create row annotation by joining with the gene-to-class lookup
      x <- proteins %>%
        dplyr::select(Gene.Symbol, Class) %>%
        dplyr::distinct()
      Class <- data.frame(
        gene = rownames(merged_mat),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::inner_join(x, by = c("gene" = "Gene.Symbol")) %>%
        tibble::column_to_rownames("gene")

      # Clean up cancer subtype labels
      annotations$Cancer_subtype <- gsub("_", "-", annotations$Cancer_subtype)
      annotations$Cancer_subtype[annotations$Cancer_subtype == "Serous-Adenocarcinoma"] <- "Serous Adenocarcinoma"
      annotations$Cancer_subtype[annotations$Cancer_subtype == "Terminal-Respiratory-Unit"] <- "Terminal Respiratory Unit"
      annotations$Cancer_subtype[annotations$Cancer_subtype == "IDH-mutant"] <- "IDH mutant"

      # Load color definitions from file.
      colors <- read.csv("./data/colors.csv")
      present_classes <- unique(as.character(Class$Class))
      present_cohorts <- unique(as.character(annotations$Cohort))
      present_subtypes <- unique(as.character(annotations$Cancer_subtype))

      colors.cohort <- setNames(
        object = colors$color[match(present_cohorts, colors$id)],
        nm = present_cohorts
      )
      colors.subtypes <- setNames(
        object = colors$color[match(present_subtypes, colors$id)],
        nm = present_subtypes
      )
      colors.class <- setNames(
        object = colors$color[match(present_classes, colors$id)],
        nm = present_classes
      )

      colors.cohort[is.na(colors.cohort)] <- "#000000"
      colors.subtypes[is.na(colors.subtypes)] <- "#000000"
      colors.class[is.na(colors.class)] <- "#000000"

      # Build the heatmap using the display matrix.
      ht <- ComplexHeatmap::Heatmap(
        display_mat,
        name = "LogFC",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = show_names,
        show_column_names = FALSE,
        split = Class$Class, # group rows by UPS Class
        row_gap = grid::unit(0, "mm"),
        row_title = NULL,
        column_title = "Differentially Expressed Proteins (DEPs) in cancer subtypes",
        column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          df = annotations[, c("Cohort", "Cancer_subtype")],
          col = list(
            Cohort = colors.cohort,
            Cancer_subtype = colors.subtypes
          ),
          which = "column"
        ),
        left_annotation = ComplexHeatmap::HeatmapAnnotation(
          Class = Class$Class,
          col = list(
            Class = colors.class
          ),
          which = "row"
        ),
        use_raster = FALSE
      )

      list(ht = ht, mat = merged_mat, class_df = Class)
    })

    # Render the heatmap.
    output$heatmapPlot <- renderPlot(
      {
        req(heatmap_obj())
        result <- heatmap_obj()
        grid::grid.newpage()
        ComplexHeatmap::draw(result$ht, newpage = TRUE)
      },
      height = 800
    )

    # Download handler for CSV export.
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("cancer_heatmap_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(heatmap_obj()$mat)
        write.csv(as.data.frame(heatmap_obj()$mat), file)
      }
    )

    # Download handler for PDF export.
    output$download <- downloadHandler(
      filename = function() {
        paste0("heatmap_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 11, height = 8)
        ComplexHeatmap::draw(heatmap_obj()$ht, newpage = TRUE)
        dev.off()
      }
    )
  })
}
