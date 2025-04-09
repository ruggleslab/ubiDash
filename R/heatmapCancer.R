heatmapCancerUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Customize UPS Protein Heatmap"),
      selectInput(ns("classSelect"), "UPS Class",
                  choices = NULL,  # will be updated dynamically
                  multiple = TRUE),
      selectInput(ns("cohortSelect"), "Cohort",
                  choices = NULL,  # will be updated dynamically
                  multiple = TRUE),
      selectInput(ns("cancerSubtypeSelect"), "Cancer Subtype",
                  choices = NULL,  # will be updated dynamically
                  multiple = TRUE),
      # NEW: Protein search input with server‐side search.
      selectizeInput(
        ns("filter_proteins"),
        label    = "Protein(s) of interest:",
        choices  = NULL,
        selected = NULL,
        multiple = TRUE,
        options  = list(
          placeholder = "Type or select one or more protein symbols",
          server = TRUE
        )
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
      prot <- read.csv("./data/DEqMS.results_cancer_subtypes_with_tumotPurity.csv")  # adjust file path
      pn <- openxlsx::read.xlsx(
        "./data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
        sheet = "dense"
      )[, -1] %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System")
      list(prot = prot, pn = pn)
    })
    
    # Update select input choices based on loaded data
    observe({
      dat <- data_raw()
      
      updateSelectInput(session, "classSelect",
                        choices = unique(dat$pn$Class),
                        selected = unique(dat$pn$Class))
      
      updateSelectInput(session, "cohortSelect",
                        choices = unique(dat$prot$cohort),
                        selected = unique(dat$prot$cohort))
      
      cancer_subtypes <- unique(sapply(strsplit(dat$prot$contrast, "_vs_"), function(x) x[1]))
      updateSelectInput(session, "cancerSubtypeSelect",
                        choices = cancer_subtypes,
                        selected = cancer_subtypes)
      
      # Update the protein search input (using Gene.Symbol from network data)
      updateSelectizeInput(session, "filter_proteins",
                           choices = sort(unique(dat$pn$Gene.Symbol)),
                           selected = NULL)
    })
    
    # Build the heatmap object
    heatmap_obj <- reactive({
      dat <- data_raw()
      prot <- dat$prot
      pn <- dat$pn
      
      # Filter the network data (pn) by selected UPS class if provided.
      if (!is.null(input$classSelect) && length(input$classSelect) > 0) {
        pn <- pn %>% dplyr::filter(Class %in% input$classSelect)
      }
      
      # Filter the expression data (prot) by selected cohorts.
      if (!is.null(input$cohortSelect) && length(input$cohortSelect) > 0) {
        prot <- prot %>% dplyr::filter(cohort %in% input$cohortSelect)
      }
      
      # Merge network and differential expression data (one row per gene)
      proteins <- pn %>% 
        dplyr::group_by(Gene.Symbol) %>%
        dplyr::slice(1) %>%  
        dplyr::ungroup() %>%
        dplyr::left_join(prot, by = c("Gene.Symbol" = "gene"))
      
      # Build a matrix for each cohort, then merge them.
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
      
      # --- NEW: Filter Entire UPS Class Based on Protein Selection ---
      # If protein(s) of interest are selected, determine their UPS class(es)
      # and filter the merged matrix to include all genes from that class.
      if (!is.null(input$filter_proteins) && length(input$filter_proteins) > 0) {
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
      }
      
      # --- NEW: Create a display matrix that only shows the selected protein names ---
      # In the display matrix, if protein(s) are selected, blank out row names that are NOT in the selection.
      display_mat <- merged_mat
      if (!is.null(input$filter_proteins) && length(input$filter_proteins) > 0) {
        rn <- rownames(display_mat)
        rn[!(rn %in% input$filter_proteins)] <- ""
        rownames(display_mat) <- rn
        show_names <- TRUE
      } else {
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
      x <- proteins %>% dplyr::select(Gene.Symbol, Class) %>% dplyr::distinct()
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
        split = Class$Class,  # group rows by UPS Class
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
    output$heatmapPlot <- renderPlot({
      req(heatmap_obj())
      result <- heatmap_obj()
      grid::grid.newpage()
      ComplexHeatmap::draw(result$ht, newpage = TRUE)
    }, height = 800)
    
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