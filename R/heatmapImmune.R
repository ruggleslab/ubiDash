heatmapImmuneUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Customize UPS Protein Heatmap"),
      
      # Existing filters
      selectInput(ns("classSelectImmune"), "UPS Class",
                  choices = NULL,  # dynamically updated later
                  multiple = TRUE),
      selectInput(ns("cohortSelectImmune"), "Cohort",
                  choices = NULL,  # dynamically updated later
                  multiple = TRUE),
      selectInput(ns("ImmuneSubtypeSelectImmune"), "Immune Subtype",
                  choices = NULL,  # dynamically updated later
                  multiple = TRUE),
      
      # New: Protein search input (with server-side search capability)
      selectizeInput(
        ns("filter_proteins"),
        label    = "Protein(s) of interest:",
        choices  = NULL,
        selected = NULL,  # full heatmap by default
        multiple = TRUE,
        options  = list(
          placeholder = "Type or select one or more protein symbols",
          server = TRUE
        )
      )
    ),
    mainPanel(
      shinycssloaders::withSpinner(
        plotOutput(ns("heatmapPlotImmune"), height = "800px"),
        color = "#6A2798"
      ),
      downloadButton(ns("downloadImmune"), "Download PDF"),
      downloadButton(ns("downloadDataImmune"), "Download Data")
    )
  )
}
heatmapImmuneServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load raw data once as a reactive expression.
    data_raw <- reactive({
      prot <- read.csv("./data/DEqMS.results_Immune_subtype_percan_with_tumotPurity.csv")  # adjust file path if needed
      pn <- openxlsx::read.xlsx(
        "./data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
        sheet = "dense"
      )[, -1] %>%
        dplyr::filter(Branch == "Ubiquitin Proteasome System")
      list(prot = prot, pn = pn)
    })
    
    # Update input choices using the loaded data.
    observe({
      dat <- data_raw()
      
      updateSelectInput(session, "classSelectImmune",
                        choices = unique(dat$pn$Class),
                        selected = unique(dat$pn$Class))
      
      updateSelectInput(session, "cohortSelectImmune",
                        choices = unique(dat$prot$cohort),
                        selected = unique(dat$prot$cohort))
      
      immune_subtypes <- unique(sapply(strsplit(dat$prot$contrast, "_vs_"), function(x) x[1]))
      updateSelectInput(session, "ImmuneSubtypeSelectImmune",
                        choices = immune_subtypes,
                        selected = immune_subtypes)
      
      # Update the protein search input from the gene symbols in the network.
      updateSelectizeInput(session, "filter_proteins",
                           choices = sort(unique(dat$pn$Gene.Symbol)),
                           selected = NULL)
    })
    
    # Build the heatmap object.
    heatmap_obj <- reactive({
      dat <- data_raw()
      prot <- dat$prot
      pn <- dat$pn
      
      # Filter the protein network by selected UPS Class.
      if (!is.null(input$classSelectImmune) && length(input$classSelectImmune) > 0) {
        pn <- pn %>% dplyr::filter(Class %in% input$classSelectImmune)
      }
      
      # Filter the differential expression data by selected cohorts.
      if (!is.null(input$cohortSelectImmune) && length(input$cohortSelectImmune) > 0) {
        prot <- prot %>% dplyr::filter(cohort %in% input$cohortSelectImmune)
      }
      
      # Join network data (pn) with proteomics data (prot) by gene symbol.
      proteins <- pn %>% 
        dplyr::group_by(Gene.Symbol) %>%
        dplyr::slice(1) %>%  
        dplyr::ungroup() %>%
        dplyr::left_join(prot, by = c("Gene.Symbol" = "gene"))
      
      # Build a matrix for each cohort and merge into one matrix.
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
      
      # --- New Filtering Step: Select Entire UPS Class Based on Protein Selection ---
      # If one or more proteins of interest are chosen, find their corresponding UPS class(es)
      # and filter merged_mat to include all genes belonging to those classes.
      if (!is.null(input$filter_proteins) && length(input$filter_proteins) > 0) {
        # Create a lookup table for gene symbol to UPS class.
        gene_class_lookup <- proteins %>% 
          dplyr::select(Gene.Symbol, Class) %>% 
          dplyr::distinct()
        # Determine the UPS classes for the selected proteins.
        selected_classes <- unique(
          gene_class_lookup$Class[
            gene_class_lookup$Gene.Symbol %in% input$filter_proteins
          ]
        )
        
        # For each gene in the merged matrix, get its UPS class.
        gene_info <- data.frame(gene = rownames(merged_mat), stringsAsFactors = FALSE) %>%
          dplyr::inner_join(gene_class_lookup, by = c("gene" = "Gene.Symbol"))
        # Keep rows where the gene's class is among the selected classes.
        keep_rows <- gene_info$Class %in% selected_classes
        merged_mat <- merged_mat[keep_rows, , drop = FALSE]
      }
      
      # Create column annotations.
      annotations <- data.frame(
        Cohort = gsub("_.*", "", colnames(merged_mat)),
        Immune_subtype = gsub("^[^_]*_|_vs_Others$", "", colnames(merged_mat)),
        row.names = colnames(merged_mat)
      )
      
      # Apply the Immune Subtype filter if set.
      if (!is.null(input$ImmuneSubtypeSelectImmune) && length(input$ImmuneSubtypeSelectImmune) > 0) {
        extracted_subtypes <- gsub("^[^_]*_|_vs_Others$", "", colnames(merged_mat))
        keep_cols <- extracted_subtypes %in% input$ImmuneSubtypeSelectImmune
        merged_mat <- merged_mat[, keep_cols, drop = FALSE]
        annotations <- annotations[keep_cols, , drop = FALSE]
      }
      
      # Create row annotation (lookup of gene symbol to UPS class).
      x <- proteins %>% dplyr::select(Gene.Symbol, Class) %>% dplyr::distinct()
      Class <- data.frame(
        gene = rownames(merged_mat), stringsAsFactors = FALSE
      ) %>%
        dplyr::inner_join(x, by = c("gene" = "Gene.Symbol")) %>%
        tibble::column_to_rownames("gene")
      
      annotations$Immune_subtype <- gsub("_", "-", annotations$Immune_subtype)
      
      # Load color definitions.
      colors <- read.csv("./data/colors.csv")
      present_classes <- unique(as.character(Class$Class))
      present_cohorts <- unique(as.character(annotations$Cohort))
      present_subtypes <- unique(as.character(annotations$Immune_subtype))
      
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
      
      # --- NEW: Create a display matrix that only shows selected protein names ---
      display_mat <- merged_mat
      if (!is.null(input$filter_proteins) && length(input$filter_proteins) > 0) {
        # Clean rownames and selected proteins by trimming whitespace and converting to uppercase.
        rn <- toupper(trimws(rownames(display_mat)))
        sel <- toupper(trimws(input$filter_proteins))
        print(rn)
        print(sel)
        # Set row names: only retain the name if it is in the selected proteins; otherwise blank.
        rn_new <- ifelse(rn %in% sel, rownames(display_mat), "")
        rownames(display_mat) <- rn_new
        show_names <- TRUE
      } else {
        show_names <- (nrow(display_mat) < 30)
      }
      
      # Build the heatmap using the display matrix (which may have blank row names).
      ht <- ComplexHeatmap::Heatmap(
        display_mat,
        name = "LogFC",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = show_names,
        show_column_names = FALSE,
        split = Class$Class,  # group rows by UPS class
        row_gap = grid::unit(0, "mm"),
        row_title = NULL,
        column_title = "Differentially Expressed Proteins (DEPs) in Immune Subtypes",
        column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          df = annotations[, c("Cohort", "Immune_subtype")],
          col = list(
            Cohort = colors.cohort,
            Immune_subtype = colors.subtypes
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
    output$heatmapPlotImmune <- renderPlot({
      req(heatmap_obj())
      result <- heatmap_obj()
      grid::grid.newpage()
      ComplexHeatmap::draw(result$ht, newpage = TRUE)
    }, height = 800)
    
    # Download handler for CSV export.
    output$downloadDataImmune <- downloadHandler(
      filename = function() {
        paste("immune_heatmap_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(heatmap_obj()$mat)
        write.csv(as.data.frame(heatmap_obj()$mat), file)
      }
    )
    
    # Download handler for PDF export.
    output$downloadImmune <- downloadHandler(
      filename = function() {
        paste0("immune_heatmap_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 11, height = 8)
        ComplexHeatmap::draw(heatmap_obj()$ht, newpage = TRUE)
        dev.off()
      }
    )
    
  })
}