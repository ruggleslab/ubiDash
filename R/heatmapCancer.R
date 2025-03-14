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
                  multiple = TRUE)
    ),
    mainPanel(
      shinycssloaders::withSpinner(plotOutput(ns("heatmapPlot"), height = "800px"), color = "#6A2798"),
      br(),
      downloadButton(ns("download"), "Download PDF"),
      downloadButton(ns("downloadData"), "Download Data")
    )
  )
}

heatmapCancerServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load raw data once and store in a reactive so that you can update inputs
    data_raw <- reactive({
      prot <- read.csv("./data/DEqMS.results_cancer_subtypes_with_tumotPurity.csv")  # adjust file path
      pn <- openxlsx::read.xlsx("./data/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
                                sheet = "dense")[, -1]%>%
        filter(Branch == "Ubiquitin Proteasome System") 
      list(prot = prot, pn = pn)
    })
    
    
    # Update select input choices once the raw data is available
    observe({
      dat <- data_raw()
      # For UPS Class: use unique values from pn$Class
      updateSelectInput(session, "classSelect",
                        choices = unique(dat$pn$Class),
                        selected = unique(dat$pn$Class))
      
      # For Cohort: use unique values from prot$cohort
      updateSelectInput(session, "cohortSelect",
                        choices = unique(dat$prot$cohort),
                        selected = unique(dat$prot$cohort))
      
      # For Cancer Subtype: Extract correctly from contrast
      # Extract what's before "_vs_Others" but after the first underscore
      cancer_subtypes <- unique(sapply(strsplit(dat$prot$contrast, "_vs_"), function(x) x[1]))
    
      
      updateSelectInput(session, "cancerSubtypeSelect",
                        choices = cancer_subtypes,
                        selected = cancer_subtypes)

    })
    
    # Create the heatmap when the user clicks "Generate Heatmap"
    heatmap_obj <- reactive({
      dat <- data_raw()
      prot <- dat$prot
      pn <- dat$pn
      
      # Filter the protein network (pn) by selected UPS Class
      if (!is.null(input$classSelect) && length(input$classSelect) > 0) {
        pn <- pn %>% dplyr::filter(Class %in% input$classSelect)
      }
      
      # Filter the differential expression data (prot) by selected Cohort
      if (!is.null(input$cohortSelect) && length(input$cohortSelect) > 0) {
        prot <- prot %>% dplyr::filter(cohort %in% input$cohortSelect)
      }
      
      # Join pn with prot; keep one row per Gene.Symbol
      proteins <- pn %>% 
        dplyr::group_by(Gene.Symbol) %>%
        dplyr::slice(1) %>%  
        dplyr::ungroup() %>%
        dplyr::left_join(prot, by = c("Gene.Symbol" = "gene"))
      
      # Build a matrix for each cohort, then merge into one matrix
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
      
      # Create column annotations
      annotations <- data.frame(
        Cohort = gsub("_.*", "", colnames(merged_mat)),
        Cancer_subtype = gsub("^[^_]*_|_vs_Others$", "", colnames(merged_mat)),
        row.names = colnames(merged_mat)
      )
      
      # If a cancer subtype filter is set, subset the columns accordingly
      if (!is.null(input$cancerSubtypeSelect) && length(input$cancerSubtypeSelect) > 0) {
        # Filter using the same transformation we used for UI
        extracted_subtypes <- gsub("^[^_]*_|_vs_Others$", "", colnames(merged_mat))
        keep_cols <- extracted_subtypes %in% input$cancerSubtypeSelect
        merged_mat <- merged_mat[, keep_cols, drop = FALSE]
        annotations <- annotations[keep_cols, , drop = FALSE]
      }
      
      # Create a row annotation based on the protein Class
      x <- proteins %>% dplyr::select(Gene.Symbol, Class) %>% dplyr::distinct()
      
      # Only keep the classes for genes that are actually in the matrix
      Class <- data.frame(
        gene = rownames(merged_mat),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::inner_join(x, by = c("gene" = "Gene.Symbol")) %>%
        tibble::column_to_rownames("gene")
      
      # Clean up Cancer_subtype labels if needed
      annotations$Cancer_subtype <- gsub("_", "-", annotations$Cancer_subtype)
      annotations$Cancer_subtype[annotations$Cancer_subtype == "Serous-Adenocarcinoma"] <- "Serous Adenocarcinoma"
      annotations$Cancer_subtype[annotations$Cancer_subtype == "Terminal-Respiratory-Unit"] <- "Terminal Respiratory Unit"
      annotations$Cancer_subtype[annotations$Cancer_subtype == "IDH-mutant"] <- "IDH mutant"
      
      # Load color definitions (assuming you have a colors.csv file)
      colors <- read.csv("./data/colors.csv")
      
      # Only use the classes that are actually present in the filtered data
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
      
      # Handle potential NA values in color mapping
      colors.cohort[is.na(colors.cohort)] <- "#000000"  # Default to black
      colors.subtypes[is.na(colors.subtypes)] <- "#000000"
      colors.class[is.na(colors.class)] <- "#000000"
      
     
      
      # Build the heatmap using ComplexHeatmap
      ht <- ComplexHeatmap::Heatmap(
        merged_mat,
        name = "LogFC",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        split = Class$Class,  # Use Class$Class directly to match the actual data
        row_gap = grid::unit(0, "mm"),
        row_title = NULL,
        column_title = "Differentially Expressed Proteins (DEPs) in cancer subtypes",  # fixed title
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
    output$heatmapPlot <- renderPlot({
      req(heatmap_obj())
      result <- heatmap_obj()
      ht <- result$ht
      merged_mat <- result$mat
      grid::grid.newpage()
      ComplexHeatmap::draw(ht, newpage = TRUE)
    },height = 800)
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("cancer_heatmap_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(heatmap_obj()$mat)
        write.csv(as.data.frame(heatmap_obj()$mat), file)
      }
    )
    
    # Download handler to export the heatmap as a SVG
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