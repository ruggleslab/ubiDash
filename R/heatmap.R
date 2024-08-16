heatmapUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Heatmap"),
      
      p("(description paragraph)"),
      
      tabsetPanel(
        id = ns("datatype"),
        selected = "mRNA",
        tabPanel("mRNA", ),
        tabPanel("Protein", ),
      ),
      div(
        br(),
        # make sure these do not interfere
        selectizeInput(
          inputId = ns("genes"),
          choices = NULL,
          label = "Genes",
          options = list(`live-search` = TRUE),
          multiple = TRUE
        ),
        
        selectizeInput(
          inputId = ns("dataset"),
          choices = c("CPTAC", "CCLE", "PRIDE"),
          selected = "CPTAC",
          label = "Dataset",
          options = list(`live-search` = TRUE),
          multiple = FALSE
        ),
        
        conditionalPanel(
          condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
          selectizeInput(
            inputId = ns("cohort"),
            choices = cohorts,
            selected = cohorts,
            label = "Select cohorts",
            options = list(`live-search` = TRUE),
            multiple = TRUE
          ),
        ),
        conditionalPanel(
          condition = paste0("input['", ns("dataset"), "'] != 'PRIDE'"),
          selectizeInput(
            inputId = ns("subtype"),
            choices = cohorts,
            selected = cohorts,
            label = "Subtypes",
            options = list(`live-search` = TRUE),
            multiple = TRUE
          ),
        ),
        conditionalPanel(
          condition = paste0("input['", ns("dataset"), "'] == 'PRIDE' && input['", ns("datatype"), "'] == 'mRNA'"),
          selectizeInput(
            inputId = ns("pridedataset"),
            choices = c("fpkm", "tpm"),
            selected = "fpkm",
            label = "Pride Dataset",
            options = list(`live-search` = TRUE),
            multiple = FALSE
          ),
        ),
        checkboxInput(ns("clusterCol"), "Cluster columns",
                      value = TRUE
        ),
        checkboxInput(ns("clusterRow"), "Cluster rows",
                      value = TRUE
        ),
        
        # toggle clustering and separation by tumor and normal
        conditionalPanel(
          condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
          checkboxInput(ns("separateNandT"), "Separate normal and tumor",
                        value = TRUE
          ),
        ),
        actionButton(ns("generateHeatmap"), "Generate Heatmap"),
      ),
      p("Download this data"),
      downloadButton(outputId = ns("download"), label = "Download Data"),
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] != 'PRIDE'"),
        downloadButton(outputId = ns("downloadSampleInfo"), label = "Download Sample Info")
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput(ns("rowscut")) %>%
        tagAppendAttributes(class = "error"),
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")), color = "#6A2798"),
    )
  )
}

heatmapServer <-
  function(id,
           cohorts,
           #  cptac_protein,
           #  cptac_rna,
           #  ccle_protein,
           #  ccle_rna,
           #  pride_protein,
           tnlabels,
           colorsTbl,
           filteredProteinMat,
           filteredMRNAMat,
           #  filteredTnlabels,
           # defaultHeatmap,
           genesSubset) {
    moduleServer(id, function(input, output, session) {
      
      observe({
        input$dataset
        
        if (input$dataset == "CPTAC") {
          if (input$datatype == "Protein") {
            updateSelectizeInput(session,
                                 "genes",
                                 choices = cptac_protein_genes,
                                 server = TRUE
            )
          } else{
            updateSelectizeInput(session,
                                 "genes",
                                 choices = cptac_rna_genes,
                                 server = TRUE
            )
          }
          
          updateSelectizeInput(
            session,
            "subtype",
            choices = tnlabels$Cancer_subtype,
            selected = c(),
            server = TRUE
          )
        } else {
          if (input$datatype == "Protein") {
            updateSelectizeInput(session,
                                 "genes",
                                 choices = ccle_protein_genes,
                                 server = TRUE
            )
          } else{
            updateSelectizeInput(session,
                                 "genes",
                                 choices = ccle_rna_genes,
                                 server = TRUE
            )
          }
          updateSelectizeInput(
            session,
            "subtype",
            choices = ccle_samples_subtypes,
            selected = c(),
            server = TRUE
          )
        }
      })
      
      
      # default plot height
      plotHeight <- reactiveVal(956)
      data <- reactiveVal()
      sample_info <- reactiveVal()
      
      
      output$heatmap <- renderPlot({
        input$generateHeatmap
        
        ccle_samples <- mem_get_from_drive("ccle_sample_info.rds")
        
        cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
        cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
       
        # pride_protein <- 
        
        
        filteredProteinMat <- cptac_protein[genesSubset$genes, rownames(tnlabels)]
        # remove rows with over 50% NA
        filteredProteinMat <- filteredProteinMat[which(rowMeans(!is.na(filteredProteinMat)) > 0.5), ]
        
        filteredTnlabels = tnlabels[rownames(tnlabels) %in% colnames(cptac_rna),]
        filteredMRNAMat <- cptac_rna[genesSubset$genes, rownames(filteredTnlabels)]
        # remove rows with over 50% NA
        filteredMRNAMat <- filteredMRNAMat[which(rowMeans(!is.na(filteredMRNAMat)) > 0.5), ]
        
        
        instanceTnLabels = NULL
        dataset = NULL
        
        output$rowscut <- renderText({
          paste("")
        })
        
        if (isolate(input$dataset) == "CPTAC") {
          
          if (isolate(input$datatype) == "Protein") {
            legendName = "Log2(Relative Abundance)"
            dataset = cptac_protein
            instanceTnLabels = tnlabels
            
          } else {
            legendName = "Log2(TPM + 1)"
            dataset = cptac_rna
            instanceTnLabels = filteredTnlabels
          }
          
          cohortSubsetLabels = instanceTnLabels
          
          if (length(isolate(input$subtype)) != 0) {
            cohortSubsetLabels = filter(cohortSubsetLabels,
                                        Cancer_subtype %in% isolate(input$subtype))
          }
          if (length(isolate(input$cohort)) != 0) {
            cohortSubsetLabels = filter(cohortSubsetLabels,
                                        cohort_code %in% isolate(input$cohort))
          }
          
          inputGenes <- isolate(input$genes)
          if (!isTruthy(inputGenes)) {
            inputGenes = genesSubset$genes
          }
          
          inputGenesMat <-
            dataset[inputGenes, rownames(cohortSubsetLabels)]
          
          if(isolate(input$clusterCol)){
            inputGenesMat <-
              as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
          }
          
          plotHeight(nrow(inputGenesMat) * 12 + 200)
          
          if (nrow(inputGenesMat) < length(isolate(input$genes))) {
            output$rowscut = renderText({
              "Unable to include some rows."
            })
          }
          
          colors = colorsTbl$color
          names(colors) = colorsTbl$id
          
          ha <-
            HeatmapAnnotation(
              Type = cohortSubsetLabels$sample_type,
              Cohort = cohortSubsetLabels$cohort_code,
              Subtype = cohortSubsetLabels$Cancer_subtype,
              col = list(
                Type = colors,
                Subtype = colors,
                Cohort = colors
                )
            )
          
          hm <- NULL
          
          colsplit <- NULL
          if (isolate(input$separateNandT)) {
            colsplit <- cohortSubsetLabels$sample_type
          }
          
          data(inputGenesMat)
          sample_info(cohortSubsetLabels)
          
          
          tryCatch(
            {
              output$err <- renderText({
                ""
              })
              drawHeatmap(
                inputGenesMat,
                ha,
                colsplit,
                isolate(input$clusterCol),
                isolate(input$clusterRow),
                legendName
              )
            },
            error = function(e) {
              print(e)
              output$err <- renderText({
                "Warning: Not enough data to cluster columns"
              })
              drawHeatmap(inputGenesMat, ha, colsplit, FALSE, FALSE, legendName)
            }
          )
          
        } else if (isolate(input$dataset) == "CCLE") {
          if (isolate(input$datatype) == "Protein") {
            legendName <- "Log2(Relative Abundance)"
            dataset <- mem_get_from_drive("ccle_tmt_protein_expression.rds")

          } else {
            legendName <- "Log2(TPM + 1)"
            ccle_rna <- mem_get_from_drive("ccle_tpm_rna_expression.rds")
            
            dataset <- ccle_rna[rownames(ccle_rna) %in% ccle_samples$BROAD_ID, ]
            ccle_samples_sub <- ccle_samples[ccle_samples$BROAD_ID %in% rownames(ccle_rna), ]
            ccle_samples_sub$id <- rownames(ccle_samples_sub)
            rownames(ccle_samples_sub) <- ccle_samples_sub$BROAD_ID
            ids <- (ccle_samples_sub[rownames(ccle_rna), ])[, c("id", "BROAD_ID")]
            ids <- ids[rowSums(is.na(ids)) != ncol(ids), ]
            rownames(dataset) <- ids$id
            
            dataset <- t(dataset)
          }
          
          inputGenes <- isolate(input$genes)
          if (!isTruthy(inputGenes)) {
            inputGenes <- genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
          }
          
          sub.ccle <- dataset[inputGenes, ]
          
          if (nrow(sub.ccle) < length(isolate(input$genes))) {
            output$rowscut = renderText({
              "Unable to include some rows."
            })
          }
          
          sample_info(ccle_samples)
          
          
          if (length(isolate(input$subtype)) != 0) {
            sub.samples <- ccle_samples[ccle_samples$Subtype %in% isolate(input$subtype), ]
            sub.samples <- sub.samples[rownames(sub.samples) %in% colnames(sub.ccle), ]
            sample_info(sub.samples)
            sub.ccle <- sub.ccle[, rownames(sub.samples)]
          }
          
          if(isolate(input$clusterCol)){
            sub.ccle <-
              as.matrix(sub.ccle[which(rowMeans(!is.na(sub.ccle)) > 0.5), ])
          }
          
          ha <- HeatmapAnnotation(subtype = ccle_samples[colnames(sub.ccle), ]$Cancer_subtype)
          
          data(sub.ccle)
          
          
          tryCatch(
            {
              output$err <- renderText({
                ""
              })
              Heatmap(
                sub.ccle,
                top_annotation = ha,
                column_names_gp = grid::gpar(fontsize = 3),
                show_column_names = FALSE, show_row_names = TRUE,
                cluster_columns = isolate(input$clusterCol),
                cluster_rows = isolate(input$clusterRow)
              )
            },
            error = function(e) {
              print(e)
              output$err <- renderText({
                "Warning: Not enough data to cluster columns"
              })
              Heatmap(
                sub.ccle,
                top_annotation = ha,
                column_names_gp = grid::gpar(fontsize = 3),
                show_column_names = FALSE, show_row_names = TRUE,
                cluster_columns = FALSE,
                cluster_rows = FALSE
              )
            }
          )
          
        } else if (isolate(input$dataset) == "PRIDE") {
          
          if (isolate(input$datatype) == "Protein") {
            dataset <- mem_get_from_drive("pride_log2ibaq_protein_expression.rds")
            
          } else {
            
            if(isolate(input$pridedataset) == "tpm")
              dataset <- as.matrix(mem_get_from_drive("pride_tpm_rna_expression.rds"))
            
            else{
              dataset <- as.matrix(mem_get_from_drive("pride_fpkm_rna_expression.rds"))
            }
            
          }
          
          
          inputGenes <- isolate(input$genes)
          
          if (!isTruthy(inputGenes)) {
            inputGenes <- genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
          }
          
          sub.pride <- dataset[inputGenes, ]
          
          data(sub.pride)
          
          tryCatch(
            {
              output$err <- renderText({
                ""
              })
              Heatmap(
                sub.pride,
                column_names_gp = grid::gpar(fontsize = 8),
                show_column_names = FALSE, show_row_names = TRUE,
                cluster_columns = isolate(input$clusterCol),
                cluster_rows = isolate(input$clusterRow)
              )
            },
            error = function(e) {
              print(e)
              output$err <- renderText({
                "Warning: Not enough data to cluster columns"
              })
              Heatmap(
                sub.pride,
                column_names_gp = grid::gpar(fontsize = 8),
                show_column_names = FALSE, show_row_names = TRUE,
                cluster_columns = FALSE,
                cluster_rows = FALSE
              )
            }
          )
          
          
        }
      },
      height = plotHeight
      )
      
      
      # Download buttons ----
      
      output$downloadSampleInfo <- downloadHandler(
        filename = function() {
          paste("sample-info-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(sample_info(), file)
        }
      )
      output$download <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          if(isTruthy(data())){
            write.csv(data(), file)
          }
        }
      )
    })
  }
