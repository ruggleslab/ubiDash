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
            choices = NULL,
            selected = NULL,
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
        checkboxInput(ns("zscores"), "Row Z Scores",
                      value = TRUE
        ),
        checkboxInput(ns("showSubtypeLegend"), "Show Subtype Legend",
                      value = FALSE
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
    
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput(ns("rowscut")) %>%
        tagAppendAttributes(class = "error"),
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")), color = "#6A2798"),
      p("Download this data"),
      downloadButton(outputId = ns("download"), label = "Download Data"),
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] != 'PRIDE'"),
        downloadButton(outputId = ns("downloadSampleInfo"), label = "Download Sample Info")
      )
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
        
        # validate(
        #   need(!isolate(input$separateNandT) || !isolate(input$zscores), 'Can\'t separate normal and tumor when Row Z Score are enabled.'),
        # )
        
        cptac_samples <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")
       
        
        ccle_samples <- mem_get_from_drive("ccle_sample_info.rds")
        
        cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
        cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
       
       
        
        
        instanceTnLabels = NULL
        dataset = NULL
        
        output$rowscut <- renderText({
          paste("")
        })
        
        if (isolate(input$dataset) == "CPTAC") {
          
          if (isolate(input$datatype) == "Protein") {
            legendName = "Log2(Relative Abundance)"
            dataset = cptac_protein
            
            
            samples = filter(cptac_samples, datatype == "protein")
            print(nrow(samples))
            missing_rows <- !(samples$harmonized_data_id %in% colnames(dataset))
            samples <- samples[!missing_rows, ]
            print(nrow(samples))
            samples = samples[!duplicated(samples$harmonized_data_id),]
            rownames(samples) = samples$harmonized_data_id
            
            
            tnlabels <- subset(samples, type %in% c("Tumor", "Normal"), select = c("type", "cohort", "cancer_subtype"))
            instanceTnLabels = tnlabels
            
            
          } else {
            legendName = "Log2(TPM + 1)"
            dataset = cptac_rna
            
            samples = filter(cptac_samples, datatype == "rna")
            missing_rows <- !(samples$harmonized_data_id %in% colnames(dataset))
            samples <- samples[!missing_rows, ]
            samples = samples[!duplicated(samples$harmonized_data_id),]
            rownames(samples) = samples$harmonized_data_id
            
            
            tnlabels <- subset(samples, type %in% c("Tumor", "Normal"), select = c("type", "cohort", "cancer_subtype"))
            filteredTnlabels = tnlabels[rownames(tnlabels) %in% colnames(cptac_rna),]
            
            instanceTnLabels = filteredTnlabels
          }
          
      
          
            # filteredProteinMat <- cptac_protein[genesSubset$genes, rownames(tnlabels)]
            # filteredProteinMat <- filteredProteinMat[which(rowMeans(!is.na(filteredProteinMat)) > 0.5), ]
       
          # filteredMRNAMat <- cptac_rna[genesSubset$genes, rownames(filteredTnlabels)]
          # # remove rows with over 50% NA
          # filteredMRNAMat <- filteredMRNAMat[which(rowMeans(!is.na(filteredMRNAMat)) > 0.5), ]
          
          cohortSubsetLabels = instanceTnLabels
          
          if (length(isolate(input$subtype)) != 0) {
            cohortSubsetLabels = filter(cohortSubsetLabels,
                                        cancer_subtype %in% isolate(input$subtype))
          }
          if (length(isolate(input$cohort)) != 0) {
            cohortSubsetLabels = filter(cohortSubsetLabels,
                                        cohort %in% isolate(input$cohort))
          }
          
          inputGenes <- isolate(input$genes)
          if (!isTruthy(inputGenes)) {
            inputGenes = genesSubset$genes
          }
          
          # print(cohortSubsetLabels)
          # e = rownames(cohortSubsetLabels) %in% colnames(dataset)
          # print("e")
          # print(colnames(dataset)[!e])
          # print(length(colnames(dataset)))
          # 
          # rownames_in_dataset = rownames(cohortSubsetLabels) %in% colnames(dataset)
          
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
              Type = cohortSubsetLabels$type,
              Cohort = cohortSubsetLabels$cohort,
              Subtype = cohortSubsetLabels$cancer_subtype,
              col = list(
                Type = colors,
                Subtype = colors,
                Cohort = colors
                ),
              
              show_legend = c(T,T,isolate(input$showSubtypeLegend))
            )
          
          hm <- NULL
          
          if(isolate(input$zscores)){
            inputGenesMat = z_score(inputGenesMat, 1)
            legendName = "Row Z Score"
           
          }
          
          colsplit <- NULL
          if (isolate(input$separateNandT)) {
            colsplit <- cohortSubsetLabels$type
          }
          
          data(inputGenesMat)
          sample_info(cohortSubsetLabels)
      
          
          output$err <- renderText({
            ""
          })
          drawHeatmap(
            inputGenesMat,
            ha,
            colsplit,
            isolate(input$clusterCol),
            isolate(input$clusterRow),
            legendName,
            F
            # output$err
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
          
          ha <- HeatmapAnnotation(subtype = ccle_samples[colnames(sub.ccle), ]$Cancer_subtype, show_legend = isolate(input$showSubtypeLegend))
          
          
          if(isolate(input$zscores)){
            sub.ccle = z_score(sub.ccle, 1)
            legendName = "Row Z Score"
            
          }
          data(sub.ccle)
          
          
          output$err <- renderText({
            ""
          })
          
          drawHeatmap(
            sub.ccle,
            ha,
            NULL,
            isolate(input$clusterCol),
            isolate(input$clusterRow),
            legendName,
            F
            # output$err
          )
          
          
        } else if (isolate(input$dataset) == "PRIDE") {
          
          
          if (isolate(input$datatype) == "Protein") {
            legendName <- "Log2(Relative Abundance)"
            dataset <- mem_get_from_drive("pride_log2ibaq_protein_expression.rds")
            
          } else {
            legendName <- "Log2(TPM + 1)"
            
            if(isolate(input$pridedataset) == "tpm"){
              dataset <- as.matrix(mem_get_from_drive("pride_tpm_rna_expression.rds"))
            }
            
            else{
              dataset <- as.matrix(mem_get_from_drive("pride_fpkm_rna_expression.rds"))
            }
            
          }
          
          
          
          inputGenes <- isolate(input$genes)
          
          if (!isTruthy(inputGenes)) {
            inputGenes <- genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
          }
          
          sub.pride <- dataset[inputGenes, ]
          
          
          if(isolate(input$zscores)){
            sub.pride = z_score(sub.pride, 1)
            legendName = "Row Z Score"
            
          }
          data(sub.pride)
          
          
          output$err <- renderText({
            ""
          })
          drawHeatmap(
            sub.pride,
            NULL,
            NULL,
            isolate(input$clusterCol),
            isolate(input$clusterRow),
            legendName,
            T
          )
        }
      },
      height = 1000
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
