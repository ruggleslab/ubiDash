heatmapUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Heatmap"),
      tabsetPanel(
        id = ns("datatype"),
        selected = "Protein",
        tabPanel("Protein", ),
        tabPanel("mRNA", ),
      ),
      div(
        br(),
        # make sure these do not interfere
        selectizeInput(
          inputId = ns("genes"),
          choices = NULL,
          label = "Input genes",
          options = list(`live-search` = TRUE),
          multiple = TRUE
        ),
        # p("OR"),
        # fileInput(ns("inputGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
        # # actionButton("resetHeatmapGenesCSV","clear"),
        # selectInput(ns("heatmapDelims"), "Delimiter", c("tab", "comma")),
        #
        # cohort seleciton checkboxes
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
            label = "Select subtypes",
            options = list(`live-search` = TRUE),
            multiple = TRUE
          ),
        ),
<<<<<<< HEAD
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
=======
        checkboxInput(ns("clusterCol"), "Cluster columns", value = TRUE),
        
        checkboxInput(ns("clusterRow"), "Cluster rows", value = TRUE),
        
        # toggle clustering and separation by tumor and normal
        conditionalPanel(
          condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
          checkboxInput(ns("separateNandT"), "Separate normal and tumor", value = TRUE),
>>>>>>> develop
        ),
        actionButton(ns("generateHeatmap"), "Generate Heatmap"),
      ),
      p("Download this data"),
      downloadButton(outputId = ns("downloadData"), label = "Download Data"),
      downloadButton(outputId = ns("downloadAnno"), label = "Download Annotation")
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
           cptac_protein,
           cptac_rna,
           ccle_protein,
           ccle_rna,
           pride_protein,
           tnlabels,
           colorsList,
           filteredProteinMat,
           filteredMRNAMat,
           filteredTnlabels,
           defaultHeatmap,
           genesSubset) {
    moduleServer(id, function(input, output, session) {
      updateSelectizeInput(session,
        "genes",
        choices = as.list(rownames(cptac_protein)),
        server = TRUE
      )



      updateSelectizeInput(
        session,
        "subtype",
        choices = tnlabels$Cancer_subtype,
        selected = tnlabels$Cancer_subtype,
        server = TRUE
      )


      observe({
        input$dataset

        if (input$dataset == "CPTAC") {
          updateSelectizeInput(session,
            "genes",
            choices = as.list(rownames(cptac_protein)),
            server = TRUE
          )
        } else {
          updateSelectizeInput(session,
            "genes",
            choices = as.list(rownames(ccle_protein)),
            server = TRUE
          )
        }
      })
      observe({
        input$dataset

        if (input$dataset == "CPTAC") {
          updateSelectizeInput(
            session,
            "subtype",
            choices = tnlabels$Cancer_subtype,
            selected = c(),
            server = TRUE
          )
        } else {
          updateSelectizeInput(
            session,
            "subtype",
            choices = unique(ccle_samples$Subtype),
            selected = c(),
            server = TRUE
          )
        }
      })

      # default plot height: maybe automate this
      plotHeight <- reactiveVal(956)
<<<<<<< HEAD


      output$heatmap <- renderPlot(
        {
          input$generateHeatmap

          instanceTnLabels <- NULL
          dataset <- NULL




          if (isolate(input$dataset) == "CPTAC") {
            if (length(isolate(input$cohort)) == 0 |
              length(isolate(input$subtype)) == 0) {
              output$rowscut <- renderText({
                paste("Error: Please select at least one cohort and one subtype")
              })
              return()
            }

            if (isolate(input$datatype) == "Protein") {
              legendName <- "Log2(Relative Abundance)"
              dataset <- cptac_protein
              instanceTnLabels <- tnlabels
            } else {
              legendName <- "Log2(TPM + 1)"
              dataset <- cptac_rna
              instanceTnLabels <- filteredTnlabels
            }



            cohortSubsetLabels <- filter(
              instanceTnLabels,
              Cancer_subtype %in% isolate(input$subtype) &
                cohort_code %in% isolate(input$cohort)
=======
      
      
      output$heatmap <- renderPlot({
        input$generateHeatmap
        
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
          
          if(length(isolate(input$subtype)) != 0){
          cohortSubsetLabels = filter(
            instanceTnLabels,
            Cancer_subtype %in% isolate(input$subtype))
          
          }
          
          
          if(length(isolate(input$cohorts)) != 0){
            cohortSubsetLabels = filter(
              instanceTnLabels,
              cohort_code %in% isolate(input$cohort))
            
          }
          
          
          inputGenes <- isolate(input$genes)
          if (!isTruthy(inputGenes)) {
            inputGenes = genesSubset$genes
          }
          
          inputGenesMat <-
            dataset[inputGenes, rownames(cohortSubsetLabels)]
          inputGenesMat <-
            as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
          
          plotHeight(nrow(inputGenesMat) * 12 + 200)
          
          if (nrow(inputGenesMat) < length(isolate(input$genes))) {
            output$rowscut = renderText({
              "Some genes removed for NA rows."
            })
          }
          
          
          ha <-
            HeatmapAnnotation(
              Type = cohortSubsetLabels$sample_type,
              Cohort = cohortSubsetLabels$cohort_code,
              Subtype = cohortSubsetLabels$Cancer_subtype,
              col = colorList
>>>>>>> develop
            )


            inputGenes <- isolate(input$genes)
            if (!isTruthy(inputGenes)) {
              inputGenes <- genesSubset$genes
            }

            inputGenesMat <-
              dataset[inputGenes, rownames(cohortSubsetLabels)]
            inputGenesMat <-
              as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])

            plotHeight(nrow(inputGenesMat) * 12 + 200)

            if (nrow(inputGenesMat) < length(isolate(input$genes))) {
              output$rowscut <- renderText({
                "Some genes removed for NA rows."
              })
            }


            ha <-
              HeatmapAnnotation(
                Type = cohortSubsetLabels$sample_type,
                Cohort = cohortSubsetLabels$cohort_code,
                Subtype = cohortSubsetLabels$Cancer_subtype,
                col = colorList
              )
            hm <- NULL

            colsplit <- NULL
            if (isolate(input$separateNandT)) {
              colsplit <- cohortSubsetLabels$sample_type
            }


            tryCatch(
              {
                output$err <- renderText({
                  ""
                })
                return(drawHeatmap(
                  inputGenesMat,
                  ha,
                  colsplit,
                  isolate(input$clusterCol),
                  isolate(input$clusterRow),
                  legendName
                ))
              },
              error = function(e) {
                print(e)
                output$err <- renderText({
                  "Warning: Not enough data to cluster columns"
                })
                return(drawHeatmap(inputGenesMat, ha, colsplit, FALSE, FALSE, legendName))
              }
            )
          } else if (isolate(input$dataset) == "CCLE") {
            if (isolate(input$datatype) == "Protein") {
              legendName <- "Log2(Relative Abundance)"
              dataset <- ccle_protein
            } else {
              legendName <- "Log2(TPM + 1)"

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



            print(isolate(input$subtype))
            if (length(isolate(input$subtype)) != 0) {
              sub.samples <- ccle_samples[ccle_samples$Subtype %in% isolate(input$subtype), ]
              sub.samples <- sub.samples[rownames(sub.samples) %in% colnames(sub.ccle), ]
              sub.ccle <- sub.ccle[, rownames(sub.samples)]
            }

            sub.ccle <-
              as.matrix(sub.ccle[which(rowMeans(!is.na(sub.ccle)) > 0.5), ])

            print(ncol(sub.ccle))
            print(nrow(sub.ccle))


            ha <- HeatmapAnnotation(subtype = ccle_samples[colnames(sub.ccle), ]$Cancer_subtype)


            output$err <- renderText({
              ""
            })
<<<<<<< HEAD

            Heatmap(
              sub.ccle,
              top_annotation = ha,
              column_names_gp = grid::gpar(fontsize = 3),
              cluster_columns = isolate(input$clusterCol),
              cluster_rows = isolate(input$clusterRow)
            )
          } else if (isolate(input$dataset) == "PRIDE") {
            dataset <- (pride_protein)


            inputGenes <- isolate(input$genes)

            if (!isTruthy(inputGenes)) {
              inputGenes <- genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
            }

            sub.pride <- pride_protein[inputGenes, ]


            Heatmap(
              sub.pride,
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = isolate(input$clusterCol),
              cluster_rows = isolate(input$clusterRow)
            )
          }
        },
        height = plotHeight
      )



=======
            return(drawHeatmap(
              inputGenesMat,
              ha,
              colsplit,
              isolate(input$clusterCol),
              isolate(input$clusterRow),
              legendName
            ))
          }, error = function(e) {
            print(e)
            output$err <- renderText({
              "Warning: Not enough data to cluster columns"
            })
            return(drawHeatmap(inputGenesMat, ha, colsplit, FALSE, FALSE, legendName))
            
          })
        }  else if (isolate(input$dataset) == "CCLE") {
          if (isolate(input$datatype) == "Protein") {
            legendName = "Log2(Relative Abundance)"
            dataset = ccle_protein
            
          } else {
            legendName = "Log2(TPM + 1)"
            
            dataset = ccle_rna[rownames(ccle_rna) %in% ccle_samples$BROAD_ID, ]
            ccle_samples_sub = ccle_samples[ccle_samples$BROAD_ID %in% rownames(ccle_rna), ]
            ccle_samples_sub$id = rownames(ccle_samples_sub)
            rownames(ccle_samples_sub) = ccle_samples_sub$BROAD_ID
            ids = (ccle_samples_sub[rownames(ccle_rna), ])[, c("id", "BROAD_ID")]
            ids = ids[rowSums(is.na(ids)) != ncol(ids), ]
            rownames(dataset) = ids$id
            
            dataset = t(dataset)
          }
          
          inputGenes = isolate(input$genes)
          
          if (!isTruthy(inputGenes)) {
            inputGenes = genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
          }
          
          sub.ccle = dataset[inputGenes, ]
          
          
          
          print(isolate(input$subtype))
          if (length(isolate(input$subtype)) != 0) {
            sub.samples = ccle_samples[ccle_samples$Subtype %in% isolate(input$subtype), ]
            sub.samples = sub.samples[rownames(sub.samples) %in% colnames(sub.ccle), ]
            sub.ccle = sub.ccle[, rownames(sub.samples)]
          }
          
          sub.ccle <-
            as.matrix(sub.ccle[which(rowMeans(!is.na(sub.ccle)) > 0.5), ])
          
          
          ha = HeatmapAnnotation(subtype = ccle_samples[colnames(sub.ccle), ]$Cancer_subtype)
          
          
          output$err <- renderText({
            ""
          })
          
          Heatmap(
            sub.ccle,
            top_annotation = ha,
            column_names_gp = grid::gpar(fontsize = 3),
            cluster_columns = isolate(input$clusterCol),
            cluster_rows = isolate(input$clusterRow)
          )
          
          
        } else if (isolate(input$dataset) == "PRIDE") {
          dataset = (pride_protein)
          
          
          inputGenes = isolate(input$genes)
          
          if (!isTruthy(inputGenes)) {
            inputGenes = genesSubset$genes[genesSubset$genes %in% rownames(dataset)]
          }
          
          sub.pride = pride_protein[inputGenes, ]
          
          
          Heatmap(
            sub.pride,
            column_names_gp = grid::gpar(fontsize = 8),
            cluster_columns = isolate(input$clusterCol),
            cluster_rows = isolate(input$clusterRow)
          )
          
        }
        
      }, height = plotHeight)
      
      
      
>>>>>>> develop
      # Download buttons ----

      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(cptac_protein, file)
        }
      )
      output$downloadAnno <- downloadHandler(
        filename = function() {
          paste("annotation-data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(cptac_samples, file)
        }
      )
    })
  }
