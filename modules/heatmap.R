heatmapUI <- function(id, cohorts){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Heatmap"),
      
      tabsetPanel(id = "GvGpanel",selected = "Protein",
                  tabPanel("Protein",
                           br(),
                           actionButton(ns("generateProteinHeatmap"), "Generate Protein Heatmap"),
                           br(),
                  ),
                  tabPanel("mRNA",
                           br(),
                           actionButton(ns("generateMRNAHeatmap"), "Generate mRNA Heatmap"),
                           br(),
                  ),
      ),
      div(
                           br(),
        # make sure these do not interfere
        selectizeInput(inputId = ns("inputGenesTXT"),
                       choices = NULL,
                       label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
        # p("OR"),
        # fileInput(ns("inputGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
        # # actionButton("resetHeatmapGenesCSV","clear"),
        # selectInput(ns("heatmapDelims"), "Delimiter", c("tab", "comma")),               
        # 
        # cohort seleciton checkboxes
        h3("Cohorts"),
        selectizeInput(inputId = ns("cohortSelect"),
                       choices = cohorts,selected = cohorts,
                       label = "Select cohorts",options = list(`live-search` = TRUE),  multiple = TRUE),
        checkboxInput(ns("cluster"), "Cluster columns",
                      value = TRUE),
        
        # toggle clustering and separation by tumor and normal
        checkboxInput(ns("separateNandT"), "Separate normal and tumor",
                      value = TRUE),
      ),
      
      
      p("Download this data"),
      downloadButton(outputId = ns("downloadData"), label = "Download")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput(ns("rowscut"))%>% 
        tagAppendAttributes(class = 'error'),
      textOutput(ns("err"))%>% 
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")), color = "#6A2798"),
    )
  )
  
}

heatmapServer <- function(id, cohorts, cptac_protein,cptac_rna,tnlabels, 
                          colorsList, filteredProteinMat, filteredMRNAMat, 
                          filteredTnlabels, defaultHeatmap, genesSubset) {
  moduleServer(id,function(input,output,session){
    
    updateSelectizeInput(session, 'inputGenesTXT', 
                         choices = as.list(rownames(cptac_protein)), 
                         server = TRUE)
    
    # default plot height: maybe automate this 
    plotHeight <- reactiveVal(956)
    
    # matrix with only genes to plot
    getGenesMat <- function(inputGenes,cohortSubsetLabels,protein) {
      if (length(inputGenes) == 0 &
          length(isolate(isolate(input$cohortSelect))) == length(cohorts)) {
        return(NULL)
      }
      
      if (length(inputGenes) == 0) {
        # print(cohortSubsetLabels)
        if(protein){
          return(filteredProteinMat[, rownames(cohortSubsetLabels)])
        }else{
          return(filteredMRNAMat[, rownames(cohortSubsetLabels)])
          
        }
      }
      if(protein){
        return(cptac_protein[inputGenes, rownames(cohortSubsetLabels)])
      } else{
        return(cptac_rna[inputGenes, rownames(cohortSubsetLabels)])
        
      }
    }
    
    protein = reactiveVal(TRUE)
    
    observeEvent(input$generateProteinHeatmap,{
      protein(TRUE)
    })
    observeEvent(input$generateMRNAHeatmap,{
      protein(FALSE)
    })
    
    
    output$heatmap <- renderPlot(
      {
        instanceTnLabels = NULL;
        if(protein()){
          legendName = "Log2(Relative Abundance)"
          instanceTnLabels = tnlabels;
        } else {
          legendName = "Log2(TPM + 1)"
          instanceTnLabels = filteredTnlabels;
          
        }
        
        # Button will trigger rerender
        input$generateProteinHeatmap
        input$generateMRNAHeatmap
        
        
        # get cohorts selected by user
        cohortSubsetLabels = NULL
        cohortSubsetLabels =getCohortPatients(instanceTnLabels, isolate(input$cohortSelect))
        
        
        if(nrow(cohortSubsetLabels) == 0){
          output$rowscut <- renderText({
            paste("Error: Please select at least one cohort")
          })
          return()
        }
        
        inputGenes <- getInputGenes(isolate(input$inputGenesCSV),isolate(input$inputGenesTXT), delim = isolate(input$heatmapDelims))
        inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels, protein())
        # print(inputGenesMat)
        
        if (is.null(inputGenesMat)) {
          if (isolate(input$cluster) && isolate(input$separateNandT) & protein()) {
            return(defaultHeatmap)
          }
          
          if(protein()){
            inputGenesMat = filteredProteinMat
          } else{
            inputGenesMat = filteredMRNAMat
          }
          
          
        } else {
          inputGenesMat <- as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
          plotHeight(nrow(inputGenesMat) * 12 + 200)
        }
        
        cutout <- inputGenes[!(inputGenes %in% toupper(rownames(inputGenesMat)))]
        
        if (length(cutout) > 0 & (is.null(isolate(input$inputGenesCSV)) | length(inputGenesMat) > 0)) {
          output$rowscut <- renderText({
            paste("No Data:", paste(cutout, collapse = ", "))
          })
        } else if (length(cutout) > 0) {
          output$rowscut <- renderText({
            paste("No Data:", paste(cutout, collapse = ", "), ". Make sure delimiter is correct.")
          })
        } else {
          output$rowscut <- renderText({
            ""
          })
        }
        
        ha <- HeatmapAnnotation(Type = cohortSubsetLabels$sample_type, Cohort = cohortSubsetLabels$cohort_code, Subtype = cohortSubsetLabels$Cancer_subtype, col = colorList)
        hm <- NULL
        
        colsplit = NULL
        if(isolate(input$separateNandT)){
          colsplit = getCohortPatients(instanceTnLabels,isolate(input$cohortSelect))$sample_type
        }
        
        
        inputGenesMat <- as.matrix(inputGenesMat)
        
        tryCatch(
          {
            output$err <- renderText({
              ""
            })
            return(drawHeatmap(inputGenesMat, ha, colsplit, isolate(input$cluster), legendName))
          },
          error = function(e) {
            output$err <- renderText({
              "Warning: Not enough data to cluster columns"
            })
            return(drawHeatmap(inputGenesMat, ha, colsplit, FALSE,legendName))
            
          }
        )
      },
      height = plotHeight
    )
    
    # output$downloadData <- downloadHandler(
    #   filename = function() {
    #     paste("data-", Sys.Date(), ".csv", sep = "")
    #   },
    #   content = function(file) {
    #     cohortSubsetLabels <- getCohortPatients(tnlabels, isolate(input$cohortSelect))
    #     inputGenes <- getInputGenes(isolate(input$inputGenesCSV),isolate(input$inputGenesTXT), delim = isolate(input$volcanoDelims))
    #     inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels)
    #     
    #     if (length(inputGenesMat) == 0) {
    #       inputGenesMat <- filteredMat
    #     }
    #     write.csv(inputGenesMat, file)
    #   }
    # )
  })
}