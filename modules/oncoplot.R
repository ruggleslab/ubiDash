oncoplotUI <- function(id, cohorts){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      
      # gene inputs
      
      selectizeInput(inputId = ns("inputOncoGenesTXT"),
                     choices = NULL,
                     label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
      p("OR"),
      # fileInput(ns("inputOncoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      # selectInput(ns("oncoDelims"), "Delimiter", c("tab", "comma")),               
      # 
      # options
      sliderInput(ns("maxOncoGenes"), "Max Number of Genes to Display:",
                  min = 1, max = 50,
                  value = 10),
      
      checkboxInput(ns("sortByCohort"), "Sort by Cohort", value = FALSE),
      
      h3("Cohorts"),
      selectizeInput(inputId = ns("oncoCohortSelect"),
                     choices = cohorts,selected = cohorts,
                     label = "Select cohorts",options = list(`live-search` = TRUE),  multiple = TRUE),
      
      actionButton(ns("generateOnco"), "Generate Plot"),
      
      # download options
      p("Download this data"),
      downloadButton(outputId = ns("downloadOncoMutations"), label = "Download Mutations Data"),               
      # downloadButton(outputId = ns("downloadOncoCohorts"), label = "Download Cohort Data"),               
    ),
    mainPanel(
      textOutput(ns("oncoError"))%>% 
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotOutput(ns("mutationOncoplot")), color = "#6A2798")
      
    )
  )
}


oncoplotServer <- function(id,cptac_protein,cptac_samples,maf,crl_genes) {
  moduleServer(id,function(input,output,session){
     
    
    updateSelectizeInput(session, 'inputOncoGenesTXT', 
                         choices = as.list(rownames(cptac_protein)), 
                         server = TRUE)
    
    # mutations oncoplot
    output$mutationOncoplot = renderPlot({
      input$generateOnco
      
      inputGenes = isolate(input$inputOncoGenesTXT)
      # print(inputGenes)
      
      
      if(!isTruthy(inputGenes)){
        inputGenes = crl_genes
        # inputGenes = NULL
      }
      
      cohorts = maf@clinical.data[cohort_code %in% isolate(input$oncoCohortSelect)]
      filteredMaf = subsetMaf(maf, tsb = cohorts$Tumor_Sample_Barcode, genes = inputGenes)
      
      
      return(oncoplot(maf = filteredMaf,
                      top = isolate(input$maxOncoGenes),
                      # genes = inputGenes,
                      annotationColor = cohort_colors,
                      sortByAnnotation = isolate(input$sortByCohort),
                      clinicalFeatures = 'cohort_code'
      )
      )
    })
    
    output$downloadOncoMutations <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(maf, file)
      }
    )
    
  })
  
}