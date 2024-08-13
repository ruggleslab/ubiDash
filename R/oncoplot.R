oncoplotUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Oncoplot"),

      # gene inputs

      selectizeInput(
        inputId = ns("genes"),
        choices = NULL,
        label = "Genes",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      # p("OR"),
      # fileInput(ns("inputOncoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      # selectInput(ns("oncoDelims"), "Delimiter", c("tab", "comma")),
      #
      # options
      sliderInput(ns("maxOncoGenes"), "Max Number of Genes to Display:",
        min = 1, max = 50,
        value = 10
      ),
      checkboxInput(ns("sortByCohort"), "Sort by Cohort", value = FALSE),
      h3("Cohorts"),
      selectizeInput(
        inputId = ns("cohorts"),
        choices = cohorts,
        selected = cohorts,
        label = "Select cohorts",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      actionButton(ns("refresh"), "Generate Plot"),

      # download options
      p("Download this data"),
      downloadButton(outputId = ns("downloadOncoMutations"), label = "Download Mutations Data"),
      # downloadButton(outputId = ns("downloadOncoCohorts"), label = "Download Cohort Data"),
    ),
    mainPanel(
      textOutput(ns("error")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("mutationOncoplot")), color = "#6A2798")
    )
  )
}


oncoplotServer <- function(id,cptac_samples, top_crls, cohorts) {
  moduleServer(id,function(input,output,session){

    
    updateSelectizeInput(session, 'genes', 
                         choices = cptac_protein_genes, 
                         server = TRUE)
    
    filteredMaf = eventReactive(input$refresh, {
      
      validate(
        need(length(isolate(input$genes)) == 0 || length(isolate(input$genes)) > 1, 'Two genes are needed to generate plot.'),
      )
      
      maf = mem_readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")
    
      filteredMaf = maf
      if (length(isolate(input$cohorts)) != 0 && length(isolate(input$cohorts)) != length(cohorts)) {
        inputCohorts <- maf@clinical.data[cohort_code %in% isolate(input$cohorts)]
        filteredMaf = subsetMaf(maf, tsb = inputCohorts$Tumor_Sample_Barcode)
      }
      
      return(filteredMaf)
      
    },ignoreNULL = FALSE, ignoreInit=FALSE)
    
    output$mutationOncoplot <- renderPlot({

      inputGenes <- isolate(input$genes)
      
      if (!isTruthy(inputGenes)) {
        inputGenes <- top_crls[1:isolate(input$maxOncoGenes)]
      }
      
      return(oncoplot(
        maf = filteredMaf(),
        # top = isolate(input$maxOncoGenes),
        genes = inputGenes,
        annotationColor = cohort_colors,
        sortByAnnotation = isolate(input$sortByCohort),
        # drawColBar = F,
        # cBioPortal = T,
        clinicalFeatures = "cohort_code"
      ))
    })

    output$downloadOncoMutations <- downloadHandler(

      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(filteredMaf()@data, file)
      }
    )
  })
}
