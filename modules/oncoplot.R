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
      
      # checkboxInput(ns("showEmptyColumns"), "Show empty columns", value = TRUE),
      
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

source('./RScripts/OncoplotColors.R')

oncoplotServer <- function(id) {
  moduleServer(id,function(input,output,session){
    
    # this is being loaded twice which is not optimal
    # load("./data/defaults.RData")
    cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
    
    cptac_samples <- readRDS("./data/cptac_sample_info.rds")
    
    maf = read.maf("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.maf")
    patients = maf@clinical.data$Tumor_Sample_Barcode
    patients = data.frame(patients)
    colnames(patients) = "harmonized_data_id"
    
    maf@clinical.data$cohort_code=  merge(x=patients,y=cptac_samples, 
                                          by="harmonized_data_id", all.x=TRUE)$cohort_code
    
    genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
    
    
    CRLS = genes
    crl_genes = unique(CRLS$hgnc_symbol)
     
    
    updateSelectizeInput(session, 'inputOncoGenesTXT', 
                         choices = as.list(rownames(cptac_protein)), 
                         server = TRUE)
    
    # mutations oncoplot
    output$mutationOncoplot = renderPlot({
      input$generateOnco
      
      inputGenes = isolate(input$inputOncoGenesTXT)
      # print(inputGenes)
      
      if(!isTruthy(inputGenes)){
        # inputGenes = crl_genes
        inputGenes = NULL
      }
      
      
      cohort_colors = list(cohort_code = c(BRCA = "royalblue4",
                                           CCRCC = "maroon4",
                                           COAD = "goldenrod1",
                                           GBM = "aquamarine",
                                           HCC = "pink3",
                                           HNSCC = "coral4",
                                           LSCC = "slategray3",
                                           LUAD = "mediumpurple",
                                           OV = "seagreen",
                                           PDA = "lightgoldenrod2",
                                           UCEC = "orangered3"))
      
      return(oncoplot(maf = maf,
                      genes = inputGenes,
                      annotationColor = cohort_colors,
                      sortByAnnotation = TRUE,
                      clinicalFeatures = 'cohort_code',
                      top = isolate(input$maxOncoGenes)
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