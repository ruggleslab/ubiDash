oncoplotUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Oncoplot"),

      # gene inputs
      p("(description paragraph)"),
      
      sliderInput(ns("maxOncoGenes"), "Max Number of Genes to Display:",
        min = 1, max = 50,
        value = 10
      ),
      selectizeInput(
        inputId = ns("genes"),
        choices = NULL,
        label = "Custon Gene Set",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      # p("OR"),
      # fileInput(ns("inputOncoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      # selectInput(ns("oncoDelims"), "Delimiter", c("tab", "comma")),
      #
      # options
      selectizeInput(
        inputId = ns("cohorts"),
        choices = cohorts,
        selected = cohorts,
        label = "Select cohorts",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      checkboxInput(ns("sortByCohort"), "Sort by Cohort", value = FALSE),
      actionButton(ns("refresh"), "Generate Plot"),

      # download options
     
      # downloadButton(outputId = ns("downloadOncoCohorts"), label = "Download Cohort Data"),
    ),
    mainPanel(
      textOutput(ns("error")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("mutationOncoplot")), color = "#6A2798"),
       p("Download this data"),
      downloadButton(outputId = ns("downloadOncoMutations"), label = "Download Mutations Data"),
    )
  )
}


oncoplotServer <- function(id, top_crls, cohorts, colorsTbl) {
  moduleServer(id,function(input,output,session){

    
    updateSelectizeInput(session, 'genes', 
                         choices = cptac_protein_genes, 
                         server = TRUE)

    inputCohorts = reactiveVal()
    
    filteredMaf = eventReactive(input$refresh, {
      
      validate(
        need(length(isolate(input$genes)) == 0 || length(isolate(input$genes)) > 1, 'Two genes are needed to generate plot.'),
      )
      
      maf = mem_readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")

      inputCohorts(cohorts)
       
      filteredMaf = maf
      if (length(isolate(input$cohorts)) != 0 && length(isolate(input$cohorts)) != length(cohorts)) {
        inputCohorts(maf@clinical.data[cohort_code %in% isolate(input$cohorts)])
        filteredMaf = subsetMaf(maf, tsb = inputCohorts()$Tumor_Sample_Barcode)
      }
      
      return(filteredMaf)
      
    },ignoreNULL = FALSE, ignoreInit=FALSE)
    
    output$mutationOncoplot <- renderPlot({
      
      # cptac_samples <- mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")

      inputGenes <- isolate(input$genes)
      
      if (!isTruthy(inputGenes)) {
        inputGenes <- top_crls[1:isolate(input$maxOncoGenes)]
      }

      colors = colorsTbl$color
      names(colors) =  colorsTbl$id
    
      colors = colors[names(colors) %in% inputCohorts()]
      
      # vc_colors = c("#EF476F", "#F78C6B", "#FFD166", "#83D483", "#06D6A0", "#0CB0A9", "#118AB2","#073B4C", "#000000")
      vc_colors = c("#0AA0BF","#F6114A", "#F36E98", "#78B177", "#9862A2", "#F05006", "#25998F","#FCA00C", "#000000")
      names(vc_colors) = c(
        'Missense_Mutation',
        'Frame_Shift_Del',
        'Frame_Shift_Ins',
        'Nonsense_Mutation',
        'Splice_Site',
        'In_Frame_Ins',
        'In_Frame_Del',
        'Nonstop_Mutation',
        'Multi_Hit'
      )

      return(oncoplot(
        maf = filteredMaf(),
        # top = isolate(input$maxOncoGenes),
        genes = inputGenes,
        annotationColor = list(cohort_code = colors),
        sortByAnnotation = isolate(input$sortByCohort),
        gene_mar = 8,
        # barcode_mar = 20,
        fontSize  = 1,
        annotationFontSize = 1.5,
        legendFontSize = 1.5,
        colors = vc_colors,
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
