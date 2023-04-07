oncoplotUI <- function(id, cohorts){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      
      # gene inputs
      
      selectizeInput(inputId = ns("inputOncoGenesTXT"),
                     choices = NULL,
                     label =  "Input genes",options = list(`live-search` = TRUE),  multiple = TRUE),
      p("OR"),
      fileInput(ns("inputOncoGenesCSV"), "Upload tab delimited CSV or TSV", multiple = FALSE, accept = c('.csv','text/csv','text/comma-separated-values,text/plain'), width = NULL),
      selectInput(ns("oncoDelims"), "Delimiter", c("tab", "comma")),               
      
      # options
      sliderInput(ns("maxOncoGenes"), "Max Number of Genes to Display:",
                  min = 1, max = 50,
                  value = 10),
      h3("Cohorts"),
      selectizeInput(inputId = ns("oncoCohortSelect"),
                     choices = cohorts,selected = cohorts,
                     label = "Select cohorts",options = list(`live-search` = TRUE),  multiple = TRUE),
      
      actionButton(ns("generateOnco"), "Generate Plot"),
      
      # download options
      p("Download this data"),
      downloadButton(outputId = ns("downloadOncoMutations"), label = "Download Mutations Data"),               
      downloadButton(outputId = ns("downloadOncoCohorts"), label = "Download Cohort Data"),               
    ),
    mainPanel(
      textOutput(ns("oncoError"))%>% 
        tagAppendAttributes(class = 'error'),
      # loading icons
      shinycssloaders::withSpinner(plotOutput(ns("cohortOncoplot")), color = "#6A2798"),
      shinycssloaders::withSpinner(plotOutput(ns("mutationOncoplot")), color = "#6A2798")
    )
  )
}

source('./RScripts/OncoplotColors.R')

oncoplotServer <- function(id) {
  moduleServer(id,function(input,output,session){
    
    
    # this is being loaded twice which is not optimal
    load("./data/defaults.RData")
    cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
    
    genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
    
    maf <- readRDS("./data/oncoplot_maf.rds")
    
    # total patients
    nPatients = length(unique(maf$Tumor_Sample_Barcode))
    
    CRLS = genes
    crl_genes = unique(CRLS$hgnc_symbol)
    print("loaded oncoplot vars")
    
    ctbl = reactiveVal(NULL)
    mtbl = reactiveVal(NULL)
    renderedGenes = reactiveVal(NULL)
    mutationsreactive = reactiveVal(NULL)
    anno = reactiveVal(NULL)
    modifiedOncoInput = reactiveVal(FALSE)
    
    
    updateSelectizeInput(session, 'inputOncoGenesTXT', 
                         choices = as.list(rownames(cptac_protein)), 
                         server = TRUE)
    
    observe({
      input$generateOnco
      
      # see if the settings are at their defaults
      modifiedOncoInput(!(
        is.null(isolate(input$inputOncoGenesTXT)) &
          is.null(isolate(input$inputOncoGenesCSV)) &
          isolate(input$maxOncoGenes) == 10 &
          length(isolate(input$oncoCohortSelect)) == 11
      ))
      
      if(modifiedOncoInput()){
        
        # inputGenes = isolate(input$inputOncoGenesTXT)
        #inputGenes <- getInputGenes(isolate(input$inputOncoGenesCSV),isolate(input$inputOncoGenesTXT), delim = isolate(input$oncoDelims))
        inputGenes <- getInputGenes(isolate(input$inputOncoGenesCSV),isolate(input$inputOncoGenesTXT), delim = isolate(input$oncoDelims))
        
        if(length(inputGenes) == 0){
          inputGenes = crl_genes
        }
        
        # maf with only most frequently mutated genes
        smallFM = filter(maf, toupper(Hugo_Symbol) %in% inputGenes & cohort_code %in% isolate(input$oncoCohortSelect))
        
        uniqueGenes = unique(smallFM$Hugo_Symbol)
        
        mutationsTbl = matrix(ncol = length(unique(smallFM$Tumor_Sample_Barcode)), nrow = length(uniqueGenes))
        rownames(mutationsTbl) = uniqueGenes
        colnames(mutationsTbl) = unique(smallFM$Tumor_Sample_Barcode)
        cohortsTbl = mutationsTbl
        
        for(i in seq_along(rownames(smallFM))){
          id = smallFM[i,Tumor_Sample_Barcode]
          # create tables to be plotted
          mutationsTbl[smallFM[i,Hugo_Symbol],id] = levels(smallFM$Variant_Classification)[smallFM[i,Variant_Classification]]
          cohortsTbl[smallFM[i,Hugo_Symbol],id] = smallFM[i,cohort_code]
          
        }
        
        # set NA values to empty string
        mutationsTbl[is.na(mutationsTbl)] <- ""
        cohortsTbl[is.na(cohortsTbl)] <- ""
        
        cohortsTbl = as.data.frame(cohortsTbl)
        mutationsTbl = as.data.frame(mutationsTbl)
        
        # sort both tables the same way
        sortedTbl = as.data.frame(cohortsTbl) %>% arrange(rowSums(cohortsTbl == ""))
        topGenes = NULL
        if(nrow(sortedTbl) > isolate(input$maxOncoGenes)){
          topGenes = rownames(sortedTbl[1:isolate(input$maxOncoGenes),])
        }else {
          topGenes = rownames(sortedTbl)
        }
        
        cohortsTbl = filter(cohortsTbl, rownames(cohortsTbl) %in% topGenes)
        mutationsTbl = filter(mutationsTbl, rownames(mutationsTbl) %in% topGenes)
        
        # get percents
        pcts = matrix(ncol = 2, nrow = nrow(cohortsTbl))
        rownames(pcts) = rownames(cohortsTbl)
        pcts[,1] = rowSums(cohortsTbl != "")
        pcts[,2] = round((pcts[,1] / 1069) * 100, digits = 1)
        
        if(nrow(pcts) > 1)
          pcts = pcts[order(pcts[,2],decreasing=TRUE),]
        
        pcts[,2] = paste0(pcts[,2], "%")
        mutationsTbl = mutationsTbl[order(pcts[,2],decreasing=TRUE),]
        anno(pcts[,2])
        
        mutationsTbl = mutationsTbl[rownames(pcts),]
        cohortsTbl = cohortsTbl[rownames(pcts),]
        
        mutations = levels(smallFM$Variant_Classification)
        mutationsreactive(mutations)
        ctbl(cohortsTbl)
        mtbl(mutationsTbl)
        renderedGenes(topGenes)
        
      }
      
    })
    # mutations oncoplot
    output$mutationOncoplot = renderPlot({
      input$generateOnco
      
      if(!modifiedOncoInput()) return(mutationsOnco)
      
      mutation_colors = get_mutation_colors(mutationsreactive())
      mutation_alter_fun = get_mutation_alter_fun(mutation_colors)
      
      
      tryCatch(
        {
          output$oncoError = renderText("")
          
          oncoPrint(
            mtbl(), alter_fun = mutation_alter_fun, col = mutation_colors, 
            row_names_side = "left", 
            show_pct = FALSE,
            
            row_order = renderedGenes(),
            
            remove_empty_columns = TRUE,
            
            heatmap_legend_param = list(title = "Alterations"),
            
            right_annotation = rowAnnotation(pct = anno_text(anno()), barplot = anno_oncoprint_barplot())
          )
        },
        error = function(e){
          output$oncoError = renderText("No Data")
          return()
        }
      )
    })
    # cohorts oncoplot
    output$cohortOncoplot = renderPlot({
      input$generateOnco
      
      if(!modifiedOncoInput()) return(cohortsOnco)
      
      cohort_colors = get_cohort_colors()
      cohort_alter_fun = get_cohort_alter_fun(cohort_colors)
      
      tryCatch({
        oncoPrint(
          ctbl(), alter_fun = cohort_alter_fun, col = cohort_colors, 
          row_names_side = "left", 
          show_pct = FALSE,
          
          row_order = renderedGenes(),
          
          remove_empty_columns = TRUE,
          
          heatmap_legend_param = list(title = "Cohorts"),
          
          right_annotation = rowAnnotation(pct = anno_text(anno()), barplot = anno_oncoprint_barplot())
          
          
        )},
        error = function(e){
          return()
        }
      )
    })
    
    output$downloadOncoMutations <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(mtbl(), file)
      }
    )
    output$downloadOncoCohorts <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(ctbl(), file)
      }
    )
    
  })
  
}