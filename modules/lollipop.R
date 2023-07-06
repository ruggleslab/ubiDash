lollipopUI <- function(id, cohortSheets){
  ns <- NS(id)
  
  sidebarLayout(
    
    sidebarPanel(
      selectizeInput(inputId = ns("inputGene"),
                     choices = NULL,
                     multiple = FALSE,
                     label =  "Input gene",options = list(`live-search` = TRUE),),
    ),
    
    mainPanel(
      textOutput(ns("err"))%>% 
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotOutput(ns("lollipop")), color = "#6A2798"),
    )
  )
}

lollipopServer<- function(id) {
  moduleServer(id,function(input,output,session){
    
    maf = read.maf("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.maf")
    maf@gene.summary$Hugo_Symbol = gsub('"',"", maf@gene.summary$Hugo_Symbol)
    
    
    updateSelectizeInput(session, 'inputGene', 
                         choices = as.list(maf@gene.summary$Hugo_Symbol), 
                         server = TRUE, selected = NULL,)
    
    
    output$lollipop = renderPlot({
      inputGene = input$inputGene
      
      if(!isTruthy(inputGene)){ 
        inputGene = "VHL"
      }
      
      return (lollipopPlot(
        maf = maf,
        gene = inputGene,
        AACol = 'Protein_Change',
        showMutationRate = TRUE,
        # labelPos = 882
      ))
    })
  })
}
