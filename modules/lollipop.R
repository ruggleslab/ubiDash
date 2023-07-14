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

lollipopServer<- function(id,maf) {
  moduleServer(id,function(input,output,session){
    
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
