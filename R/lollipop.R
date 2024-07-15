lollipopUI <- function(id){
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Lollipop Plot"),
      selectizeInput(
        inputId = ns("inputGene"),
        choices = NULL,
        multiple = FALSE,
        label = "Input gene",
        options = list(`live-search` = TRUE),
      ),
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("lollipop")),
        color = "#6A2798"
      ),
    )
  )
}

lollipopServer <- function(id) {
  moduleServer(id,function(input,output,session){

    maf = mem_readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")

    
    updateSelectizeInput(session, 'inputGene', 
                         choices = as.list(maf@gene.summary$Hugo_Symbol), 
                         server = TRUE, selected = "VHL")
      
    output$lollipop = renderPlot({
      inputGene = input$inputGene
      
      return (lollipopPlot(
        maf = maf,
        gene = inputGene,
        AACol = "Protein_Change",
        showMutationRate = TRUE,
        roundedRect = FALSE,
        # labelPos = 882
      ))
    })
  })
}
