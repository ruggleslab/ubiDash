lollipopUI <- function(id){
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Lollipop Plot"),
      p("(description paragraph)"),
      
      selectizeInput(
        inputId = ns("gene"),
        choices = NULL,
        multiple = FALSE,
        label = "Gene",
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


    
    updateSelectizeInput(session, 'gene', 
                         choices = maf_genes, 
                         server = TRUE, selected = "VHL")
      
    output$lollipop = renderPlot({
      
      maf = mem_readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")
      
      gene = input$gene
      
      return (lollipopPlot(
        maf = maf,
        gene = gene,
        AACol = "Protein_Change",
        showMutationRate = TRUE,
        roundedRect = FALSE,
        # labelPos = 882
      ))
    })
  })
}
