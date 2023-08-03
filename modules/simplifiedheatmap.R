simplifiedHeatmapUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      
      selectizeInput(
        inputId = ns("inputGeneset"),
        choices = c("E1", "E2", "E3", "CRLs"),
        multiple = FALSE,
        selected = "E2",
        label =  "Geneset",
        options = list(`live-search` = TRUE),
      ),
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")),
                                   color = "#6A2798"),
    )
  )
}


simplifiedHeatmapServer <-
  function(id, prt.n, tertiles, mTert, my.colors) {
    moduleServer(id, function(input, output, session) {
      
      height = reactiveVal(600)
      
      
      output$heatmap = renderPlot({
        selected_geneset = input$inputGeneset
        height(1000)
        
        input = read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
        
        if (selected_geneset == "E1") {
          input <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 2) ### update
            height(300)
        }
        if (selected_geneset == "E2") {
          input <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 3) ### update
            height(600)
        }
        if (selected_geneset == "E3") {
          input <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 4) ### update
            height(2500)
        }
        
        
        ### update
        genes = input$hgnc_symbol
        geneset = selected_geneset  ### update
        column_title = glue("Human Healthy Tissues-{geneset}")
        
        # sub-setting mTert
        sub.dat <- mTert[rownames(mTert) %in% genes,]
        
        
        # e3 heatmap has more rows than can be displayed on the screen
        ht <- Heatmap(
          sub.dat,
          column_title = column_title,
          na_col = "grey",
          col = my.colors,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Protein"),
          rect_gp = gpar(col = "white"),
          column_names_rot = 45,
          border = T,
        )
        # ht
        
        # pdf(glue("../output/{geneset}-profile_human-healthy-tissue_cptac-dataset.pdf"))
     
        draw(ht)
        
      }, height = height)
      # dev.off()
      
      # save corresponding abundance table
      # sub.dat.table <- prt.n[rownames(prt.n) %in% genes, ]
      # write.csv(
      #   sub.dat.table,
      #   file = glue(
      #     "../output/{geneset}-profile_human-healthy-tissue_cptac-dataset.csv"
      #   ),
      #   row.names = T
      # )
      
    })
  }
