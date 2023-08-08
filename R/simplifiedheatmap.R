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
        label = "Geneset",
        options = list(`live-search` = TRUE),
      ),
      conditionalPanel(
       condition = paste0("input['", ns("inputGeneset"), "'] == 'CRLs'"),
        selectizeInput(
          inputId = ns("inputSubGeneset"),
          choices = c("all", "CRL1", "CRL2-5", "CRL3", "CRL4"),
          multiple = FALSE,
          selected = "all",
          label = "E3 Complex",
          options = list(`live-search` = TRUE),
        ),
       # ns = ns
      ),

      p("Download this data"),
      downloadButton(outputId = ns("downloadMtert"), label = "Tertile data"),
      downloadButton(outputId = ns("downloadPrtn"), label = "Original Data")
      
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")),
                                   color = "#6A2798"),
    )
  )
}


simplifiedHeatmapServer <-
  function(id, prt.n, tertiles, mTert, my.colors) {
    moduleServer(id, function(input, output, session) {
      height <- reactiveVal(600)
      
      
      output$heatmap <- renderPlot({
        input$inputSubGeneset
        input$inputGeneset
        
        selected_geneset <- input$inputGeneset
        height(1000)
        geneset_df <-
          read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
        
        
        if (selected_geneset == "CRLs") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
          height(1000)
          
          if(input$inputSubGeneset != "all"){
            geneset_df = geneset_df[geneset_df$E3_Complex == input$inputSubGeneset,]
            height(800)
          }
          
        }
        if (selected_geneset == "E1") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 2)
          height(300)
        }
        if (selected_geneset == "E2") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 3)
          height(600)
        }
        if (selected_geneset == "E3") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 4)
          height(2500)
        }
        
        
        ### update
        genes <- geneset_df$hgnc_symbol
        geneset <- selected_geneset ### update
        column_title <- glue("Human Healthy Tissues-{geneset}")
        
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
     
        
        draw(ht)
      },
      height = height)
     
      output$downloadMtert =  downloadHandler(
        filename = function() {
          paste("mtert-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(mTert, file)
        }
      )
      output$downloadPrtn =  downloadHandler(
        filename = function() {
          paste("prtn-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(prt.n, file)
        }
      )
    })
    
    
  }
