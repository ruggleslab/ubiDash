simplifiedHeatmapUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      p("description paragraph here"),
      selectizeInput(
        inputId = ns("inputGeneset"),
        choices = c("E1", "E2", "E3", "CRLs"),
        multiple = FALSE,
        selected = "E2",
        label = "Geneset",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("tissueType"),
        choices = c("Normal", "Tumor"),
        multiple = FALSE,
        selected = "Normal",
        label = "Tissue Type",
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
      actionButton(ns("generate"), label = "Generate plot"),
      
      
      p("Download this data"),
      downloadButton(outputId = ns("downloadMtert"), label = "Tertile Data"),
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
  function(id,
           tertiles,
           prt.global.AI,
           labels,
           my.colors) {
    moduleServer(id, function(input, output, session) {
      height <- reactiveVal(600)
      fontSize <- reactiveVal(10)
      
      mTert_reactive = reactiveVal(NULL)
      prt.n_reactive = reactiveVal(NULL)
      
      
      output$heatmap <- renderPlot({
        input$generate
        
        
        selected_geneset <- isolate(input$inputGeneset)
        height(1000)
        fontSize(5)
        geneset_df <-
          read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
        
        
        if (selected_geneset == "CRLs") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
          height(1000)
          fontSize(5)
          
          if (isolate(input$inputSubGeneset) != "all") {
            geneset_df = geneset_df[geneset_df$E3_Complex == isolate(input$inputSubGeneset),]
            height(800)
            fontSize(8)
          }
          
        }
        if (selected_geneset == "E1") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 2)
          height(300)
          fontSize(10)
        }
        if (selected_geneset == "E2") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 3)
          height(600)
          fontSize(10)
        }
        if (selected_geneset == "E3") {
          geneset_df <-
            read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx", sheet = 4)
          height(2500)
          fontSize(4)
        }
        
        ### update
        genes <- geneset_df$hgnc_symbol
        geneset <- selected_geneset ### update
        
        column_title <-
          glue("Human {isolate(input$tissueType)} Tissues-{geneset}")
        
        prt = NULL
        
        
       
          prt = prt.global.AI
     
        
        prt <- prt[rownames(prt) %in% genes,]
        
        
        
        selected_sample_type = isolate(input$tissueType)
        
        subset <-
          subset(
            labels,
            Dataset == "CPTAC" &
              Sample_type == selected_sample_type,
            select = c(Sample_id, Tissue_type)
          )
        subset <- subset[subset$Sample_id %in% colnames(prt),]
        
        prt <- prt[, subset$Sample_id]
        prt.n <-
          aggregate_by(input = prt, label = subset$Tissue_type)
        
        print(head(prt.n))
        
        # Apply function to each column of the matrix
        mTert <- apply(prt.n, 2, tertiles)
        row.names(mTert) <- row.names(prt.n)
        
        mTert_reactive(mTert)
        prt.n_reactive(prt.n)
        
        # e3 heatmap has more rows than can be displayed on the screen
        ht <- Heatmap(
          mTert,
          column_title = column_title,
          na_col = "grey",
          col = my.colors,
          row_names_gp = gpar(fontsize = fontSize()),
          column_names_gp = gpar(fontsize = 16),
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
          write.csv(mTert_reactive(), file)
        }
      )
      output$downloadPrtn =  downloadHandler(
        filename = function() {
          paste("prtn-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(prt.n_reactive(), file)
        }
      )
    })
    
    
  }
