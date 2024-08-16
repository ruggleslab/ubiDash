simplifiedHeatmapUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      p("(description paragraph)"),
      selectizeInput(
        inputId = ns("geneset"),
        choices = NULL,
        multiple = FALSE,
        selected = "E2",
        label = "Geneset",
        options = list(`live-search` = TRUE),
      ),
      # selectizeInput(
      #   inputId = ns("tissueType"),
      #   choices = c("normal"),
      #   multiple = FALSE,
      #   selected = "Normal",
      #   label = "Tissue Type",
      #   options = list(`live-search` = TRUE),
      # ),
      selectizeInput(
        inputId = ns("datatype"),
        choices = c("rna fpkm", "rna tpm", "pro ibaq"),
        multiple = FALSE,
        selected = "rna",
        label = "Data Type",
        options = list(`live-search` = TRUE),
      ),
      
      actionButton(ns("generate"), label = "Generate plot"),
      p("Download this data"),
      downloadButton(outputId = ns("download"), label = "Download Data"),
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("heatmap")),
                                   color = "#6A2798"
      ),
    )
  )
}


simplifiedHeatmapServer <-
  function(id,
           tertiles,
           # prt.global.AI,
           labels,
           my.colors) {
    moduleServer(id, function(input, output, session) {
      
      
      ups_genes <- read.xlsx("./data//Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06-tgr.xlsx",
                             sheet = "dense") %>% 
        filter(Branch == "Ubiquitin Proteasome System") %>% 
        select(Gene.Symbol, Class, Group, Type, Subtype)
      
      # below let user decide whether to split by Class, Group or Subtype
      gns <- split(ups_genes$Gene.Symbol, ups_genes[, "Class"]) 
      
      col_fun = colorRamp2(c(0, 50, 100), c("#EFECEA", "#86A6DF", "#2C4E80"))
      
      updateSelectizeInput(
        session,
        "geneset",
        choices = names(gns),
        # selected = c(),
        server = TRUE
      )
      
      missing = reactiveVal()
      data_type = reactiveVal()
      crl_num = reactiveVal()
      crl_total = reactiveVal()
      crl = reactiveVal()
      
      
      dataset = eventReactive(input$generate, {
        
        input_datatype = isolate(input$datatype)
        if(!isTruthy(input_datatype)){
          input_datatype = "rna_fpkm"
        }
        
        i = isolate(input$geneset)
        if(!isTruthy(i)){
          i = "DUBs and UBL demodifiers"
        }
        
        
        crl <- names(gns[i])
        crl_total <- length(gns[[i]])
        
        if(input_datatype == "rna fpkm"){
          dataset = mem_get_from_drive("pride_fpkm_rna_expression.rds")
          data_type(c("RNA","FPKM"))
          
        } else if(input_datatype == "rna tpm"){
          dataset = mem_get_from_drive("pride_tpm_rna_expression.rds") 
          data_type(c("RNA","TPM"))
          
        } else{
          dataset = mem_get_from_drive("pride_log2ibaq_protein_expression.rds")
          dataset = 2^dataset
          data_type(c("PROTEIN","IBAQ"))
        }
        
        
        crl_num(nrow(dataset))
        
        crl <- names(gns[i])
        
        crl_total <- length(gns[[i]])
        
        
        dataset <- dataset[rownames(dataset) %in% gns[[i]], ]
        crl_num <- nrow(dataset)
        
        dataset_percent <- as.data.frame(apply(t(dataset), 2, function(x) x / max(x) * 100))
        dataset_percent <- as.data.frame(t(dataset_percent))
        
        dataset_percent = na.omit(dataset_percent)
        
        missing(setdiff(gns[[i]], rownames(dataset_percent)))
        
        return(dataset_percent)
      },ignoreNULL = FALSE, ignoreInit=FALSE)
      
      
      output$heatmap = renderPlot({
        
        req(dataset())
        
        type <- "normal"
        
        
        ht <- draw(Heatmap(as.matrix(dataset()),
                           name = "ht",
                           col = col_fun,
                           column_title = glue("{crl()} {data_type()[1]} profiles of {type} human tissues\n
                          N (SR) = {crl_total()} of {crl_num()}"),
                           cluster_columns = F,
                           row_names_gp = gpar(fontsize = 9),
                           border_gp = gpar(col = "black"),
                           rect_gp = gpar(col = "white"),
                           heatmap_legend_param = list(title = glue("{data_type()[2]}-{data_type()[1]}\nLevels, %")
                                                       # legend_direction = "horizontal",
                                                       # title_position = "topleft"
                           ),
                           column_names_side =  "top",
                           column_names_rot = 45
        ))
        
          
        ht
       
      }, height = 1000)
      
      output$download <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          req(dataset())
          write.csv(dataset(), file)
        }
      )
    })
  }
