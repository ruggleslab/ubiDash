simplifiedHeatmapUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Summarized Profile"),
      p("description paragraph here"),
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
      downloadButton(outputId = ns("downloadMtert"), label = "Tertile Data"),
      downloadButton(outputId = ns("downloadPrtn"), label = "Original Data")
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
      
      output$heatmap = renderPlot({
        
        i = isolate(input$geneset)
        if(!isTruthy(i)){
          i = "DUBs and UBL demodifiers"
        }
        
        type <- isolate(input$tissueType)
        crl <- names(gns[i])
        crl_total <- length(gns[[i]])
        
        
        input$generate
        
        if(isolate(input$datatype) == "rna fpkm"){
          
          rna_fpkm <- mem_get_from_drive("pride_fpkm_rna_expression.rds")
          rna_fpkm <- data.frame(rna_fpkm[rownames(rna_fpkm) %in% gns[[i]],])
          
          crl_num <- nrow(rna_fpkm)
          
          data_type <- "fpkm"
          
          rna_percent <- as.data.frame(apply(t(rna_fpkm), 2, function(x) x / max(x) * 100))
          rna_percent <- as.data.frame(t(rna_percent))
          
         
          
          missing <- setdiff(gns[[i]], rownames(rna_percent))
          

          ht <- draw(Heatmap(as.matrix(rna_percent),
                             name = "ht",
                             col = col_fun,
                             column_title = glue("{crl} RNA profiles of {type} human tissues\n
                            N (SR) = {crl_num} of {crl_total}"),
                             cluster_columns = F,
                             row_names_gp = gpar(fontsize = 9),
                             border_gp = gpar(col = "black"),
                             rect_gp = gpar(col = "white"),
                             heatmap_legend_param = list(title = "FPKM-RNA\nLevels, %"
                                                         # legend_direction = "horizontal",
                                                         # title_position = "topleft"
                             ),
                             column_names_side =  "top",
                             column_names_rot = 45
          ))
          
          return(ht)
        } else if(isolate(input$datatype) == "rna tpm"){
          
          rna_tpm <- mem_get_from_drive("pride_tpm_rna_expression.rds")                
          
          rna_tpm <- rna_tpm[rownames(rna_tpm) %in% gns[[i]],]
          
          crl_num <- nrow(rna_tpm)
          
          
          data_type <- "tpm"
          
          rna_percent <- as.data.frame(apply(t(rna_tpm), 2, function(x) x / max(x) * 100))
          rna_percent <- as.data.frame(t(rna_percent))
          
          missing <- setdiff(gns[[i]], rownames(rna_percent))
          
          ht <- draw(Heatmap(as.matrix(rna_percent),
                             name = "ht",
                             col = col_fun,
                             column_title = glue("{crl} RNA profiles of {type} human tissues\n
                              N (SR) = {crl_num} of {crl_total}"),
                             cluster_columns = F,
                             row_names_gp = gpar(fontsize = 9),
                             border_gp = gpar(col = "black"),
                             rect_gp = gpar(col = "white"),
                             heatmap_legend_param = list(title = "TPM-RNA\nLevels, %"
                                                         # legend_direction = "horizontal",
                                                         # title_position = "topleft"
                             ),
                             column_names_side =  "top",
                             column_names_rot = 45
          ))
          
          return(ht)
          
          
        } else {
          
          i = 1
          type <- "normal"
          crl <- names(gns[i])
          crl_total <- length(gns[[i]])
          
          pro_ibaq <- mem_get_from_drive("pride_log2ibaq_protein_expression.rds")
          
          pro_ibaq <- pro_ibaq[rownames(pro_ibaq) %in% gns[[i]], ]
          crl_num <- nrow(pro_ibaq)
          
          data_type <- "ibaq"
          
          pro_percent <- as.data.frame(apply(t(pro_ibaq), 2, function(x) x / max(x) * 100))
          pro_percent <- as.data.frame(t(pro_percent))
          
          pro_percent = na.omit(pro_percent)
          
          missing <- setdiff(gns[[i]], rownames(pro_percent))
          
          ht <- draw(Heatmap(as.matrix(pro_percent),
                              col = col_fun,
                              column_title = glue("{crl} Protein profiles of {type} human tissues\n
                              N (SR) = {crl_num} of {crl_total}"),
                              cluster_columns = F,
                              # cluster_rows = F,
                              row_names_gp = gpar(fontsize = 9),
                              border_gp = gpar(col = "black"),
                              rect_gp = gpar(col = "white"),
                              heatmap_legend_param = list(title = "Protein\nLevels, %"
                                                          # legend_direction = "horizontal",
                                                          # title_position = "topleft"
                              ),
                              column_names_side =  "top",
                              column_names_rot = 45
          ))
          
          return(ht)
        }
        
      }, height = 1000)
      
      
    })
  }
