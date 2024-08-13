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
        choices = c("normal"),
        multiple = FALSE,
        selected = "Normal",
        label = "Tissue Type",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("datatype"),
        choices = c("rna", "pro"),
        multiple = FALSE,
        selected = "rna",
        label = "Data Type",
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
      
      
      output$heatmap = renderPlot({
        
        
        input$generate
        
        if(isolate(input$datatype) == "rna"){
          # RNA----
          # rna <- read.xlsx("../../datasets/pride/Table_EV2_RNAseq_across_tissues.xlsx",
          #                  sheet = 3, check.names = F)
          
          # for (i in seq_along(gns)) {
          i = 1
          
          type <- isolate(input$tissueType)
          crl <- names(gns[i])
          crl_total <- length(gns[[i]])
          
          # FPKM values
          # rna_fpkm <- data.frame(rna[, 2:31]) %>%
          #   group_by(Gene.name) %>%
          #   summarise_all(mean) %>%
          #   column_to_rownames(var = "Gene.name")
          # 
          # saveRDS(rna_fpkm, file = "~/Google Drive/My Drive/Ruggles Lab Team Folder/projects/tania/ubiDash/data/shiny_accesible_data/pride_fpkm_rna_expression.rds")
          
          rna_fpkm <- mem_get_from_drive("pride_fpkm_rna_expression.rds")
          # rna_fpkm <- data.frame(rna[rna$Gene.name %in% gns[[i]], 2:31], row.names = 1)
          
          crl_num <- nrow(rna_fpkm)
          
          data_type <- "fpkm"
          
          rna_percent <- as.data.frame(apply(t(rna_fpkm), 2, function(x) x / max(x) * 100))
          rna_percent <- as.data.frame(t(rna_percent))
          # write.csv(rna_percent, file = glue("./output/{date}_{crl}_rna_{data_type}_percent_levels.csv"))
          
          missing <- setdiff(gns[[i]], rownames(rna_percent))
          # write(missing, file = glue("./output/{date}_{crl}_rna_{data_type}_percent_levels_missing_genes.txt"))
          
          library(circlize)
          col_fun = colorRamp2(c(0, 50, 100), c("#EFECEA", "#86A6DF", "#2C4E80"))
          
          # svg(filename = glue("./output/{date}_{crl}_rna_{data_type}_percent_levels.svg"), width = 8, height = 12)
          # ht <- draw(Heatmap(as.matrix(rna_percent),
          #                    name = "ht",
          #                    col = col_fun,
          #                    column_title = glue("{crl} RNA profiles of {type} human tissues\n
          #                             N (SR) = {crl_num} of {crl_total}"),
          #                    cluster_columns = F,
          #                    row_names_gp = gpar(fontsize = 9),
          #                    border_gp = gpar(col = "black"),
          #                    rect_gp = gpar(col = "white"),
          #                    heatmap_legend_param = list(title = "FPKM-RNA\nLevels, %"
          #                                                # legend_direction = "horizontal",
          #                                                # title_position = "topleft"
          #                    ),
          #                    column_names_side =  "top",
          #                    column_names_rot = 45
          # ))
          # 
          # dev.off()
          # ht
          
          # TPM values
          # rna_tpm <- as.data.frame(apply(rna[,3:31], 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)) %>% 
          #   mutate(.before=1, Gene.name = rna$Gene.name)
          # 
          # rna_tpm <- rna_tpm %>%
          #   group_by(Gene.name) %>%
          #   summarise_all(mean) %>%
          #   column_to_rownames(var = "Gene.name")
          # 
          # saveRDS(rna_tpm, file = "~/Google Drive/My Drive/Ruggles Lab Team Folder/projects/tania/ubiDash/data/shiny_accesible_data/pride_tpm_rna_expression.rds")
          
          rna_tpm <- mem_get_from_drive("pride_tpm_rna_expression.rds")                
          
          rna_tpm <- rna_tpm[rownames(rna_tpm) %in% gns[[i]],]
          
          
          data_type <- "tpm"
          
          rna_percent <- as.data.frame(apply(t(rna_tpm), 2, function(x) x / max(x) * 100))
          rna_percent <- as.data.frame(t(rna_percent))
          # write.csv(rna_percent, file = glue("./output/{date}_{crl}_rna_{data_type}_percent_levels.csv"))
          
          # missing <- setdiff(gns[[i]], rownames(rna_percent))
          # write(missing, file = glue("./output/{date}_{crl}_rna_{data_type}_levels_missing_genes.txt"))
          
          # svg(filename = glue("./output/{date}_{crl}_rna_{data_type}_percent_levels.svg"), width = 8, height = 12)
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
          
          # dev.off()
          # ht
          # break()
          # }
        } else {
          
          # Protein----
          # pro <- read.xlsx("../../datasets/pride/Table_EV1_protein_iBAQ_across_tissues.xlsx",
          #                  sheet = 4, check.names = F)
          
          # for (i in seq_along(gns)) {
          
          i = 1
          type <- isolate(input$tissueType)
          crl <- names(gns[i])
          crl_total <- length(gns[[i]])
          
          # pro_ibaq <- data.frame(pro[, 2:31]) %>%
          #   group_by(Gene.name) %>%
          #   summarise_all(mean) %>%
          #   column_to_rownames(var = "Gene.name")
          # 
          # saveRDS(pro_ibaq, file = "~/Google Drive/My Drive/Ruggles Lab Team Folder/projects/tania/ubiDash/data/shiny_accesible_data/pride_ibaq_protein_expression.rds")
          
          pro_ibaq <- mem_get_from_drive("pride_log2ibaq_protein_expression.rds")
          
          # pro_ibaq <- pro[rownames(pro) %in% gns[[i]], ]
          
          
          
          crl_num <- nrow(pro_ibaq)
          
          data_type <- "ibaq"
          
          pro_percent <- as.data.frame(apply(t(pro_ibaq), 2, function(x) x / max(x) * 100))
          pro_percent <- as.data.frame(t(pro_percent))
          # write.csv(pro_percent, file = glue("./output/{date}_{crl}_protein_{data_type}_percent_levels.csv"))
          
          missing <- setdiff(gns[[i]], rownames(pro_percent))
          # write(missing, file = glue("./output/{date}_{crl}_protein_{data_type}_percent_levels_missing_genes.txt"))
          col_fun = colorRamp2(c(0, 50, 100), c("#EFECEA", "#86A6DF", "#2C4E80"))
          
          
          # svg(filename = glue("./output/{date}_{crl}_protein_{data_type}_percent_levels.svg"), width = 8, height = 8)
          ht <- draw(Heatmap(as.matrix(pro_percent),
                             col = col_fun,
                             column_title = glue("{crl} Protein profiles of {type} human tissues\n
                              N (SR) = {crl_num} of {crl_total}"),
                             cluster_columns = F,
                             cluster_rows = F,
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
          # dev.off()
          
          ht
        }
        
      }, height = 1000)
      
      # }
      # height <- reactiveVal(600)
      # fontSize <- reactiveVal(10)
      # 
      # mTert_reactive <- reactiveVal(NULL)
      # prt.n_reactive <- reactiveVal(NULL)
      # 
      # 
      # output$heatmap <- renderPlot({
      #     
      #     input$generate
      #   
      #     prt.global.AI <- mem_prt.global.AI()
      # 
      # 
      #     selected_geneset <- isolate(input$inputGeneset)
      #     height(1000)
      #     fontSize(5)
      #     geneset_df <-
      #       read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
      # 
      # 
      #     if (selected_geneset == "CRLs") {
      #       geneset_df <-
      #         read_excel(path = "./data/simplifiedheatmap/AllCRLs.v6.xlsx")
      #       height(1000)
      #       fontSize(5)
      # 
      #       if (isolate(input$inputSubGeneset) != "all") {
      #         geneset_df <- geneset_df[geneset_df$E3_Complex == isolate(input$inputSubGeneset), ]
      #         height(800)
      #         fontSize(8)
      #       }
      #     }
      #     if (selected_geneset == "E1") {
      #       geneset_df <-
      #         read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx",
      #                    sheet = 2)
      #       height(300)
      #       fontSize(10)
      #     }
      #     if (selected_geneset == "E2") {
      #       geneset_df <-
      #         read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx",
      #                    sheet = 3)
      #       height(600)
      #       fontSize(10)
      #     }
      #     if (selected_geneset == "E3") {
      #       geneset_df <-
      #         read_excel(path = "./data/simplifiedheatmap/UPS_genes.xlsx",
      #                    sheet = 4)
      #       height(2500)
      #       fontSize(4)
      #     }
      # 
      #     ### update
      #     genes <- geneset_df$hgnc_symbol
      #     geneset <- selected_geneset ### update
      # 
      #     column_title <-
      #       glue("Human {isolate(input$tissueType)} Tissues-{geneset}")
      # 
      #  
      #     prt <- prt.global.AI
      # 
      # 
      #     prt <- prt[rownames(prt) %in% genes, ]
      # 
      # 
      # 
      #     selected_sample_type <- isolate(input$tissueType)
      # 
      #     subset <-
      #       subset(
      #         labels,
      #         Dataset == "CPTAC" &
      #           Sample_type == selected_sample_type,
      #         select = c(Sample_id, Tissue_type)
      #       )
      #     subset <- subset[subset$Sample_id %in% colnames(prt), ]
      # 
      #     prt <- prt[, subset$Sample_id]
      #     prt.n <-
      #       aggregate_by(input = prt, label = subset$Tissue_type)
      # 
      #     print(head(prt.n))
      # 
      #     # Apply function to each column of the matrix
      #     mTert <- apply(prt.n, 2, tertiles)
      #     row.names(mTert) <- row.names(prt.n)
      # 
      #     mTert_reactive(mTert)
      #     prt.n_reactive(prt.n)
      # 
      #     # e3 heatmap has more rows than can be displayed on the screen
      #     ht <- Heatmap(
      #       mTert,
      #       column_title = column_title,
      #       na_col = "grey",
      #       col = my.colors,
      #       row_names_gp = gpar(fontsize = fontSize()),
      #       column_names_gp = gpar(fontsize = 16),
      #       heatmap_legend_param = list(title = "Protein"),
      #       rect_gp = gpar(col = "white"),
      #       column_names_rot = 45,
      #       border = T,
      #     )
      # 
      # 
      #     draw(ht)
      #   },
      #   height = height
      # )
      
      # output$downloadMtert <- downloadHandler(
      #   filename = function() {
      #     paste("mtert-", Sys.Date(), ".csv", sep = "")
      #   },
      #   content = function(file) {
      #     write.csv(mTert_reactive(), file)
      #   }
      # )
      # output$downloadPrtn <- downloadHandler(
      #   filename = function() {
      #     paste("prtn-", Sys.Date(), ".csv", sep = "")
      #   },
      #   content = function(file) {
      #     write.csv(prt.n_reactive(), file)
      #   }
      # )
    })
  }
