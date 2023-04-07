violinUI <- function(id, cohortSheets){
  ns <- NS(id)
  
  shinycssloaders::withSpinner(plotlyOutput(ns("violin")), color = "#6A2798")
}


violinServer<- function(id) {
  moduleServer(id,function(input,output,session){
    
    ccle_rna <- readRDS(glue("./data/ccle_tpm_rna_expression.rds"))
    
    # Boxplot for RNA values
    # labeling instructions
    p <- c("GMCL1", "USP28", "TP53BP1", "APC","RB1", "TP53", "MYC", "CCND1")
    POIs <- ccle_rna[, p, drop = F]
    POIs$DepMap_ID <- rownames(POIs)
    
    lab <- data.frame(genes = colnames(POIs[1:8]),
                      driver = c(rep("Unknown", 3), # "GMCL1"   "USP28"   "TP53BP1"
                                 rep("Tumor Suppressor", 3), # "APC"     "RB1"  "TP53" 
                                 rep("Oncogene", 2))) # "MYC"     "CCND1"  
    
    # constructing violin.input table for ggplot in the long format
    violin.input <- reshape2::melt(POIs, id.vars="DepMap_ID") %>% 
      mutate(driver = rep("Unknown", nrow(.))) 
    
    violin.input$driver <- ifelse(violin.input$variable %in% c("MYC","CCND1"),
                                  "Oncogene", violin.input$driver) 
    
    violin.input$driver <- ifelse(violin.input$variable %in% c("APC","RB1","TP53"),
                                  "Tumor Suppressor", violin.input$driver)
    
    
    output$violin = renderPlotly({
      violin = ggplot(violin.input, aes(x = variable, y = value, fill = driver)) + 
        
        geom_violin(show.legend = T) + # use geom_boxplot for boxplots format
        geom_jitter(size = 0.1, show.legend = F) + # resizing individual dots
        theme_minimal() +
        scale_fill_brewer(palette = "Set1") + # setting colors
        
        # re-naming x and y axis
        labs(x = NULL, y = "mRNA expression, log2(TPM +1)", 
             title = paste("CCLE cell lines, n =", nrow(POIs))
        )
      ggplotly(violin)
    })
    
  })
}
