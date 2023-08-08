violinUI <- function(id, cohortSheets){
  ns <- NS(id)
  
  shinycssloaders::withSpinner(plotlyOutput(ns("violin")), color = "#6A2798")
}


violinServer<- function(id, ccle_rna, violin.input, POIS) {
  moduleServer(id,function(input,output,session){
    
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
