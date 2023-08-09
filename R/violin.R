violinUI <- function(id, cohortSheets) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Violin Plot"),
      
      selectizeInput(
        inputId = ns("genes"),
        choices = NULL,
        multiple = T,
        selected = c("VHL"),
        label =  "Input genes",
        options = list(`live-search` = TRUE),
      ),
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotlyOutput(ns("violin")),
                                   color = "#6A2798"),
    )
  )
}


violinServer <- function(id, ccle_rna, cptac_protein,cptac_samples) {
  moduleServer(id, function(input, output, session) {
    
    updateSelectizeInput(
      session,
      'genes',
      choices = as.list(maf@gene.summary$Hugo_Symbol),
      server = TRUE,
      selected = "VHL"
    )
    
    
    output$violin = renderPlotly({
      dataset = t(cptac_protein)
      
      p <- input$genes
      
      
      POIs <- dataset[, p, drop = F]
      POIs = as.data.frame(POIs)
      POIs$sample_id <- rownames(POIs)
      
      # lab <- data.frame(genes = colnames(POIs[1:2]))
      
      
      violin.input <- reshape2::melt(POIs, id.vars = "sample_id") %>%
        mutate(driver = rep("Unknown", nrow(.)))
      
      violin.input = merge(
        x = violin.input,
        y = cptac_samples,
        by.x = "sample_id",
        by.y = "harmonized_data_id"
      )[c(
        "sample_id",
        "variable",
        "value",
        "cohort_code",
        "Sample_type",
        "Cancer_subtype"
      )]
      
      
      violin = ggplot(violin.input, aes(x = variable, y = value)) +
        
        
        geom_violin(show.legend = T) + # use geom_boxplot for boxplots format
        geom_jitter(aes(color = cohort_code),
                    size = 0.1,
                    show.legend = F,
        ) + # resizing individual dots
        theme_minimal() +
        scale_fill_brewer(palette = "Set1") + # setting colors
        
        # re-naming x and y axis
        labs(
          x = NULL,
          y = "mRNA expression, log2(TPM +1)",
          title = paste("CPTAC Proteins, n =", nrow(POIs))
        )
      ggplotly(violin)
    })
    
  })
}
