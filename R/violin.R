violinUI <- function(id, cohortSheets) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Violin Plot"),
      tabsetPanel(
        id = ns("datatype"),
        selected = "protein",
        
        tabPanel(
          "protein",
          selectizeInput(
            inputId = ns("proteinGenes"),
            choices = NULL,
            multiple = T,
            selected = c("VHL"),
            label =  "Input genes",
            options = list(`live-search` = TRUE),
          ),
          
        ),
        
        tabPanel(
          "mRNA",
          selectizeInput(
            inputId = ns("mRNAGenes"),
            choices = NULL,
            multiple = T,
            selected = c("BRWD3"),
            label =  "Input genes",
            options = list(`live-search` = TRUE),
          ),
          
        )
      ),
      selectizeInput(
        inputId = ns("dataset"),
        choices = c("CPTAC"),
        multiple = F,
        selected = "CPTAC",
        label =  "Dataset",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("colorby"),
        choices = c("Cancer Type", "Sample Type"),
        multiple = F,
        selected = "Sample Type",
        label =  "Color by",
        options = list(`live-search` = TRUE),
      ),
      actionButton(ns("generate"), label = "Generate plot")
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotlyOutput(ns("violin")),
                                   color = "#6A2798"),
    )
  )
}


violinServer <-
  function(id,
           cptac_rna,
           cptac_protein,
           cptac_samples,
           cohort_colors) {
    moduleServer(id, function(input, output, session) {
      updateSelectizeInput(
        session,
        'mRNAGenes',
        choices = rownames(cptac_protein),
        server = TRUE,
        selected = "BRWD3"
      )
      
      updateSelectizeInput(
        session,
        'proteinGenes',
        choices = rownames(cptac_protein),
        server = TRUE,
        selected = "VHL"
      )
      
      output$violin = renderPlotly({
        input$generate
        
        
        p = NULL
        dataset = NULL
        
        
        if (isolate(input$datatype) == "protein") {
          dataset = t(cptac_protein)
          p <- isolate(input$proteinGenes)
        } else{
          # there appear to be duplicate columns
          dataset = t(cptac_rna[, !duplicated(colnames(cptac_rna))])
          p <- isolate(input$mRNAGenes)
          
        }
        
        
        POIs <- dataset[, p, drop = F]
        POIs = as.data.frame(POIs)
        POIs$sample_id <- rownames(POIs)
        
        # lab <- data.frame(genes = colnames(POIs[1:2]))
        
        
        violin.input <-
          reshape2::melt(POIs, id.vars = "sample_id") %>%
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
          geom_jitter(aes(color = if(isolate(input$colorby) == "Cancer Type") cohort_code else Sample_type),
                      size = 0.1,
                      show.legend = F,) + # resizing individual dots
          theme_minimal() +
          
          # re-naming x and y axis
          labs(
            x = NULL,
            y = "mRNA expression, log2(TPM +1)",
            title = paste("CPTAC", isolate(input$datatype), ", n =", nrow(POIs))
          )
        ggplotly(violin)
      })
      
    })
  }
