violinUI <- function(id, cohorts) {
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
      
      selectizeInput(inputId = ns("cohorts"),
                     choices = cohorts,
                     selected = cohorts,
                     label = "Select cohorts",
                     options = list(`live-search` = TRUE),  
                     multiple = TRUE),
      
      actionButton(ns("generate"), label = "Generate plot")
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotlyOutput(ns("violin"),  height="1000px"),
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
        
        violin.input = violin.input[violin.input$cohort_code %in% isolate(input$cohorts),]
      
        
        ylab = NULL
        
        if(isolate(input$datatype) == "protein"){
          ylab = paste(isolate(input$datatype), "expression, log2(Relative Intensity)")
        } else {
          ylab = paste(isolate(input$datatype), "expression, log2(TPM + 1)")
        }
        
        violin = ggplot(violin.input, aes(x = Sample_type, y = value, fill = Sample_type)) +
          
          geom_boxplot(
            show.legend = FALSE
            ) + # use geom_boxplot for boxplots format
          geom_jitter(aes(color = Sample_type,
                          text=
                            paste0(
                              "Sample ID: ", sample_id, "\n",
                              "Sample Type: ", Sample_type, "\n",
                              "Cohort Code: ", cohort_code, "\n"
                            )),
                      size = 0.1,
                      show.legend = F,) + # resizing individual dots
          facet_wrap(~variable)+
          theme_prism() +
          
          # re-naming x and y axis
          labs(
            x = NULL,
            y = ylab,
            title = paste("CPTAC", isolate(input$datatype), ", n =", nrow(POIs))
          ) + 
          
          facet_grid(rows = vars(violin.input$cohort_code), cols = vars(violin.input$variable))
        
        ggplotly(violin, tooltip="text")
      })
        
    })
  }
