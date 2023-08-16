violinUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Descriptive Statistics"),
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
          conditionalPanel(
            condition = paste0("input['", ns("dataset"), "'] != 'PRIDE'"),
            selectizeInput(
              inputId = ns("mRNAGenes"),
              choices = NULL,
              multiple = T,
              selected = c("BRWD3"),
              label =  "Input genes",
              options = list(`live-search` = TRUE),
            )
          ),
          conditionalPanel(condition = paste0(
            "input['", ns("dataset"), "'] == 'PRIDE'"
          ),
          p("no data"))
          
        )
      ),
      selectizeInput(
        inputId = ns("dataset"),
        choices = c("CPTAC", "CCLE", "PRIDE"),
        multiple = F,
        selected = "CPTAC",
        label =  "Dataset",
        options = list(`live-search` = TRUE),
      ),
      
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
        selectizeInput(
          inputId = ns("colorby"),
          choices = c("Cancer Type", "Sample Type"),
          multiple = F,
          selected = "Sample Type",
          label =  "Color by",
          options = list(`live-search` = TRUE),
        ),
        
        selectizeInput(
          inputId = ns("cohorts"),
          choices = cohorts,
          selected = cohorts,
          label = "Select cohorts",
          options = list(`live-search` = TRUE),
          multiple = TRUE
        ),
      ),
      conditionalPanel(
        condition = paste0("!(input['", ns("dataset"), "'] == 'PRIDE' && input['", ns("datatype"), "'] == 'mRNA')"),
        actionButton(ns("generate"), label = "Generate plot")
      )
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = 'error'),
      shinycssloaders::withSpinner(plotlyOutput(ns("violin"),  height =
                                                  "1000px"),
                                   color = "#6A2798"),
    )
  )
}


violinServer <-
  function(id,
           cptac_rna,
           cptac_protein,
           cptac_samples,
           cohort_colors,
           pride_protein) {
    moduleServer(id, function(input, output, session) {
      updateSelectizeInput(
        session,
        'mRNAGenes',
        choices = rownames(cptac_rna),
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
      
      observe({
        input$dataset
        
        if (input$dataset == "CCLE") {
          updateSelectizeInput(
            session,
            'mRNAGenes',
            choices = colnames(ccle_rna),
            server = TRUE,
            selected = "TSPAN6"
          )
          
          updateSelectizeInput(
            session,
            'proteinGenes',
            choices = rownames(ccle_protein),
            server = TRUE,
            selected = "A1CF"
          )
        }
        else if (input$dataset == "CPTAC") {
          updateSelectizeInput(
            session,
            'mRNAGenes',
            choices = rownames(cptac_rna),
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
        }
        else if (input$dataset == "PRIDE") {
          updateSelectizeInput(
            session,
            'proteinGenes',
            choices = rownames(pride_protein),
            server = TRUE,
            selected = "A1BG"
          )
        }
        
      })
      
      
      output$violin = renderPlotly({
        input$generate
        
        
        p = NULL
        dataset = NULL
        samples = NULL
        
        if (isolate(input$dataset) == "CPTAC") {
          if (isolate(input$datatype) == "protein") {
            dataset = t(cptac_protein)
            p <- isolate(input$proteinGenes)
          } else{
            # there appear to be duplicate columns
            dataset = t(cptac_rna[, !duplicated(colnames(cptac_rna))])
            p <- isolate(input$mRNAGenes)
          }
        } else if (isolate(input$dataset) == "CCLE") {
          if (isolate(input$datatype) == "protein") {
            dataset = t(ccle_protein)
            p <- isolate(input$proteinGenes)
          } else{
            # there appear to be duplicate columns
            dataset = ccle_rna
            p <- isolate(input$mRNAGenes)
          }
          
        } else if (isolate(input$dataset) == "PRIDE") {
          dataset = t(pride_protein)
          p <- isolate(input$proteinGenes)
          
        }
        
        if (is.null(p)) {
          p = "VHL"
        }
        if (is.null(dataset)) {
          dataset =  t(cptac_protein)
        }
        
        
        POIs <- dataset[, p, drop = F]
        POIs = as.data.frame(POIs)
        POIs$sample_id <- rownames(POIs)
        
        
        violin.input <-
          reshape2::melt(POIs, id.vars = "sample_id")
        # %>%
        #   mutate(driver = rep("Unknown", nrow(.)))
        
        if (isolate(input$dataset) == "CPTAC") {
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
          
          violin.input = violin.input[!is.na(violin.input$Sample_type),]
          violin.input = violin.input[violin.input$cohort_code %in% isolate(input$cohorts),]
          
        } else if (isolate(input$dataset == "CCLE")) {
          print(head(violin.input))
          
          violin.input = merge(
            x = violin.input,
            y = ccle_samples,
            by.x = "sample_id",
            by.y =  (if (isolate(input$datatype) == "protein")
              "row.names"
              else
                "parent_depmap_id")
            # all.x = TRUE
          )[c("Subtype",
              "variable",
              "value")]
        }
        
        
        
        ylab = NULL
        
        if (isolate(input$datatype) == "protein") {
          ylab = paste(isolate(input$datatype),
                       "expression, log2(Relative Intensity)")
        } else {
          ylab = paste(isolate(input$datatype), "expression, log2(TPM + 1)")
        }
        
        violin = NULL
        
        
        if (isolate(input$dataset) == "CPTAC") {
          violin = ggplot(violin.input,
                          aes(
                            x = cohort_code,
                            y = value,
                            fill =
                              Sample_type
                          ))
          
        } else if (isolate(input$dataset) == "CCLE") {
          violin = ggplot(violin.input,
                          aes(x =  Subtype, y = value), fill = "cyan")
        }
        
        
        if (isolate(input$dataset) != "PRIDE") {
          violin = violin +
            geom_boxplot(show.legend = FALSE,
                         position = position_dodge())  # use geom_boxplot for boxplots format
          
          
          if (isolate(input$dataset) == "CPTAC") {
            violin = violin + geom_point(
              position = (position_jitterdodge(dodge.width = 0.9)),
              alpha = .5,
              size = .75,
              stroke = .1
            )
          } else if (isolate(input$dataset) == "CCLE") {
            violin = violin + geom_point(
              aes(
                text = paste0("Subtype: ",
                              Subtype,
                              "\n",
                              "Value: ",
                              value,
                              "\n")
              ),
              position = (position_jitter()),
              alpha = .5,
              size = .75,
              stroke = .1
            )
          }
          
          
          violin = violin +
            theme_prism() +
            
            # re-naming x and y axis
            labs(
              x = NULL,
              y = ylab,
              title = paste(
                isolate(input$dataset),
                isolate(input$datatype),
                ", n =",
                nrow(POIs)
              )
            ) +
            facet_grid(rows = vars(violin.input$variable)) +
            theme(axis.text.x = element_text(
              angle = 60,
              vjust = 0.5,
              hjust = 1,
              size = 5
            ))
          
        } else {
          violin = ggplot(violin.input,
                          aes(x =  0, y = value), fill = "cyan") +
            geom_boxplot(show.legend = FALSE,
                         position = position_dodge()) + # use geom_boxplot for boxplots format
            
            geom_point(
              aes(text =
                    paste0("Type: ",
                           sample_id)),
              position = (position_jitter())
              ,
              alpha = .5,
              size = .75,
              stroke = .1
              
            ) +
            theme_prism() +
            labs(title = paste("PRIDE n =",
                               nrow(POIs))) +
            facet_grid(rows = vars(violin.input$variable))
        }
        
        ggplotly(violin, tooltip = "text") %>%
          # https://stackoverflow.com/questions/58586314/keeping-dodge-position-in-boxplot-passed-to-plotly
          layout(boxmode = "group")
      })
      
    })
  }
