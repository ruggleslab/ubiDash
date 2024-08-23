violinUI <- function(id, cohorts) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Descriptive Statistics"),
      p("(description paragraph)"),
      
      tabsetPanel(
        id = ns("datatype"),
        selected = "mRNA",
        tabPanel(
          "mRNA",
          selectizeInput(
            inputId = ns("mRNAGene"),
            choices =  c("KBTBD11"),
            multiple = F,
            selected = c("KBTBD11"),
            label = "Gene",
            options = list(`live-search` = TRUE),
          ),
        ),
        tabPanel(
          "protein",
          selectizeInput(
            inputId = ns("proteinGene"),
            choices =  c("KBTBD11"),
            multiple = F,
            selected = c("KBTBD11"),
            label = "Gene",
            options = list(`live-search` = TRUE),
          ),
        ),
      ),
      
      selectizeInput(
        inputId = ns("dataset"),
        choices = c("CPTAC", "CCLE"),
        multiple = F,
        selected = "CPTAC",
        label = "Dataset",
        options = list(`live-search` = TRUE),
      ),
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
        selectizeInput(
          inputId = ns("cohorts"),
          choices = cohorts,
          selected = cohorts,
          label = "Cohorts",
          options = list(`live-search` = TRUE),
          multiple = TRUE
        ),
      ),
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] == 'CPTAC'"),
        checkboxInput(ns("independentY"), "Independent Y Axes",
                      value = FALSE
        ),
      ),
      conditionalPanel(
        condition = paste0("input['", ns("dataset"), "'] == 'CCLE'"),
        selectizeInput(
          inputId = ns("subtype"),
          choices = c(),
          label = "Select subtype",
          options = list(`live-search` = TRUE),
          multiple = TRUE
        ),
      ),
      actionButton(ns("generate"), label = "Generate plot"),
      
      
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(
        plotOutput(ns("boxplot"),
                     height =
                       "1000px"
        ),
        color = "#6A2798"
      ),
      # uiOutput(ns("plots_ui")),
      p("Download this data"),
      downloadButton(outputId = ns("download"), label = "Download")
    )
  )
}


violinServer <-
  function(id, colorsTbl
           # cptac_rna,
           #  cptac_protein,
           # cptac_samples,
  ) {
    moduleServer(id, function(input, output, session) {
      ns <- session$ns  # Access the namespace function
      
      updateSelectizeInput(
        session,
        "mRNAGene",
        choices = cptac_rna_genes,
        server = TRUE,
        selected = "KBTBD11"
      )
      
      updateSelectizeInput(
        session,
        "subtype",
        choices = ccle_samples_subtypes,
        server = TRUE,
      )
      
      updateSelectizeInput(
        session,
        "proteinGene",
        choices = cptac_protein_genes,
        server = TRUE,
        selected = "KBTBD11"
      )
      
      observe({
        input$dataset
        
        if (input$dataset == "CCLE") {
          updateSelectizeInput(
            session,
            "mRNAGene",
            choices = ccle_rna_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
          
          updateSelectizeInput(
            session,
            "proteinGene",
            choices = ccle_protein_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
        } else if (input$dataset == "CPTAC") {
          updateSelectizeInput(
            session,
            "mRNAGene",
            choices = cptac_rna_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
          
          updateSelectizeInput(
            session,
            "proteinGene",
            choices = cptac_protein_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
        }
      })
      
      POIs_reactive = reactiveVal(NULL)
      
      violin.input = eventReactive(input$generate, {
        p <- NULL
        dataset <- NULL
        samples <- NULL
        
        
        
        if (isolate(input$dataset) == "CPTAC") {
          if (isolate(input$datatype) == "protein") {
            
            cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
            dataset <- t(cptac_protein)
            
            
            p <- isolate(input$proteinGene)
          } else {
            # there appear to be duplicate columns
            cptac_rna <- mem_get_from_drive("cptac_transformed_rna_expression_v1-1.rds")
            dataset <- t(cptac_rna)
            # dataset <- t(cptac_rna[, !duplicated(colnames(cptac_rna))])
            p <- isolate(input$mRNAGene)
          }
        } else if (isolate(input$dataset) == "CCLE") {
          if (isolate(input$datatype) == "protein") {
            ccle_protein <- mem_get_from_drive("ccle_tmt_protein_expression.rds")
            dataset <- t(ccle_protein)
            p <- isolate(input$proteinGene)
          } else {
            # there appear to be duplicate columns
            ccle_rna <- mem_get_from_drive("ccle_tpm_rna_expression.rds")
            dataset <- ccle_rna
            p <- isolate(input$mRNAGene)
          }
        }
        
        
        # if (is.null(dataset)) {
        #   cptac_protein <- mem_get_from_drive("cptac_tmt-log2-AI-v1.1_protein_expression.rds")
        #   dataset <- t(cptac_protein)
        # }
        
        
        dataset = dataset[!duplicated(rownames(dataset)),]
        POIs <- dataset[, p, drop = F]
        POIs <- as.data.frame(POIs)
        POIs$sample_id <- rownames(POIs)
        
        POIs_reactive(POIs)
        
        
        
        # print(POIs$sample_id)
        
        
        violin.input <-
          reshape2::melt(POIs, id.vars = "sample_id")
        
        # %>%
        #   mutate(driver = rep("Unknown", nrow(.)))
        
        
        if (isolate(input$dataset) == "CPTAC") {
          cptac_samples = mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds")
          cptac_samples = filter(cptac_samples, datatype == ifelse(isolate(input$datatype) == "mRNA", "rna", "protein"))
          cptac_samples = filter(cptac_samples, type == "Normal" | type == "Tumor")
          
          
          
          violin.input <- merge(
            x = violin.input,
            y = cptac_samples,
            by.x = "sample_id",
            by.y = "harmonized_data_id"
          )[c(
            "sample_id",
            "variable",
            "value",
            "cohort",
            "type",
            "cancer_subtype"
          )]
          
          
          
          violin.input <- violin.input[!is.na(violin.input$type), ]
          if(length(isolate(input$cohorts)) > 0){
            violin.input <- violin.input[violin.input$cohort %in% isolate(input$cohorts), ]
          }
          
          
        } else if (isolate(input$dataset == "CCLE")) {
          
          ccle_samples <- mem_get_from_drive("ccle_sample_info.rds")
          
          violin.input <- merge(
            x = violin.input,
            y = ccle_samples,
            by.x = "sample_id",
            by.y = (if (isolate(input$datatype) == "protein") {
              "row.names"
            } else {
              "parent_depmap_id"
            })
            # all.x = TRUE
          )[c(
            "Subtype",
            "variable",
            "value"
          )]
          if (length(isolate(input$subtype)) != 0) {
            violin.input <- violin.input[violin.input$Subtype %in% isolate(input$subtype), ]
          }
          
          violin.input$Sample_type <- "Tumor"
        }
        
        return(violin.input)
        
      },ignoreNULL = FALSE, ignoreInit=FALSE)
      
      colors  = colorsTbl$color
      names(colors) = colorsTbl$id
      
     
      
    
      
      output$boxplot = renderPlot({
        
        gene = "KBTBD11"
        
        
        if(isolate(input$datatype) == "protein"){
          gene = isolate(input$proteinGene)
        } else{
          gene = isolate(input$mRNAGene)
        }
        
        if(!isTruthy(gene)){
          gene = "KBTBD11"
        }
        
        ylab <- NULL
        
        if (isolate(input$datatype) == "protein") {
          ylab <- paste(
            isolate(input$datatype),
            "expression, log2(Intensity)"
          )
        } else {
          ylab <- paste(isolate(input$datatype), "expression, log2(TPM + 1)")
        }
        
        violin <- NULL
        
        sub_violin.input = filter(violin.input(), variable == gene)
        
        
        if (isolate(input$dataset) == "CPTAC") {
          violin <- ggplot(
            sub_violin.input,
            aes(
              x = type,
              y = value,
              fill =
                type
            )
            
          ) +
            stat_compare_means(method = "wilcox.test", size = 3) +  # Add Wilcoxon test comparison
            
            geom_boxplot(
              show.legend = FALSE,
              position = position_dodge()
            ) 
        } else if (isolate(input$dataset) == "CCLE") {
          violin <- ggplot(
            sub_violin.input,
            aes(x = Subtype, y = value)
          ) +
            # stat_compare_means(method = "wilcox.test") +  # Add Wilcoxon test comparison
            
            geom_boxplot(
              show.legend = FALSE,
              position = position_dodge(),
              fill = colors[["Tumor"]]
            ) 
          
        }
        
        # print(colors)
        
        if (isolate(input$dataset) == "CPTAC") {
          violin <- violin + geom_point(
            position = (position_jitterdodge(dodge.width = 0.9)),
            alpha = .5,
            size = 2,
            stroke = 1,
            mapping = aes(fill = type),
            pch = 21,
            color = "black"
          ) +
            theme_minimal() +
            theme(
              # strip.text.x = element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            scale_fill_manual(values = colors) +
            scale_color_manual(values = colors) +
            facet_wrap(facets = vars(sub_violin.input$cohort), scales = ifelse(isolate(input$independentY), "free", "free_x"), axes = "all_x") 
          
            # ggh4x::facet_grid2(cols = vars(sub_violin.input$cohort), rows = NULL, scales = ifelse(isolate(input$independentY), "free", "fixed"), independent = ifelse(isolate(input$independentY), "y", "none")) 
        } else if (isolate(input$dataset) == "CCLE") {
          violin <- violin + geom_point(
            aes(text = paste0(
              "Subtype: ",
              Subtype,
              "\n",
              "Value: ",
              value,
              "\n"
            )),
            position = (position_jitter()),
            alpha = .5,
            size = .75,
            stroke = .1
          ) +
            theme_minimal() +
            
            theme(axis.text.x = element_text(size = 8, angle = 90))
          
          # +
          #   facet_grid(cols = vars(sub_violin.input$Subtype), rows = NULL) 
        }
        
        violin <- violin +
          # theme_prism() +
          
          # re-naming x and y axis
          labs(
            x = NULL,
            y = ylab,
            title = paste(
              isolate(input$dataset),
              isolate(input$datatype),
              ", n =",
              nrow(sub_violin.input)
            )
          ) 
          # theme(axis.text.x = element_text(
          #   angle = 60,
          #   vjust = 0.5,
          #   hjust = 1,
          #   size = 5
          # )) +
        
        violin
        
        # return(ggplotly(violin, tooltip = "text")
        #        %>%
        #          # https://stackoverflow.com/questions/58586314/keeping-dodge-position-in-boxplot-passed-to-plotly
        #          layout(boxmode = "group", height = 800))
        # # }))
      }, height = 800)
      
      output$download <- downloadHandler(
        filename = function() {
          if (isolate(input$datatype) == "protein") {
            paste("descriptive-statistics-", isolate(input$dataset), "-protein-",isolate(input$proteinGenes), "-", Sys.Date(), ".csv", sep = "")
          } else {
            paste("descriptive-statistics-", isolate(input$dataset), "-mRNA-",isolate(input$mRNAGenes), "-", Sys.Date(), ".csv", sep = "")
          }
        },
        content = function(file) {
          req(violin.input())
          write.csv(violin.input(), file)
        }
      )
      
      
    })
    
    
  }
