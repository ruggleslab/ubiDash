violinUI <- function(id, cohorts) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Descriptive Statistics"),
      tabsetPanel(
        id = ns("datatype"),
        selected = "mRNA",
        tabPanel(
          "mRNA",
          selectizeInput(
            inputId = ns("mRNAGenes"),
            choices = NULL,
            multiple = T,
            selected = c("KBTBD11"),
            label = "Genes",
            options = list(`live-search` = TRUE),
          ),
        ),
        tabPanel(
          "protein",
          selectizeInput(
            inputId = ns("proteinGenes"),
            choices = NULL,
            multiple = T,
            selected = c("KBTBD11"),
            label = "Genes",
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
          label = "Select cohorts",
          options = list(`live-search` = TRUE),
          multiple = TRUE
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
      actionButton(ns("generate"), label = "Generate plot")
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(
        plotlyOutput(ns("violin"),
          height =
            "1000px"
        ),
        color = "#6A2798"
      ),
    )
  )
}


violinServer <-
  function(id,
           # cptac_rna,
          #  cptac_protein,
           # cptac_samples,
           cohort_colors) {
    moduleServer(id, function(input, output, session) {
      

      updateSelectizeInput(
        session,
        "mRNAGenes",
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
        "proteinGenes",
        choices = cptac_protein_genes,
        server = TRUE,
        selected = "KBTBD11"
      )

      observe({
        input$dataset

        if (input$dataset == "CCLE") {
          updateSelectizeInput(
            session,
            "mRNAGenes",
            choices = ccle_rna_genes,
            server = TRUE,
            selected = "KBTBD11"
          )

          updateSelectizeInput(
            session,
            "proteinGenes",
            choices = ccle_protein_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
        } else if (input$dataset == "CPTAC") {
          updateSelectizeInput(
            session,
            "mRNAGenes",
            choices = cptac_rna_genes,
            server = TRUE,
            selected = "KBTBD11"
          )

          updateSelectizeInput(
            session,
            "proteinGenes",
            choices = cptac_protein_genes,
            server = TRUE,
            selected = "KBTBD11"
          )
        }
      })


      output$violin <- renderPlotly({
        input$generate
        
        validate(
          need(isolate(input$datatype) != "protein" || length(isolate(input$proteinGenes)) >= 1, 'Please select at least one gene.'),
          need(isolate(input$datatype) != "mRNA" || length(isolate(input$mRNAGenes)) >= 1, 'Please select at least one gene.'),
        )
       

        p <- NULL
        dataset <- NULL
        samples <- NULL

        if (isolate(input$dataset) == "CPTAC") {
          if (isolate(input$datatype) == "protein") {
           
             cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
             dataset <- t(cptac_protein)
             
             
             p <- isolate(input$proteinGenes)
          } else {
            # there appear to be duplicate columns
            cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
            dataset <- t(cptac_rna)
            # dataset <- t(cptac_rna[, !duplicated(colnames(cptac_rna))])
            p <- isolate(input$mRNAGenes)
          }
        } else if (isolate(input$dataset) == "CCLE") {
          if (isolate(input$datatype) == "protein") {
            ccle_protein <- mem_get_from_drive("ccle_tmt_protein_expression.rds")
            dataset <- t(ccle_protein)
            p <- isolate(input$proteinGenes)
          } else {
            # there appear to be duplicate columns
            ccle_rna <- mem_get_from_drive("ccle_tpm_rna_expression.rds")
            dataset <- ccle_rna
            p <- isolate(input$mRNAGenes)
          }
        }

       
        if (is.null(dataset)) {
          cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
          dataset <- t(cptac_protein)
        }


        POIs <- dataset[, p, drop = F]
        POIs <- as.data.frame(POIs)
        POIs$sample_id <- rownames(POIs)


        violin.input <-
          reshape2::melt(POIs, id.vars = "sample_id")
        # %>%
        #   mutate(driver = rep("Unknown", nrow(.)))

        if (isolate(input$dataset) == "CPTAC") {
          violin.input <- merge(
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

          violin.input <- violin.input[!is.na(violin.input$Sample_type), ]
          if(length(isolate(input$cohorts)) > 0){
            violin.input <- violin.input[violin.input$cohort_code %in% isolate(input$cohorts), ]
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



        ylab <- NULL

        if (isolate(input$datatype) == "protein") {
          ylab <- paste(
            isolate(input$datatype),
            "expression, log2(Relative Intensity)"
          )
        } else {
          ylab <- paste(isolate(input$datatype), "expression, log2(TPM + 1)")
        }

        violin <- NULL


        if (isolate(input$dataset) == "CPTAC") {
          violin <- ggplot(
            violin.input,
            aes(
              x = cohort_code,
              y = value,
              fill =
                Sample_type
            )
          ) +
            geom_boxplot(
              show.legend = FALSE,
              position = position_dodge()
            )
        } else if (isolate(input$dataset) == "CCLE") {
          violin <- ggplot(
            violin.input,
            aes(x = Subtype, y = value)
          ) +
            geom_boxplot(
              show.legend = FALSE,
              position = position_dodge(),
              fill = "#FFBD42"
            )
        }




        if (isolate(input$dataset) == "CPTAC") {
          violin <- violin + geom_point(
            position = (position_jitterdodge(dodge.width = 0.9)),
            alpha = .5,
            size = .75,
            stroke = .1
          ) +
            scale_fill_manual(values = c(
              "Tumor" = "#FFBD42",
              "Normal" = "#BA5CF4"
            ))
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
          )
        }


        violin <- violin +
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
          )) +
          theme(axis.text.x = element_text(size = 8))



        ggplotly(violin, tooltip = "text") %>%
          # https://stackoverflow.com/questions/58586314/keeping-dodge-position-in-boxplot-passed-to-plotly
          layout(boxmode = "group")
      })
    })
  }
