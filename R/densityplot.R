densityUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Density Plot"),
      selectizeInput(
        inputId = ns("gene"),
        choices = NULL,
        multiple = FALSE,
        label = "Gene",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("datatype"),
        choices = c("rna", "protein"),
        multiple = FALSE,
        label = "Data Type",
        options = list(`live-search` = TRUE),
      ),
      actionButton(ns('refresh'), "Generate Plot"),
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("plot"), ), color = "#6A2798"),
      p("Download this data"),
      downloadButton(outputId = ns("download"), label = "Download")
    )
  )
}


densityServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    updateSelectizeInput(
      session,
      'gene',
      choices = maf_genes,
      server = TRUE,
      selected = "KBTBD11"
    )
    
    sample = eventReactive(input$refresh, {
      
      validate(
        need(isolate(input$gene), 'Please select a gene.'),
        need(isolate(input$datatype), 'Please select a data type'),
      )
      
      
      ups_gn <- isolate(input$gene)
      data_type <- isolate(input$datatype)
      xlab <- ifelse(data_type == "rna", "Log2(TPM+1)", "Log2(Intensity)")
      
      # data loading
      # prot <- readRDS(file = "../../ubiDash/data/cptac_tmt-log2-AI-v1.1_protein_expression.rds")
      # sample <- readRDS("../../ubiDash/data/cptac_samples_info_across_datatypes_v1.1.rds")
      
      dataset = NULL
      if (data_type == 'rna') {
        dataset <- as.matrix(mem_get_from_drive("duplicates_averaged_cptac_rna.rds"))
      } else{
        dataset <- mem_get_from_drive("cptac_tmt-log2-AI_protein_expression.rds")
      }
      
      sample <- readRDS("./data/cptac_samples_info_across_datatypes_v1.1.rds")
      
      # parameters
      sample <- subset(
        sample,
        datatype == data_type &
          harmonized_data_id %in% colnames(dataset) &
          include == "TRUE"
      ) %>% # samples that passed QC
        select(c(harmonized_data_id, cohort, type, cancer_subtype)) %>%
        unique()
      
      sample[[ups_gn]] <- as.numeric(dataset[ups_gn, sample$harmonized_data_id])
      
      # reorder cohorts based on the median expression of the gene of interest
      sample <- sample %>%
        mutate(cohort = fct_reorder(
          cohort,!!sym(ups_gn),
          .fun = median,
          .na_rm = TRUE
        ))
      
      return(sample)
      
    },ignoreNULL = FALSE, ignoreInit=FALSE)
    
    output$plot = renderPlot({
      # NOTES on data source
      # sample_v1.1 <- openxlsx::read.xlsx(xlsxFile = "../EXP169/output/cptac_samples_info_across_datatypes_v1.1.xlsx")
      # saveRDS(sample_v1.1, file = "../../ubiDash/data/cptac_samples_info_across_datatypes_v1.1.rds")
      # input$refresh
      # variables
      ups_gn <- isolate(input$gene)
      data_type <- isolate(input$datatype)
      xlab <- ifelse(data_type == "rna", "Log2(TPM+1)", "Log2(Intensity)")
      
      # calculating  sample sizes and making plot annotations
      annotations <- sample() %>%
        group_by(cohort, type) %>%
        summarise(n = n(), # calculating sample sizes
                  label_x = median(!!sym(ups_gn), na.rm = TRUE)) %>% # calculating medians
        ungroup() %>%
        mutate(label = paste0("N = ", n))
      
      # Define the x-axis range for the plot
      x_range <- range(sample()[[ups_gn]], na.rm = TRUE)
      
      # Ensure label positions are within the x-axis range
      annotations <- annotations %>%
        mutate(label_x = pmax(pmin(label_x, x_range[2]), x_range[1]))
      
      
      # visualization
      ggplot(sample()) +
        geom_density_ridges(
          aes(
            x = !!sym(ups_gn),
            y = cohort,
            fill = type
          ),
          alpha = 0.5,
          bandwidth = 0.3,
          scale = 1,
          quantiles = 2,
          quantile_lines = T
        ) +
        
        geom_text_repel(
          data = annotations,
          # geom_text(data = annotations,
          aes(
            x = label_x,
            y = cohort,
            label = label,
            color = type
          ),
          # position = position_nudge(y = 0.8)) + # use with geom_text()
          nudge_y = 0.7,
          segment.color = NA
        ) + # use with geom_text_repel
        
        labs(
          x = glue("{ups_gn}\n{xlab}"),
          y = NULL,
          title = glue::glue("{ups_gn} Global {data_type} expression")
        ) +
        theme_prism(border = T) +
        scale_fill_manual(values = c("Normal" = "orange", "Tumor" = "purple")) +
        scale_color_manual(values = c("Normal" = "orange", "Tumor" = "purple"),
                           guide = F)
    }, height=600)
    
    output$download <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(sample())
        write.csv(sample(), file)
      }
    )
    
  })
}