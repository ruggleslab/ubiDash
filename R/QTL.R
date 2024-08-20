QTLUI <- function(id, cohorts){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Quantitative Trait Locus (QTL) Analysis"),
      p("(description paragraph)"),
      selectizeInput(
        inputId = ns("mutgene"),
        choices = NULL,
        multiple = FALSE,
        label = "Mutant Gene",
        options = list(`live-search` = TRUE),
      ),
      
      selectizeInput(
        inputId = ns("gene"),
        choices = NULL,
        multiple = FALSE,
        label = "UPS Gene",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("datatype"),
        choices = c("rna","pro"),
        multiple = FALSE,
        label = "Data Type",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("cohorts"),
        choices = cohorts,
        selected = cohorts,
        label = "Cohorts",
        options = list(`live-search` = TRUE),
        multiple = TRUE
      ),
      h5("Volcano Plot"),
      selectizeInput(
        inputId = ns("plottype"),
        choices = c("Mutation Centric", "UPS Centric"),
        selected = "UPS Centric", 
        multiple = FALSE,
        label = "Plot Type",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("genelist"),
        choices = NULL,
        selected = "Mutation Centric", 
        multiple = FALSE,
        label = "Gene List",
        options = list(`live-search` = TRUE),
      ),
      textInput(ns("pvalcutoff"), "P value Cutoff:",
                placeholder = "5e-8",value = "5e-8"),
      sliderInput(ns("n_genes"), "Minimum Number of Samples per Mutation Group",
                  min = 0, max = 1000, value = 1, step = 1),
      
      actionButton(ns('refresh'), "Generate Plot"),
   
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotlyOutput(ns("boxplot")),
                                   color = "#6A2798"
      ),
      shinycssloaders::withSpinner(plotlyOutput(ns("volcano")),
                                   color = "#6A2798"
      ),
         p("Download this data"),
      downloadButton(outputId = ns("downloadDat"), label = "Download Data Table"),
      
    )
  )
}

QTLServer <- function(id, genelists, colorsTbl, cohorts) {
  moduleServer(id,function(input,output,session){
    
    mut = mem_get_from_drive("cptac_v1.1_mut.rds")
    
    
    observe({
      input$datatype
      
      if(input$datatype == "pro"){
        updateSelectizeInput(session, 'gene', 
                             choices = cptac_protein_genes, 
                             server = TRUE, selected = "ABL1")
      }
      else{
        updateSelectizeInput(session, 'gene', 
                             choices = cptac_rna_genes, 
                             server = TRUE, selected = "ABL1")
      }
    })
    
    updateSelectizeInput(session, 'gene', 
                         choices = cptac_rna_genes, 
                         server = TRUE, selected = "ABL1")
    
    updateSelectizeInput(session, 'genelist', 
                         choices = names(genelists), 
                         server = TRUE, selected = "hCD")
    
    
    updateSelectizeInput(session, 'mutgene', 
                         choices = colnames(mut), 
                         server = TRUE, selected = "VHL.mut")
    
    output$boxplot = renderPlotly({
      ### Pancan boxplots ###
      
      input$refresh
      
      # Retrieve command-line arguments
      cptac_samples <- as.data.frame(mem_get_from_drive("cptac_samples_info_across_datatypes_v1.1.rds"))

      cptac_samples = cptac_samples[cptac_samples$type == "Tumor",]
      
      cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
      # mut <- mem_get_from_drive("cptac_v1.1_mut.rds")
      
      load("./data/dat~cptac_v1.1_list_object_mut_cnv_rna_pro.RData")
      
      # Assign command-line arguments to variables
      selected_datatype <- isolate(input$datatype)
      gn <- isolate(input$gene)
      m = isolate(input$mutgene)
      
      if(!isTruthy(gn)){
        gn = "ABL1"
      }
      if(!isTruthy(m)){
        m = "VHL.mut"
      }
      
      colors <- colorsTbl
      colors <- setNames(object = colors$color, # colors
                         nm = colors$id) # names
      
      
      selected_cohorts = input$cohorts
      if(!isTruthy(selected_cohorts)){
        selected_cohorts = cohorts
      }
      
      samples = filter(cptac_samples, cohort %in% cohorts & datatype == ifelse(selected_datatype == "rna", "rna", "protein"))

      # print(samples)
      # print(1)

      # samples <- dat[["samples"]] %>% 
      samples <- samples %>% 
        select(c(harmonized_data_id, cohort))
      
      # print(samples)
      # print(2)
      
      cohorts <- levels(as.factor(cptac_samples$cohort))
      
      if(selected_datatype == "pro"){
        E3 <- t(cptac_protein)
      } else{
        E3 <- t(cptac_rna)
      }
      

      samples = samples[!duplicated(samples$harmonized_data_id),]
      rownames(samples) = samples$harmonized_data_id
      sampls <- samples$harmonized_data_id
      sampls = sampls[sampls %in% rownames(E3) & sampls %in% rownames(mut)]
      
      
      # 
      df <- as.data.frame(E3[sampls, gn, drop=F]) %>%
        mutate(harmonized_data_id = rownames(.),
               cohort = samples[rownames(.), "cohort"]) %>%
        left_join(.,
                  data.frame(harmonized_data_id = rownames(mut[sampls, m, drop=F]),
                             mut[sampls, m, drop=F]),
                  join_by("harmonized_data_id"))
      
      # To specify WT / Mutant boxplot order
      df[[m]] <- factor(df[[m]], levels = c("WT", "Mutant"))
      
      
      ggplot(df, aes(x = !!sym(m), y = !!sym(gn))) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = cohort),
                    position = position_jitter(0.2),
                    alpha = 0.6,
                    size = 1.2) +
        scale_color_manual(values = colors[cohorts]) + # defined at startup
        theme_minimal() +
        stat_compare_means(method = "wilcox.test") +  # Add Wilcoxon test comparison
        labs(x = m, y = gn,
             title = glue("CPTAC Pancan {selected_datatype}"))
      
    })


    df = eventReactive(input$refresh, {

      dat = get_from_drive("CPTAC_Pancan_UPS_hCD_hTF_hChR_e_pQTLs.rds")

      i = isolate(input$datatype)

      if(isolate(input$plottype) == "Mutation Centric"){

        gn = gsub(".mut", "", isolate(input$mutgene))
        if(!isTruthy(gn)){
          gn = "VHL"
        }

        df <- dat[[i]] %>%
          filter(n_mut > n &
                   is.numeric(LFC_median) &
                   !is.nan(LFC_median) &
                   !is.infinite(LFC_median) &
                   mutation == gn) %>%
          mutate(Label = paste(ups_gene, mutation, sep = "-")) %>%
          mutate(Class = as.factor(Class))

        return(df)

      } else {


        gn = gsub(".mut", "", isolate(input$gene))
        if(!isTruthy(gn)){
          gn = "ABL1"
        }

        df <-
          filter(dat[[i]],
                 ups_gene == gn) %>%
          mutate(Label = paste(ups_gene, mutation, sep = "-")) %>%
          mutate(Class = as.factor(Class))

        return(df)

      }


    },ignoreNULL = FALSE, ignoreInit=FALSE)


    output$volcano = renderPlotly({

      #  colorsTbl <- read.csv("./data/colors.csv")

      n <- isolate(input$n_genes)
      i = isolate(input$datatype)
      # p_adj = 5e-8
      p_adj = as.numeric(gsub('\\*10\\^', 'e', isolate(input$pvalcutoff)))
      if(!isTruthy(p_adj)){
        p_adj = 5e-8
      }

      if(isolate(input$plottype) == "Mutation Centric"){

        colors <- setNames(object =  colorsTbl[ colorsTbl$id %in% levels(df$Class),]$color, # colors
                           nm = levels(df$Class)) # names


        gn = gsub(".mut", "", isolate(input$mutgene))
        if(!isTruthy(gn)){
          gn = "VHL"
        }

        p = ggplot(df(), aes(x = LFC_median, y = -log10(fdr), label = ups_gene, color = Class)) +
          geom_point(alpha = 0.5, size = 2.5) +
          theme_prism(border = TRUE) +
          scale_color_manual(values = colors, name = "UPS Class") +
          labs(x = glue("log2({gn}-mut/wt)"), y = "-log10(p-adj)",
               title = glue("QTL-{i}: 10 CPTAC Cohorts"),
               # subtitle = glue("Integration of adjusted p-values"),
               caption = glue("Filtered samples with n_mut > {n}\nDashed lines: p.adj = {p_adj}")) +
          geom_hline(yintercept = -log10(p_adj), linetype = "dashed", color = "black") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "black") # +

        return(p)

      } else {

        selected_genelist = genelists[[isolate(input$genelist)]]
        if(!isTruthy(selected_genelist)){
          selected_genelist = genelists[["hCD"]]
        }

        gn = gsub(".mut", "", isolate(input$gene))
        if(!isTruthy(gn)){
          gn = "ABL1"
        }


        p= ggplot(df(), aes(x = LFC_median, y = -log10(fdr), label = mutation, color = mutation %in% selected_genelist)) + #, color = Class)) +
          geom_point(alpha = 0.5, size = 1.2) +
          scale_color_manual(name=glue("{isolate(input$genelist)} Gene"),  values = c("TRUE" = "#44BBFF", "FALSE" = "#55555555"))+
          theme_prism(border = TRUE) +
          # scale_color_manual(values = colors, name = "UPS Class") +
          labs(x = glue("log2(mut/wt)\n{gn}-LFC"), y = "-log10(p-adj)",
               title = glue("QTL-{i}: 10 CPTAC Cohorts"),
               # subtitle = glue("Integration of adjusted p-values"),
               caption = glue("Filtered samples with n_mut > {n}\nDashed lines: p.adj = {p_adj}")) +
          geom_hline(yintercept = -log10(p_adj), linetype = "dashed", color = "black") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "black") # +

        return(p)
      }

    })
    
    output$downloadDat <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        
        write.csv(df(), file)
      }
    )
    
  })
}
