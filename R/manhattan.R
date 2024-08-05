manhattanUI <- function(id){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Manhattan Plot"),
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
        label = "Gene",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("datatype"),
        choices = c("rna","pro"),
        multiple = FALSE,
        label = "Data Type",
        options = list(`live-search` = TRUE),
      ),
      actionButton(ns('refresh'), "Generate Plot")
      
    ),
    
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      shinycssloaders::withSpinner(plotOutput(ns("manhattan")),
                                   color = "#6A2798"
      ),
      shinycssloaders::withSpinner(plotOutput(ns("boxplot")),
                                   color = "#6A2798"
      ),
    )
  )
}

manhattanServer <- function(id, dat) {
  moduleServer(id,function(input,output,session){
    
    mem_rdata_from_drive("dat~cptac_v1.1_list_object_mut_cnv_rna_pro.RData")
    
    observe({
      input$datatype
      
      updateSelectizeInput(session, 'gene', 
                           choices = colnames(dat[[input$datatype]]), 
                           server = TRUE, selected = "ABL1.rna")
    })
    
    updateSelectizeInput(session, 'gene', 
                         choices = colnames(dat[["rna"]]), 
                         server = TRUE, selected = "ABL1.rna")
    
    
    updateSelectizeInput(session, 'mutgene', 
                         choices = colnames(dat[["mut"]]), 
                         server = TRUE, selected = "VHL.mut")
    
    output$boxplot = renderPlot({
      ### Pancan boxplots ###
      
      input$refresh
      
      # "dat" object listing: samples, mut, cnv, rna, pro
      
      # Retrieve command-line arguments
      args <- commandArgs(trailingOnly = TRUE)
      
      # Assign command-line arguments to variables
      datatype <- isolate(input$datatype)
      gn <- isolate(input$gene)
      m = isolate(input$mutgene)
      
      colors <- read.csv("./data/colors.csv")
      colors <- setNames(object = colors$color, # colors
                         nm = colors$id) # names
      
      samples <- dat[["samples"]] %>% 
        select(c(harmonized_data_id, cohort))
      
      cohorts <- levels(as.factor(dat$samples$cohort))
      sampls <- samples %>% rownames()
      
      
      E3 <- dat[[datatype]]
      mut <- dat[["mut"]]
      
      
      df <- as.data.frame(E3[row.names(samples), gn, drop=F]) %>%
        mutate(harmonized_data_id = rownames(.),
               cohort = samples[row.names(.), "cohort"]) %>%
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
             title = glue("CPTAC Pancan {datatype}"))
      
    })
    
    output$manhattan = renderPlot({
      # 
      # dat <- list(rna = mem_get_from_drive("./data/CPTAC_Pancan_UPS_genome-wide_eQTLs.rds"),
      #             pro = mem_get_from_drive("./data/CPTAC_Pancan_UPS_genome-wide_pQTLs.rds"))
      # 
      # # filter table for n_mut > n
      # n = 1
      # p_adj = 5e-8
      # 
      # # for (i in names(dat)) { # 2 datatypes
      # 
      # # Manhattan plot: All cohorts together
      # 
      # type = "rna"
      # print("1")
      # # input data type
      # df <- dat[[type]] %>% 
      #   filter(n_mut > n & 
      #            is.numeric(LFC_median) & 
      #            !is.nan(LFC_median) & 
      #            !is.infinite(LFC_median)) %>% 
      #   mutate(Class = as.factor(Class)) %>% 
      #   mutate(Class_ups_gene = paste(Class, ups_gene, sep = "_"))
      # 
      # colors <- setNames(object = brewer.pal(length(levels(df$Class)), "Set1"), # colors
      #                    nm = levels(df$Class)) # names
      # print("2")  
      # sig_vals <- df %>% 
      #   mutate(Label = paste(ups_gene, mutation, sep = "-")) %>% 
      #   # arrange(fdr) %>% # ordered min to max
      #   filter(fdr < p_adj)
      # print("3")  
      # p = ggplot(df, aes(x = Class_ups_gene, y = log10fdr, color = Class)) +
      #   geom_point(alpha = 0.5, size = 1.2) +
      #   theme_prism(border = TRUE) +
      #   scale_color_manual(values = colors, name = "UPS Class") +
      #   labs(x = "UPS genes", y = "-log10(p-adj)", 
      #        title = glue("Manhattan Plot: Quantitative Trait Locus"),
      #        subtitle = glue("10 CPTAC Cohorts: {type} "),
      #        caption = glue("Filtered samples with n_mut > {n}. Dashed line: p.adj = {p_adj}.")) +
      #   geom_hline(yintercept = -log10(p_adj), linetype = "dashed", color = "black") +
      #   ggrepel::geom_text_repel(data = sig_vals, aes(label = Label),
      #                            vjust = -1, size = 3, color = "black") +
      #   theme(axis.text.x = element_blank(),
      #         axis.ticks.x = element_blank())
      # 
      # 
      # print("4")
      
      # return(p)
    })
    
    
  
  })
}
