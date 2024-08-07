QTLUI <- function(id){
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Quantitative Trait Locus (QTL) Analysis"),
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
      # shinycssloaders::withSpinner(plotOutput(ns("manhattan")),
      #                              color = "#6A2798"
      # ),
      shinycssloaders::withSpinner(plotlyOutput(ns("boxplot")),
                                   color = "#6A2798"
      ),
    )
  )
}

QTLServer <- function(id) {
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
    
    
    updateSelectizeInput(session, 'mutgene', 
                         choices = colnames(mut), 
                         server = TRUE, selected = "VHL.mut")
    
    output$boxplot = renderPlotly({
      ### Pancan boxplots ###
      
      input$refresh
      
      # Retrieve command-line arguments
      cptac_samples <- as.data.frame(mem_get_from_drive("cptac_sample_info_with_subtypes.rds"))
      cptac_samples = cptac_samples[cptac_samples$sample_type == "T",]
      
      cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
      cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
      mut <- mem_get_from_drive("cptac_v1.1_mut.rds")
      
      load("./data/dat~cptac_v1.1_list_object_mut_cnv_rna_pro.RData")
      
      # Assign command-line arguments to variables
      datatype <- isolate(input$datatype)
      gn <- isolate(input$gene)
      m = isolate(input$mutgene)
      
      if(!isTruthy(gn)){
        gn = "ABL1"
      }
      if(!isTruthy(m)){
        m = "VHL.mut"
      }
      
      colors <- read.csv("./data/colors.csv")
      colors <- setNames(object = colors$color, # colors
                         nm = colors$id) # names
      
      # samples <- dat[["samples"]] %>% 
      samples <- cptac_samples %>% 
        select(c(harmonized_data_id, cohort_code))
      
      cohorts <- levels(as.factor(cptac_samples$cohort))
      
      if(input$datatype == "pro"){
        E3 <- t(cptac_protein)
      } else{
        E3 <- t(cptac_rna)
      }
      
      
      sampls <- samples %>% rownames()
      sampls = sampls[sampls %in% rownames(E3) & sampls %in% rownames(mut)]
      
      # print(rownames(samples))
      # print(E3)
      # 
      # print(rownames(samples) %in% colnames(E3))
      # 
      # E3[rownames(samples),]
      # sampls
      # E3[sampls, gn]
      # 
      # rownames(E3)
      
      
      # print(E3[gn, row.names(samples), drop=F])
      # 
      # rownames(samples) %in% E3
      # 
      # print(E3[row.names(samples),gn])
      
      # mut[sampls, m, drop=F]
      
      # 
      # gn = "A1BG"
      # m = "VHL.mut"
      
      # print(gn)
      # print(m)
      # 
      df <- as.data.frame(E3[sampls, gn, drop=F]) %>%
        mutate(harmonized_data_id = rownames(.),
               cohort = samples[row.names(.), "cohort_code"]) %>%
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
    
    # output$manhattan = renderPlot({
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
    # })
    
    
  
  })
}
