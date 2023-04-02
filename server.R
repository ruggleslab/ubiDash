# library(shiny)
# library(ComplexHeatmap)
# library(readxl)
# library(tidyverse)
# library(RColorBrewer)
# library(plotly)

packagelist = c("shiny", "ComplexHeatmap", "readxl", "tidyverse", "glue", "stringr", "reshape2", "RColorBrewer", "plotly")

lapply(packagelist, function(x) {
  suppressMessages(require(x, character.only = T,
                           quietly=T, 
                           warn.conflicts = F))
}
)

load("./data/defaults.RData")
# loading default objects below:
# genesSubset
# defaultHeatmap
# cohortsOnco
# mutationsOnco

source("./RScripts/SplitInputGenes.R")

# heatmap variables ----

source('./RScripts/HeatmapFunctions.R')

# ccle dataset
ccle_rna <- readRDS(glue("./data/ccle_tpm_rna_expression.rds"))
ccle_protein <- readRDS("./data/ccle_tmt_protein_expression.rds")
ccle_samples <- readRDS("./data/ccle_sample_info.rds")

# cptac dataset
cptac_rna <- readRDS("./data/cptac_tpm_rna_expression.rds")
cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
cptac_samples <- readRDS("./data/cptac_sample_info.rds")


# cohort and sample type of all normal and tumor samples
# labels <- subset(cptac_samples, sample_type %in% c("T", "N","G"), select = c("sample_type", "cohort_code"))
tnlabels <- subset(cptac_samples, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code"))
# get list of cohorts
cohorts <- unique(tnlabels$cohort_code)


# cptac_protein <- cptac_protein[, ]
filteredMat <- cptac_protein[genesSubset$genes, rownames(tnlabels)]

# remove rows with over 50% NA
filteredMat <- filteredMat[which(rowMeans(!is.na(filteredMat)) > 0.5), ]

colorsTbl = read.csv("./data/colors.csv")
colors = colorsTbl$color[1:11]
colorList <- list(
  Type = setNames(colorsTbl$color[12:13], nm = c("T", "N")),
  Cohort = setNames(colors, nm = cohorts)
)



# volcano plot variables ----

filepath <- "./data/DE_protein.xlsx"
genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
unique_genes <- unique(genes$hgnc_symbol)

# oncoplot variables ----
source('./RScripts/OncoplotColors.R')

# cohortsOnco
maf <- readRDS("./data/oncoplot_maf.rds")

# total patients
nPatients = length(unique(maf$Tumor_Sample_Barcode))

CRLS = genes
crl_genes = unique(CRLS$hgnc_symbol)

# correlation variables ----
plot_tab <- readRDS("./data/cptac_multi_omic_correlations.rds")
source("./RScripts/GenerateCorrPlot.R")
lines = list(geom_hline(yintercept=.5, color = "red") ,
             geom_hline(yintercept=-.5,  color = "red") ,
             geom_vline(xintercept=-.5,color = "red") ,
             geom_vline(xintercept=--.5, color = "red") ,
             geom_hline(yintercept=0, size = 1, color = "black") ,
             geom_vline(xintercept=-0, size = 1, color = "black"))



# boxplot variables ----
  
  # Boxplot for RNA values
  # labeling instructions
  p <- c("GMCL1", "USP28", "TP53BP1", "APC","RB1", "TP53", "MYC", "CCND1")
  POIs <- ccle_rna[, p, drop = F]
  POIs$DepMap_ID <- rownames(POIs)
  
  lab <- data.frame(genes = colnames(POIs[1:8]),
                    driver = c(rep("Unknown", 3), # "GMCL1"   "USP28"   "TP53BP1"
                               rep("Tumor Suppressor", 3), # "APC"     "RB1"  "TP53" 
                               rep("Oncogene", 2))) # "MYC"     "CCND1"  
  
  # constructing violin.input table for ggplot in the long format
  violin.input <- reshape2::melt(POIs, id.vars="DepMap_ID") %>% 
    mutate(driver = rep("Unknown", nrow(.))) 
  
  violin.input$driver <- ifelse(violin.input$variable %in% c("MYC","CCND1"),
                         "Oncogene", violin.input$driver) 
  
  violin.input$driver <- ifelse(violin.input$variable %in% c("APC","RB1","TP53"),
                         "Tumor Suppressor", violin.input$driver)


  # server ----

shinyServer(function(input, output, session) {

  
  # oncoplot tab ----
  # reactive vals are able to be shared across both oncoplots
  ctbl = reactiveVal(NULL)
  mtbl = reactiveVal(NULL)
  renderedGenes = reactiveVal(NULL)
  mutationsreactive = reactiveVal(NULL)
  anno = reactiveVal(NULL)
  modifiedOncoInput = reactiveVal(FALSE)
  
  
  updateSelectizeInput(session, 'inputOncoGenesTXT', 
                       choices = as.list(rownames(cptac_protein)), 
                       server = TRUE)
  
  observe({
    input$generateOnco
  
    # see if the settings are at their defaults
    modifiedOncoInput(!(
                      is.null(isolate(input$inputOncoGenesTXT)) &
                      is.null(isolate(input$inputOncoGenesCSV)) &
                      isolate(input$maxOncoGenes) == 10 &
                      length(isolate(input$oncoCohortSelect)) == 11
                      ))
    
    if(modifiedOncoInput()){
    
      # inputGenes = isolate(input$inputOncoGenesTXT)
      #inputGenes <- getInputGenes(isolate(input$inputOncoGenesCSV),isolate(input$inputOncoGenesTXT), delim = isolate(input$oncoDelims))
      inputGenes <- getInputGenes(isolate(input$inputOncoGenesCSV),isolate(input$inputOncoGenesTXT), delim = isolate(input$oncoDelims))

      if(length(inputGenes) == 0){
        inputGenes = crl_genes
      }
      
      # maf with only most frequently mutated genes
      smallFM = filter(maf, toupper(Hugo_Symbol) %in% inputGenes & cohort_code %in% isolate(input$oncoCohortSelect))
      
      uniqueGenes = unique(smallFM$Hugo_Symbol)
      
      mutationsTbl = matrix(ncol = length(unique(smallFM$Tumor_Sample_Barcode)), nrow = length(uniqueGenes))
      rownames(mutationsTbl) = uniqueGenes
      colnames(mutationsTbl) = unique(smallFM$Tumor_Sample_Barcode)
      cohortsTbl = mutationsTbl
      
      for(i in seq_along(rownames(smallFM))){
        id = smallFM[i,Tumor_Sample_Barcode]
        # create tables to be plotted
        mutationsTbl[smallFM[i,Hugo_Symbol],id] = levels(smallFM$Variant_Classification)[smallFM[i,Variant_Classification]]
        cohortsTbl[smallFM[i,Hugo_Symbol],id] = smallFM[i,cohort_code]
        
      }
      
      # set NA values to empty string
      mutationsTbl[is.na(mutationsTbl)] <- ""
      cohortsTbl[is.na(cohortsTbl)] <- ""
      
      cohortsTbl = as.data.frame(cohortsTbl)
      mutationsTbl = as.data.frame(mutationsTbl)
      
      # sort both tables the same way
      sortedTbl = as.data.frame(cohortsTbl) %>% arrange(rowSums(cohortsTbl == ""))
      topGenes = NULL
      if(nrow(sortedTbl) > isolate(input$maxOncoGenes)){
        topGenes = rownames(sortedTbl[1:isolate(input$maxOncoGenes),])
      }else {
        topGenes = rownames(sortedTbl)
      }
  
      cohortsTbl = filter(cohortsTbl, rownames(cohortsTbl) %in% topGenes)
      mutationsTbl = filter(mutationsTbl, rownames(mutationsTbl) %in% topGenes)
      
      # get percents
      pcts = matrix(ncol = 2, nrow = nrow(cohortsTbl))
      rownames(pcts) = rownames(cohortsTbl)
      pcts[,1] = rowSums(cohortsTbl != "")
      pcts[,2] = round((pcts[,1] / 1069) * 100, digits = 1)
      
      if(nrow(pcts) > 1)
        pcts = pcts[order(pcts[,2],decreasing=TRUE),]
      
      pcts[,2] = paste0(pcts[,2], "%")
      mutationsTbl = mutationsTbl[order(pcts[,2],decreasing=TRUE),]
      anno(pcts[,2])
      
      mutationsTbl = mutationsTbl[rownames(pcts),]
      cohortsTbl = cohortsTbl[rownames(pcts),]
      
      mutations = levels(smallFM$Variant_Classification)
      mutationsreactive(mutations)
      ctbl(cohortsTbl)
      mtbl(mutationsTbl)
      renderedGenes(topGenes)
      
    }
    
  })
  # mutations oncoplot
  output$mutationOncoplot = renderPlot({
    input$generateOnco
    
    if(!modifiedOncoInput()) return(mutationsOnco)
  
    mutation_colors = get_mutation_colors(mutationsreactive())
    mutation_alter_fun = get_mutation_alter_fun(mutation_colors)
    
    
    tryCatch(
      {
        output$oncoError = renderText("")
        
        oncoPrint(
        mtbl(), alter_fun = mutation_alter_fun, col = mutation_colors, 
        row_names_side = "left", 
        show_pct = FALSE,
        
        row_order = renderedGenes(),
        
        remove_empty_columns = TRUE,
        
        heatmap_legend_param = list(title = "Alterations"),

        right_annotation = rowAnnotation(pct = anno_text(anno()), barplot = anno_oncoprint_barplot())
      )
        },
      error = function(e){
        output$oncoError = renderText("No Data")
        return()
      }
    )
  })
  # cohorts oncoplot
  output$cohortOncoplot = renderPlot({
    input$generateOnco
    
    if(!modifiedOncoInput()) return(cohortsOnco)
    
    cohort_colors = get_cohort_colors()
    cohort_alter_fun = get_cohort_alter_fun(cohort_colors)
    
    tryCatch({
    oncoPrint(
      ctbl(), alter_fun = cohort_alter_fun, col = cohort_colors, 
      row_names_side = "left", 
      show_pct = FALSE,
      
      row_order = renderedGenes(),
      
      remove_empty_columns = TRUE,
      
      heatmap_legend_param = list(title = "Cohorts"),
      
      right_annotation = rowAnnotation(pct = anno_text(anno()), barplot = anno_oncoprint_barplot())
      
      
    )},
    error = function(e){
      return()
    }
    )
  })
  
  output$downloadOncoMutations <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(mtbl(), file)
    }
  )
  output$downloadOncoCohorts <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(ctbl(), file)
    }
  )
  

  # heatmap tab ----
  
  updateSelectizeInput(session, 'inputGenesTXT', 
                       choices = as.list(rownames(cptac_protein)), 
                       server = TRUE)

  # add download heatmap button
  
  # default plot height maybe automate this 
  plotHeight <- reactiveVal(956)

  # matrix with only genes to plot
  getGenesMat <- function(inputGenes,cohortSubsetLabels) {
    if (length(inputGenes) == 0 &
        length(isolate(isolate(input$cohortSelect))) == length(cohorts)) {
      return(NULL)
    }
    
    if (length(inputGenes) == 0) {
      print(cohortSubsetLabels)
      return(filteredMat[, rownames(cohortSubsetLabels)])
    }
    return(cptac_protein[inputGenes, rownames(cohortSubsetLabels)])
  }
  
  
  output$heatmap <- renderPlot(
    {
      
      # Button will trigger rerender
      input$generateHeatmap
      
      # get cohorts selected by user
      cohortSubsetLabels <- getCohortPatients(tnlabels, isolate(input$cohortSelect))
     
      if(nrow(cohortSubsetLabels) == 0){
        output$rowscut <- renderText({
          paste("Error: Please select at least one cohort")
        })
        return()
      }
      
      inputGenes <- getInputGenes(isolate(input$inputGenesCSV),isolate(input$inputGenesTXT), delim = isolate(input$heatmapDelims))
      inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels)
      # print(inputGenesMat)

      if (is.null(inputGenesMat)) {

        if (isolate(input$cluster)) {
          return(defaultHeatmap)
        }
        inputGenesMat <- filteredMat
      } else {
        inputGenesMat <- as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
        plotHeight(nrow(inputGenesMat) * 12 + 200)
      }

      cutout <- inputGenes[!(inputGenes %in% toupper(rownames(inputGenesMat)))]
      if (length(cutout) > 0 & (is.null(isolate(input$inputGenesCSV)) | length(inputGenesMat) > 0)) {
        output$rowscut <- renderText({
          paste("No Data:", paste(cutout, collapse = ", "))
        })
      } else if (length(cutout) > 0) {
        output$rowscut <- renderText({
          paste("No Data:", paste(cutout, collapse = ", "), ". Make sure delimiter is correct.")
        })
      } else {
        output$rowscut <- renderText({
          ""
        })
      }

      ha <- HeatmapAnnotation(Type = cohortSubsetLabels$sample_type, Cohort = cohortSubsetLabels$cohort_code, col = colorList)
      hm <- NULL
      
      colsplit = NULL
      if(isolate(input$separateNandT)){
        colsplit = getCohortPatients(tnlabels,isolate(input$cohortSelect))$sample_type
      }
      

      inputGenesMat <- as.matrix(inputGenesMat)
      
      tryCatch(
        {
          output$err <- renderText({
            ""
          })
          return(drawHeatmap(inputGenesMat, ha, colsplit, isolate(input$cluster)))
        },
        error = function(e) {
          output$err <- renderText({
            "Warning: Not enough data to cluster columns"
          })
          return(drawHeatmap(inputGenesMat, ha, colsplit, FALSE))
          
        }
      )
    },
    height = plotHeight
  )

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      cohortSubsetLabels <- getCohortPatients(tnlabels, isolate(input$cohortSelect))
      inputGenes <- getInputGenes(isolate(input$inputGenesCSV),isolate(input$inputGenesTXT), delim = isolate(input$volcanoDelims))
      inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels)
      
      if (length(inputGenesMat) == 0) {
        inputGenesMat <- filteredMat
      }
      write.csv(inputGenesMat, file)
    }
  )

  # output$downloadHeatmap <- downloadHandler(
  #   filename = function() {
  #     paste("heatmap-", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
  #     ggsave(file, plot = plotInput(), device = device)
  #   }
  # )
  
  # volcano plot tab ----
  updateSelectizeInput(session, 'inputVolcanoGenesTXT', 
                       choices = as.list(rownames(cptac_protein)), 
                       server = TRUE)
  
  output$distPlot <- renderPlotly({
    # allow user to save plots

    input$generateVolcano
    
    inputGenes <- getInputGenes(isolate(input$inputVolcanoGenesCSV),isolate(input$inputVolcanoGenesTXT), delim = isolate(input$volcanoDelims))
    if(length(inputGenes) == 0){
      inputGenes = c("")
    }
    # cutoffs
    logFCCutoff <- isolate(input$logFCCutoff)
    PValCutoff <- -log10(as.double(isolate(input$PValCutoff)))

    # prepare dataframe
    data <- read_excel(path = filepath, sheet = isolate(input$cohort))
    data$POI <- "Other"
    data$POI[data$gene %in% unique_genes] = "CRL"
    data$POI[toupper(data$gene) %in% inputGenes] = "Potential Substrate"
    data$log.scal.pval <- -log10(data$sca.adj.pval)

    # get significant genes
    labeledPoints <- filter(data, (gene %in% unique_genes & (abs(logFC) > logFCCutoff & abs(log.scal.pval) > PValCutoff)) | toupper(gene) %in% inputGenes)

    p <- ggplot(data, aes(x = logFC, y = log.scal.pval, color = POI)) +
      ylab("-Log10(adjusted p-val)") +
      xlab("Log2(Tumor/Normal)") +
      labs(title = isolate(input$cohort)) +
      geom_point(shape = 16, alpha = .5, mapping = aes(text = paste(
        "Gene: ", gene, "\n",
        "POI: ", POI, "\n",
        "LogFC: ", logFC, "\n",
        "P Value: ", log.scal.pval, "\n",
        sep = ""
      ))) +
      scale_color_manual(breaks = c("CRL", "Potential Substrate", "Other"), values = c("#FF0000","#0000FF","#888888")) +
      geom_vline(xintercept = logFCCutoff, linetype = "solid", color = "red") +
      geom_vline(xintercept = -logFCCutoff, linetype = "solid", color = "red") +
      geom_hline(yintercept = PValCutoff, linetype = "solid", color = "red")

    if (nrow(labeledPoints) > 0) {
      p <- p + geom_text(
        hjust = 1, vjust = 1, data = labeledPoints,
        mapping = aes(label = gene, fontface = "bold"),
        color = "#000000"
      )
    }

    ggplotly(p, tooltip = "text") %>% layout(dragmode = FALSE)
  })

  output$downloadVolcanoData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(read_excel(path = filepath, sheet = input$cohort), file)
    }
  )
  # correlation tab ----
  
  warning = reactiveVal("")
  
  # proteome tumor - proteome normal
  output$pt_pn = renderPlotly({
    input$generateCorrelation
    output$corrErr1 <- renderText({
      warning()
    })
    tryCatch(
    {
      p = generateCorrPlot(isolate(input$CorrelationCutoff), plot_tab, "pr_tumor.cor", "pr_normal.cor", 
                           isolate(input$inputCorrGenesTXT), warning)
      ggplotly(p, tooltip = "text")
    },
      error = function(e) {
        print(e)
        output$corrErr1 <- renderText({
          "Error: No data available for input genes"
      })
        return()
    })
    
  })
  # proteome tunor - transcriptome tunor
  output$pt_rt = renderPlotly({
    input$generateCorrelation
    output$corrErr2 <- renderText({
      warning()
    })
    tryCatch(
      {
        p = generateCorrPlot(isolate(input$CorrelationCutoff), plot_tab, "pr_tumor.cor","rn_tumor.cor", 
                             isolate(input$inputCorrGenesTXT), warning)
        ggplotly(p, tooltip = "text")
      },
      error = function(e) {
        output$corrErr2 <- renderText({
          "Error: No data available for input genes"
        })
        return()
      })
    
  })
  
  # transcriptome tumor - transcriptome normal
  output$rt_rn = renderPlotly({
    input$generateCorrelation
    output$corrErr3 <- renderText({
      warning()
    })
    tryCatch(
      {
        p = generateCorrPlot(isolate(input$CorrelationCutoff), plot_tab, "rn_tumor.cor", "rn_normal.cor", 
                             isolate(input$inputCorrGenesTXT), warning)
        ggplotly(p, tooltip = "text")
      },
      error = function(e) {
        output$corrErr3 <- renderText({
          "Error: No data available for input genes"
        })
        return()
      })
    
  })

  output$downloadCorrelationData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(plot_tab, file)
    }
  )
  # boxplot tab ----
  output$violin = renderPlotly({
  violin = ggplot(violin.input, aes(x = variable, y = value, fill = driver)) + 
    
    geom_violin(show.legend = T) + # use geom_boxplot for boxplots format
    geom_jitter(size = 0.1, show.legend = F) + # resizing individual dots
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") + # setting colors
    
    # re-naming x and y axis
    labs(x = NULL, y = "mRNA expression, log2(TPM +1)", 
         title = paste("CCLE cell lines, n =", nrow(POIs))
    )
  ggplotly(violin)
  })
  
  # protein specific correlation ----
  source(file = "./RScripts/corrAnalysis.R")
  
  #Gene vs Gene correlation RNA
  updateSelectizeInput(session, 'inputGvgCorrRnaNameA', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
  updateSelectizeInput(session, 'inputGvgCorrRnaNameB', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
  
  resGvgCorrRna <- eventReactive(input$generateGvgCorrRna,{singleCorr(geneNameA = input$inputGvgCorrRnaNameA,
                                                                      geneNameB = input$inputGvgCorrRnaNameB,
                                                                      cohort = input$inputGvgCorrRnaCohort,
                                                                      tissueType = input$inputGvgCorrRnaStype, 
                                                                      cptac_protein = NULL,
                                                                      cptac_rna = cptac_rna,
                                                                      cptac_samples = cptac_samples)})
  
  output$outputPlotGvgCorrRna <- renderPlot({resGvgCorrRna()$singCorrPlot})
  output$outputTabGvgCorrRna <- renderTable({resGvgCorrRna()$singCorrTab})#,rownames = T, colnames = F)
  
  #Protein vs Protein correlation RNA
  updateSelectizeInput(session, 'inputGvgCorrProteinNameA', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
  updateSelectizeInput(session, 'inputGvgCorrProteinNameB', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
  
  resGvgCorrProtein <- eventReactive(input$generateGvgCorrProtein,{singleCorr(geneNameA = input$inputGvgCorrProteinNameA,
                                                                      geneNameB = input$inputGvgCorrProteinNameB,
                                                                      cohort = input$inputGvgCorrProteinCohort,
                                                                      tissueType = input$inputGvgCorrProteinStype, 
                                                                      cptac_protein = cptac_protein,
                                                                      cptac_rna = NULL, 
                                                                      cptac_samples = cptac_samples)})
  
  output$outputPlotGvgCorrProtein <- renderPlot({resGvgCorrProtein()$singCorrPlot})
  output$outputTabGvgCorrProtein <- renderTable({resGvgCorrProtein()$singCorrTab})#,rownames = T, colnames = F)
  
  #Protein vs Gene correlation
  updateSelectizeInput(session, 'inputGvgCorrBothNameA', choices = as.list(intersect(rownames(cptac_protein),rownames(cptac_rna))), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type Protein name'))
  updateSelectizeInput(session, 'inputGvgCorrBothNameB', choices = as.list(intersect(rownames(cptac_protein),rownames(cptac_rna))), server = TRUE, selected = "A2M", options = list(placeholder = 'Type gene name'))
  
  resGvgCorrBoth <- eventReactive(input$generateGvgCorrBoth,{singleCorr(geneNameA = input$inputGvgCorrBothNameA,
                                                                        geneNameB = input$inputGvgCorrBothNameB,
                                                                        cohort = input$inputGvgCorrBothCohort,
                                                                        tissueType = input$inputGvgCorrBothStype, 
                                                                        cptac_protein = cptac_protein,
                                                                        cptac_rna = cptac_rna,
                                                                        cptac_samples = cptac_samples)})
  
  output$outputPlotGvgCorrBoth <- renderPlot({resGvgCorrBoth()$singCorrPlot})
  output$outputTabGvgCorrBoth <- renderTable({resGvgCorrBoth()$singCorrTab})#,rownames = T, colnames = F)
  
  #Gene vs All correlation for proteins
  updateSelectizeInput(session, 'inputGvACorrProteinGeneName', choices = as.list(rownames(cptac_protein)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
  resGvACorrProtein <- eventReactive(input$generateGvACorrProtein,{multiCorr(geneName = input$inputGvACorrProteinGeneName,
                                                                 cutOff =  input$inputGvACorrProteinR, 
                                                                 cohort = input$inputGvACorrProteinCohort,
                                                                 tissueType =input$inputGvACorrProteinType,
                                                                 cptac_protein = cptac_protein,
                                                                 cptac_rna = NULL,
                                                                 cptac_samples = cptac_samples
                                                                 )
  })
  
  output$outputTabGvACorrProtein <- DT::renderDataTable(server = FALSE,{resGvACorrProtein()$result},rownames = F,extensions = 'Buttons', options = list(
    dom = 'frtBip',#Bfrtip '<<t>B>' B - buttons, f - search panel r - processing display element t - the table itself i - number of items p - paging
    buttons = list(
      list(extend = "excel", text = "Download results", filename = paste0(input$inputGvACorrRnaGeneName, "_ProteinCorrelations"),
           exportOptions = list(
             modifier = list(page = "all")
           )))))
  
  output$outputPlotGvACorrProtein <- renderPlot({resGvACorrProtein()$corrPlot})
  
  #Gene vs All correlation for RNA
  updateSelectizeInput(session, 'inputGvACorrRnaGeneName', choices = as.list(rownames(cptac_rna)), server = TRUE, selected = "TNNC1", options = list(placeholder = 'Type gene name'))
  resGvACorrRna <- eventReactive(input$generateGvACorrRna,{multiCorr(geneName = input$inputGvACorrRnaGeneName,
                                                                       cutOff =  input$inputGvACorrRnaR,
                                                                       tissueType =input$inputGvACorrRnaType, 
                                                                       cohort = input$inputGvACorrRnaCohort,
                                                                       cptac_protein = NULL,
                                                                       cptac_rna = cptac_rna,
                                                                       cptac_samples = cptac_samples
  )
  })
  
  output$outputTabGvACorrRna <- DT::renderDataTable(server = FALSE,{resGvACorrRna()$result},rownames = F,extensions = 'Buttons', options = list(
    dom = 'frtBip',#Bfrtip '<<t>B>' B - buttons, f - search panel r - processing display element t - the table itself i - number of items p - paging
    buttons = list(
      list(extend = "excel", text = "Download results", filename = paste0(input$inputGvACorrRnaGeneName, "_RNAcorrelations"),
           exportOptions = list(
             modifier = list(page = "all")
           )))))
  
  output$outputPlotGvACorrRna <- renderPlot({resGvACorrRna()$corrPlot})
  
  #Gene panel correlation for RNA
  updateSelectizeInput(session, 'inputGvpCorrRnaGeneName', choices = as.list(rownames(cptac_rna)), server = TRUE,  options = list(placeholder = 'Type gene name'))
  
  resGvpCorrRna <- eventReactive(input$generateGvpCorrRna,{panelCorr(geneNamesRNA = input$inputGvpCorrRnaGeneName,
                                                                      geneNamesProtein = NULL, 
                                                                      tissueType = input$inputGvpCorrRnaStype, 
                                                                      cohort = input$inputGvpCorrRnaCohort, 
                                                                      cptac_protein = NULL, 
                                                                      cptac_rna = cptac_rna, 
                                                                      cptac_samples = cptac_samples)})
  
  output$outputPlotGvpCorrRna <- renderPlot({resGvpCorrRna()$corrPlot})
  
  #Gene panel correlation for Protein
  updateSelectizeInput(session, 'inputGvpCorrProteinGeneName', choices = as.list(rownames(cptac_protein)), server = TRUE, options = list(placeholder = 'Type gene name'))
  
  resGvpCorrProtein <- eventReactive(input$generateGvpCorrProtein,{panelCorr(geneNamesProtein = input$inputGvpCorrProteinGeneName,
                                                                     geneNamesRNA = NULL, 
                                                                     tissueType = input$inputGvpCorrProteinStype, 
                                                                     cohort = input$inputGvpCorrProteinCohort, 
                                                                     cptac_protein = cptac_protein, 
                                                                     cptac_rna = NULL, 
                                                                     cptac_samples = cptac_samples)})
  
  output$outputPlotGvpCorrProtein <- renderPlot({resGvpCorrProtein()$corrPlot})
  
  #Gene panel correlation for Protein/RNA
  updateSelectizeInput(session, 'inputGvpCorrBothGeneNameA', choices = as.list(rownames(cptac_rna)), server = TRUE, options = list(placeholder = 'Type gene name'))
  updateSelectizeInput(session, 'inputGvpCorrBothGeneNameB', choices = as.list(rownames(cptac_protein)), server = TRUE, options = list(placeholder = 'Type protein name'))
  
  resGvpCorrBoth <- eventReactive(input$generateGvpCorrBoth,{panelCorr(geneNamesRNA = input$inputGvpCorrBothGeneNameA,
                                                                     geneNamesProtein = input$inputGvpCorrBothGeneNameB, 
                                                                     tissueType = input$inputGvpCorrBothStype, 
                                                                     cohort = input$inputGvpCorrBothCohort, 
                                                                     cptac_protein = cptac_protein, 
                                                                     cptac_rna = cptac_rna, 
                                                                     cptac_samples = cptac_samples)})
  
  output$outputPlotGvpCorrBoth <- renderPlot({resGvpCorrBoth()$corrPlot})
  
  
  
})
