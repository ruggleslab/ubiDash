library(shiny)
library(ComplexHeatmap)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(plotly)

# TODO ----
# add cohort select to oncoplot
# fix download button on heatmap
# add download data button ton oncoplot
# Fix issue when no cohorts are selected


# heatmap variables ----

source('./RScripts/HeatmapFunctions.R')

load("./data/genesSubset.Rdata")
load("./data/20220214_ANA044_Proteome_Pancan_mat.combo_observed_norm.Rdata")
load("./data/sample_info.Rdata")
load("./data/defaultHeatmap.RData")


labels <- subset(samples, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code"))
cohorts <- unique(labels$cohort_code)

mat <- as.data.frame(mat)
largeMat <- mat[, rownames(labels)]
filteredMat <- largeMat[genesSubset$genes, ]

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

load("./data/sample_info.Rdata")
load("./data/maf.Rdata")
nPatients = length(unique(maf$Tumor_Sample_Barcode))

CRLS = genes
crl_genes = unique(CRLS$hgnc_symbol)

# correlation variables ----
plot_tab = read.csv("./data/plotTabv4.csv")
source("./RScripts/GenerateCorrPlot.R")
lines = list(geom_hline(yintercept=.5, color = "red") ,
             geom_hline(yintercept=-.5,  color = "red") ,
             geom_vline(xintercept=-.5,color = "red") ,
             geom_vline(xintercept=--.5, color = "red") ,
             geom_hline(yintercept=0, size = 1, color = "black") ,
             geom_vline(xintercept=-0, size = 1, color = "black"))



# server ----

shinyServer(function(input, output) {
  
  
  # oncoplot tab ----
  ctbl = reactiveVal(NULL)
  mtbl = reactiveVal(NULL)
  renderedGenes = reactiveVal(NULL)
  mutationsreactive = reactiveVal(NULL)
  anno = reactiveVal(NULL)
  
  observe({
    input$generateOnco
    
    inputGenes = isolate(input$inputOncoGenesTXT)
    if(inputGenes == ""){
      inputGenes = crl_genes
    } else {
      inputGenes = gsub(",", " ", inputGenes)
      inputGenes = gsub("  ", " ", inputGenes)
      inputGenes = unlist(strsplit(inputGenes, " "))
    }
    
    smallFM = filter(maf, Hugo_Symbol %in% inputGenes & cohort_code %in% isolate(input$oncoCohortSelect))
    
    uniqueGenes = unique(smallFM$Hugo_Symbol)
    
    mutationsTbl = matrix(ncol = length(unique(smallFM$Tumor_Sample_Barcode)), nrow = length(uniqueGenes))
    rownames(mutationsTbl) = uniqueGenes
    colnames(mutationsTbl) = unique(smallFM$Tumor_Sample_Barcode)
    cohortsTbl = mutationsTbl
    
    for(i in seq_along(rownames(smallFM))){
      id = smallFM[i,Tumor_Sample_Barcode]
      mutationsTbl[smallFM[i,Hugo_Symbol],id] = levels(smallFM$Variant_Classification)[smallFM[i,Variant_Classification]]
      cohortsTbl[smallFM[i,Hugo_Symbol],id] = smallFM[i,cohort_code]
      
    }
    
    mutationsTbl[is.na(mutationsTbl)] <- ""
    cohortsTbl[is.na(cohortsTbl)] <- ""
    
    cohortsTbl = as.data.frame(cohortsTbl)
    mutationsTbl = as.data.frame(mutationsTbl)
    
    sortedTbl = as.data.frame(cohortsTbl) %>% arrange(rowSums(cohortsTbl == ""))
    topGenes = NULL
    if(nrow(sortedTbl) > isolate(input$maxOncoGenes)){
      topGenes = rownames(sortedTbl[1:isolate(input$maxOncoGenes),])
    }else {
      topGenes = rownames(sortedTbl)
    }

    cohortsTbl = filter(cohortsTbl, rownames(cohortsTbl) %in% topGenes)
    mutationsTbl = filter(mutationsTbl, rownames(mutationsTbl) %in% topGenes)
    
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
    
    
    
  })
  output$mutationOncoplot = renderPlot({
    input$generateOnco
    
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

        right_annotation = rowAnnotation(pct = anno_text(anno()))
      )
        },
      error = function(e){
        output$oncoError = renderText("No Data")
        return()
      }
    )
  })
  output$cohortOncoplot = renderPlot({
    input$generateOnco
    
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
      
      right_annotation = rowAnnotation(pct = anno_text(anno()))
      
      
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

  # add download heatmap button
  
  plotHeight <- reactiveVal(956)
  
  getGenesMat <- function(inputGenes,cohortSubsetLabels) {
    if (length(inputGenes) == 0 &
        length(isolate(isolate(input$cohortSelect))) == length(cohorts)) {
      return(NULL)
    }
    
    if (length(inputGenes) == 0) {
      return(filteredMat[, rownames(cohortSubsetLabels)])
    }
    return(largeMat[inputGenes, rownames(cohortSubsetLabels)])
  }

  output$heatmap <- renderPlot(
    {
      
      # Button will trigger rerender
      input$generateHeatmap
      
      
      cohortSubsetLabels <- getCohortPatients(labels, isolate(input$cohortSelect))
     
      if(nrow(cohortSubsetLabels) == 0){
        output$rowscut <- renderText({
          paste("Error: Please select at least one cohort")
        })
        return()
      }
      
      inputGenes <- getInputGenes(input$inputGenesCSV,isolate(input$inputGenesTXT))
      inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels)

      if (is.null(inputGenesMat)) {

        if (isolate(input$cluster)) {
          return(defaultHeatmap)
        }
        inputGenesMat <- filteredMat
      } else {
        inputGenesMat <- as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
        plotHeight(nrow(inputGenesMat) * 12 + 200)
      }

      cutout <- inputGenes[!(inputGenes %in% rownames(inputGenesMat))]
      if (length(cutout) > 0) {
        output$rowscut <- renderText({
          paste("No Data:", paste(cutout, collapse = ", "))
        })
      } else {
        output$rowscut <- renderText({
          ""
        })
      }

      ha <- HeatmapAnnotation(Type = cohortSubsetLabels$sample_type, Cohort = cohortSubsetLabels$cohort_code, col = colorList)
      hm <- NULL
      colsplit = getCohortPatients(labels,isolate(input$cohortSelect))$sample_type

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
      cohortSubsetLabels <- getCohortPatients(labels, isolate(input$cohortSelect))
      inputGenes <- getInputGenes(input$inputGenesCSV,isolate(input$inputGenesTXT))
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
  output$distPlot <- renderPlotly({
    # allow user to save plots

    input$generateVolcano

    logFCCutoff <- isolate(input$logFCCutoff)
    PValCutoff <- -log10(as.double(isolate(input$PValCutoff)))

    data <- read_excel(path = filepath, sheet = isolate(input$cohort))
    data$CRL <- data$gene %in% unique_genes
    data$log.scal.pval <- -log10(data$sca.adj.pval)

    labeledPoints <- filter(data, gene %in% unique_genes & (abs(logFC) > logFCCutoff & abs(log.scal.pval) > PValCutoff))

    p <- ggplot(data, aes(x = logFC, y = log.scal.pval, color = CRL)) +
      ylab("-Log10(adjusted p-val)") +
      xlab("Log2(Tumor/Normal)") +
      labs(title = isolate(input$cohort)) +
      geom_point(shape = 16, alpha = .5, mapping = aes(text = paste(
        "Gene: ", gene, "\n",
        "CRL: ", CRL, "\n",
        "LogFC: ", logFC, "\n",
        "P Value: ", log.scal.pval, "\n",
        sep = ""
      ))) +
      scale_color_manual(breaks = c(TRUE, FALSE), values = c("#FF0000", "#888888"), labels = c("CRL", "Non-CRL")) +
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
  output$rt_rn = renderPlotly({
    input$generateCorrelation
    output$corrErr2 <- renderText({
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
})
