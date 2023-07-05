# library(shiny)
# library(ComplexHeatmap)
# library(readxl)
# library(tidyverse)
# library(RColorBrewer)
# library(plotly)



packagelist = c("shiny", "ComplexHeatmap", "openxlsx", "readxl", "tidyverse", "glue", "stringr", "reshape2", "RColorBrewer", "plotly")

lapply(packagelist, function(x) {
  suppressMessages(require(x, character.only = T,
                           quietly=T, 
                           warn.conflicts = F))
})

load("./data/defaults.RData")
# loading default objects below:
# genesSubset
# defaultHeatmap
# cohortsOnco
# mutationsOnco

source("./RScripts/SplitInputGenes.R")

# heatmap variables ----
source("./modules/heatmap.R")

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
source("./modules/volcano.R")

# oncoplot variables ----
source('./modules/oncoplot.R')


# correlation variables ----
source("./modules/correlations.R")
# boxplot variables ----
source("./modules/violin.R")


# server ----

shinyServer(function(input, output, session) {
  
  
  # oncoplot tab ----
  oncoplotServer("oncoTab")
  
  # heatmap tab ----
  
  heatmapServer("heatmapTab")
  
  # volcano plot tab ----
  volcanoServer("volcanoTab")
  
  # correlation tab (----
  correlationsServer("corrTab")
  
  # boxplot tab ----
  violinServer("violinTab")
  
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
  output$outputTabGvgCorrRna <- renderTable({resGvgCorrRna()$singCorrTab},rownames = T, colnames = F)
  
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
  output$outputTabGvgCorrProtein <- renderTable({resGvgCorrProtein()$singCorrTab},rownames = T, colnames = F)
  
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
  output$outputTabGvgCorrBoth <- renderTable({resGvgCorrBoth()$singCorrTab},rownames = T, colnames = F)
  
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
      list(extend = "excel", text = "Download results", filename = paste0(input$inputGvACorrProteinGeneName, "_ProteinCorrelations"),
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
  
  output$dataDownLoadGvpCorrRna <- renderUI({
    req(input$generateGvpCorrRna)
    downloadButton("downLoadGvpCorrRna", label = "Download correlation data") })
  
  output$downLoadGvpCorrRna <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                               content = function(file) {
                                                 write.xlsx(resGvpCorrRna()$corTab, file,row.names = TRUE)})
  
  
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
  
  output$dataDownLoadGvpCorrProtein <- renderUI({
    req(input$generateGvpCorrProtein)
    downloadButton("downLoadGvpCorrProtein", label = "Download correlation data") })
  
  output$downLoadGvpCorrProtein <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                                   content = function(file) {
                                                     write.xlsx(resGvpCorrProtein()$corTab, file,row.names = TRUE)})
  
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
  
  output$dataDownLoadGvpCorrBoth <- renderUI({
    req(input$generateGvpCorrBoth)
    downloadButton("downLoadGvpCorrBoth", label = "Download correlation data") })
  
  output$downLoadGvpCorrBoth <- downloadHandler(filename = function() {paste0("ubiDashCorrData",".xlsx")},
                                                content = function(file) {
                                                  write.xlsx(resGvpCorrBoth()$corTab, file,row.names = TRUE)})
  
})
