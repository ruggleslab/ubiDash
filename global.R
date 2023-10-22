library("glue")
library("readxl")
library("tidyverse")

time = as.numeric(Sys.time())
startTime = time

skip_load_data = F

if(!skip_load_data){
  # Datasets ----
  
  # ccle dataset
  ccle_rna <- readRDS("./data/ccle_tpm_rna_expression.rds")
  ccle_protein <- readRDS("./data/ccle_tmt_protein_expression.rds")
  ccle_samples <- readRDS("./data/ccle_sample_info.rds")
  
  # cptac dataset
  # cptac_rna <- readRDS("./data/cptac_tpm_rna_expression.rds")
  # cptac_rna <- readRDS("./data/cptac_log-tpm_rna_expression.rds")
  cptac_rna <- readRDS("./data/duplicates_averaged_cptac_rna.rds")
  cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
  cptac_samples <- readRDS("./data/cptac_sample_info_with_subtypes.rds")
  
  
  #pride dataset
  pride_protein <- readRDS("./data/pride_lfq_protein_expression.rds")
  
  # for correlations tab
  plot_tab <- readRDS("./data/cptac_multi_omic_correlations.rds")
  
  
  currentTime =  as.numeric(Sys.time())
  print(glue("read datasets in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  # maf for lollipop and oncoplot tabs 
  maf = readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")
  
  currentTime =  as.numeric(Sys.time())
  print(glue("read maf in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  
  # Information used by various modules ----
  
  # cohort  and sample type of all normal and tumor samples
  tnlabels <- subset(cptac_samples, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code", "Cancer_subtype"))
  # get list of cohorts
  cohorts <- unique(tnlabels$cohort_code)
  
  genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
  crl_genes = unique(genes$hgnc_symbol)
  
  getInputGenes <- function(CSV, TXT, delim = "tab") {
    delims = c(
      tab = "\t",
      comma = ","
    )
    inFile <- CSV
    if (!is.null(inFile)) {
      return(unlist(read.table(file = inFile$datapath, header = FALSE, sep = delims[delim])))
    }
    return(TXT)
  }
  
  
  getInputGenes <- function(CSV, TXT, delim = "tab") {
    delims = c(
      tab = "\t",
      comma = ","
    )
    inFile <- CSV
    if (!is.null(inFile)) {
      return(unlist(read.table(file = inFile$datapath, header = FALSE, sep = delims[delim])))
    }
    return(TXT)
  }
  
  
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Initialized global functions in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  
  # Tab-specific variables ----
  
  # Oncoplot ----
  
  cohort_colors = list(cohort_code = c(BRCA = "royalblue4",
                                       CCRCC = "maroon4",
                                       COAD = "goldenrod1",
                                       GBM = "aquamarine",
                                       HCC = "pink3",
                                       HNSCC = "coral4",
                                       LSCC = "slategray3",
                                       LUAD = "mediumpurple",
                                       OV = "seagreen",
                                       PDA = "lightgoldenrod2",
                                       UCEC = "orangered3"))
  
  top_crls = readRDS("./data/topcrls.rds")
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Prepared oncoplot data in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  
  # Correlations ----
  source("./scripts/GenerateCorrPlot.R")
  
  lines = list(geom_hline(yintercept=.5, color = "red") ,
               geom_hline(yintercept=-.5,  color = "red") ,
               geom_vline(xintercept=-.5,color = "red") ,
               geom_vline(xintercept=--.5, color = "red") ,
               geom_hline(yintercept=0, size = 1, color = "black") ,
               geom_vline(xintercept=-0, size = 1, color = "black"))
  
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Prepared correlations data in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  # Heatmap ----
  load("./data/defaults.RData")
  
  filteredProteinMat <- cptac_protein[genesSubset$genes, rownames(tnlabels)]
  # remove rows with over 50% NA
  filteredProteinMat <- filteredProteinMat[which(rowMeans(!is.na(filteredProteinMat)) > 0.5), ]
  
  filteredTnlabels = tnlabels[rownames(tnlabels) %in% colnames(cptac_rna),]
  filteredMRNAMat <- cptac_rna[genesSubset$genes, rownames(filteredTnlabels)]
  # remove rows with over 50% NA
  filteredMRNAMat <- filteredMRNAMat[which(rowMeans(!is.na(filteredMRNAMat)) > 0.5), ]
  
  
  colorsTbl = read.csv("./data/colors.csv")
  colors = colorsTbl$color[1:11]
  colorList <- list(
    Type = setNames(colorsTbl$color[12:13], nm = c("T", "N")),
    Cohort = setNames(colors, nm = cohorts)
  )
  
  # getCohortPatients <- function(labels, cohortSelect) {
  #   patients <- filter(labels, cohort_code %in% cohortSelect)
  #   return(patients)
  # }
  
  
  drawHeatmap <- function(inputGenesMat, annotation,colsplit, clusterCol, clusterRow, legendName){
    hm <- Heatmap(inputGenesMat,
                  name = legendName,
                  show_column_names = FALSE, show_row_names = TRUE,
                  top_annotation = annotation, heatmap_legend_param = list(direction = "horizontal"),
                  column_split = colsplit, cluster_columns = clusterCol, cluster_rows = clusterRow
    )
    return(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right"))
  }
  
  
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Prepared heatmap data in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  # Drug Sensitivity ----
  drug <- read.csv("./data/drugsensitivity/drug_longdf.csv", row.names = 1)
  #drug pathway and targets annotation
  drugpath <- read.csv("./data/drugsensitivity/drugpath.csv", row.names = 1)
  
  # exp <- read.delim("./data/drugsensitivity/Cell_line_RMA_proc_basalExp.txt", header = T)
  exp <- readRDS("./data/drugsensitivity/Cell_line_RMA_proc_basalExp.rds")
  cell = read.csv("./data/drugsensitivity/MatchedCellAnno.csv", row.names = 1)
  
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Prepared drug sensitivity data in {(currentTime- time) * 1000} ms"))
  time = currentTime
  
  
  #Simplified Heatmap ----
  
  source("./data/simplifiedheatmap/tgr_misc.R")
  # prt.n = readRDS("./data/simplifiedheatmap/prt_n.rds")
  
  # prt.global.RI <-
  #   readRDS(file = "./data/cptac_tmt-log2-RI_protein_expression.rds")
  # prt.global.RI <- data_norm(prt.global.RI)

  prt.global.AI <-
    readRDS(file = "./data/cptac_tmt-log2-AI_protein_expression.rds")
  prt.global.AI <- data_norm(prt.global.AI, col_center = T)
  
  labels <- read.csv("./data/simplifiedheatmap//sub_types.csv")
  
  
  tertiles <- function(x) {
    vTert <- quantile(x, c(0:3 / 3), na.rm = TRUE)
    cut(
      x,
      vTert,
      include.lowest = TRUE,
      labels = c("Low", "Medium", "High")
    )
  }
  
  
  # Apply function to each column of the matrix
  # mTert <- apply(prt.n, 2, tertiles)
  # row.names(mTert) <- row.names(prt.n)
  
  
  my.colors <-
    setNames(c("#070A52", "#6C9BCF", "#AFD3E2"),
             # colors to use
             nm = c("High", "Medium", "Low")) # for each category
  
  currentTime =  as.numeric(Sys.time())
  print(glue("Prepared simplified heatmap data in {(currentTime- time) * 1000} ms"))
  time = currentTime
}

print(glue("Total startup time was {(time-startTime) * 1000} ms"))
rm(time,currentTime,startTime)


