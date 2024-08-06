library("glue")
library("readxl")
library("tidyverse")
library("googledrive")
library("memoise")

# loading R packages
source("./scripts/load.packages.R")
pks_lst <- readLines("./www/packages.txt")
load.packages(pks_lst)
skip_load_data = exists("skip_load_data") && get("skip_load_data")

if(!skip_load_data){
time = as.numeric(Sys.time())
startTime = time

# drive config ----

options(
  gargle_oauth_email = "ruggleslab.shinyseq.backend@gmail.com",
  gargle_oauth_cache = ".secrets/" ## location within the shiny structure to store the auth token
)

## Getting the auth token (only needs to be done once, so if you have the file you don't need to do this)
#drive_auth()

## List all of the files that you have available to download
# drive_find()
# options(shiny.autoreload = TRUE)


get_from_drive <- function(filename){
  print(filename)
  ## download the file in current wd and save the object so you know the location of the download
  ## This function can fail so wrap it in an error object to rerun it if it fails
  file_data <- drive_download(filename, overwrite = TRUE)
  
  ## Loading the Rds object based on the local path that it's stored in
  data <- readRDS(file_data$local_path)
  
  ## Deleting the file once it's loaded in to save on storage space
  if (file.exists(file_data$local_path)) {
    #Delete file if it exists
    file.remove(file_data$local_path)
  }
  
  return(data)
}

# rdata_from_drive <- function(filename){
#   ## download the file in current wd and save the object so you know the location of the download
#   ## This function can fail so wrap it in an error object to rerun it if it fails
#   file_data <- drive_download(filename, overwrite = TRUE)
# 
#   ## Loading the Rds object based on the local path that it's stored in
#   load(file_data$local_path)
# 
#   ## Deleting the file once it's loaded in to save on storage space
#   if (file.exists(file_data$local_path)) {
#     #Delete file if it exists
#     file.remove(file_data$local_path)
#   }
# 
# }


mem_get_from_drive <- memoise(get_from_drive)
# mem_rdata_from_drive <- memoise(rdata_from_drive)



mem_readRDS <- memoise(readRDS)

# Datasets ----

# ccle dataset
# ccle_rna <- mem_get_from_drive("ccle_tpm_rna_expression.rds")
# ccle_protein <- mem_get_from_drive("ccle_tmt_protein_expression.rds")
# ccle_samples <- mem_get_from_drive("ccle_sample_info.rds")

# cptac dataset
# cptac_rna <- mem_get_from_drive("duplicates_averaged_cptac_rna.rds")
# cptac_protein <- mem_get_from_drive("cptac_tmt-log2-RI_protein_expression.rds")
cptac_samples <- mem_get_from_drive("cptac_sample_info_with_subtypes.rds")

load("./data/searchlists.RData")


#pride dataset
# pride_protein <- mem_get_from_drive("pride_lfq_protein_expression.rds")


# for correlations tab
# plot_tab <- mem_get_from_drive("cptac_multi_omic_correlations.rds")


currentTime =  as.numeric(Sys.time())
print(glue("read datasets in {(currentTime- time) * 1000} ms"))
time = currentTime

# # maf for lollipop and oncoplot tabs 
# maf = readRDS("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds")

# currentTime =  as.numeric(Sys.time())
# print(glue("read maf in {(currentTime- time) * 1000} ms"))
# time = currentTime


# Information used by various modules ----

# cohort  and sample type of all normal and tumor samples
tnlabels <- subset(cptac_samples, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code", "Cancer_subtype"))
# get list of cohorts
cohorts <- unique(tnlabels$cohort_code)

genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
crl_genes = unique(genes$hgnc_symbol)

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
# load("./data/defaults.RData")




colorsTbl = read.csv("./data/colors.csv")
colors = colorsTbl$color[1:11]
colorList <- list(
  Type = setNames(colorsTbl$color[12:13], nm = c("T", "N")),
  Cohort = setNames(colors, nm = cohorts)
)


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

exp <- readRDS("./data/drugsensitivity/Cell_line_RMA_proc_basalExp.rds")
cell = read.csv("./data/drugsensitivity/MatchedCellAnno.csv", row.names = 1)


currentTime =  as.numeric(Sys.time())
print(glue("Prepared drug sensitivity data in {(currentTime- time) * 1000} ms"))
time = currentTime


#Simplified Heatmap ----

source("./data/simplifiedheatmap/tgr_misc.R")


get_prt.global.AI <- function(){
  data = get_from_drive("cptac_tmt-log2-AI_protein_expression.rds")
  data = data_norm(data, col_center = T)
  return(data)
}

mem_prt.global.AI <- memoise(get_prt.global.AI)


labels <- read.csv("./data/simplifiedheatmap/sub_types.csv")


tertiles <- function(x) {
  vTert <- quantile(x, c(0:3 / 3), na.rm = TRUE)
  cut(
    x,
    vTert,
    include.lowest = TRUE,
    labels = c("Low", "Medium", "High")
  )
}


currentTime =  as.numeric(Sys.time())
print(glue("Prepared simplified heatmap data in {(currentTime- time) * 1000} ms"))
time = currentTime


my.colors <-
  setNames(c("#070A52", "#6C9BCF", "#AFD3E2"),
           # colors to use
           nm = c("High", "Medium", "Low")) # for each category

# volcano ----
DE_protein_cohorts <- excel_sheets("./data/DE_protein.xlsx")
DE_protein <- lapply(DE_protein_cohorts, function(cohort) read_excel("./data/DE_protein.xlsx", sheet = cohort))
DE_protein <- lapply(DE_protein, as.data.frame)
names(DE_protein) <- DE_protein_cohorts


currentTime =  as.numeric(Sys.time())
print(glue("Prepared volcano plot data in {(currentTime- time) * 1000} ms"))
time = currentTime
# dat <- list(rna = mem_get_from_drive("CPTAC_Pancan_UPS_genome-wide_eQTLs.rds"),
#             pro = mem_get_from_drive("CPTAC_Pancan_UPS_genome-wide_pQTLs.rds"))

print(glue("Total startup time was {(time-startTime) * 1000} ms"))
rm(time, currentTime, startTime)

skip_load_data = T
}
