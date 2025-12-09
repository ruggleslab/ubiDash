library("glue")
library("readxl")
library("tidyverse")
library("googledrive")
library("memoise")
library("plotly")
library("ComplexUpset")
library(scales) # for scales::alpha
library(ggpubr) # for stat_cor
library(ggprism) # for theme_prism
library(patchwork) # for combining plots
library(jsonlite)

tab_summaries <- fromJSON("./www/tab_summaries.json")


options(rsconnect.max.bundle.size = 3145728000)

# loading R packages
source("./scripts/load.packages.R")
pks_lst <- readLines("./www/packages.txt")
load.packages(pks_lst)

# Source all R files from R/ subfolders (Shiny only auto-sources root R/ files)
r_files <- list.files("./R", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (f in r_files) {
  source(f)
}

skip_load_data <- exists("skip_load_data") && get("skip_load_data")

if (!skip_load_data) {
  time <- as.numeric(Sys.time())
  startTime <- time

  # drive config ----

  options(
    gargle_oauth_email = "ruggleslab.shinyseq.backend@gmail.com",
    gargle_oauth_cache = ".secrets/" ## location within the shiny structure to store the auth token
  )

  ## Getting the auth token (only needs to be done once, so if you have the file you don't need to do this)
  # drive_auth()

  ## List all of the files that you have available to download
  # drive_find()
  # options(shiny.autoreload = TRUE)


  get_from_drive <- function(filename) {
    # print(filename)
    ## download the file in current wd and save the object so you know the location of the download
    ## This function can fail so wrap it in an error object to rerun it if it fails
    suppressMessages(
      file_data <- drive_download(filename, overwrite = TRUE)
    )

    ## Loading the Rds object based on the local path that it's stored in
    data <- readRDS(file_data$local_path)

    ## Deleting the file once it's loaded in to save on storage space
    if (file.exists(file_data$local_path)) {
      # Delete file if it exists
      file.remove(file_data$local_path)
    }

    return(data)
  }


  mem_get_from_drive <- memoise(get_from_drive)
  # mem_rdata_from_drive <- memoise(rdata_from_drive)


  mem_readRDS <- memoise(readRDS)


  load("./data/searchlists.RData")


  cohorts1.1 <- c(
    "BRCA", "COAD", "HGSC", "LUAD", "ccRCC", "UCEC", "PDAC", "LSCC", "GBM", "HNSCC"
    # ,"MB","HCC"
  )
  cohorts <- c("BRCA", "COAD", "HGSC", "LUAD", "CCRCC", "UCEC", "PDA", "LSCC", "GBM", "HNSCC", "OV")


  currentTime <- as.numeric(Sys.time())
  print(glue("read datasets in {(currentTime- time) * 1000} ms"))
  time <- currentTime


  genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
  crl_genes <- unique(genes$hgnc_symbol)

  currentTime <- as.numeric(Sys.time())
  print(glue("Initialized global functions in {(currentTime- time) * 1000} ms"))
  time <- currentTime

  genelists <- readRDS("./data/combined_genelists.rds")


  e3_ligase_list <- readRDS("./data/UbiDash_NYU+WEHI_E3_ligase_list.rds")

  # Top 15 UPS genes by mutation frequency (pre-computed to avoid MAF subsetting issues)
  # See EXP084_UPS_mutations analysis for methodology
  top_ups_genes <- c(
    "KMT2D", "KMT2C", "VHL", "HERC2", "BIRC6", "UBR4", "HUWE1", "HECW1",
    "RNF213", "MAGEL2", "CREBBP", "HERC1", "CMYA5", "HECTD4", "FBXW7"
  )

  # Top 15 E3 ligase genes by mutation frequency (pre-computed to avoid MAF subsetting issues)
  # See EXP084_UPS_mutations analysis for methodology
  top_e3_genes <- c(
    "VHL", "HERC2", "BIRC6", "UBR4", "HUWE1", "HECW1", "RNF213", "HERC1",
    "HECTD4", "FBXW7", "MYCBP2", "RNF43", "BRWD3", "MAP3K1", "KLHL1"
  )

  # Use top E3 genes as default for backwards compatibility
  top_crls <- top_e3_genes

  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared oncoplot data in {(currentTime- time) * 1000} ms"))
  time <- currentTime


  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared correlations data in {(currentTime- time) * 1000} ms"))
  time <- currentTime

  z_score <- function(data_matrix, which = c(1, 2)) { # 1 rows, 2 columns
    # Remove rows with all NA values
    data_matrix <- data_matrix[!apply(is.na(data_matrix), which, all), , drop = FALSE]
    # data_matrix <- data_matrix[!apply(is.nan(data_matrix), which, all), , drop = FALSE]

    # Calculate the row means and standard deviations
    row_means <- rowMeans(data_matrix, na.rm = TRUE)
    row_sds <- apply(data_matrix, which, sd, na.rm = TRUE)

    # Calculate Z-scores for each row
    z_scores <- (data_matrix - row_means) / row_sds

    return(z_scores)
  }

  colorsTbl <- read.csv("./data/colors.csv")

  drawHeatmap <- function(inputGenesMat, annotation, colsplit, clusterCol, clusterRow, legendName, column_names) {
    hm <- Heatmap(inputGenesMat,
      name = legendName,
      show_column_names = column_names, show_row_names = TRUE,
      top_annotation = annotation, heatmap_legend_param = list(direction = "horizontal"),
      column_split = colsplit, cluster_columns = clusterCol, cluster_rows = clusterRow,
      column_names_rot = 45,
      column_names_gp = grid::gpar(fontsize = 13.2),
      column_names_side = "top"
    )

    return(draw(hm,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "right"
    ))
  }


  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared heatmap data in {(currentTime- time) * 1000} ms"))
  time <- currentTime

  # Drug Sensitivity ----
  drug <- read.csv("./data/drugsensitivity/drug_longdf.csv", row.names = 1)
  # drug pathway and targets annotation
  drugpath <- read.csv("./data/drugsensitivity/drugpath.csv", row.names = 1)

  exp <- readRDS("./data/drugsensitivity/Cell_line_RMA_proc_basalExp.rds")
  cell <- read.csv("./data/drugsensitivity/MatchedCellAnno.csv", row.names = 1)


  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared drug sensitivity data in {(currentTime- time) * 1000} ms"))
  time <- currentTime


  # Simplified Heatmap ----


  # get_prt.global.AI <- function(){
  #   data = get_from_drive("cptac_tmt-log2-AI_protein_expression.rds")
  #   data = data_norm(data, col_center = T)
  #   return(data)
  # }

  # mem_prt.global.AI <- memoise(get_prt.global.AI)


  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared simplified heatmap data in {(currentTime- time) * 1000} ms"))
  time <- currentTime


  my.colors <-
    setNames(c("#070A52", "#6C9BCF", "#AFD3E2"),
      # colors to use
      nm = c("High", "Medium", "Low")
    ) # for each category

  # volcano ----
  # Retrieve all sheet names from the Excel file
  DE_protein_cohorts <- excel_sheets("./data/Table S2A. DEqMS_Percancer_CPTAC_imputed_proteome_matchedTN.xlsx")
  # Exclude the sheet named "long table"
  DE_protein_cohorts <- DE_protein_cohorts[DE_protein_cohorts != "long-table"]
  # Read only the sheets you want into a list of data frames
  DE_protein <- lapply(DE_protein_cohorts, function(cohort) {
    read_excel("./data/Table S2A. DEqMS_Percancer_CPTAC_imputed_proteome_matchedTN.xlsx", sheet = cohort)
  })
  # Convert each element to a data.frame (if needed)
  DE_protein <- lapply(DE_protein, as.data.frame)
  # Name the list elements after the corresponding cohorts
  names(DE_protein) <- DE_protein_cohorts

  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared volcano plot data in {(currentTime- time) * 1000} ms"))
  time <- currentTime

  # correlations ----


  currentTime <- as.numeric(Sys.time())
  print(glue("Prepared volcano plot data in {(currentTime- time) * 1000} ms"))
  time <- currentTime

  print(glue("Total startup time was {(time-startTime) * 1000} ms"))
  rm(time, currentTime, startTime)

  skip_load_data <- T
}
