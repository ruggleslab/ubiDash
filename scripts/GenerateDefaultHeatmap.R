load("./data/defaults.RData")
source("./scripts/HeatmapFunctions.R")

# cptac dataset
cptac_rna <- readRDS("./data/cptac_tpm_rna_expression.rds")
cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
cptac_samples <- readRDS("./data/cptac_sample_info.rds")
cptac_samples_subtypes <- readRDS("./data/cptac_sample_info_with_subtypes.rds")


# cohort and sample type of all normal and tumor samples
# labels <- subset(cptac_samples, sample_type %in% c("T", "N","G"), select = c("sample_type", "cohort_code"))
tnlabels <- subset(cptac_samples_subtypes, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code", "Cancer_subtype"))
# get list of cohorts
cohorts <- unique(tnlabels$cohort_code)


# cptac_protein <- cptac_protein[, ]
filteredMat <- cptac_protein[genesSubset$genes, rownames(tnlabels)]

# remove rows with over 50% NA
filteredMat <- filteredMat[which(rowMeans(!is.na(filteredMat)) > 0.5), ]

colorsTbl <- read.csv("./data/colors.csv")
colors <- colorsTbl$color[1:11]
colorList <- list(
  Type = setNames(colorsTbl$color[12:13], nm = c("T", "N")),
  Cohort = setNames(colors, nm = cohorts)
)


# matrix with only genes to plot
getGenesMat <- function(inputGenes, cohortSubsetLabels) {
  return(filteredMat[, rownames(cohortSubsetLabels)])
}


# get cohorts selected by user
cohortSubsetLabels <- getCohortPatients(tnlabels, cohorts)

<<<<<<< HEAD
inputGenesMat <- getGenesMat(inputGenes,cohortSubsetLabels)
=======
# inputGenes <- getInputGenes(c(),C(), ",")
inputGenesMat <- getGenesMat(inputGenes, cohortSubsetLabels)
>>>>>>> 35c8a31a35d5b21dd589cc26d1ba85476d8ae864
# print(inputGenesMat)

if (is.null(inputGenesMat)) {
  if (isolate(input$cluster) && isolate(input$separateNandT)) {
    return(defaultHeatmap)
  }
  inputGenesMat <- filteredMat
} else {
  inputGenesMat <- as.matrix(inputGenesMat[which(rowMeans(!is.na(inputGenesMat)) > 0.5), ])
}

# cutout <- inputGenes[!(inputGenes %in% toupper(rownames(inputGenesMat)))]


ha <- HeatmapAnnotation(Type = cohortSubsetLabels$sample_type, Cohort = cohortSubsetLabels$cohort_code, Subtype = cohortSubsetLabels$Cancer_subtype, col = colorList)
hm <- NULL

colsplit <- getCohortPatients(tnlabels, cohorts)$sample_type

inputGenesMat <- as.matrix(inputGenesMat)

h <- drawHeatmap(inputGenesMat, ha, colsplit, TRUE)


heatmap <- Heatmap(inputGenesMat,
  name = "Log2(Relative Abundance)",
  show_column_names = FALSE, show_row_names = TRUE,
  top_annotation = ha, heatmap_legend_param = list(direction = "horizontal"),
  column_split = colsplit, cluster_columns = TRUE
)


heatmap <- draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "right")

defaultHeatmap <- heatmap
save(cohortsOnco, defaultHeatmap, genesSubset, mutationsOnco, file = "./data/defaults.RData")
