source('./RScripts/subset.R')

load("./data/20220214_ANA044_Proteome_Pancan_mat.combo_observed_norm.Rdata")
load("./data/sample_info.Rdata")
labels = subset(samples, sample_type %in% c("T", "N"), select = c("sample_type", "cohort_code"))

mat = as.data.frame(mat)
largeMat = mat[, rownames(labels)]
filteredMat = largeMat[genesSubset$genes,]

filteredMat <- filteredMat[which(rowMeans(!is.na(filteredMat)) > 0.5),]

typeCols = list(Type  = c("T" = "purple", "N" = "orange"))

colors = c("royalblue4",
           "maroon",
           "goldenrod1",
           "aquamarine",
           "pink3",
           "coral4",
           "slategray3",
           "mediumpurple",
           "seagreen",
           "lightgoldenrod2",
           "orangered3")

colorList = list(Type = setNames(c("purple", "orange"), nm = c("T", "N")),
                 Cohort = setNames(colors, nm = unique(samples$cohort_code)))


ha = HeatmapAnnotation(Type = labels$sample_type, Cohort = labels$cohort_code, col = colorList)
# cohortAnnotation = HeatmapAnnotation()

defaultHeatmap = Heatmap(filteredMat, name = "Log2(Relative Abundance)",
            show_column_names = FALSE, show_row_names = TRUE,
            top_annotation = ha,heatmap_legend_param = list(direction = "horizontal"),
            column_split = labels$sample_type)

defaultHeatmap = draw(defaultHeatmap, heatmap_legend_side="bottom", annotation_legend_side="right")

save(defaultHeatmap, file = "./data/defaultHeatmap.RData")

