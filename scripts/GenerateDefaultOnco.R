maf <- readRDS("./data/oncoplot_maf.rds")

genes <- read_excel(path = "./data/AllCRLs.v5.xlsx", sheet = 1)
crl_genes = unique(genes$hgnc_symbol)

inputGenes = crl_genes

smallFM = filter(maf, Hugo_Symbol %in% inputGenes)

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
if(nrow(sortedTbl) > 10){
  topGenes = rownames(sortedTbl[1:10,])
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

mutationsTbl = mutationsTbl[rownames(pcts),]
cohortsTbl = cohortsTbl[rownames(pcts),]

mutations = levels(smallFM$Variant_Classification)

source("./scripts/OncoplotColors.R")

mutation_colors = get_mutation_colors(mutations)
mutation_alter_fun = get_mutation_alter_fun(mutation_colors)
    
cohort_colors = get_cohort_colors()
cohort_alter_fun = get_cohort_alter_fun(cohort_colors)

mutationsOnco = oncoPrint(
  mutationsTbl, alter_fun = mutation_alter_fun, col = mutation_colors, 
  row_names_side = "left", 
  show_pct = FALSE,
  
  row_order = topGenes,
  
  # remove_empty_columns = TRUE,
  
  heatmap_legend_param = list(title = "Alterations"),
  
  right_annotation = rowAnnotation(pct = anno_text(pcts[,2]), barplot = anno_oncoprint_barplot())
)
cohortsOnco = oncoPrint(
  cohortsTbl, alter_fun = cohort_alter_fun, col = cohort_colors, 
  row_names_side = "left", 
  show_pct = FALSE,
  
  row_order = topGenes,
  
  # remove_empty_columns = TRUE,
  
  heatmap_legend_param = list(title = "Alterations"),
  
  right_annotation = rowAnnotation(pct = anno_text(pcts[,2]), barplot = anno_oncoprint_barplot())
)

save(cohortsOnco, mutationsOnco, file = "./data/defaultOncoplot.RData") 
