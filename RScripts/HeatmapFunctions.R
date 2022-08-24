getCohortPatients <- function(labels, cohortSelect) {
  patients <- filter(labels, cohort_code %in% cohortSelect)
  return(patients)
}

getInputGenes <- function(CSV, TXT) {
  inFile <- CSV
  if (!is.null(inFile)) {
    return(unlist(read.csv(file = inFile$datapath, header = FALSE)))
  }
  TXT <- gsub(",", " ", TXT)
  TXT <- gsub("  ", " ", TXT)
  return(unlist(strsplit(TXT, " ")))
}

drawHeatmap <- function(inputGenesMat, annotation,colsplit, cluster){
  hm <- Heatmap(inputGenesMat,
                name = "Log2(Relative Abundance)",
                show_column_names = FALSE, show_row_names = TRUE,
                top_annotation = annotation, heatmap_legend_param = list(direction = "horizontal"),
                column_split = colsplit, cluster_columns = cluster
  )
  return(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right"))
}