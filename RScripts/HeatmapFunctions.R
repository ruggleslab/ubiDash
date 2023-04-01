getCohortPatients <- function(labels, cohortSelect) {
  patients <- filter(labels, cohort_code %in% cohortSelect)
  return(patients)
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

drawHeatmap <- function(inputGenesMat, annotation,colsplit, cluster){
  hm <- Heatmap(inputGenesMat,
                name = "Log2(Relative Abundance)",
                show_column_names = FALSE, show_row_names = TRUE,
                top_annotation = annotation, heatmap_legend_param = list(direction = "horizontal"),
                column_split = colsplit, cluster_columns = cluster
  )
  return(draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right"))
}