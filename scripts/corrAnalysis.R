panelCorr <- function(geneNamesRNA, geneNamesProtein, tissueType, cohort, cptac_protein, cptac_rna, cptac_samples) {
  if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    # cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort %in% cohort, ]

    tumNames <- as.character(
      intersect(
        cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id,
        intersect(colnames(cptac_protein), colnames(cptac_rna))
      )
    )
    normNames <- as.character(intersect(
      cptac_samples_filt[cptac_samples_filt$type == "Normal", ]$harmonized_data_id,
      intersect(colnames(cptac_protein), colnames(cptac_rna))
    ))

    if (tissueType == "Tumor") {
      usedTissue <- tumNames
    } else if (tissueType == "Normal") {
      usedTissue <- normNames
    }

    if (length(geneNamesRNA) < 5 | length(geneNamesProtein) < 5) {
      stop("Minimum number of input genes/proteins is five.")
    }

    # Filter for genes
    selectedGenesProt <- cptac_protein[rownames(cptac_protein) %in% geneNamesProtein, colnames(cptac_protein) %in% usedTissue]
    selectedGenesRna <- cptac_rna[rownames(cptac_rna) %in% geneNamesRNA, colnames(cptac_rna) %in% usedTissue]

    protNAs <- rowSums(is.na(selectedGenesProt))
    protsToRemove <- ncol(selectedGenesProt) - protNAs
    protsToRemove <- names(protsToRemove[protsToRemove < 10])
    selectedGenesProt <- selectedGenesProt[!(rownames(selectedGenesProt) %in% protsToRemove), ]

    # keep common samples and genes
    commonRows <- intersect(rownames(selectedGenesProt), rownames(selectedGenesRna))
    commonCols <- intersect(colnames(selectedGenesProt), colnames(selectedGenesRna))
    selectedGenesProt <- selectedGenesProt[commonRows, commonCols]
    selectedGenesRna <- selectedGenesRna[commonRows, commonCols]

    if (ncol(selectedGenesProt) < 10 | ncol(selectedGenesRna) < 10) {
      stop("One or more genes have missing values, therefore decreasing the number of samples available for analysis. Please, adjust input parameters (e.g., add cohorts and/or additional genes).")
    }

    cor.mat <- cor(x = t(selectedGenesProt), y = t(selectedGenesRna), use = "pairwise.complete.obs", method = "spearman")
    p.mat <- cor.mtest(cor.mat, method = "spearman")

    flattenCorrMatrix <- function(cormat, pmat) {
      data.frame(
        row = rownames(cormat)[row(cormat)],
        column = rownames(cormat)[col(cormat)],
        rho = c(cormat),
        p = c(pmat)
      )
    }

    df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p + 1e-7)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- ifelse(df$rho > 0, "Positive", "Negative")

    corrPlot <- ggballoonplot(
      data = df,
      size = "log_pval",
      facet.by = NULL,
      size.range = c(1, 10),
      shape = 21,
      color = ifelse(df$qval < 0.01, "black", "transparent"),
      color.range = c(0, 5),
      fill = "rho",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE
    ) +
      gradient_fill(c("purple", "white", "red"))

    # Remove axis labels (i.e. the columns label)
    corrPlot <- corrPlot + labs(x = NULL, y = NULL) +
      aes(text = paste("Adjusted p-value:", qval))

    ggplotly(corrPlot, tooltip = "text")
  } else if (!is.null(cptac_protein) | !is.null(cptac_rna)) { # for RNA or Protein separately
    varToAssign <- list(cptac_protein, cptac_rna)
    inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein), !is.null(cptac_rna)))])
    rm(varToAssign)
    geneName <- c(geneNamesProtein, geneNamesRNA)
    if (length(geneName) < 5) {
      stop("Minimum number of input genes/proteins is five.")
    }

    colnames(inPut) <- gsub(pattern = ".?^X", replacement = "", x = colnames(inPut))
    colnames(inPut) <- gsub(pattern = "\\.", replacement = "-", x = colnames(inPut))

    # cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort %in% cohort, ]

    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id, colnames(inPut)))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal", ]$harmonized_data_id, colnames(inPut)))

    if (tissueType == "Tumor") {
      usedTissue <- tumNames
    } else if (tissueType == "Normal") {
      usedTissue <- normNames
    }

    selectedGenes <- inPut[rownames(inPut) %in% geneName, colnames(inPut) %in% usedTissue]
    if (ncol(selectedGenes) < 10) {
      stop("One or more genes have missing values, therefore decreasing the number of samples available for analysis. Please, adjust input parameters (e.g., add cohorts and/or additional genes).")
    }

    cor.mat <- cor(t(selectedGenes), use = "pairwise.complete.obs", method = "spearman")
    p.mat <- cor.mtest(mat = cor.mat, method = "spearman")

    flattenCorrMatrix <- function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        rho = cormat[ut],
        p = pmat[ut]
      )
    }
    df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p + 1e-7)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- ifelse(df$rho > 0, "Positive", "Negative")

    corrPlot <- ggballoonplot(
      data = df,
      size = "log_pval",
      facet.by = NULL,
      size.range = c(1, 10),
      shape = 21,
      color = ifelse(df$qval < 0.01, "black", "transparent"),
      color.range = c(0, 5),
      fill = "rho",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE
    ) +
      theme_minimal() +
      gradient_fill(c("purple", "white", "red"))

    # Remove axis labels (columns labels) from the plot
    corrPlot <- corrPlot + labs(x = NULL, y = NULL) +
      aes(text = paste("Adjusted p-value:", qval))

    ggplotly(corrPlot, tooltip = "text")
  }

  res <- list(
    corrPlot = corrPlot,
    corTab = list(cor = cor.mat, pValue = p.mat$p),
    usedCohorts = paste0(unique(cptac_samples_filt$cohort), collapse = "_")
  )
  return(res)
}
