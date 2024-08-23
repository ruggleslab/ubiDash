# library(corrplot)
# setwd("E:/ubiDash/")
# #cohort <- c("GBM","COAD", "CCRCC", "HCC", "HNSCC", "LSCC", "OV", "PDA")
# cohort <- c("BRCA")
# tissueType <- "Tumor"
# cutOff <- c(-0.05,0.05)
#
# cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
# cptac_rna <- readRDS("./data/cptac_tpm_rna_expression.rds")
# #cptac_rna <- NULL
# cptac_samples <- readRDS("./data/cptac_sample_info.rds")
# #panel
# geneNamesRNA <- c("A1BG", "A2M", "A2ML1", "AAAS", "A1CF", "AACS")#rownames(cptac_rna)[1:17]
# geneNamesProtein <- c("A1BG", "A2M", "A2ML1", "AAAS", "A1CF", "AACS")#rownames(cptac_protein)[1:20]
# #GvA
# geneName <- "A1BG"
#
# #single
# geneNameA <- "A1BG"
# geneNameB <- "A2M"
#
# asd <- singleCorr(geneNameA = "A1BG",geneNameB = "A2M",tissueType = "Tumor",cohort = c("BRCA","COAD", "CCRCC"),
#                   cptac_protein = NULL,
#                   cptac_rna = cptac_rna,
#                   cptac_samples = cptac_samples)
# asd$singCorrTab

singleCorr <- function(geneNameA, geneNameB, tissueType, cohort, cptac_protein, cptac_rna, cptac_samples, colorsTbl) {
  xlab = "Log2(Intensity)"
  ylab = "Log2(TPM + 1)"
  geneTypeA = "protein"
  geneTypeB = "rna"
  if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    # cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohor %in% cohort, ]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id, intersect(colnames(cptac_protein), colnames(cptac_rna))))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal" | cptac_samples_filt$sample_type == "G", ]$harmonized_data_id, intersect(colnames(cptac_protein), colnames(cptac_rna))))
    
    if (tissueType == "Tumor") {
      usedTissue <- tumNames
    } else if (tissueType == "Normal") {
      usedTissue <- normNames
    }
    cptac_samples_filt <- cptac_samples_filt[cptac_samples_filt$harmonized_data_id %in% usedTissue, ]
    
    # Filter for genes
    selectedGeneA <- as.data.frame(cptac_protein)[rownames(cptac_protein) == geneNameA, colnames(cptac_protein) %in% usedTissue]
    selectedGeneB <- as.data.frame(cptac_rna)[rownames(cptac_rna) == geneNameB, colnames(cptac_rna) %in% usedTissue]
    # remove NAs
    selectedGeneA <- selectedGeneA[, !is.na(selectedGeneA)]
    selectedGeneB <- selectedGeneB[, !is.na(selectedGeneB)]
    # keep common samples
    selectedGeneA <- selectedGeneA[, intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    selectedGeneB <- selectedGeneB[, intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
  } else if (!is.null(cptac_protein) | !is.null(cptac_rna)) { # for RNA or Protein separately
    
    if(!is.null(cptac_protein)){
      ylab = "Log2(Intensity)"
      geneTypeB = "protein"
    }
    if(!is.null(cptac_rna)){
      xlab = "Log2(TPM + 1)"
      geneTypeA = "rna"
    }
    
    varToAssign <- list(cptac_protein, cptac_rna)
    inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein), !is.null(cptac_rna)))])
    rm(varToAssign)
    
    colnames(inPut) <- gsub(pattern = ".?^X", replacement = "", x = colnames(inPut))
    colnames(inPut) <- gsub(pattern = "\\.", replacement = "-", x = colnames(inPut))
    # cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort %in% cohort, ]

    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id, colnames(inPut)))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal" | cptac_samples_filt$sample_type == "G", ]$harmonized_data_id, colnames(inPut)))
    
    
    if (tissueType == "Tumor") {
      usedTissue <- tumNames
    } else if (tissueType == "Normal") {
      usedTissue <- normNames
    }
    
    cptac_samples_filt <- cptac_samples_filt[cptac_samples_filt$harmonized_data_id %in% usedTissue, ]
    
    # Filter for genes
    selectedGeneA <- inPut[rownames(inPut) == geneNameA, colnames(inPut) %in% usedTissue]
    selectedGeneB <- inPut[rownames(inPut) == geneNameB, colnames(inPut) %in% usedTissue]
    # remove NAs
    selectedGeneA <- selectedGeneA[, !is.na(selectedGeneA)]
    selectedGeneB <- selectedGeneB[, !is.na(selectedGeneB)]
    # keep common samples
    selectedGeneA <- selectedGeneA[, intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    selectedGeneB <- selectedGeneB[, intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
  }
  
  if (length(selectedGeneA) < 10 | length(selectedGeneB) < 10) {
    stop("One or more genes have missing values, therefore decreasing the number of samples available for analysis. Please, adjust input parameters (e.g., add cohorts and/or additional genes).")
  }
  
  selectedGeneA <- as.data.frame(t(selectedGeneA))
  selectedGeneB <- as.data.frame(t(selectedGeneB))
  
  selectedGeneA_Cohort <- cptac_samples_filt[cptac_samples_filt$harmonized_data_id %in% rownames(selectedGeneA) & cptac_samples_filt$datatype == geneTypeA, ]$cohort
  selectedGeneB_Cohort <- cptac_samples_filt[cptac_samples_filt$harmonized_data_id %in% rownames(selectedGeneB) & cptac_samples_filt$datatype == geneTypeB, ]$cohort
  
  selectedGeneA_Cohort = selectedGeneA_Cohort[!duplicated(selectedGeneA_Cohort)]
  selectedGeneB_Cohort = selectedGeneB_Cohort[!duplicated(selectedGeneB_Cohort)]
  
  selectedGeneA$Cohort = selectedGeneA_Cohort
  selectedGeneB$Cohort = selectedGeneB_Cohort
  
  colnames(selectedGeneA) <- c("ExpValueA", "Cohort")
  colnames(selectedGeneB) <- c("ExpValueB", "Cohort")
  
  corrSum <- cor.test(x = as.numeric(selectedGeneA$ExpValueA), y = as.numeric(selectedGeneB$ExpValueB), method = "spearman")
  
  # visualization
  if (is.null(cptac_rna)) {
    dat <- data.frame("geneAExp" = as.numeric(selectedGeneA$ExpValueA), "geneBExp" = as.numeric(selectedGeneB$ExpValueB)) # Protein
  } else if (is.null(cptac_protein)) {
    dat <- data.frame("geneAExp" = log2(as.numeric(selectedGeneA$ExpValueA)), "geneBExp" = log2(as.numeric(selectedGeneB$ExpValueB))) # RNA
  } else if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    dat <- data.frame("geneAExp" = as.numeric(selectedGeneA$ExpValueA), "geneBExp" = log2(as.numeric(selectedGeneB$ExpValueB))) # Both
  }
  
  if (identical(rownames(selectedGeneA), rownames(selectedGeneB))) {
    dat$Cohort <- selectedGeneA$Cohort
  }
  
  
  
  # singCorrTab <- data.frame("N" = nrow(dat),
  #                           "Spearman's Rho" = round(x = corrSum$estimate,digits = 2),
  #                           "pValue" = format(x = corrSum$p.value,digits = 2, scientific = TRUE), row.names = NULL)
  
  pValue <- format(x = corrSum$p.value, digits = 2, scientific = TRUE)
  rho <- round(x = corrSum$estimate, digits = 2)
  singCorrTab <- t(data.frame(
    "Cohorts" = names(table(dat$Cohort)),
    "N" = as.numeric(unlist(table(dat$Cohort)))
  ))
  
  xmax <- max(dat$geneAExp)
  ymax <- max(dat$geneBExp)
  
  colors  = colorsTbl$color
  names(colors) = colorsTbl$id
  
  scatPlot <- ggplot(dat, aes(x = geneAExp, y = geneBExp, color = Cohort)) +
    geom_point(size = 2, alpha = .6) +
    scale_color_manual(values = colors) +
    geom_smooth(method = "lm", size = 2) +
    scale_fill_manual(values = colors) +
    # geom_smooth(aes(x = geneAExp, y = geneBExp,color = NULL),method = "lm", size = 2) +
    xlab(glue("{xlab} {geneNameA}")) +
    ylab(glue("{ylab} {geneNameB}")) +
    annotate(geom = "text", x = xmax, y = ymax, label = paste0("p = ", pValue, "\n", "R = ", rho)) +
    # coord_cartesian(clip = "off")+
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = NA, colour = "grey50"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, vjust = 2),
      axis.title.x = element_text(size = 12, vjust = 2)
    )
  # scatPlot
  
  
  res <- list(singCorrPlot = scatPlot, singCorrTab = singCorrTab)
  return(res)
}


multiCorr <- function(geneName, tissueType, cutOff, cohort, cptac_protein, cptac_rna, cptac_samples) {
  varToAssign <- list(cptac_protein, cptac_rna)
  inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein), !is.null(cptac_rna)))])
  rm(varToAssign)
  
  colnames(inPut) <- gsub(pattern = ".?^X", replacement = "", x = colnames(inPut))
  colnames(inPut) <- gsub(pattern = "\\.", replacement = "-", x = colnames(inPut))
  
  # cohort filter
  cptac_samples_filt <- cptac_samples[cptac_samples$cohor %in% cohort, ]
  
  tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id, colnames(inPut)))
  normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal" | cptac_samples_filt$sample_type == "G", ]$harmonized_data_id, colnames(inPut)))
  
  
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
  } else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  # make the selected gene to a numeric vector and perform correlation test
  selectedGene <- inPut[rownames(inPut) == geneName, colnames(inPut) %in% usedTissue]
  selectedGene <- selectedGene[, !is.na(selectedGene)]
  rawForAn <- inPut[, colnames(inPut) %in% colnames(selectedGene)]
  rawForAn <- as.data.frame(na.omit(rawForAn))
  selectedGene <- as.numeric(selectedGene)
  if (length(selectedGene) < 10) {
    stop("One or more genes have missing values, therefore decreasing the number of samples available for analysis. Please, adjust input parameters (e.g., add cohorts and/or additional genes).")
  }

  # correlation analysis
  
  corResult <- suppressWarnings(apply(rawForAn, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4, 3)]))

  # store the results as a dataframe and convert the list to a dataframe
  result <- as.data.frame(matrix(nrow = nrow(rawForAn), ncol = 3))
  colnames(result) <- c("Protein", "Spearman's Rho", "p-value")

  result[, 1] <- names(corResult)
  result[, 2] <- round(x = sapply(corResult, "[[", 1), digits = 2)
  result[, 3] <- format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
  
  result <- result[result$`Spearman's Rho` >= min(cutOff) & result$`Spearman's Rho` <= max(cutOff), ]
  result <- result[order(result$`Spearman's Rho`, decreasing = T), ]

  plotTab <- result
  plotTab$From <- geneName
  plotTab <- plotTab[, c(4, 1, 2, 3)]
  colnames(plotTab) <- c("From", "Genes", "R_value", "P_value")
  plotTab$R_value <- as.numeric(plotTab$R_value)
  plotTab$P_value <- as.numeric(plotTab$P_value)
  

  corrPlot <- ggplot(data = plotTab, aes(x = R_value, y = -log10(P_value), label = Genes)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel(size = 3, max.overlaps = 20) +
    # geom_vline(xintercept=c(cutOff), col="black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), col = "black", linetype = "dashed") +
    ylab("-log10(p-value)") +
    xlab("Rho Value") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 12, vjust = 2, face = "bold"),
      axis.title.y = element_text(size = 12, vjust = 2, face = "bold")
    ) +
    guides(col = guide_legend(nrow = 2, byrow = TRUE))
  
  
  # corrPlot
  res <- list(result = result, corrPlot = corrPlot, usedCohorts = paste0(unique(cptac_samples_filt$cohort), collapse = "_"))
  return(res)
}


panelCorr <- function(geneNamesRNA, geneNamesProtein, tissueType, cohort, cptac_protein, cptac_rna, cptac_samples) {
  if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    # cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohor %in% cohort, ]
    
    tumNames <- as.character(
      intersect(
        cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id,
        intersect(colnames(cptac_protein), colnames(cptac_rna))
      )
    )
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal" | cptac_samples_filt$sample_type == "G", ]$harmonized_data_id, intersect(colnames(cptac_protein), colnames(cptac_rna))))
    
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
    selectedGenesProt <- selectedGenesProt[intersect(rownames(selectedGenesProt), rownames(selectedGenesRna)), intersect(colnames(selectedGenesProt), colnames(selectedGenesRna))]
    selectedGenesRna <- selectedGenesRna[intersect(rownames(selectedGenesProt), rownames(selectedGenesRna)), intersect(colnames(selectedGenesProt), colnames(selectedGenesRna))]
    
    if (ncol(selectedGenesProt) < 10 | ncol(selectedGenesRna) < 10) {
      stop("One or more genes have missing values, therefore decreasing the number of samples available for analysis. Please, adjust input parameters (e.g., add cohorts and/or additional genes).")
    }
    
    
    cor.mat <- cor(x = t(selectedGenesProt), y = t(selectedGenesRna), use = "pairwise.complete.obs", method = "spearman")
    p.mat <- cor.mtest(cor.mat, method = "spearman")
    
    
    # asd <- outer(t(na.omit(selectedGenesProt)), t(selectedGenesRna), Vectorize(function(a, b) cor.test(a, b)$p.value))
    
    flattenCorrMatrix <- function(cormat, pmat) {
      ut <- cormat
      data.frame(
        row = rownames(cormat)[row(cormat)],
        column = rownames(cormat)[col(cormat)],
        cor = c(cormat),
        p = c(pmat)
      )
    }
    # testing p values
    # table(is.na(selectedGenesProt[rownames(selectedGenesProt) == "AARS2",]))
    # testCorr <- cor.test(as.numeric(selectedGenesProt[rownames(selectedGenesProt) == "AARS2",]),
    #                      as.numeric(selectedGenesRna[rownames(selectedGenesRna) == "A1BG",]),method = "spearman")
    
    
    df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p + 0.0000001)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1, 2)]
    
    corrPlot <- ggballoonplot(
      data = df,
      x = "row",
      y = "column",
      size = "log_pval",
      facet.by = NULL,
      size.range = c(1, 10),
      shape = 21,
      color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
      color.range = c(0, 5),
      fill = "cor",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE
    ) +
      gradient_fill(c("purple", "white", "red"))
    
    # corrPlot
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
    cptac_samples_filt <- cptac_samples[cptac_samples$cohor %in% cohort, ]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Tumor", ]$harmonized_data_id, colnames(inPut)))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$type == "Normal" | cptac_samples_filt$sample_type == "G", ]$harmonized_data_id, colnames(inPut)))
    
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
    
    # corrPlot <- corrplot(cor.mat, p.mat = p.mat[[1]], sig.level = 0.01, tl.col = "black",type="lower",diag = FALSE)
    flattenCorrMatrix <- function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor = (cormat)[ut],
        p = pmat[ut]
      )
    }
    df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p + 0.0000001)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1, 2)]
    
    
    corrPlot <- ggballoonplot(
      data = df,
      x = "row", # pride
      y = "column", # cptac_n
      size = "log_pval",
      facet.by = NULL,
      size.range = c(1, 10),
      shape = 21,
      color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
      color.range = c(0, 5),
      fill = "cor",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE
    ) +
      theme_minimal() +
      
      # scale_fill_gradient(low = "purple",high = "red") +
      gradient_fill(c("purple", "white", "red"))
    
    # corrPlot
  }
  
  res <- list(corrPlot = corrPlot, corTab = list(cor = cor.mat, pValue = p.mat$p), usedCohorts = paste0(unique(cptac_samples_filt$cohort), collapse = "_"))
  return(res)
}
#
# asd <- panelCorr(geneNamesRNA = c("A1BG", "A2M", "AAAS", "AACS", "AARD"),
#                  geneNamesProtein = c("A1BG", "A2M", "AAAS", "AACS", "AARD"),
#                  tissueType = "Tumor", cohort = "COAD", cptac_protein = cptac_protein, cptac_rna = cptac_protein, cptac_samples = cptac_samples)
# asd$usedCohorts
