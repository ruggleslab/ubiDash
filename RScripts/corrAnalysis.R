# library(corrplot)
# cohort <- c("BRCA","COAD", "CCRCC")
# tissueType <- "Tumor"
# cutOff <- c(-0.05,0.05)
# 
# cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
# cptac_rna <- readRDS("./data/cptac_tpm_rna_expression.rds")
# #cptac_rna <- NULL
# cptac_samples <- readRDS("./data/cptac_sample_info.rds")
# #panel
# geneNamesRNA <- rownames(cptac_rna)[1:13]
# geneNamesProtein <- rownames(cptac_protein)[1:13]
# 
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

singleCorr <- function(geneNameA, geneNameB, tissueType, cohort, cptac_protein, cptac_rna, cptac_samples){
  
  if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    #cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort_code %in% cohort,]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "T",5], intersect(colnames(cptac_protein),colnames(cptac_rna))))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "N" | cptac_samples_filt$sample_type == "G" ,5], intersect(colnames(cptac_protein),colnames(cptac_rna))))
    
    if (tissueType == "Tumor") {
      usedTissue <- tumNames
      
    }else if (tissueType == "Normal") {
      usedTissue <- normNames
    }
    cptac_samples_filt <- cptac_samples_filt[rownames(cptac_samples_filt) %in% usedTissue,]
    
    #Filter for genes
    selectedGeneA <- cptac_protein[rownames(cptac_protein) == geneNameA,colnames(cptac_protein) %in% usedTissue]
    selectedGeneB <- cptac_rna[rownames(cptac_rna) == geneNameB,colnames(cptac_rna) %in% usedTissue]
    #remove NAs
    selectedGeneA <- selectedGeneA[,!is.na(selectedGeneA)]
    selectedGeneB <- selectedGeneB[,!is.na(selectedGeneB)]
    #keep common samples
    selectedGeneA <- selectedGeneA[,intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    selectedGeneB <- selectedGeneB[,intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    
  } else if (!is.null(cptac_protein) | !is.null(cptac_rna)) {#for RNA or Protein separately
    varToAssign <- list(cptac_protein,cptac_rna)
    inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein),!is.null(cptac_rna)))])
    rm(varToAssign)
    
    colnames(inPut) <- gsub(pattern = ".?^X",replacement = "",x = colnames(inPut))
    colnames(inPut) <- gsub(pattern = "\\.", replacement = "-",x = colnames(inPut))
    #cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort_code %in% cohort,]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "T",5], colnames(inPut)))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "N" | cptac_samples_filt$sample_type == "G" ,5], colnames(inPut)))
    
   
    if (tissueType == "Tumor") {
      usedTissue <- tumNames
      
    }else if (tissueType == "Normal") {
      usedTissue <- normNames
    }
    
    cptac_samples_filt <- cptac_samples_filt[rownames(cptac_samples_filt) %in% usedTissue,]
    
    #Filter for genes
    selectedGeneA <- inPut[rownames(inPut) == geneNameA,colnames(inPut) %in% usedTissue]
    selectedGeneB <- inPut[rownames(inPut) == geneNameB,colnames(inPut) %in% usedTissue]
    #remove NAs
    selectedGeneA <- selectedGeneA[,!is.na(selectedGeneA)]
    selectedGeneB <- selectedGeneB[,!is.na(selectedGeneB)]
    #keep common samples
    selectedGeneA <- selectedGeneA[,intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    selectedGeneB <- selectedGeneB[,intersect(colnames(selectedGeneA), colnames(selectedGeneB))]
    
  }
  
  selectedGeneA <- as.data.frame(t(selectedGeneA))
  selectedGeneB <- as.data.frame(t(selectedGeneB))
  
  selectedGeneA$Cohort <- cptac_samples_filt[rownames(cptac_samples_filt) %in% rownames(selectedGeneA),1]
  selectedGeneB$Cohort <- cptac_samples_filt[rownames(cptac_samples_filt) %in% rownames(selectedGeneB),1]
  colnames(selectedGeneA) <- c("ExpValueA", "Cohort")
  colnames(selectedGeneB) <- c("ExpValueB", "Cohort")
  
  corrSum <- cor.test(x = as.numeric(selectedGeneA$ExpValueA),y = as.numeric(selectedGeneB$ExpValueB), method = "spearman")
  
  #visualization
  if (is.null(cptac_rna)) {
    dat <- data.frame("geneAExp" = as.numeric(selectedGeneA$ExpValueA), "geneBExp" = as.numeric(selectedGeneB$ExpValueB))#Protein
  } else if (is.null(cptac_protein)) {
    dat <- data.frame("geneAExp" = log2(as.numeric(selectedGeneA$ExpValueA)), "geneBExp" = log2(as.numeric(selectedGeneB$ExpValueB)))#RNA
  } else if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    dat <- data.frame("geneAExp" = as.numeric(selectedGeneA$ExpValueA), "geneBExp" = log2(as.numeric(selectedGeneB$ExpValueB)))#Both
  } 
  
  if (identical(rownames(selectedGeneA), rownames(selectedGeneB))) {
    dat$Cohort <- selectedGeneA$Cohort
  }
   
  
  
  # singCorrTab <- data.frame("N" = nrow(dat),
  #                           "Spearman's Rho" = round(x = corrSum$estimate,digits = 2),
  #                           "pValue" = format(x = corrSum$p.value,digits = 2, scientific = TRUE), row.names = NULL)
   
  pValue <- format(x = corrSum$p.value,digits = 2, scientific = TRUE)
  rho <- round(x = corrSum$estimate,digits = 2)
  singCorrTab <- t(data.frame("Cohorts" = names(table(dat$Cohort)),
                               "N" = as.numeric(unlist(table(dat$Cohort)))))

  xmax <- max(dat$geneAExp)
  ymax <- max(dat$geneBExp)
  
  library(ggplot2)
  scatPlot <- ggplot(dat, aes(x = geneAExp, y = geneBExp, color = Cohort)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", size = 2) +
    #geom_smooth(aes(x = geneAExp, y = geneBExp,color = NULL),method = "lm", size = 2) +
    xlab(paste0("Log2 expression of \n", geneNameA)) + 
    ylab(paste0("Log2 expression of \n", geneNameB)) + 
    annotate(geom = "text",x = xmax,y = ymax, label = paste0( "p = ",pValue, "\n", "R = ",rho))+
    #coord_cartesian(clip = "off")+
    theme(legend.position = "right", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2),
          axis.title.x = element_text(size = 26, face = "bold",  vjust = 2))
  scatPlot
  
  res <- list(singCorrPlot = scatPlot, singCorrTab = singCorrTab)
  return(res)
  
}


multiCorr <- function(geneName, tissueType, cutOff, cohort, cptac_protein, cptac_rna, cptac_samples){
  varToAssign <- list(cptac_protein,cptac_rna)
  inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein),!is.null(cptac_rna)))])
  rm(varToAssign)
  
  colnames(inPut) <- gsub(pattern = ".?^X",replacement = "",x = colnames(inPut))
  colnames(inPut) <- gsub(pattern = "\\.", replacement = "-",x = colnames(inPut))
  
  #cohort filter
  cptac_samples_filt <- cptac_samples[cptac_samples$cohort_code %in% cohort,]
  
  tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "T",5], colnames(inPut)))
  normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "N" | cptac_samples_filt$sample_type == "G" ,5], colnames(inPut)))
  
  
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
    
  }else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  #make the selected gene to a numeric vector and perform correlation test
  selectedGene <- inPut[rownames(inPut) == geneName,colnames(inPut) %in% usedTissue]
  selectedGene <- selectedGene[,!is.na(selectedGene)]
  rawForAn <- inPut[,colnames(inPut) %in% colnames(selectedGene)]
  rawForAn <- as.data.frame(na.omit(rawForAn)) 
  selectedGene <- as.numeric(selectedGene)
  #correlation analysis
  corResult <- apply(rawForAn, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4,3)])
  
  #store the results as a dataframe and convert the list to a dataframe
  result = as.data.frame(matrix(nrow = nrow(rawForAn), ncol = 3))
  colnames(result) <- c("Protein", "Spearman's Rho", "p-value")
  
  result[, 1] = names(corResult)
  result[, 2] = round(x = sapply(corResult, "[[", 1), digits = 2)
  result[, 3] = format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
  
  result <- result[result$`Spearman's Rho` >= min(cutOff) & result$`Spearman's Rho` <=max(cutOff),]
  result <- result[order(result$`Spearman's Rho` ,decreasing = T),]
  
  plotTab <- result
  plotTab$From <- geneName
  plotTab <- plotTab[,c(4,1,2,3)]
  colnames(plotTab) <-c("From", "Genes", "R_value", "P_value")
  plotTab$R_value <- as.numeric(plotTab$R_value)
  plotTab$P_value <- as.numeric(plotTab$P_value)  
  
  library(ggrepel)
  library(ggplot2)
  corrPlot <- ggplot(data=plotTab, aes(x=R_value, y=-log10(P_value), label = Genes)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(size = 3, max.overlaps = 20) +
    #geom_vline(xintercept=c(cutOff), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed") +
    ylab("-log10(p-value)")+
    xlab("Rho Value")+
    theme(legend.title =  element_blank(), 
          legend.text = element_text(size=25),
          legend.position="bottom",
          axis.text.x = element_text(size = 25, colour = "black"),
          axis.text.y = element_text(size = 25, colour = "black"),
          axis.title.x = element_text(size = 25,  vjust = 2, face = "bold"),
          axis.title.y = element_text(size = 25,  vjust = 2, face = "bold"))+
    guides(col=guide_legend(nrow=2,byrow=TRUE))
  
 
  corrPlot
  res <- list(result = result, corrPlot = corrPlot)
  return(res)
}


panelCorr <- function(geneNamesRNA,geneNamesProtein, tissueType, cohort, cptac_protein, cptac_rna, cptac_samples){

  if (!is.null(cptac_protein) & !is.null(cptac_rna)) {
    #cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort_code %in% cohort,]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "T",5], intersect(colnames(cptac_protein),colnames(cptac_rna))))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "N" | cptac_samples_filt$sample_type == "G" ,5], intersect(colnames(cptac_protein),colnames(cptac_rna))))

    if (tissueType == "Tumor") {
      usedTissue <- tumNames

    }else if (tissueType == "Normal") {
      usedTissue <- normNames
    }
    
    #Filter for genes
    selectedGenesProt <- cptac_protein[rownames(cptac_protein) %in% geneNamesProtein,colnames(cptac_protein) %in% usedTissue]
    selectedGenesRna <- cptac_rna[rownames(cptac_rna) %in% geneNamesRNA,colnames(cptac_rna) %in% usedTissue]
    
    #keep common samples and genes
    selectedGenesProt <- selectedGenesProt[intersect(rownames(selectedGenesProt), rownames(selectedGenesRna)),intersect(colnames(selectedGenesProt), colnames(selectedGenesRna))]
    selectedGenesRna <- selectedGenesRna[intersect(rownames(selectedGenesProt), rownames(selectedGenesRna)),intersect(colnames(selectedGenesProt), colnames(selectedGenesRna))]
    library(corrplot)
    cor.mat <- cor(t(selectedGenesProt), t(selectedGenesRna), use = "pairwise.complete.obs", method = "spearman")
    p.mat <- cor.mtest(mat = cor.mat, method = "spearman")
    
    #corrPlot <- corrplot(cor.mat, p.mat = p.mat[[1]], sig.level = 0.01, tl.col = "black")
    #title(ylab = substitute(paste(bold("RNA"))), main = "Protein",cex.lab = 2, cex.main = 2)
    flattenCorrMatrix <- function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
      )
    }
    
    df <- flattenCorrMatrix(cormat = cor.mat,pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p)
    df$log_pval_test <- -log10(df$p + 0.000001)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1,2)]
    
    library(ggpubr)
    corrPlot <- ggballoonplot(
      data =df,
      x = "row", 
      y = "column",
      size = "cor",
      facet.by = NULL,
      size.range = c(1,10),
      shape = 21,
      color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
      color.range = c(0,5),
      fill = "log_pval",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE) +
      gradient_fill(c("purple","white","red"))
      
    #corrPlot
    
  
    
  } else if (!is.null(cptac_protein) | !is.null(cptac_rna)) {#for RNA or Protein separately
    varToAssign <- list(cptac_protein,cptac_rna)
    inPut <- as.data.frame(varToAssign[which(c(!is.null(cptac_protein),!is.null(cptac_rna)))])
    rm(varToAssign)
    geneName <- c(geneNamesProtein, geneNamesRNA)
    
    colnames(inPut) <- gsub(pattern = ".?^X",replacement = "",x = colnames(inPut))
    colnames(inPut) <- gsub(pattern = "\\.", replacement = "-",x = colnames(inPut))

    #cohort filter
    cptac_samples_filt <- cptac_samples[cptac_samples$cohort_code %in% cohort,]
    
    tumNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "T",5], colnames(inPut)))
    normNames <- as.character(intersect(cptac_samples_filt[cptac_samples_filt$sample_type == "N" | cptac_samples_filt$sample_type == "G" ,5], colnames(inPut)))
    
    if (tissueType == "Tumor") {
      usedTissue <- tumNames
      
    }else if (tissueType == "Normal") {
      usedTissue <- normNames
    }

    selectedGenes <- inPut[rownames(inPut) %in% geneName,colnames(inPut) %in% usedTissue]
    
    cor.mat <- cor(t(selectedGenes), use = "pairwise.complete.obs", method = "spearman")
    p.mat <- cor.mtest(mat = cor.mat, method = "spearman")
    
    #corrPlot <- corrplot(cor.mat, p.mat = p.mat[[1]], sig.level = 0.01, tl.col = "black",type="lower",diag = FALSE)
    flattenCorrMatrix <- function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
      )
    }
    df <- flattenCorrMatrix(cormat = cor.mat,pmat = p.mat[[1]])
    df$log_pval <- -log10(df$p)
    df$log_pval_test <- -log10(df$p + 0.000001)
    df$qval <- p.adjust(df$p, method = "BH")
    df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1,2)]
    
    library(ggpubr)
    corrPlot <- ggballoonplot(
      data =df,
      x = "row", #pride
      y = "column", #cptac_n
      size = "cor",
      facet.by = NULL,
      size.range = c(1, 10),
      shape = 21,
      color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
      color.range = c(0,5),
      fill = "log_pval",
      font.label = list(size = 12, color = "black"),
      rotate.x.text = TRUE) +
      #scale_fill_gradient(low = "purple",high = "red") +
      gradient_fill(c("purple","white","red"))
    
    #corrPlot
    
    
    

  }
  
  res <- list(corrPlot = corrPlot, corTab = list(cor = cor.mat, pValue = p.mat$p))
  return(res)

}
