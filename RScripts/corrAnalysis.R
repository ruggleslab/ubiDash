# geneName <- "TNNC1"
# tissueType <- "Tumor"
# cutOff <- -0.01
# 
# cptac_protein <- readRDS("./data/cptac_tmt-log2-RI_protein_expression.rds")
# cptac_samples <- readRDS("./data/cptac_sample_info.rds")


multiCorr <- function(geneName, tissueType, cutOff,cptac_protein, cptac_samples){
  #opening clinical table
  # load(file = "data/cptac_sample_info.RData")
  # cptac_samples <- samples
  # rm(samples) 
  
  #opening expression table
  #load(file = "data/cptac_tmt_protein_expression.RData")
  #cptac_protein <- mat
  #rm(mat)
  
  tumNames <- as.character(intersect(cptac_samples[cptac_samples$sample_type == "T",5], colnames(cptac_protein)))
  normNames <- as.character(intersect(cptac_samples[cptac_samples$sample_type == "N" | cptac_samples$sample_type == "G" ,5], colnames(cptac_protein)))
  #cptac_samples$sample_type %in% c(“N”, “G”)
  
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
    
  }else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  #make the selected gene to a numeric vector and perform correlation test
  selectedGene <- cptac_protein[rownames(cptac_protein) == geneName,colnames(cptac_protein) %in% usedTissue]
  selectedGene <- selectedGene[,!is.na(selectedGene)]
  rawForAn <- cptac_protein[,colnames(cptac_protein) %in% colnames(selectedGene)]
  rawForAn <- as.data.frame(na.omit(rawForAn)) 
  selectedGene <- as.numeric(selectedGene)
  
  corResult <- apply(rawForAn, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4,3)])
  
  #store the results as a dataframe and convert the list to a dataframe
  result = as.data.frame(matrix(nrow = nrow(rawForAn), ncol = 3))
  colnames(result) <- c("Protein", "Spearman's Rho", "p-value")
  
  result[, 1] = names(corResult)
  result[, 2] = round(x = sapply(corResult, "[[", 1), digits = 2)
  result[, 3] = format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
  
  #order the results based on R and remove the values below cutoff
  if ((cutOff < 0) == TRUE) {
    result <- result[result$`Spearman's Rho` < 0,]
    result <- subset(result, result$`Spearman's Rho` <= cutOff)
    result <- result[order(result$`Spearman's Rho` ,decreasing = F),]
    #result$`Spearman's Rho` <- abs(result$`Spearman's Rho`)
    
  } else {
    result <- subset(result, result$`Spearman's Rho` >= cutOff)
    result <- result[order(result$`Spearman's Rho` ,decreasing = T),]
  }
  
  library(ggplot2)
  library(ggrepel)
  
  plotTab <- result
  plotTab$From <- geneName
  plotTab <- plotTab[,c(4,1,2,3)]
  colnames(plotTab) <-c("From", "Genes", "R_value", "P_value")
  plotTab$R_value <- as.numeric(plotTab$R_value)
  plotTab$P_value <- as.numeric(plotTab$P_value)  
    
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
  
 
  #corrPlot
  res <- list(result = result, corrPlot = corrPlot)
  return(res)
}
