generateCorrPlot = function(threshold, plot_tab,yaxis, inputGenes, warning){
  
  warning("")
  
  inputGenes = gsub(",", " ", inputGenes)
  inputGenes = gsub("  ", " ", inputGenes)
  inputGenes = unlist(strsplit(inputGenes, " "))
  
  plot_tab = filter(plot_tab, abs(pr_normal.cor) > threshold & abs(pr_tumor.cor) > threshold & abs(rn_tumor.cor) > threshold)
 
  if(length(inputGenes) > 0){
    plot_tab = filter(plot_tab,row %in% inputGenes)
    if(length(unique(plot_tab$row)) != length(inputGenes)){
      warning("Data unavailable for some genes")
    }
  }
  
   text = paste("<b>",plot_tab$gene_pairs, "</b> \n",
               "pr_normal.cor: ", plot_tab$pr_normal.cor, "\n",
               "pr_normal.p: ", plot_tab$pr_normal.p, "\n",
               "pr_tumor.cor: ", plot_tab$pr_tumor.cor, "\n",
               "pr_tumor.p: ", plot_tab$pr_tumor.p, "\n",
               "rn_tumor.cor: ", plot_tab$rn_tumor.cor, "\n",
               "rn_tumor.p: ", plot_tab$rn_tumor.p, "\n")
  
  lines = list(geom_hline(yintercept=threshold, color = "red") ,
               geom_hline(yintercept=-threshold,  color = "red") ,
               geom_vline(xintercept=-threshold,color = "red") ,
               geom_vline(xintercept=threshold, color = "red") ,
               geom_hline(yintercept=0, size = 1, color = "black") ,
               geom_vline(xintercept=-0, size = 1, color = "black"))
  
  return(suppressWarnings(
    ggplot(plot_tab, aes(y = plot_tab[,yaxis], x = pr_tumor.cor)) + 
      ylab(yaxis) +
      geom_point(aes(text = text)) +
      lines
  ))
}