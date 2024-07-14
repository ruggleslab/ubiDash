generateCorrPlot <- function(threshold, plot_tab, xaxis, yaxis, inputGenes, warning) {
  warning("")

  inputGenes <- splitInputGenes(inputGenes)

  plot_tab <- filter(plot_tab, abs(plot_tab[, xaxis]) > threshold &
    abs(plot_tab[, yaxis]) > threshold)

  if (length(inputGenes) > 0) {
    plot_tab <- filter(plot_tab, toupper(row) %in% inputGenes)
    if (length(unique(plot_tab$row)) != length(inputGenes)) {
      warning("Data unavailable for some genes")
    }
  }

  text <- paste(
    "<b>", plot_tab$gene_pairs, "</b> \n",
    "pr_normal.cor: ", round(plot_tab$pr_normal.cor, 2), "\n",
    "pr_tumor.cor: ", round(plot_tab$pr_tumor.cor, 2), "\n",
    "rn_tumor.cor: ", round(plot_tab$rn_tumor.cor, 2), "\n",
    "rn_normal.cor: ", round(plot_tab$rn_normal.cor, 2), "\n"
  )

  lines <- list(
    geom_hline(yintercept = threshold, color = "gray"),
    geom_hline(yintercept = -threshold, color = "gray"),
    geom_vline(xintercept = -threshold, color = "gray"),
    geom_vline(xintercept = threshold, color = "gray"),
    geom_segment(aes(x = threshold, xend = 1.3, y = threshold, yend = threshold), color = "red"),
    geom_segment(aes(x = threshold, xend = threshold, y = threshold, yend = 1.3), color = "red"),
    geom_segment(aes(x = -threshold, xend = -1.3, y = -threshold, yend = -threshold), color = "blue"),
    geom_segment(aes(x = -threshold, xend = -threshold, y = -threshold, yend = -1.3), color = "blue"),
    geom_segment(aes(x = -threshold, xend = -1.3, y = threshold, yend = threshold), color = "purple"),
    geom_segment(aes(x = -threshold, xend = -threshold, y = threshold, yend = 1.3), color = "purple"),
    geom_segment(aes(x = threshold, xend = 1.3, y = -threshold, yend = -threshold), color = "green"),
    geom_segment(aes(x = threshold, xend = threshold, y = -threshold, yend = -1.3), color = "green"),
    geom_hline(yintercept = 0, size = 1, color = "black"),
    geom_vline(xintercept = -0, size = 1, color = "black")
  )

  return(suppressWarnings(
    ggplot(plot_tab, aes(y = plot_tab[, yaxis], x = plot_tab[, xaxis])) +
      ylab(yaxis) +
      xlab(xaxis) +
      geom_point(aes(text = text), alpha = 0.35) +
      lines +
      coord_cartesian(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05))
  ))
}
