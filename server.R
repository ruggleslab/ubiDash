packagelist <- c("shiny", "DT", "ComplexHeatmap", "openxlsx", "readxl", "tidyverse", "glue", "stringr", "reshape2", "RColorBrewer", "plotly", "maftools", "circlize", "ggrepel", "corrplot", "ggpubr")

lapply(packagelist, function(x) {
  suppressMessages(require(x,
    character.only = T,
    quietly = T,
    warn.conflicts = F
  ))
})

load("./data/defaults.RData")
# loading default objects below:
# genesSubset
# defaultHeatmap
# cohortsOnco
# mutationsOnco

source("./scripts/SplitInputGenes.R")


# server ----

shinyServer(function(input, output, session) {
  # heatmap tab ----
  heatmapServer(
    "heatmapTab", cohorts, cptac_protein, cptac_rna, tnlabels,
    colorsList, filteredProteinMat, filteredMRNAMat,
    filteredTnlabels, defaultHeatmap, genesSubset
  )

  # volcano plot tab ----
  volcanoServer("volcanoTab", cptac_protein, crl_genes)

  # correlation tab (----
  correlationsServer("corrTab", plot_tab, lines)

  # boxplot tab ----
  violinServer("violinTab", ccle_rna, violin.input, POIS)

  # oncoplot tab ----
  oncoplotServer("oncoTab", cptac_protein, cptac_samples, maf, crl_genes, top_crls, cohorts)

  # lollipop tab
  lollipopServer("lollipopTab", maf)

  corrAnalysisServer("corrAnalysisTab", cohorts, cptac_protein, cptac_rna, cptac_samples)

  drugResponseServer("drugResponseTab", drug, drugpath, exp, cell, crl_genes)

  simplifiedHeatmapServer("simplifiedHeatmapTab", prt.n, tertiles, mTert, my.colors)
})
