packagelist <- c("shiny", "DT", "ggprism", "ComplexHeatmap", "openxlsx", "readxl", "tidyverse", "glue", "stringr", "reshape2", "RColorBrewer", "plotly", "maftools", "circlize", "ggrepel", "corrplot", "ggpubr")

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
    "heatmapTab", cohorts,
    # cptac_protein,
    # cptac_rna,
    # ccle_protein,
    # ccle_rna,
    # pride_protein,
    tnlabels,
    colorsList,
    filteredProteinMat,
    filteredMRNAMat,
    # filteredTnlabels,
    defaultHeatmap,
    genesSubset
  )

  # volcano plot tab ----
  volcanoServer("volcanoTab", crl_genes, DE_protein) # cptac_protein

  # correlation tab ----
  correlationsServer("corrTab", lines) # plot_tab

  # boxplot tab ----
  violinServer("violinTab",cptac_samples,  cohort_colors) # cptac_protein, cptac_rna

  # oncoplot tab ----
  oncoplotServer("oncoTab", cptac_samples, crl_genes, top_crls, cohorts) # cptac_protein, maf

  # lollipop tab
  lollipopServer("lollipopTab") # maf

  corrAnalysisServer("corrAnalysisTab", cohorts, cptac_samples) # cptac_protein, cptac_rna

  drugResponseServer("drugResponseTab", drug, drugpath, exp, cell, crl_genes)

  simplifiedHeatmapServer("simplifiedHeatmapTab",  tertiles,prt.global.AI, labels, my.colors)
})
