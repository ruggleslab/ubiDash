load("./data/defaults.RData")
# loading default objects below:
# genesSubset
# defaultHeatmap
# cohortsOnco
# mutationsOnco

# source("./scripts/SplitInputGenes.R")


source('./R/corranalysis.R')

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
    colorsTbl,
    filteredProteinMat,
    filteredMRNAMat,
    # filteredTnlabels,
    # defaultHeatmap,
    genesSubset
  )

  # volcano plot tab ----
  volcanoServer("volcanoTab", crl_genes, DE_protein) # cptac_protein

  # correlation tab ----
  # correlationsServer("corrTab", lines, E3_genes) # plot_tab

  # boxplot tab ----
  violinServer("violinTab") # cptac_protein, cptac_rna, cptac_samples

  # oncoplot tab ----
  oncoplotServer("oncoTab", cptac_samples, top_crls, cohorts, colorsTbl) # cptac_protein, maf

  # lollipop tab
  lollipopServer("lollipopTab") # maf
  
  QTLServer("QTLTab", genelists, colorsTbl)
  
    densityServer("densityTab") # maf

  corrAnalysisServer("corrAnalysisTab", cohorts, cptac_samples) # cptac_protein, cptac_rna

  drugResponseServer("drugResponseTab",
                     drug,
                     drugpath,
                     exp,
                     cell,
                     crl_genes)

  simplifiedHeatmapServer("simplifiedHeatmapTab",
                          tertiles,
                          # prt.global.AI,
                          labels,
                          my.colors)
  
  dataDownloadServer("downloadTab")
  
})
