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
    "heatmapTab",
    cohorts,
    cptac_protein,
    cptac_rna,
    ccle_protein,
    ccle_rna,
    pride_protein,
    tnlabels,
    colorsList,
    filteredProteinMat,
    filteredMRNAMat,
    filteredTnlabels,
    defaultHeatmap,
    genesSubset
  )

  # volcano plot tab ----
  volcanoServer("volcanoTab",
                cptac_protein,
                crl_genes)

  # correlation tab ----
  correlationsServer("corrTab",
                     plot_tab,
                     lines)

  # boxplot tab ----
  violinServer("violinTab",
               cptac_rna,
               cptac_protein,
               cptac_samples,
               cohort_colors)

  # oncoplot tab ----
  oncoplotServer("oncoTab",
                 cptac_protein,
                 cptac_samples,
                 maf,
                 crl_genes,
                 top_crls,
                 cohorts)

  # lollipop tab
  lollipopServer("lollipopTab",
                 maf)

  corrAnalysisServer("corrAnalysisTab",
                     cohorts,
                     cptac_protein,
                     cptac_rna,
                     cptac_samples)

  drugResponseServer("drugResponseTab",
                     drug,
                     drugpath,
                     exp,
                     cell,
                     crl_genes)

  simplifiedHeatmapServer("simplifiedHeatmapTab",
                          tertiles,
                          prt.global.AI,
                          labels,
                          my.colors)
})
