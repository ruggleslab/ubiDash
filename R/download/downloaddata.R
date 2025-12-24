dataDownloadUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Download Data"),
    p(tab_summaries$`Download Data`),
    tags$div(
      class = "download-grid",
      tags$div(
        class = "download-block",
        h4("CPTAC Protein Expression"),
        tags$small(a("PMID: 37582339", href = "https://pubmed.ncbi.nlm.nih.gov/37582339/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_cptac_tmt"), "CPTAC_TMT_Protein_Expression.rds")),
          tags$li(downloadLink(ns("dl_mat_prot"), "CPTAC_Protein_Matrix.rds")),
          tags$li(downloadLink(ns("dl_ubidash_proteome_clean"), "CPTAC_PanCan_Proteome_Clean.rds")),
          tags$li(downloadLink(ns("dl_cptac_proteome_rds"), "CPTAC_proteome.rds")),
          tags$li(downloadLink(ns("dl_cptac_proteome_scaled_rds"), "CPTAC_proteome_scaled.rds (pre-processed)"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("CPTAC Sample Metadata"),
        tags$small(a("PMID: 37582339", href = "https://pubmed.ncbi.nlm.nih.gov/37582339/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_cptac_samples"), "CPTAC_Sample_Info.rds")),
          tags$li(downloadLink(ns("dl_cptac_plus_samples"), "CPTAC_Extended_Sample_Info.rds")),
          tags$li(downloadLink(ns("dl_ubidash_sample_mapping"), "CPTAC_Sample_ID_Mapping.rds"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("CPTAC Mutations & QTL"),
        tags$small(a("PMID: 37582339", href = "https://pubmed.ncbi.nlm.nih.gov/37582339/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_pancan_maf"), "CPTAC_PanCan_Mutations.rds")),
          tags$li(downloadLink(ns("dl_cptac_pqtls"), "CPTAC_pQTLs.rds"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("CPTAC Differential Expression"),
        tags$small(a("PMID: 37582339", href = "https://pubmed.ncbi.nlm.nih.gov/37582339/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_deqms_percancer"), "Table S2A. DEqMS_Percancer_CPTAC_imputed_proteome_matchedTN.xlsx")),
          tags$li(downloadLink(ns("dl_deqms_cancer_subtypes"), "DEqMS_Cancer_Subtypes.csv")),
          tags$li(downloadLink(ns("dl_deqms_immune"), "DEqMS_Immune_Subtypes.csv")),
          tags$li(downloadLink(ns("dl_deqms_matrix"), "DEqMS_Signed_Padj_Matrix.csv"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("CCLE Protein Expression"),
        tags$small(a("PMID: 35839778", href = "https://pubmed.ncbi.nlm.nih.gov/35839778/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_ccle_lfq_protein"), "CCLE_Protein_Expression.rds")),
          tags$li(downloadLink(ns("dl_ccle_lfq_sample"), "CCLE_Sample_Info.rds"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("Drug Sensitivity (CCLE)"),
        tags$small(a("DepMap Portal", href = "https://depmap.org/portal/data_page/?tab=allData", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_drug_longdf"), "Drug_Sensitivity_Data.csv")),
          tags$li(downloadLink(ns("dl_drugpath"), "Drug_Pathways.csv")),
          tags$li(downloadLink(ns("dl_cell_line_exp"), "Cell_Line_Basal_Expression.rds")),
          tags$li(downloadLink(ns("dl_matched_cell"), "Matched_Cell_Annotations.csv")),
          tags$li(downloadLink(ns("dl_repurposing"), "Drug_Repurposing_Screen_24Q2.rds")),
          tags$li(downloadLink(ns("dl_mutations_matrix"), "Somatic_Mutations_Hotspots.rds"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("Harmonized & Correlation Data"),
        tags$small("Multi-cohort harmonized datasets and covariance matrices"),
        tags$ul(
          tags$li(downloadLink(ns("dl_harmonized_proteome"), "Harmonized_Proteome_Data.rds")),
          tags$li(downloadLink(ns("dl_ubidash_sample_mapping2"), "Harmonized_Sample_Mapping.rds")),
          tags$li(downloadLink(ns("dl_ups_df"), "UPS_Protein_Annotations.rds")),
          tags$li(downloadLink(ns("dl_pancan_covariance"), "PanCan_Protein_Correlations.rds"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("LSCC Ubiquitylation"),
        tags$small(a("PMID: 34534465", href = "https://pubmed.ncbi.nlm.nih.gov/34534465/", target = "_blank")),
        tags$ul(
          tags$li(downloadLink(ns("dl_lscc_prot"), "LSCC_Protein_Expression.csv")),
          tags$li(downloadLink(ns("dl_lscc_ub_prot"), "LSCC_Ubiquitylation_Sites.csv"))
        )
      ),
      tags$div(
        class = "download-block",
        h4("Reference & Gene Lists"),
        tags$small("Supporting annotations and gene lists used throughout the app"),
        tags$ul(
          tags$li(downloadLink(ns("dl_biogrid"), "BioGRID_Interactions.rds")),
          tags$li(downloadLink(ns("dl_combined_genelists"), "Combined_Gene_Lists.rds")),
          tags$li(downloadLink(ns("dl_e3_ligase"), "E3_Ligase_List.rds")),
          tags$li(downloadLink(ns("dl_ups_network"), "Human_Proteostasis_Network.rds")),
          tags$li(downloadLink(ns("dl_colors"), "Color_Scheme.csv"))
        )
      )
    )
  )
}

dataDownloadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # CPTAC Protein Expression
    output$dl_cptac_tmt <- downloadHandler(
      filename = "CPTAC_TMT_Protein_Expression.rds",
      content = function(file) {
        file.copy("./data/cptac_tmt-log2-AI-v1.1_protein_expression.rds", file)
      }
    )

    output$dl_mat_prot <- downloadHandler(
      filename = "CPTAC_Protein_Matrix.rds",
      content = function(file) {
        file.copy("./data/mat_prot.rds", file)
      }
    )

    output$dl_ubidash_proteome_clean <- downloadHandler(
      filename = "CPTAC_PanCan_Proteome_Clean.rds",
      content = function(file) {
        file.copy("./data/UbiDash_CPTAC_PanCan_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed_clean.rds", file)
      }
    )


    output$dl_cptac_proteome_rds <- downloadHandler(
      filename = "CPTAC_proteome.rds",
      content = function(file) {
        file.copy("./data/CPTAC_proteome.rds", file)
      }
    )

    output$dl_cptac_proteome_scaled_rds <- downloadHandler(
      filename = "CPTAC_proteome_scaled.rds",
      content = function(file) {
        file.copy("./data/CPTAC_proteome_scaled.rds", file)
      }
    )

    # CPTAC Sample Metadata
    output$dl_cptac_samples <- downloadHandler(
      filename = "CPTAC_Sample_Info.rds",
      content = function(file) {
        file.copy("./data/cptac_samples_info_across_datatypes_v1.1.rds", file)
      }
    )

    output$dl_cptac_plus_samples <- downloadHandler(
      filename = "CPTAC_Extended_Sample_Info.rds",
      content = function(file) {
        file.copy("./data/cptac_plus_samples_info_across_datatypes_v1.1.rds", file)
      }
    )

    output$dl_ubidash_sample_mapping <- downloadHandler(
      filename = "CPTAC_Sample_ID_Mapping.rds",
      content = function(file) {
        file.copy("./data/UbiDash_CPTAC_PanCan_sample_ID_mapping_across_datatypes_v1.1.rds", file)
      }
    )

    # CPTAC Mutations & QTL
    output$dl_pancan_maf <- downloadHandler(
      filename = "CPTAC_PanCan_Mutations.rds",
      content = function(file) {
        file.copy("./data/PanCan_Union_Maf_Broad_WashU_v1.1_TGR.rds", file)
      }
    )

    output$dl_cptac_pqtls <- downloadHandler(
      filename = "CPTAC_pQTLs.rds",
      content = function(file) {
        file.copy("./data/CPTAC_Pancan_pQTLs_with_class.rds", file)
      }
    )

    # CPTAC Differential Expression
    output$dl_deqms_percancer <- downloadHandler(
      filename = "Table S2A. DEqMS_Percancer_CPTAC_imputed_proteome_matchedTN.xlsx",
      content = function(file) {
        file.copy("./data/Table S2A. DEqMS_Percancer_CPTAC_imputed_proteome_matchedTN.xlsx", file)
      }
    )

    output$dl_deqms_cancer_subtypes <- downloadHandler(
      filename = "DEqMS_Cancer_Subtypes.csv",
      content = function(file) {
        file.copy("./data/DEqMS.results_cancer_subtypes_with_tumotPurity.csv", file)
      }
    )

    output$dl_deqms_immune <- downloadHandler(
      filename = "DEqMS_Immune_Subtypes.csv",
      content = function(file) {
        file.copy("./data/DEqMS.results_Immune_subtype_percan_with_tumotPurity.csv", file)
      }
    )

    output$dl_deqms_matrix <- downloadHandler(
      filename = "DEqMS_Signed_Padj_Matrix.csv",
      content = function(file) {
        file.copy("./data/EXP193_DEqMS_CPTAC_imputed_proteome_signed_p.adj_matrix.csv", file)
      }
    )

    # CCLE Protein Expression
    output$dl_ccle_lfq_protein <- downloadHandler(
      filename = "CCLE_Protein_Expression.rds",
      content = function(file) {
        file.copy("./data/ccle_lfq_protein_expression.rds", file)
      }
    )

    output$dl_ccle_lfq_sample <- downloadHandler(
      filename = "CCLE_Sample_Info.rds",
      content = function(file) {
        file.copy("./data/ccle_lfq_sample_info.rds", file)
      }
    )

    # Drug Sensitivity
    output$dl_drug_longdf <- downloadHandler(
      filename = "Drug_Sensitivity_Data.csv",
      content = function(file) {
        file.copy("./data/drugsensitivity/drug_longdf.csv", file)
      }
    )

    output$dl_drugpath <- downloadHandler(
      filename = "Drug_Pathways.csv",
      content = function(file) {
        file.copy("./data/drugsensitivity/drugpath.csv", file)
      }
    )

    output$dl_cell_line_exp <- downloadHandler(
      filename = "Cell_Line_Basal_Expression.rds",
      content = function(file) {
        file.copy("./data/drugsensitivity/Cell_line_RMA_proc_basalExp.rds", file)
      }
    )

    output$dl_matched_cell <- downloadHandler(
      filename = "Matched_Cell_Annotations.csv",
      content = function(file) {
        file.copy("./data/drugsensitivity/MatchedCellAnno.csv", file)
      }
    )

    output$dl_repurposing <- downloadHandler(
      filename = "Drug_Repurposing_Screen_24Q2.rds",
      content = function(file) {
        file.copy("./data/Repurposing_Public_24Q2_LFC_COLLAPSED_tgr.rds", file)
      }
    )

    output$dl_mutations_matrix <- downloadHandler(
      filename = "Somatic_Mutations_Hotspots.rds",
      content = function(file) {
        file.copy("./data/OmicsSomaticMutationsMatrixHotspot.rds", file)
      }
    )

    # Harmonized & Correlation Data
    output$dl_harmonized_proteome <- downloadHandler(
      filename = "Harmonized_Proteome_Data.rds",
      content = function(file) {
        file.copy("./data/UbiDash_harmonized_proteome_data.rds", file)
      }
    )

    output$dl_ubidash_sample_mapping2 <- downloadHandler(
      filename = "Harmonized_Sample_Mapping.rds",
      content = function(file) {
        file.copy("./data/UbiDash_sample_mapping.rds", file)
      }
    )

    output$dl_ups_df <- downloadHandler(
      filename = "UPS_Protein_Annotations.rds",
      content = function(file) {
        file.copy("./data/ups_df.rds", file)
      }
    )

    output$dl_pancan_covariance <- downloadHandler(
      filename = "PanCan_Protein_Correlations.rds",
      content = function(file) {
        file.copy("./data/Pancan_protein~protein_spearman_correlations_ups_variance.rds", file)
      }
    )

    # LSCC Ubiquitylation
    output$dl_lscc_prot <- downloadHandler(
      filename = "LSCC_Protein_Expression.csv",
      content = function(file) {
        file.copy("./data/LSCC_prot_former_Combine_PanCan_CPTAC_Proteome-gene_V1.1_UMich_GENCODE34_Sinai_imputed.csv", file)
      }
    )

    output$dl_lscc_ub_prot <- downloadHandler(
      filename = "LSCC_Ubiquitylation_Sites.csv",
      content = function(file) {
        file.copy("./data/LSCC_ub_prot_former_1-s2.0-S0092867421008576-mmc3.csv", file)
      }
    )

    # Reference & Gene Lists
    output$dl_biogrid <- downloadHandler(
      filename = "BioGRID_Interactions.rds",
      content = function(file) {
        file.copy("./data/BIOGRID-ALL-4.4.241.tab3.rds", file)
      }
    )

    output$dl_combined_genelists <- downloadHandler(
      filename = "Combined_Gene_Lists.rds",
      content = function(file) {
        file.copy("./data/combined_genelists.rds", file)
      }
    )

    output$dl_e3_ligase <- downloadHandler(
      filename = "E3_Ligase_List.rds",
      content = function(file) {
        file.copy("./data/UbiDash_NYU+WEHI_E3_ligase_list.rds", file)
      }
    )

    output$dl_colors <- downloadHandler(
      filename = "Color_Scheme.csv",
      content = function(file) {
        file.copy("./data/colors.csv", file)
      }
    )

    output$dl_ups_network <- downloadHandler(
      filename = "Human_Proteostasis_Network.rds",
      content = function(file) {
        file.copy("./data/Human_Proteostasis_Network.rds", file)
      }
    )
  })
}
