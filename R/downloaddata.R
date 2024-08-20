dataDownloadUI <- function(id){
  ns <- NS(id)

  card(
    card_header("Download Data"),
    tags$table(
      tags$tr(
        tags$th(a(icon("download"),"PanCan_Union_Maf_Broad_WashU_v1.1.maf",href="https://drive.google.com/uc?export=download&id=1RBPKXfEKPW5bR8pY5jT_EWTzpnJ1IkqZ")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC_Pancan_UPS_hCD_hTF_hChR_eQTLs.csv",href="https://drive.google.com/uc?export=download&id=1zkWxDRWMyp9UVYIVpZ3qMCHZTQRq5mXb")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC_Pancan_UPS_hCD_hTF_hChR_pQTLs.csv",href="https://drive.google.com/uc?export=download&id=10n0T-6IR5_UJuoE1MKKi1x0BVeLmKOx3")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"pride_lfq_protein_expression.csv",href="https://drive.google.com/uc?export=download&id=1rOl6RJ1z4tBZaGINNqfwh6-NsOPL3ePQ")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"duplicates_averaged_cptac_rna.csv",href="https://drive.google.com/uc?export=download&id=1VfXWEt4-YdAMnk8fG8dZT68-mBuPla9j")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_v1.1_mut.csv",href="https://drive.google.com/uc?export=download&id=1kXjqsFY-ABHtSpFC0lbGGVNIXWdcmylM")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_tmt-log2-RI_protein_expression.csv",href="https://drive.google.com/uc?export=download&id=1IvdObz_HZa0pwExDqWq-HXcI6QvRHQwQ")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_tmt-log2-AI_protein_expression.csv",href="https://drive.google.com/uc?export=download&id=1N5ILXzjiS1lK9fgsAI9W9XHQKKPJtqSI")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_multi_omic_correlations.csv",href="https://drive.google.com/uc?export=download&id=1-cs3oKmLaWxkWDZfYr3bBtrDv-Xki_Ig")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_sample_info_with_subtypes.csv",href="https://drive.google.com/uc?export=download&id=1gRysoGtYpQvcEwS7y2VeoMqZbCanPqm1")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"cptac_samples_info_across_datatypes_v1.1.rds",href="https://drive.google.com/uc?export=download&id=16Wf6gpvVQUBb_yau-MfzEFb6MscvNKIx")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"ccle_tpm_rna_expression.csv",href="https://drive.google.com/uc?export=download&id=1PydHmcj3HSCAWDI9TjrVGApfoOMstd0J")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"ccle_tmt_protein_expression.csv",href="https://drive.google.com/uc?export=download&id=1U32-VxNUx9T59lSg0zyU2lbjk33ZDFYN")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"ccle_sample_info.csv",href="https://drive.google.com/uc?export=download&id=1tu-OWS07EYACjHnbvhNlq6CDiskDgxlq")),
      )
    )
  )
}

dataDownloadServer <- function(id) {
  moduleServer(id,function(input,output,session){

  })
}
