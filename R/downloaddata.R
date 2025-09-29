dataDownloadUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Download Data"),
    p(tab_summaries$`Download Data`),
    tags$table(
      tags$tr(
        tags$th(tags$strong("PRIDE"), a("PMID: 30777892", href = "https://pubmed.ncbi.nlm.nih.gov/30777892/")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "PRIDE TPM RNA Expression", href = "https://drive.google.com/uc?export=download&id=1x5ZGeGv0urNVx9lRqSkHN9pMt8lF4nxD")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "PRIDE FKPM RNA Expression", href = "https://drive.google.com/uc?export=download&id=1wm40Ww8GzWT2Z0DzRU4NkA-BAtxVEti2")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "PRIDE log2(IBAQ) Protein Expression", href = "https://drive.google.com/uc?export=download&id=1Vi82Vuqyd6eB7MzvIWIFnaadiiROoC1_")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "PRIDE LFQ Protein Expression", href = "https://drive.google.com/uc?export=download&id=1rOl6RJ1z4tBZaGINNqfwh6-NsOPL3ePQ")),
      ),
      tags$tr(
        tags$th(tags$strong("CPTAC"), a("PMID: 37582339", href = "https://pubmed.ncbi.nlm.nih.gov/37582339/")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC TMT log2(AI) Protein Expression", href = "https://drive.google.com/uc?export=download&id=1s54KghEG7pvVY_Gws0Dx1s9Kp3SmVqET")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC log2(TPM + 1) RNA Expression", href = "https://drive.google.com/uc?export=download&id=1aDmO7hwc3gYmheIMHqvnd_IovFnhSUez")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC TPM RNA Expression", href = "https://drive.google.com/uc?export=download&id=1PiOPhZsjJo73CdeSEWfH28i1yX96vme_")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC Sample Info", href = "https://drive.google.com/uc?export=download&id=1kCtS4Ix-L5ZhN7it5T4tGZeQPTeReVtx")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC Pancan UPS eQTLs", href = "https://drive.google.com/uc?export=download&id=1zkWxDRWMyp9UVYIVpZ3qMCHZTQRq5mXb")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC Pancan UPS pQTLs", href = "https://drive.google.com/uc?export=download&id=10n0T-6IR5_UJuoE1MKKi1x0BVeLmKOx3")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC Mutations Matrix", href = "https://drive.google.com/uc?export=download&id=1kXjqsFY-ABHtSpFC0lbGGVNIXWdcmylM")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CPTAC Pancan Mutations Maf", href = "https://drive.google.com/uc?export=download&id=1RBPKXfEKPW5bR8pY5jT_EWTzpnJ1IkqZ")),
      ),
      tags$tr(
        tags$th(tags$strong("CCLE TPM"), a("DepMap", href = "https://depmap.org/portal/data_page/?tab=allData")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CCLE TPM RNA Expression", href = "https://drive.google.com/uc?export=download&id=1PydHmcj3HSCAWDI9TjrVGApfoOMstd0J")),
      ),
      tags$tr(
        tags$th(tags$strong("CCLE TMT"), a("PMID: 31978347", href = "https://pubmed.ncbi.nlm.nih.gov/31978347/")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CCLE TMT Protein Expression", href = "https://drive.google.com/uc?export=download&id=1U32-VxNUx9T59lSg0zyU2lbjk33ZDFYN")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CCLE TMT Sample Info", href = "https://drive.google.com/uc?export=download&id=1tu-OWS07EYACjHnbvhNlq6CDiskDgxlq")),
      ),
      tags$tr(
        tags$th(tags$strong("CCLE LFQ"), a("PMID: 35839778", href = "https://pubmed.ncbi.nlm.nih.gov/35839778/")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CCLE LFQ Sample Info", href = "https://drive.google.com/uc?export=download&id=1LFv8WxNCnE98fWa1ZkZiaGBIdDJSMDBK")),
      ),
      tags$tr(
        tags$th(a("•", icon("download"), "CCLE LFQ Protein Expression", href = "https://drive.google.com/uc?export=download&id=1Bsi8IAGgh733LDvYH2arn5dwcUVQrPyv")),
      ),
    )
  )
}

dataDownloadServer <- function(id) {
  moduleServer(id, function(input, output, session) {

  })
}
