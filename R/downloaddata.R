dataDownloadUI <- function(id){
  ns <- NS(id)

  card(
    card_header("Download Data"),
    tags$table(
      
      tags$tr(
        tags$th(tags$strong("PRIDE Dataset")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"PRIDE TPM RNA Expression", href="1x5ZGeGv0urNVx9lRqSkHN9pMt8lF4nxD")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"PRIDE FKPM RNA Expression", href="1wm40Ww8GzWT2Z0DzRU4NkA-BAtxVEti2")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"PRIDE log2(IBAQ) Protein Expression", href="1Vi82Vuqyd6eB7MzvIWIFnaadiiROoC1_")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"PRIDE LFQ Protein Expression", href="1rOl6RJ1z4tBZaGINNqfwh6-NsOPL3ePQ")),
      ),
      tags$tr(
        tags$th(tags$strong("CPTAC Dataset (v1.1)")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC TMT log2(AI) Protein Expression", href="1s54KghEG7pvVY_Gws0Dx1s9Kp3SmVqET")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC log2(TPM + 1) RNA Expression", href="1aDmO7hwc3gYmheIMHqvnd_IovFnhSUez")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC TPM RNA Expression", href="1PiOPhZsjJo73CdeSEWfH28i1yX96vme_")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC Sample Info", href="1kCtS4Ix-L5ZhN7it5T4tGZeQPTeReVtx")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC Pancan UPS eQTLs", href="1zkWxDRWMyp9UVYIVpZ3qMCHZTQRq5mXb")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC Pancan UPS pQTLs", href="10n0T-6IR5_UJuoE1MKKi1x0BVeLmKOx3")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC Mutations Matrix", href="1kXjqsFY-ABHtSpFC0lbGGVNIXWdcmylM")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CPTAC Pancan Mutations Maf", href="1RBPKXfEKPW5bR8pY5jT_EWTzpnJ1IkqZ")),
      ),
      tags$tr(
        tags$th(tags$strong("CCLE Dataset")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CCLE TPM RNA Expression", href="1PydHmcj3HSCAWDI9TjrVGApfoOMstd0J")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CCLE TMT Protein Expression", href="1U32-VxNUx9T59lSg0zyU2lbjk33ZDFYN")),
      ),
      tags$tr(
        tags$th(a(icon("download"),"CCLE Sample Info", href="1tu-OWS07EYACjHnbvhNlq6CDiskDgxlq")),
      ),
    )
  )
}

dataDownloadServer <- function(id) {
  moduleServer(id,function(input,output,session){

  })
}
