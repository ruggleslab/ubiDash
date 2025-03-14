# make plot square
cohorts1.1 = c("BRCA","COAD","HGSC","LUAD","ccRCC","UCEC","PDAC","LSCC","GBM","HNSCC"
               # ,"MB","HCC"
)
cohorts = c("BRCA","COAD","HGSC","LUAD","CCRCC","UCEC","PDA","LSCC","GBM","HNSCC", "OV")
abstract <- paste(read.delim("./www/assets/abstract.txt"))

ui <- shinyUI(
  fluidPage(
    titlePanel(
      windowTitle = "UbiDash",
      title = tags$head(tags$link(
        rel = "icon",
        href = "/favicon.png",
        type = "image/png"
      ))
    ),
    includeCSS("./www/style.css"),
    navbarPage(
      id = "mainNav",
      img(src = "assets/ubiDash_logo.png", width = "120"),
      theme = bs_theme(version = 5, bootswatch = "pulse"),
      tabPanel("Home",
          tags$div(
          h1("Integrated multi-omic analysis of Cullin-RING E3 ligases in cancer"),
          p(abstract),
          img(src = "assets/OverviewSchematic.png", align = "right") %>%
            tagAppendAttributes(class = "home-img")
        ) %>%
          tagAppendAttributes(class = "title-page"),
        tags$iframe(src="https://docs.google.com/forms/d/e/1FAIpQLScjQcuKfMg_6HrA6aPIWoE2oH0S6pyt6a-spvUZr-SMmR0Kpg/viewform?embedded=true", width="640", height="449", frameborder="0", marginheight="0", marginwidth="0")
      ),
      navbarMenu(
        "Expression Profile",
        tabPanel("Heatmap",
                 heatmapHarmonizedUI("heatmapHarmonizedTab")
        ),
        tabPanel("Boxplot",
                 boxplotUI("boxplotTab")
        ),
        tabPanel("Density Plots",
                 densityUI("densityTab")
        )
      ),
      navbarMenu(
        "Differential Expression Profile",
        # volcano plot tab ----
        tabPanel("Volcano Plots",
                 volcanoUI("volcanoTab", DE_protein_cohorts)
        ),
        tabPanel("Heatmap Cancer Subtypes",
                 heatmapCancerUI("HeatmapCancerTab")
        ),
        tabPanel("Heatmap Immune Subtypes",
                 heatmapImmuneUI("HeatmapImmuneTab")
        )
      ),
      navbarMenu(
        "Correlations",
        tabPanel("Pancan Gene vs Gene",
                 scatterPlotUI("scatterPlotTab")
        ),
        tabPanel("Pancan Gene vs all",
                 correlationGeneAllUI("correlationGeneAllTab")
        ),
        tabPanel("Pancan Gene Panel",
                 correlationGeneUI("correlationGeneTab", cohorts1.1)
        ),
        tabPanel("LSCC Ubiquityl",
                 LSCCUbiquitylUI("LSCCUbiquitylTab")
        ),
        tabPanel("LSCC Upset",
                 LSCCUpsetUI("LSCCUpsetTab")
        )
      ), 
      navbarMenu(
        "Mutations (QTL)",
        tabPanel("Oncoplot",
          oncoplotUI("oncoTab", cohorts)
        ),
        tabPanel("Lollipop",
          lollipopUI("lollipopTab")
        ),
        tabPanel("QTL Analysis",
          QTLUI("QTLTab", cohorts1.1)
        )
      ),
      navbarMenu(
        "Functional Analysis (Drug association)",
        tabPanel("Descriptive Statistics",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          violinUI("violinTab", cohorts1.1)
        ),
        tabPanel("Drug Sensitivity",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          drugResponseUI("drugResponseTab")
        ),
      ),
      tabPanel("Download Data",
               theme = bs_theme(version = 5, bootswatch = "pulse"),
      dataDownloadUI("downloadTab"),
      ),
    ),
    
    tags$div() %>%
      tagAppendAttributes(class = "spacer"),
    tags$footer(
      class="d-flex flex-column flex-lg-row gap-4 gap-lg-2",
      img(src = "assets/NYU-G-SoM_logo-final_White_RGB_300ppi.png", width = "200",style = "height: auto"),
      tags$div(
          class = "d-flex flex-column flex-sm-row justify-content-center text-align-center gap-2",
          style ="white-space: nowrap",
      tags$div(
        class = "d-flex flex-column align-items-center align-items-lg-end",
        p(
          strong("Contact:"), "Ruggleslab@gmail.com"
        ),
        p(
          strong("Citation:"),
          tags$a(
            href = "citation-link-here", "González-Robles et al. (coming soon)",
            style = "color: white"
          )
        ),
      ),
      tags$div(
        class = "flex-column justify-content-center align-items-center d-none d-sm-flex",
      p(
        "|"
      ),
      ),
      tags$div(
        class = "d-flex flex-column align-items-center align-items-lg-start",
        p(
          tags$a(
            href = "https://www.paganolab.org", "The Pagano Lab",
            style = "color: white"
          )
        ),
        p(
          tags$a(
            href = "https://www.ruggleslab.org", "The Ruggles Lab",
            style = "color: white"
          )
        )
      ),
      ),
      img(src = "assets/HHMI-Gilliam-Logo-cropped.png", width = "256",style = "height: auto"),
    ),
    tags$script(src = "navbar.js"),
  )
)

