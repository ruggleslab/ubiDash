# make plot square

cohorts <- c("BRCA", "CCRCC", "COAD", "GBM", "HCC", "HNSCC", "LSCC", "LUAD", "OV", "PDA", "UCEC")
abstract <- paste(read.delim("./www/assets/abstract.txt"))

shinyUI(
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


      img(src = "assets/ubiDash_logo.png", width = "120"),
      theme = bs_theme(version = 5, bootswatch = "pulse"),
      # home page tab ----
      tabPanel("Home",
        theme = bs_theme(version = 5, bootswatch = "pulse"),
        tags$div(
          h1("Integrated multi-omic analysis of Cullin-RING E3 ligases in cancer"),
          p(abstract),
          img(src = "assets/OverviewSchematic.png", align = "right") %>%
            tagAppendAttributes(class = "home-img")
        ) %>%
          tagAppendAttributes(class = "title-page")
      ),
      navbarMenu(
        "Mutations",

        # oncoplot tab ----
        tabPanel("Oncoplot",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          oncoplotUI("oncoTab", cohorts)
        ),
        # lollipop tab ----
        tabPanel("Lollipop",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          lollipopUI("lollipopTab")
        ),
        tabPanel("QTL Analysis",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          QTLUI("QTLTab", cohorts1.1)
        )
      ),
      navbarMenu(
        "Expression Profile",
        # heatmap tab ----
        tabPanel("Heatmap",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          heatmapUI("heatmapTab", cohorts1.1)
        ),
        tabPanel("Summarized Profile",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          simplifiedHeatmapUI("simplifiedHeatmapTab")
        ),
        # volcano plot tab ----
        tabPanel("Volcano Plots",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          volcanoUI("volcanoTab", DE_protein_cohorts)
        ),
        tabPanel("Density Plots",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          densityUI("densityTab")
        ),
      ),
      navbarMenu(
        "Functional Analysis",

        # correlation tab ----
        # tabPanel("Correlations",
        #   theme = bs_theme(version = 5, bootswatch = "pulse"),
        #   correlationsUI("corrTab")
        # ),
        # violin plot tap ----
        # need to let users input genes of interest
        tabPanel("Descriptive Statistics",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          violinUI("violinTab", cohorts1.1)
        ),
        tabPanel("Drug Sensitivity",
          theme = bs_theme(version = 5, bootswatch = "pulse"),
          drugResponseUI("drugResponseTab")
        ),
      ),
      corrAnalysisUI("corrAnalysisTab", cohorts1.1),
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
