library(shiny)
library(shinycssloaders)
library(plotly)
library(readxl)
library(bslib)
library(DT)
library(corrplot)
library(openxlsx)

source('./modules/oncoplot.R')
source('./modules/heatmap.R')
source('./modules/volcano.R')
source('./modules/correlations.R')
source('./modules/violin.R')
source('./modules/lollipop.R')
source('./modules/corranalysis.R')


# make plot square

filepath = "./data/DE_protein.xlsx"
cohorts = c("BRCA","CCRCC","COAD","GBM","HCC","HNSCC","LSCC","LUAD","OV","PDA","UCEC")
cohortSheets <- excel_sheets(path = filepath)
abstract =  paste(read.delim("./www/assets/abstract.txt"))

shinyUI(
  fluidPage(
    includeCSS("./www/style.css"),
    
    navbarPage(
      
      img(src="assets/ubiDash_logo.png", width="120"), theme = bs_theme(version = 5, bootswatch = "pulse"),
      # home page tab ----
      tabPanel("Home", theme = bs_theme(version = 5, bootswatch = "flatly"),
               tags$div(
                 h1("Integrated multi-omic analysis of Cullin-RING E3 ligases in cancer"),
                 p(abstract),
                 img(src='assets/OverviewSchematic.png', align = "right")%>% 
                   tagAppendAttributes(class = 'home-img')
               )  %>%  
                 tagAppendAttributes(class = 'title-page')
               
      ),
      navbarMenu("Mutations",
                 
                 # oncoplot tab ----
                 tabPanel("Oncoplot", theme = bs_theme(version = 5, bootswatch = "flatly"),
                          oncoplotUI("oncoTab", cohorts)
                 ),
                 # lollipop tab ----
                 tabPanel("Lollipop", theme = bs_theme(version = 5, bootswatch = "flatly"),
                          lollipopUI("lollipopTab", cohorts)
                 )
      ),
      navbarMenu("Expression Profile",
                 # heatmap tab ----
                 tabPanel("Heatmap", theme = bs_theme(version = 5, bootswatch = "flatly"),
                          heatmapUI("heatmapTab", cohorts)
                 ),
      ),
      # volcano plot tab ----
      tabPanel("Volcano Plots", theme = bs_theme(version = 5, bootswatch = "pulse"),
               volcanoUI("volcanoTab", cohortSheets)
      ),
      
      # correlation tab ----
      tabPanel("Correlations", theme = bs_theme(version = 5, bootswatch = "pulse"),
               correlationsUI("corrTab")
               
      ),
      # violin plot tap ----
      # need to let users input genes of interest
      tabPanel("Violin Plot", theme = bs_theme(version = 5, bootswatch = "pulse"),
               
               violinUI("violinTab")
      ),
      corrAnalysisUI("corrAnalysisTab", cohorts),
    ),
    tags$footer(
      img(src="assets/NYU-G-SoM_logo-final_White_RGB_300ppi.png", width="200"),
      tags$div(
        p("Footer placeholder text"),
        p("Footer placeholder text")
      ),
      img(src="assets/HHMI-Gilliam-Logo-cropped.png", width="256"),
    ),
    
    
    tags$script(src="navbar.js"),
    
  )
)
