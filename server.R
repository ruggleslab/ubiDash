shinyServer(function(input, output, session) {
  

  
  # Observe which navbar tab is selected
  observeEvent(input$mainNav, {
    
    # For the "Heatmap" tab:
    if (input$mainNav == "Heatmap") {
      heatmapHarmonizedServer("heatmapHarmonizedTab")
    }
    
    # For the "Boxplot" tab:
    else if (input$mainNav == "Boxplot") {
      boxplotServer("boxplotTab")
    }
    
    # For the "Volcano Plots" tab:
    else if (input$mainNav == "Volcano Plots") {
      volcanoServer("volcanoTab", DE_protein)
    } 
    else if (input$mainNav == "Heatmap Cancer Subtypes") {
      heatmapCancerServer("HeatmapCancerTab")
    }
    else if (input$mainNav == "Heatmap Immune Subtypes") {
      heatmapImmuneServer("HeatmapImmuneTab")
    }
    
    # Similarly, add conditions for other tabs/modules as needed:
    else if (input$mainNav == "Density Plots") {
      # load file if required, then initialize the module
      densityServer("densityTab")
    }
    else if (input$mainNav == "Pancan Gene vs Gene") {
      scatterPlotServer("scatterPlotTab")
    }
    else if (input$mainNav == "Pancan Gene vs all") {
      correlationGeneAllServer("correlationGeneAllTab")
    }
    else if (input$mainNav == "Pancan Gene Panel") {
      correlationGeneServer("correlationGeneTab", cohorts1.1, colorsTbl)
    }
    else if (input$mainNav == "LSCC Ubiquityl") {
      LSCCUbiquitylServer("LSCCUbiquitylTab")
    }
    else if (input$mainNav == "LSCC Upset") {
      LSCCUpsetServer("LSCCUpsetTab")
    }
    else if (input$mainNav == "Oncoplot") {
      oncoplotServer("oncoTab", top_crls, cohorts, colorsTbl)
    }
    else if (input$mainNav == "Lollipop") {
      lollipopServer("lollipopTab")
    }
    else if (input$mainNav == "QTL Analysis") {
      QTLServer("QTLTab", genelists, colorsTbl, cohorts1.1)
    }
    else if (input$mainNav == "Drug Sensitivity") {
      drugResponseServer("drugResponseTab", drug, drugpath, exp, cell, genelists$UPS)
    }
    else if (input$mainNav == "Download Data") {
      dataDownloadServer("downloadTab")
    }
  })
  
})