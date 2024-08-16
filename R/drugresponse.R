drugResponseUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      h3("Drug Response"),
      p("(description paragraph)"),
      selectizeInput(
        inputId = ns("inputGene"),
        choices = NULL,
        multiple = FALSE,
        label = "Gene",
        options = list(`live-search` = TRUE),
      ),
      selectizeInput(
        inputId = ns("inputDrugs"),
        choices = NULL,
        multiple = TRUE,
        label = "Drugs (Optional)",
        options = list(`live-search` = TRUE),
      ),
      actionButton(ns("refresh"), "Generate Plot"),
      
      p("Download this data"),
      downloadButton(outputId = ns("selected_gene_drugmat"), label = "selected_gene_drugmat"),
      downloadButton(outputId = ns("selected_gene_all_wilcoxon"), label = "selected_gene_all_wilcoxon"),
      downloadButton(outputId = ns("selected_gene_sigdrug_wilcoxon"), label = "selected_gene_sigdrug_wilcoxon"),
      downloadButton(outputId = ns("selected_gene_subdrugmat"), label = "selected_gene_subdrugmat"),
      downloadButton(outputId = ns("selected_gene_corr"), label = "selected_gene_corr")
    ),
    mainPanel(
      textOutput(ns("err")) %>%
        tagAppendAttributes(class = "error"),
      p("Data from TCGA"),
      shinycssloaders::withSpinner(plotOutput(ns("plot1")),
        color = "#6A2798"
      ),
      p("Data from CCLE"),
      shinycssloaders::withSpinner(plotOutput(ns("plot2")),
        color = "#6A2798"
      ),
      p("Data from CCLE"),
      shinycssloaders::withSpinner(plotlyOutput(ns("plot3")),
        color = "#6A2798"
      ),
      DTOutput(ns("plot3data")),
    )
  )
}

drugResponseServer <- function(id, drug, drugpath, exp, cell, crl_genes) {
  moduleServer(id, function(input, output, session) {
    # initial data loading ----
    # read data frames generated from Processing.R
    # selected_gene expression df
    # drug ic50 longerdf

    updateSelectizeInput(session, "inputGene",
      choices = as.list(crl_genes),
      server = TRUE, selected = "VHL"
    )

    updateSelectizeInput(session, "inputDrugs",
                         choices = unique(drug$DRUG_NAME),
                         server = TRUE,
    )

    selected_gene_drugmat <- reactiveVal()
    selected_gene_all_wilcoxon <- reactiveVal()
    selected_gene_sigdrug_wilcoxon <- reactiveVal()
    selected_gene_subdrugmat <- reactiveVal()
    selected_gene_corr <- reactiveVal()
    
    selected_gene = eventReactive(input$refresh, {
      # generate data for selected gene ----
      selected_gene_name <- isolate(input$inputGene)

      if (!isTruthy(selected_gene_name)) {
        selected_gene_name <- "VHL"
      }

      print(selected_gene_name)
      print(selected_gene)

      selected_gene_exp <- exp[exp$GENE_SYMBOLS == selected_gene_name, ]
      selected_gene_exp <- t(selected_gene_exp[, -c(1:2)]) # transpose
      selected_gene_exp <- as.data.frame(selected_gene_exp)


      selected_gene_exp$Data <- rownames(selected_gene_exp) # add column of sample IDs


      colnames(selected_gene_exp)[1] <- "Exp"
      rownames(selected_gene_exp) <- NULL

      # match sample annotations with expression table

      selected_gene = merge(selected_gene_exp, cell, by = "Data")
      # take average of duplicated cell lines
      dupcl <-
        selected_gene$CELL_LINE_NAME[which(duplicated(selected_gene$CELL_LINE_NAME))] # duplicated cell line names
      dup <-
        selected_gene[which(selected_gene$CELL_LINE_NAME %in% dupcl), ] # duplicated df
      selected_gene =
        selected_gene[-c(which(selected_gene$CELL_LINE_NAME %in% dupcl)), ] # remove duplicated
      dup <-
        dup[order(dup$CELL_LINE_NAME), ] # order by cell line names
      dup_ <- dup[c(T, F), ] # keep only one result for each cell line

      for (i in c(1:(dim(dup)[1] / 2))) {
        # calculate average of each pair and fill the column
        m <- (dup$Exp[2 * i - 1] + dup$Exp[2 * i]) / 2
        dup_$Exp[i] <- m
      }
      # add duplicated cell lines back to df
      selected_gene = rbind(selected_gene, dup_)
      selected_gene =
        selected_gene[order(selected_gene$Exp, decreasing = T), ] # reorder expression
      
      return(selected_gene)
    },ignoreNULL = FALSE, ignoreInit=FALSE)

      # 1. expression by type ----
      output$plot1 <- renderPlot({
        selected_gene()
        
        
        if(!isTruthy(selected_gene())){
          return(NULL)
        }
        
        print("plot1")

        ### generate expression by type

        typedf <-
          selected_gene()[, c(2, 4)] # keep expression and tcga annotation
        return(
          ggplot(typedf, aes(x = reorder(TCGA_DESC, -Exp, median), y = Exp)) +
            geom_boxplot() +
            geom_point() +
            theme(axis.text.x = element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            )) +
            xlab(paste0(isolate(input$gene), "expression ordered by median"))
          # ggsave("selected_gene/selected_gene_by_tcga.pdf", width = 20, height = 15, units = "cm")
        )
      })

      # 2. z-score heatmap for top/bot ----
      output$plot2 <- renderPlot({
        
        if(!isTruthy(selected_gene())){
          return(NULL)
        }
        
        print("plot2")
        
        selected_gene_copy = selected_gene()
        
        # top 10% and bottom 10% selected_gene
        selected_gene_copy$type <- "-"
        # subset selected_gene top/bot
        selected_gene_copy$type[c(1:100)] <- "Top"
        selected_gene_copy$type[c((nrow(selected_gene_copy) - 99):nrow(selected_gene_copy))] <-
          "Bot"
        top <- selected_gene_copy[selected_gene_copy$type == "Top", ]
        bot <- selected_gene_copy[selected_gene_copy$type == "Bot", ]
        # get top/bot cell line names
        topCL <- top$CELL_LINE_NAME
        botCL <- bot$CELL_LINE_NAME
        # subset top and bottom df
        drugt <- drug[drug$CELL_LINE_NAME %in% topCL, ]
        length(unique(drugt$CELL_LINE_NAME)) # 100 cell lines in top
        drugb <- drug[drug$CELL_LINE_NAME %in% botCL, ]
        length(unique(drugb$CELL_LINE_NAME)) # 100 cell lines in bot

        # IC50 data for heatmap
        drugdf <-
          rbind(drugt, drugb)[, c(1, 2, 5)] # keep only drugname, cell line name, ic50
        drugdf <- as.data.frame(drugdf)
        drugmat <- drugdf %>% pivot_wider(
          names_from = DRUG_NAME,
          values_from = LN_IC50,
          values_fn = list(LN_IC50 = mean)
        )
        drugmat <- as.data.frame(drugmat)
        rownames(drugmat) <- drugmat$CELL_LINE_NAME
        drugmat <- drugmat[, -1]
        drugmat <- as.matrix(drugmat)
        df <- as.data.frame(drugmat)
        # 200 cell lines + 542 drugs

        selected_gene_drugmat(df)
        # write.csv(df, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_drugmat.csv", row.names = T)

        # expression data for heatmap
        selected_genemat <- rbind(top[, c(2, 3, 5)], bot[, c(2, 3, 5)])
        selected_genemat$type <- as.factor(selected_genemat$type)

        # pathway annotation
        drugpath$PATHWAY_NAME <- as.factor(drugpath$PATHWAY_NAME)

        # generate df for wilcoxon test
        df_longer <- cbind(rownames(df), df)
        colnames(df_longer)[1] <- "CELL_LINE_NAME"
        df_longer <-
          df_longer %>% pivot_longer(!CELL_LINE_NAME,
            names_to = "drug",
            values_to = "ln(IC50)"
          )
        colnames(df_longer)[2] <- "DRUG_NAME"
        df_longer$Type <- "Top"
        df_longer$Type[which(df_longer$CELL_LINE_NAME %in% botCL)] <-
          "Bot"
        df_longer$CELL_LINE_NAME <-
          as.factor(df_longer$CELL_LINE_NAME)
        df_longer$Type <- as.factor((df_longer$Type))
        # combine drug and type info
        df_longer$label <-
          paste(df_longer$DRUG_NAME, "_", df_longer$Type)
        df_longer <- df_longer[, c(1, 3, 5)]
        # testdf for wilcoxon test,108400 obs. of 2 vars
        testdf <- df_longer[, c(2, 3)]
        testdf$label <- as.factor(testdf$label)
        testdf <- testdf[order(testdf$label), ]

        # get p values
        n <- c(1:542) # for 542 drugs in 200 cell lines
        p <- c()
        e <- c()
        for (i in n)
        {
          wil <-
            testdf[c((200 * i - 199):(200 * i)), ] # subset each cell lines # CHANGE 200 to cell line (top+bot) and 199 to (top+bot-1)
          t_i <-
            wilcox.test(`ln(IC50)` ~ label,
              data = wil,
              conf.int = TRUE
            ) # wilcoxon test within each cell line top/bot
          p <- append(p, t_i$p.value)
          e <- append(e, t_i$estimate)
        }

        # get corresponding drug lists
        druglist <- unique(testdf$label)
        druglist <- druglist[c(T, F)]
        druglist <- gsub(" _.*", "", druglist)

        test <- cbind(druglist, p, e) # test result df
        test <- as.data.frame(test)
        colnames(test)[2] <- "Wilcoxon.test.p.value"
        # p-value is the probability of getting a test statistic as large or
        # larger assuming both distributions are the same
        colnames(test)[3] <- "Wilcoxon.test.estimate"
        test$Wilcoxon.test.estimate <-
          as.numeric(test$Wilcoxon.test.estimate)
        test$Wilcoxon.test.p.value <-
          as.numeric(test$Wilcoxon.test.p.value)

        # significant drugs e>0.5(95) +p(86)
        topdrug <-
          test[which(abs(test$Wilcoxon.test.estimate) > 0.5), ]
        topdrug <-
          topdrug[which(topdrug$Wilcoxon.test.p.value < 0.05), ]
        # add annotations to topdrug
        drugindex <- match(topdrug$druglist, drugpath$DRUG_NAME)
        druganno <- drugpath$PATHWAY_NAME[drugindex]
        drug_sub <- drugpath[drugindex, ]
        topdrug$pathway <- druganno

        selected_gene_all_wilcoxon(test)
        selected_gene_sigdrug_wilcoxon(topdrug)
        # write.csv(test, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_all_wilcoxon.csv")
        # write.csv(topdrug, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_sigdrug_wilcoxon_0.5+p.csv")

        topdrug <-
          topdrug[order(topdrug$Wilcoxon.test.estimate, decreasing = F), ] # order by estimate
        topdrug_names <- topdrug$druglist # get top drugnames

        # subset ic50 data for top drugs
        drugmat_sub <-
          drugmat[, match(topdrug_names, colnames(drugmat))]

        # scale to z-score sig drug plot
        drugmat_sub <- scale(drugmat_sub, center = T, scale = T)
        drugmat_sub[is.nan(drugmat_sub)] <- NA
        # 200 cell lines + topdrugs (86)

        selected_gene_subdrugmat(drugmat_sub)
        # write.csv(drugmat_sub, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_subdrugmat.csv", row.names = T)

        ## generate heatmap
        # pdf(file="/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_sigdrugs_0.5+p_scaled.pdf",
        #     width=15,height=6) #CHANGE if heatmap size changed
        # Heatmap
        col_path <- c(
          "#bdd7e7",
          "#6baed6",
          "#3182bd",
          "#08519c",
          "#bae4b3",
          "#74c476",
          "#31a354",
          "#006d2c",
          "#fdbe85",
          "#fd8d3c",
          "#e6550d",
          "#a63603",
          "#cbc9e2",
          "#9e9ac8",
          "#756bb1",
          "#54278f",
          "#fbb4b9",
          "#f768a1",
          "#c51b8a",
          "#7a0177",
          "#cccccc",
          "#969696",
          "#636363",
          "#252525"
        ) # pathway colors
        names(col_path) <- levels(drugpath$PATHWAY_NAME)

        col_selected_gene <-
          colorRamp2(c(0, 15), c("white", "#3C2692")) # expression scale and color
        col_sub <-
          c("Bot" = "#EEC9E5", "Top" = "#7C4D79") # group label and color
        col_fun <-
          colorRamp2(c(-5, 0, 5), c("#3F2D91", "#E5E5E5", "red")) # ic50 scale and color
        # pathway annotation of drugs
        column_t <- HeatmapAnnotation(
          path = drug_sub$PATHWAY_NAME,
          annotation_name_gp = gpar(fontsize = 8),
          annotation_legend_param = list(
            title_position = "lefttop-rot",
            title = "drugpath",
            direction = "horizontal",
            ncol = 6
          ),
          col = list(path = col_path)
        )
        # Top/bot groups row annotation
        row_sub <- rowAnnotation(
          group = selected_genemat$type,
          annotation_legend_param = list(
            at = c("Top", "Bot"),
            labels = c("Top", "Bot"),
            title_position = "lefttop-rot",
            title = "subgroup"
          ),
          col = list(group = col_sub),
          annotation_name_gp = gpar(fontsize = 5.5),
          annotation_name_rot = 0,
          annotation_name_side = "bottom"
        )
        # selected_gene expression annotation
        row_exp <-
          rowAnnotation(
            selected_gene = selected_genemat$Exp,
            annotation_name_gp = gpar(fontsize = 5.5),
            annotation_legend_param = list(
              title_position = "lefttop-rot",
              title =
                paste0(isolate(input$gene), "exp")
            ),
            col = list(selected_gene = col_selected_gene),
            annotation_name_rot = 0,
            annotation_name_side = "bottom"
          )
        # main heatmap
        ht_IC50 <-
          Heatmap(
            drugmat_sub,
            name = "ln(IC50)",
            col = col_fun,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 6),
            heatmap_legend_param = list(
              title = "ln(IC50)",
              title_position = "lefttop-rot"
            ),
            column_names_rot = 90,
            heatmap_width = unit(25, "cm"),
            # CHANGE if needed
            heatmap_height = unit(10, "cm"),
            # CHANGE if needed and also change pdf size
            row_names_side = "right",
            show_row_names = F,
            # hastag this if cell line names are needed
            column_title = "Z-scored significant drugs between Top/Bot",
            column_title_side = "bottom",
            cluster_columns = F,
            cluster_rows = F,
            left_annotation = row_sub,
            right_annotation = row_exp,
            row_split = factor(selected_genemat$type, levels = c("Top", "Bot")),
            top_annotation = column_t
          )
        return(draw(
          ht_IC50,
          heatmap_legend_side = "left",
          annotation_legend_side = "top"
          # legend_grouping = "original"
        ))

        # dev.off()
      })
      # 3. volcano plots of spearman correlation among all drugs and cell lines ----
      
      
      
      output$plot3 <- renderPlotly({
        
        if(!isTruthy(selected_gene())){
          return(NULL)
        }
        
        print("plot3")
        # spearman correlation
        drugall <- drug[, c(1, 2, 5)] # keep cell line, ic50 and drugname
        # wider df
        drugall <- drugall %>% pivot_wider(
          names_from = DRUG_NAME,
          values_from = LN_IC50,
          values_fn = list(LN_IC50 = mean)
        )
        drugall <- as.data.frame(drugall)
        selected_gene_all <- selected_gene()[, c(2, 3)]
        colnames(selected_gene_all)[2] <- "CELL_LINE_NAME"
        corrdf <-
          merge(selected_gene_all, drugall, by = "CELL_LINE_NAME") # merge expression data and ic50data

        # get spearman correlation stats
        n <- c(1:542)
        p <- c() # empty vector for p values
        r <- c() # empty vector for estimates
        for (i in n)
        {
          t_i <-
            cor.test(corrdf[, 2], corrdf[, i + 2], method = "spearman") # spearman correlation between exp and ic50
          p <- append(p, t_i$p.value)
          r <- append(r, t_i$estimate)
        }

        # generate dataframe for correlation results
        druglist <- colnames(corrdf)[c(3:544)] # drug names
        testdf <-
          cbind(druglist, r, p) # combine results of correlation analysis
        testdf <- as.data.frame(testdf)
        colnames(testdf)[2] <- "corr.coef"
        colnames(testdf)[3] <- "p.value"
        testdf$p.value <- as.numeric(testdf$p.value)
        testdf$corr.coef <- as.numeric(testdf$corr.coef)
        testdf$p.adjust <- p.adjust(testdf$p.value, method = "BH")
        testdf$log10.padj <- (-1 * log10(testdf$p.adjust))
        testdf$label <- "Not significant"
        testdf$label[which(testdf$p.adjust < 0.05)] <- "<0.05"
        testdf$group <- "-"
        testdf$group[which(testdf$log10.padj > 5)] <-
          "Top" # CHANGE manual threshold for label on plot, change if needed
        testdf$sig <- "No"
        testdf$sig[which(testdf$p.adjust < 0.05)] <- "Yes"

        # add drug annotaiton to result table
        colnames(testdf)[1] <- "DRUG_NAME"
        testdf_ <- merge(drugpath, testdf)
        testdf_ <-
          testdf_[order(testdf_$corr.coef, decreasing = F), ] # order result table by correlation estimates

        selected_gene_corr(testdf_)
        # write.csv(testdf_, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_corr.csv",row.names = F)

        # generate volcano plot_met
        # pdf(file="/Users/staceyvann/Desktop/Lab/Drug Sensitivity/selected_gene/selected_gene_Corr_vln.pdf",width=11,height=9)
        # vln <-
        print(head(testdf))
        text <- paste(
          "<b>", "Drug Name: ", testdf$DRUG_NAME, "</b>", "\n",
          "corr coeff: ", testdf$corr.coef, "\n",
          "log10 padj: ", testdf$log10.padj, "\n"
        )
        
        output$plot3data <- renderDT({
          testdf
        })
        
        plt <-
          ggplot(testdf, aes(x = corr.coef, y = log10.padj)) +
          geom_point(aes(color = label, text = text), size = 3) +
          xlim(-0.4, 0.4) + # x and y range
          ylim(0, 10) +
          xlab("Spearman Correlation Coefficient") +
          ylab("-log<sub>10</sub> (Q-value)") +
          ggtitle("Spearman correlation test between selected_gene expression and drug ln(IC50)") +
          theme_bw() +
          theme(panel.background = element_blank()) +
          theme(panel.border = element_rect(
            colour = "black",
            fill = NA,
            size = 1
          )) +
          theme(text = element_text(size = 10)) +
          theme(axis.text = element_text(color = "black", size = 10)) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line()
          ) +
          theme(legend.position = "right") +
          labs(color = "padj") +
          geom_vline(xintercept = 0, linetype = "dotted") +
          # geom_text_repel(
          #   #top drugs label
          #   data = subset(testdf, group == "Top"),
          #   aes(label = testdf$DRUG_NAME[which(testdf$group == "Top")]),
          #   size = 4,
          #   box.padding = unit(0.35, "lines"),
          #   point.padding = unit(0.5, "lines"),
          #   max.overlaps = Inf
          # ) +
          # change text according to x range
          annotate(
            "text",
            x = -0.25,
            y = 0,
            label = c("<b>Sensitive</b>"),
            parse = T,
            size = 5
          ) +
          annotate(
            "text",
            x = 0.25,
            y = 0,
            label = c("<b>Resistant</b>"),
            parse = T,
            size = 5
          )
        
        labeledPoints <- filter(testdf, (DRUG_NAME %in% isolate(input$inputDrugs)))
        print(labeledPoints)
        print("DRUGS")
        
        if (nrow(labeledPoints) > 0) {
          plt <- plt + geom_text(
            hjust = 1, vjust = 1, data = labeledPoints,
            mapping = aes(label = DRUG_NAME, fontface = "bold"),
            color = "#000000"
          )
        }

        # print(vln)
        # dev.off()
        return(ggplotly(plt, tooltip = "text"))
      })
     
    
    
    # downloads ----
    output$selected_gene_drugmat <- downloadHandler(
      filename = function() {
        paste("selected_gene_drugmat-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(selected_gene_drugmat(), file)
      }
    )

    output$selected_gene_all_wilcoxon <- downloadHandler(
      filename = function() {
        paste("selected_gene_all_wilcoxon-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(selected_gene_all_wilcoxon(), file)
      }
    )

    output$selected_gene_sigdrug_wilcoxon <- downloadHandler(
      filename = function() {
        paste("selected_gene_sigdrug_wilcoxon-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(selected_gene_sigdrug_wilcoxon(), file)
      }
    )

    output$selected_gene_subdrugmat <- downloadHandler(
      filename = function() {
        paste("selected_gene_subdrugmat-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(selected_gene_subdrugmat(), file)
      }
    )

    output$selected_gene_corr <- downloadHandler(
      filename = function() {
        paste("selected_gene_corr-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(selected_gene_corr(), file)
      }
    )
  })
}
