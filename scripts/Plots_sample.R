# sample code for generating plots
# contents:
# 1. expression by type
# 2. z-score heatmap for top/bot
# 3. volcano plots of spearman correlation among all drugs and cell lines

# 47 occurrence for each target gene

setwd("/Users/staceyvann/Desktop/Lab/Drug Sensitivity")
library(tidyr)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)

# read data frames generated from Processing.R
# KBTBD11 expression df
KBTBD11 <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/KBTBD11.csv", row.names = 1)
# drug ic50 longerdf
drug <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drug_longdf.csv", row.names = 1)
# drug pathway and targets annotation
drugpath <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drugpath.csv", row.names = 1)

# 1. expression by type
### generate expression by type
typedf <- KBTBD11[, c(2, 4)] # keep expression and tcga annotation
ggplot(typedf, aes(x = reorder(TCGA_DESC, -Exp, median), y = Exp)) +
  geom_boxplot() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("KBTBD11 expression ordered by median")
ggsave("KBTBD11/KBTBD11_by_tcga.pdf", width = 20, height = 15, units = "cm")

# 2. z-score heatmap for top/bot
# top 10% and bottom 10% KBTBD11
KBTBD11$type <- "-"
# subset KBTBD11 top/bot
KBTBD11$type[c(1:100)] <- "Top"
KBTBD11$type[c(829:928)] <- "Bot"
top <- KBTBD11[KBTBD11$type == "Top", ]
bot <- KBTBD11[KBTBD11$type == "Bot", ]
# get top/bot cell line names
topCL <- top$CELL_LINE_NAME
botCL <- bot$CELL_LINE_NAME
# subset top and bottom df
drugt <- drug[drug$CELL_LINE_NAME %in% topCL, ]
length(unique(drugt$CELL_LINE_NAME)) # 100 cell lines in top
drugb <- drug[drug$CELL_LINE_NAME %in% botCL, ]
length(unique(drugb$CELL_LINE_NAME)) # 100 cell lines in bot

# IC50 data for heatmap
drugdf <- rbind(drugt, drugb)[, c(1, 2, 5)] # keep only drugname, cell line name, ic50
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
write.csv(df, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_drugmat.csv", row.names = T)

# expression data for heatmap
KBTBD11mat <- rbind(top[, c(2, 3, 5)], bot[, c(2, 3, 5)])
KBTBD11mat$type <- as.factor(KBTBD11mat$type)

# pathway annotation
drugpath$PATHWAY_NAME <- as.factor(drugpath$PATHWAY_NAME)

# generate df for wilcoxon test
df_longer <- cbind(rownames(df), df)
colnames(df_longer)[1] <- "CELL_LINE_NAME"
df_longer <- df_longer %>% pivot_longer(!CELL_LINE_NAME, names_to = "drug", values_to = "ln(IC50)")
colnames(df_longer)[2] <- "DRUG_NAME"
df_longer$Type <- "Top"
df_longer$Type[which(df_longer$CELL_LINE_NAME %in% botCL)] <- "Bot"
df_longer$CELL_LINE_NAME <- as.factor(df_longer$CELL_LINE_NAME)
df_longer$Type <- as.factor((df_longer$Type))
# combine drug and type info
df_longer$label <- paste(df_longer$DRUG_NAME, "_", df_longer$Type)
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
  wil <- testdf[c((200 * i - 199):(200 * i)), ] # subset each cell lines # CHANGE 200 to cell line (top+bot) and 199 to (top+bot-1)
  t_i <- wilcox.test(`ln(IC50)` ~ label, data = wil, conf.int = TRUE) # wilcoxon test within each cell line top/bot
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
test$Wilcoxon.test.estimate <- as.numeric(test$Wilcoxon.test.estimate)
test$Wilcoxon.test.p.value <- as.numeric(test$Wilcoxon.test.p.value)

# significant drugs e>0.5(95) +p(86)
topdrug <- test[which(abs(test$Wilcoxon.test.estimate) > 0.5), ]
topdrug <- topdrug[which(topdrug$Wilcoxon.test.p.value < 0.05), ]
# add annotations to topdrug
drugindex <- match(topdrug$druglist, drugpath$DRUG_NAME)
druganno <- drugpath$PATHWAY_NAME[drugindex]
drug_sub <- drugpath[drugindex, ]
topdrug$pathway <- druganno

write.csv(test, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_all_wilcoxon.csv")
write.csv(topdrug, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_sigdrug_wilcoxon_0.5+p.csv")

topdrug <- topdrug[order(topdrug$Wilcoxon.test.estimate, decreasing = F), ] # order by estimate
topdrug_names <- topdrug$druglist # get top drugnames

# subset ic50 data for top drugs
drugmat_sub <- drugmat[, match(topdrug_names, colnames(drugmat))]

# scale to z-score sig drug plot
drugmat_sub <- scale(drugmat_sub, center = T, scale = T)
drugmat_sub[is.nan(drugmat_sub)] <- NA
# 200 cell lines + topdrugs (86)
write.csv(drugmat_sub, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_subdrugmat.csv", row.names = T)

## generate heatmap
pdf(
  file = "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_sigdrugs_0.5+p_scaled.pdf",
  width = 15, height = 6
) # CHANGE if heatmap size changed
# Heatmap
col_path <- c(
  "#bdd7e7", "#6baed6", "#3182bd", "#08519c",
  "#bae4b3", "#74c476", "#31a354", "#006d2c",
  "#fdbe85", "#fd8d3c", "#e6550d", "#a63603",
  "#cbc9e2", "#9e9ac8", "#756bb1", "#54278f",
  "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177",
  "#cccccc", "#969696", "#636363", "#252525"
) # pathway colors
names(col_path) <- levels(drugpath$PATHWAY_NAME)

col_KBTBD11 <- colorRamp2(c(0, 15), c("white", "#3C2692")) # expression scale and color
col_sub <- c("Bot" = "#EEC9E5", "Top" = "#7C4D79") # group label and color
col_fun <- colorRamp2(c(-5, 0, 5), c("#3F2D91", "#E5E5E5", "red")) # ic50 scale and color
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
  group = KBTBD11mat$type,
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
# KBTBD11 expression annotation
row_exp <- rowAnnotation(
  KBTBD11 = KBTBD11mat$Exp,
  annotation_name_gp = gpar(fontsize = 5.5),
  annotation_legend_param = list(
    title_position = "lefttop-rot",
    title = "KBTBD11 exp"
  ),
  col = list(KBTBD11 = col_KBTBD11),
  annotation_name_rot = 0,
  annotation_name_side = "bottom"
)
# main heatmap
ht_IC50 <- Heatmap(drugmat_sub,
  name = "ln(IC50)", col = col_fun,
  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "ln(IC50)",
    title_position = "lefttop-rot"
  ),
  column_names_rot = 90, heatmap_width = unit(25, "cm"), # CHANGE if needed
  heatmap_height = unit(10, "cm"), # CHANGE if needed and also change pdf size
  row_names_side = "right",
  show_row_names = F, # hastag this if cell line names are needed
  column_title = "Z-scored significant drugs between Top/Bot",
  column_title_side = "bottom",
  cluster_columns = F, cluster_rows = F,
  left_annotation = row_sub,
  right_annotation = row_exp,
  row_split = factor(KBTBD11mat$type, levels = c("Top", "Bot")),
  top_annotation = column_t
)
draw(ht_IC50,
  heatmap_legend_side = "left",
  annotation_legend_side = "top"
  # legend_grouping = "original"
)

dev.off()




# 3. volcano plots of spearman correlation among all drugs and cell lines
# spearman correlation
drugall <- drug[, c(1, 2, 5)] # keep cell line, ic50 and drugname
# wider df
drugall <- drugall %>% pivot_wider(
  names_from = DRUG_NAME,
  values_from = LN_IC50,
  values_fn = list(LN_IC50 = mean)
)
drugall <- as.data.frame(drugall)
KBTBD11_all <- KBTBD11[, c(2, 3)]
colnames(KBTBD11_all)[2] <- "CELL_LINE_NAME"
corrdf <- merge(KBTBD11_all, drugall, by = "CELL_LINE_NAME") # merge expression data and ic50data

# get spearman correlation stats
n <- c(1:542)
p <- c() # empty vector for p values
r <- c() # empty vector for estimates
for (i in n)
{
  t_i <- cor.test(corrdf[, 2], corrdf[, i + 2], method = "spearman") # spearman correlation between exp and ic50
  p <- append(p, t_i$p.value)
  r <- append(r, t_i$estimate)
}

# generate dataframe for correlation results
druglist <- colnames(corrdf)[c(3:544)] # drug names
testdf <- cbind(druglist, r, p) # combine results of correlation analysis
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
testdf$group[which(testdf$log10.padj > 5)] <- "Top" # CHANGE manual threshold for label on plot, change if needed
testdf$sig <- "No"
testdf$sig[which(testdf$p.adjust < 0.05)] <- "Yes"

# add drug annotaiton to result table
colnames(testdf)[1] <- "DRUG_NAME"
testdf_ <- merge(drugpath, testdf)
testdf_ <- testdf_[order(testdf_$corr.coef, decreasing = F), ] # order result table by correlation estimates

write.csv(testdf_, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_corr.csv", row.names = F)

# generate volcano plot_met
pdf(file = "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/KBTBD11/KBTBD11_Corr_vln.pdf", width = 11, height = 9)
vln <- ggplot(testdf, aes(x = corr.coef, y = log10.padj)) +
  geom_point(aes(color = label), size = 3) +
  xlim(-0.4, 0.4) + # x and y range
  ylim(0, 10) +
  xlab("Spearman Correlation Coefficient") +
  ylab(expression(paste(-log[10], "(Q-value)"))) +
  ggtitle("Spearman correlation test between KBTBD11 expression and drug ln(IC50)") +
  theme_bw() +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(text = element_text(size = 10)) +
  theme(axis.text = element_text(color = "black", size = 10)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line()
  ) +
  theme(legend.position = "right") +
  labs(color = "padj") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_text_repel( # top drugs label
    data = subset(testdf, group == "Top"),
    aes(label = testdf$DRUG_NAME[which(testdf$group == "Top")]),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    max.overlaps = Inf
  ) +
  # change text according to x range
  annotate("text", x = -0.25, y = 0, label = c("bold(Sensitive)"), parse = T, size = 5) +
  annotate("text", x = 0.25, y = 0, label = c("bold(Resistant)"), parse = T, size = 5)

# print(vln)
dev.off()


# save the parameters for next time adjustment!
