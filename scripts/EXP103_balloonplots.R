# title: "CRL correlation comparison between PRIDE-CPTAC-CCLE datasets"
# author: "Tania Gonzalez"
# date: "20221201"

# ENVIRONMENT SETUP----
source("../../scripts/envs/tgr_misc.R")
source("../../scripts/envs/env-parameterSetup.R")
analysis.id <- "EXP103"

# change folder name here
dir.create(glue("./output/{date}"), recursive = T) ########
out.path <- glue("./output/{date}") ########

# DATA----
# labeling files
label_row <- read.csv(file = "../EXP102/input/label_row.csv", row.names = 1)
label_col <- read.csv(file = "../EXP102/input/label_col.csv")
colors_tab <- read.csv(file = "../EXP102/input/colors.csv", row.names = 1)
colors <- setNames(colors_tab$group_colors, colors_tab$group)

# crl lists
crls <- rownames(subset(label_row, E3_Complex != "Controls"))
# crls <- rownames(label_row)
# crls <- unique(read_excel(path = "./input/AllCRLs.v5.xlsx")[["hgnc_symbol"]])
# e3s <- unique(read_excel(path = "./input/20220607_EXP083_All_E3_Ligases.xlsx")[["hgnc_symbol"]])

# datasets
pride <- readRDS("../../ubiDash/data/pride_lfq_protein_expression.rds")
# cptac <- readRDS("../../ubiDash/data/cptac_tmt-log2-AI_protein_expression.rds")
cptac <- readRDS("../../ubiDash/data/cptac_tmt-log2-RI_protein_expression.rds")
ccle <- readRDS("../../ubiDash/data/ccle_tmt_protein_expression.rds")
# ccle_lfq <- readRDS("../../ubiDash/data/ccle_lfq_protein_expression.rds")

# data pre-processing
pride_norm <- data_norm(pride, row_center = T, col_center = T); 
cptac_norm <- cptac # if using RI format
# cptac_norm <- data_norm(cptac, row_center = T, col_center = T); 
ccle_norm <- data_norm(ccle, row_center = F,  col_center = T)

# pride_n
pride_n <- pride_norm[rownames(pride_norm) %in% crls,]

# cptac_n
normals <- subset(label_col, Dataset == "CPTAC" & Sample_type == "Normal",
                  select = c(Sample_id, Tissue_type))
cptac_n <- aggregate_by(input = cptac_norm[rownames(cptac_norm) %in% crls, 
                                           normals$Sample_id],
                        label = normals$Tissue_type)

# cptac_t
tumors <- subset(label_col, Dataset == "CPTAC" & Sample_type == "Tumor",
                 select = c(Sample_id, Tissue_type, Cancer_type))
cptac_t <- aggregate_by(input = cptac_norm[rownames(cptac_norm) %in% crls, 
                                           tumors$Sample_id],
                        label = tumors$Tissue_type)

# ccle_t
ccle_lab <- subset(label_col, Dataset == "CCLE_TMT" ,
                   select = c(Sample_id, Tissue_type, Cancer_type))
ccle_t <- aggregate_by(input = ccle_norm[rownames(ccle_norm) %in% crls, 
                                         ccle_lab$Sample_id], 
                       label = ccle_lab$Tissue_type)


# correlation test aggregating by tissues
# pride_n; cptac_n; cptac_t; ccle_t
lst <- list(rownames(pride_n), rownames(cptac_n), rownames(cptac_t), rownames(ccle_t))
common.proteins <- Reduce(intersect, lst) # 163 common proteins

lst <- list(colnames(pride_n), colnames(cptac_n), colnames(cptac_t), colnames(ccle_t))
common.tissues <-  Reduce(intersect, lst)

pride_n <- pride_n[common.proteins, common.tissues]
cptac_n <- cptac_n[common.proteins, common.tissues]
cptac_t <- cptac_t[common.proteins, common.tissues]
ccle_t <- ccle_t[common.proteins, common.tissues]


# corelation on normals----
# cor.mat <- cor(pride_n, cptac_n, use = "pairwise.complete.obs", method = "spearman")
# p.mat <- cor.mtest(mat = cor.mat, method = "spearman")

my.dat <- cbind(pride_n, cptac_n) 
colnames(my.dat) <- make.unique(colnames(my.dat))
my.dat <- as.tibble(my.dat, rownames = NA)

rows <- common.tissues 
cols <- sub("$", ".1", common.tissues)

cor.mat <- data.frame(cor_mat(my.dat, method = "spearman"), row.names = "rowname")
p.mat <- data.frame(cor_pmat(my.dat, method = "spearman"), row.names = "rowname")

cor.mat <- cor.mat[rows, cols]; p.mat <- p.mat[rows, cols]

df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat, diag = T)
df$log_pval <- -log10(df$p)
df$log_pval_test <- -log10(df$p + 0.000001)
df$qval <- p.adjust(df$p, method = "BH")
df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1,2)]

ball.plot <- ggballoonplot(
  df,
  x = "row", #pride
  y = "column", #cptac_n
  size = "log_pval_test",
  facet.by = NULL,
  size.range = c(1, 10),
  shape = 21,
  color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
  fill = "Spearman",
  font.label = list(size = 12, color = "black"),
  rotate.x.text = TRUE) 
#F67280
ggpar(ball.plot, palette = c("#8675A9", "#C06C84"))
ggsave(filename = new_unique(name = glue("{out.path}/correlogram_p-val_pride_n__cptac_n"), suffix = ".svg"),
       device = "svg",
       width = 6.35,
       height = 5.35,
       units = "in")


# corelation on normals vs tumor----
# cor.mat <- cor(cptac_t, cptac_n, use = "pairwise.complete.obs", method = "spearman")
# p.mat <- cor.mtest(mat = cor.mat, method = "spearman")

my.dat <- cbind(cptac_t, cptac_n) 
colnames(my.dat) <- make.unique(colnames(my.dat))
my.dat <- as.tibble(my.dat, rownames = NA)

rows <- common.tissues 
cols <- sub("$", ".1", common.tissues)

cor.mat <- data.frame(cor_mat(my.dat, method = "spearman"), row.names = "rowname")
p.mat <- data.frame(cor_pmat(my.dat, method = "spearman"), row.names = "rowname")

cor.mat <- cor.mat[rows, cols]; p.mat <- p.mat[rows, cols]

df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat, diag = T)
df$log_pval <- -log10(df$p)
df$log_pval_test <- -log10(df$p + 0.000001)
df$qval <- p.adjust(df$p, method = "BH")
df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1,2)]

ball.plot <- ggballoonplot(
  df,
  x = "row", #pride
  y = "column", #cptac_n
  size = "log_pval_test",
  facet.by = NULL,
  size.range = c(1, 10),
  shape = 21,
  color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
  fill = "Spearman",
  font.label = list(size = 12, color = "black"),
  rotate.x.text = TRUE) 
# c("#389393", "#F5A25D")
ggpar(ball.plot, palette = c("#8675A9", "#C06C84"))
ggsave(filename = new_unique(name = glue("{out.path}/correlogram_p-val_cptac_t__cptac_n"), suffix = ".svg"),
       device = "svg",
       width = 6.35,
       height = 5.35,
       units = "in")



# corelation on tumors----
# cor.mat <- cor(cptac_t, ccle_t, use = "pairwise.complete.obs", method = "spearman")
# p.mat <- cor.mtest(mat = cor.mat, method = "spearman")

my.dat <- cbind(cptac_t, ccle_t) 
colnames(my.dat) <- make.unique(colnames(my.dat))
my.dat <- as.tibble(my.dat, rownames = NA)

rows <- common.tissues 
cols <- sub("$", ".1", common.tissues)

cor.mat <- data.frame(cor_mat(my.dat, method = "spearman"), row.names = "rowname")
p.mat <- data.frame(cor_pmat(my.dat, method = "spearman"), row.names = "rowname")

cor.mat <- cor.mat[rows, cols]; p.mat <- p.mat[rows, cols]

df <- flattenCorrMatrix(cormat = cor.mat, pmat = p.mat, diag = T)
df$log_pval <- -log10(df$p)
df$log_pval_test <- -log10(df$p + 0.000001)
df$qval <- p.adjust(df$p, method = "BH")
df$Spearman <- c("Positive", "Negative")[ifelse(df$cor > 0, 1,2)]

ball.plot <- ggballoonplot(
  df,
  x = "row", #pride
  y = "column", #cptac_n
  size = "log_pval_test",
  facet.by = NULL,
  size.range = c(1, 10),
  shape = 21,
  color = c("black", "transparent")[ifelse(df$qval < 0.01, 1, 2)],
  fill = "Spearman",
  font.label = list(size = 12, color = "black"),
  rotate.x.text = TRUE) 

ggpar(ball.plot, palette = c("#8675A9", "#C06C84"))
ggsave(filename = new_unique(name = glue("{out.path}/correlogram_p-val_cptac_t__ccle_t"), suffix = ".svg"),
       device = "svg",
       width = 6.35,
       height = 5.35,
       units = "in")

# correlations grouping by cancer type----
# # cptac_n
# normals <- subset(label_col, Dataset == "CPTAC" & Sample_type == "Normal",
#                   select = c(Sample_id, Cancer_type))
# cptac_n <- aggregate_by(input = cptac_norm[rownames(cptac_norm) %in% crls, 
#                                            normals$Sample_id],
#                         label = normals$Cancer_type)
# 
# # cptac_t
# tumors <- subset(label_col, Dataset == "CPTAC" & Sample_type == "Tumor",
#                  select = c(Sample_id, Tissue_type, Cancer_type))
# cptac_t <- aggregate_by(input = cptac_norm[rownames(cptac_norm) %in% crls, 
#                                            tumors$Sample_id],
#                         label = tumors$Cancer_type)
# 
# # ccle_t
# ccle_lab <- subset(label_col, Dataset == "CCLE_TMT" ,
#                    select = c(Sample_id, Tissue_type, Cancer_type))
# ccle_t <- aggregate_by(input = ccle_norm[rownames(ccle_norm) %in% crls, 
#                                          ccle_lab$Sample_id], 
#                        label = ccle_lab$Cancer_type)
# 
# 
# # correlations
# # cptac_n; cptac_t; ccle_t
# 
# lst <- list(rownames(cptac_n), rownames(cptac_t), rownames(ccle_t))
# common.proteins <- Reduce(intersect, lst) # 183 common proteins
# 
# lst <- list(colnames(cptac_n), colnames(cptac_t), colnames(ccle_t))
# common.tissues <-  Reduce(intersect, lst)
# 
# cptac_n <- cptac_n[common.proteins, common.tissues]
# cptac_t <- cptac_t[common.proteins, common.tissues]
# ccle_t <- ccle_t[common.proteins, common.tissues]

### END OF CODE ###


