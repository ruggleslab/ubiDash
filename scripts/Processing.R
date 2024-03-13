# code for df cleaning and annotation

setwd("/Users/staceyvann/Desktop/Lab/Drug Sensitivity")
library(readxl)

###############################################################################################
# 1. read exp,drug,sample info
exp <- read.delim("Cell_line_RMA_proc_basalExp.txt", header = T)
drug1 <- read_xlsx("GDSC1_fitted_dose_response_24Jul22.xlsx")
drug2 <- read_xlsx("GDSC2_fitted_dose_response_24Jul22.xlsx")
sample <- read.delim("sample.txt", header = T)

###############################################################################################
# 2. merge drug data and create drug matrix
# merge drug data = longer df
length(unique(drug1$CELL_LINE_NAME)) # 970 cell lines
length(unique(drug1$DRUG_NAME)) # 378 drugs
length(unique(drug2$CELL_LINE_NAME)) # 969 cell lines
length(unique(drug2$DRUG_NAME)) # 286 drugs
sum(unique(drug1$CELL_LINE_NAME) %in% drug2$CELL_LINE_NAME) # 961 in 2
sum(unique(drug1$DRUG_NAME) %in% drug2$DRUG_NAME) # 122 in 2

drug <- rbind(drug1, drug2) # merge to a longer df
drug <- drug[, c(5, 9, 10, 11, 16, 7)] # select feature cell line name, TCGS_class, drugname, target, pathway, ln(ic50) and z-score
write.csv(drug, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drug_longdf.csv", row.names = T)
# Using drug longer df:
# drug <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drug_longdf.csv", row.names = 1)

###############################################################################################
# 3. generate drug ic50 matrix
drugdf <- drug[, c(1, 2, 5)] # select cell line, drugname, ic50
# wider df
drugmat <- drugdf %>% pivot_wider(
  names_from = DRUG_NAME,
  values_from = LN_IC50,
  values_fn = list(LN_IC50 = mean)
)
# change row names
drugmat <- as.data.frame(drugmat)
rownames(drugmat) <- drugmat$CELL_LINE_NAME
drugmat <- drugmat[, -1]
# ic50 matrix and df, 978 cell lines with 542 drugs
drugmat <- as.matrix(drugmat)
df <- as.data.frame(drugmat)
write.csv(df, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drugmat_all.csv", row.names = T)
# Using drugmat:
# drugmat <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drugmat_all.csv", row.names = 1)


###############################################################################################
# 4. sample and drug data merge = get annotations and sample ids
# get sample annotations
sample <- as.data.frame(sample)
# match sample id
sample$Source.Name <- gsub(".*_", "", sample$Source.Name)
sample$Source.Name <- sub("^", "DATA.", sample$Source.Name)
sample <- sample[, c(1, 3)] # select needed
colnames(sample)[1] <- "Data"
colnames(sample)[2] <- "CELL_LINE_NAME"

# match sample cell lines with drug data
drugcl <- unique(drug[, c(1, 6)]) # get cell lines and TCGA description
nomatch <- c(
  "786-o", "BCPAP", "Capan-2", "DIFI", "Hs746T", "HUTU-80", "JURKAT",
  "NB-TU-1-10", "NTERA-S-cl-D1", "Sci-1", "T-24"
) # cell line names in sample
a <- which(sample$CELL_LINE_NAME %in% nomatch)
sample$CELL_LINE_NAME[a] <- c(
  "786-0", "B-CPAP", "CAPAN-2", "DiFi", "Hs-746T", "HuTu-80", "Jurkat",
  "NB(TU)1-10", "NTERA-2-cl-D1", "SC-1", "T24"
)
# cell lines in common between sample & drug data = 933 cell lines
cell <- merge(drugcl, sample)
nomatch <- cell[which(is.na(cell$TCGA_DESC)), ] # these are with both data but no TCGA description
write.csv(cell, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/MatchedCellAnno.csv")
# Use cell line annotation
# sample <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/MatchedCellAnno.csv", row.names = 1)


###############################################################################################
# 5. drug annotation df
drugpath <- drug[, c(2, 3, 4)] # keep drugname, pathway, targets
drugpath <- unique(drugpath)
drugpath$PATHWAY_NAME <- as.factor(drugpath$PATHWAY_NAME)
which(duplicated(drugpath$DRUG_NAME) == T)
drugpath <- drugpath[-c(which(duplicated(drugpath$DRUG_NAME) == T)), ] # remove duplicated drugs
# save expression table
write.csv(drugpath, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drugpath.csv", row.names = T)
# Use drug pathway annotation:
# drugpath <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/drugpath.csv", row.names = 1)

###############################################################################################
# 6. tidy expression sample
# tidy exp data, CNDP1 as example
##### replace CNDP1 with target gene by Find/Replace######30 occurrence#
CNDP1_exp <- exp[exp$GENE_SYMBOLS == "CNDP1", ]
CNDP1_exp <- t(CNDP1_exp[, -c(1:2)]) # transpose
CNDP1_exp <- as.data.frame(CNDP1_exp)
CNDP1_exp$Data <- rownames(CNDP1_exp) # add column of sample IDs
colnames(CNDP1_exp)[1] <- "Exp"
rownames(CNDP1_exp) <- NULL

# match sample annotations with expression table
CNDP1 <- merge(CNDP1_exp, cell, by = "Data")
# take average of duplicated cell lines
dupcl <- CNDP1$CELL_LINE_NAME[which(duplicated(CNDP1$CELL_LINE_NAME))] # duplicated cell line names
dup <- CNDP1[which(CNDP1$CELL_LINE_NAME %in% dupcl), ] # duplicated df
CNDP1 <- CNDP1[-c(which(CNDP1$CELL_LINE_NAME %in% dupcl)), ] # remove duplicated
dup <- dup[order(dup$CELL_LINE_NAME), ] # order by cell line names
dup_ <- dup[c(T, F), ] # keep only one result for each cell line

for (i in c(1:(dim(dup)[1] / 2))) { # calculate average of each pair and fill the column
  m <- (dup$Exp[2 * i - 1] + dup$Exp[2 * i]) / 2
  dup_$Exp[i] <- m
}
# add duplicated cell lines back to df
CNDP1 <- rbind(CNDP1, dup_)
CNDP1 <- CNDP1[order(CNDP1$Exp, decreasing = T), ] # reorder expression
# save expression table
write.csv(CNDP1, "/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/CNDP1.csv", row.names = T)
# Use expression table: 928 cell lines
# CNDP1 <- read.csv("/Users/staceyvann/Desktop/Lab/Drug Sensitivity/DataProcessing/CNDP1.csv", row.names = 1)
