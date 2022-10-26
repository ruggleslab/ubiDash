# title: "Extracting GMCL1 table from CCLE protein data for Yuki"
# author: "Tania Gonzalez"
# date: "20220809"

# LOADING PACKAGES----
packagelist = c("readxl", "tidyverse", "glue", "stringr", "reshape2")

lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

# PARAMETERS SETUP----
date <- format(Sys.Date(), "%Y%m%d") # setting today's date in specific format: YYYYMMDD
# analysis.id <- "EXP096"

# DATA LOAD----
out.path <- paste0("./", analysis.id, "/output")
data.path <- "../datasets/depmap/ccle/"

file.name <- "CCLE_rna_expression.csv"
ccle_protein <- read.csv(file = glue("{data.path}{file.name}"), row.names = 1, check.names = F)

pr_names <- str_split(colnames(ccle_protein), pattern = " ")
names <- c()
for (i in seq_along(pr_names)) {
  pr <- unlist(pr_names[i])
  tmp <- pr[1]
  names <- c(names,tmp)
}
colnames(ccle_protein) <- names
# saveRDS(ccle_protein, file = glue("{data.path}CCLE_expression.rds")) # Only do once


file.name <- "sample_info.csv"
samples_info <- read.csv(file = glue("{data.path}{file.name}"))
x <- c("DepMap_ID", "CCLE_Name", "sample_collection_site", "primary_or_metastasis", "primary_disease", "Subtype")
samples_info_sub <- samples_info[, x]


file.name <- "Achilles_metadata.csv"
metadata <- read.csv(file = glue("{data.path}{file.name}"))
x <- c("DepMap_ID", "culture_medium")
metadata_sub <- metadata[, x]


# POIs info
p <- c("GMCL1", "USP28", "TP53BP1", "APC","RB1", "TP53", "MYC", "CCND1")
gmcl1 <- ccle_protein[, p, drop = F]
gmcl1$DepMap_ID <- rownames(gmcl1)

gmcl1 <- left_join(gmcl1, samples_info_sub, by = "DepMap_ID")
gmcl1 <- left_join(gmcl1, metadata_sub, by = "DepMap_ID")

gmcl1$median_expression <- matrixStats::rowMedians(as.matrix(ccle_protein), na.rm = F)

write_csv(gmcl1, file = glue("{out.path}/ccle_gmcl1_protein.csv"))

# Boxplot visualization
driver <- data.frame(genes = colnames(gmcl1[1:8]),
                     driver = c(rep("Unknown", 3), 
                                rep("Tumor Suppressor", 3),
                                rep("Oncogene", 2))) 
  
input <- reshape2::melt(gmcl1[1:9]) %>% 
  mutate(driver = rep("Unknown", nrow(input)))

input$driver <- ifelse(input$variable %in% c("MYC","CCND1"),
                       "Oncogene", input$driver)

input$driver <- ifelse(input$variable %in% c("APC","RB1","TP53"),
                       "Tumor Suppressor", input$driver)

ggplot(input, aes(x=variable, y= value, fill = driver)) + 
  geom_violin(show.legend = T) +
  geom_jitter(size = 0.1, show.legend = F) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "mRNA expression, log2(TPM +1)", 
       title = "CCLE cell lines, n = 1406")

ggsave(filename = glue("{out.path}/EXP096_ccle_gmcl1_mRNA.png"))
   

# END OF CODE #