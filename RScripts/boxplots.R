# title: "Extracting Proteins of Interest (POIs) table from CCLE data"
# author: "Tania Gonzalez"
# date: "20221113"

# LOADING PACKAGES----
packagelist = c("readxl", "tidyverse", "glue", "stringr", "reshape2")

lapply(packagelist, function(x) {
  suppressMessages(require(x, character.only = T,
                           quietly=T, 
                           warn.conflicts = F))
}
)


# PARAMETERS SETUP----
date <- format(Sys.Date(), "%Y%m%d") # setting today's date in specific format: YYYYMMDD
analysis.id <- "EXP096"

out.path <- paste0("./", analysis.id, "/output")
data.path <- "../../datasets/ccle"

# DATA LOAD----
ccle_rna <- readRDS(glue("{data.path}/ccle_rna_expression.rds"))
ccle_protein <- readRDS(glue("{data.path}/ccle_tmt_protein_expression.rds"))

ccle_samples_info <- read.csv(file = glue("{data.path}/ccle_sample_info.csv"))


# Boxplot for RNA values
# labeling instructions
p <- c("GMCL1", "USP28", "TP53BP1", "APC","RB1", "TP53", "MYC", "CCND1")
POIs <- ccle_rna[, p, drop = F]
POIs$DepMap_ID <- rownames(POIs)

lab <- data.frame(genes = colnames(POIs[1:8]),
                  driver = c(rep("Unknown", 3), # "GMCL1"   "USP28"   "TP53BP1"
                             rep("Tumor Suppressor", 3), # "APC"     "RB1"  "TP53" 
                             rep("Oncogene", 2))) # "MYC"     "CCND1"  

# constructing input table for ggplot in the long format
input <- reshape2::melt(POIs, id.vars="DepMap_ID") %>% 
  mutate(driver = rep("Unknown", nrow(.))) 

input$driver <- ifelse(input$variable %in% c("MYC","CCND1"),
                       "Oncogene", input$driver) 

input$driver <- ifelse(input$variable %in% c("APC","RB1","TP53"),
                       "Tumor Suppressor", input$driver)

# constructing boxplot viz
ggplot(input, aes(x = variable, y = value, fill = driver)) + 
  
  geom_violin(show.legend = T) + # use geom_boxplot for boxplots format
  geom_jitter(size = 0.1, show.legend = F) + # resizing individual dots
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + # setting colors
  
  # re-naming x and y axis
  labs(x = NULL, y = "mRNA expression, log2(TPM +1)", 
       title = paste("CCLE cell lines, n =", nrow(POIs))
  )

# saving file
ggsave(filename = glue("{out.path}/ccle_mRNA_violin_plot.png"))


# END OF CODE #