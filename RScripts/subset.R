library(tidyverse);
library(reshape2)
library(dplyr)
library(readxl)

filepath = "./data/DE_protien.xlsx";
sheets = excel_sheets(path = filepath)


#Generate Matrices----
logFCTable = NULL;
PValueTable = NULL;

for(i in sheets){
  cohortSheet = read_excel(path=filepath, sheet = i)
  tempfc = cbind(cohortSheet$gene, cohortSheet$logFC)
  colnames(tempfc) = c("gene", i)
  tempp = cbind(cohortSheet$gene, cohortSheet$sca.adj.pval)
  colnames(tempp) = c("gene", i)
  
  if(is.null(logFCTable)){
    logFCTable =as.data.frame(cbind(cohortSheet$gene, cohortSheet$logFC))
    colnames(logFCTable) = c("gene", i)
    
    PValueTable =as.data.frame(cbind(cohortSheet$gene, cohortSheet$sca.adj.pval))
    colnames(PValueTable) = c("gene", i)
  }else{
    logFCTable =full_join(logFCTable, as.data.frame(tempfc), by = "gene");
    PValueTable =full_join(PValueTable, as.data.frame(tempp), by = "gene");
  }
  logFCTable[[i]] = as.numeric(logFCTable[[i]])
}


genes <- read_excel(path="./data/AllCRLs.v4.xlsx", sheet = 1)

source("./RScripts/table.R")


tempfc = as.data.frame(logFCTable)
tempp = as.data.frame(PValueTable)

tempfc = filter(tempfc, gene %in% genes$hgnc_symbol)
tempp = filter(tempp, gene %in% genes$hgnc_symbol)

rownames(tempfc) = tempfc$gene
tempfc$gene = NULL

rownames(tempp) = tempp$gene
tempp$gene = NULL


# https://stackoverflow.com/questions/9856632/subset-rows-with-1-all-and-2-any-columns-larger-than-a-specific-value

for(i in rownames(tempfc)){
  
  if(is.na(any(as.numeric(tempp[i,]) < .05))){
    tempp = tempp[ !rownames(tempp) == i, ]
  }
  
  if(is.na(any(abs(as.numeric(tempfc[i,])) > .5))){
    tempfc = tempfc[ !rownames(tempfc) == i, ]
  }
}

tempfc = as.data.frame(rownames(tempfc))
colnames(tempfc) = "genes"

tempp  = as.data.frame(rownames(tempp))
colnames(tempp) = "genes"

genesSubset = inner_join(tempfc, tempp, by = "genes")


save(genesSubset, file =  "./data/genesSubset.RData")
rm(tempfc, tempp, cohortSheet, genes, logFCTable, PValueTable, genesSubset)


