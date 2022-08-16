if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt",force = TRUE)
BiocManager::install("SummarizedExperiment",force = TRUE)
library("biomaRt")

source("http://bioconductor.org/biocLite.R")
library(ensembldb)

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
hasProteinData(edb)
listTables(edb)


akid <- read.csv(file= "C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA-3/AKID-Ensembl.csv",sep=";")
df <- as.data.frame(akid)

df_ids<- select(edb, keys=row.names(df), keytype ="ENSEMBL",columns = "UNIPROTID") 
