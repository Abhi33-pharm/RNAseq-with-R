#Load required packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)

# Import metadata

options(timeout = 600)
options(download.file.method = "libcurl")
res <- getGEO("GSE32474", GSEMatrix = TRUE)

res
class(res)

# Metadata

metadata<-pData(phenoData(res[[1]]))

# Access the expression data

express<-exprs(res[[1]])

# Access the feature data

feature<-fData(res[[1]])

#Normalize the data

#Export the files

write.csv(metadata,"metadata.csv", row.names = F)
write.csv(express, "express.csv")
write.csv(feature, "feature.csv")









