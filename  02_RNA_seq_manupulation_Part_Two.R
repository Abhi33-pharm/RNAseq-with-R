#Load required packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)

# Import Raw Counts Data 

Counts_data<- read.csv("GSE183947_fpkm.csv")

# Get meta data

res<-getGEO(GEO = "GSE183947", GSEMatrix = T)
res
class(res)

# Metadata 

metadata<- pData(phenoData(res[[1]]))

# Subset metadata

metadata_subset<- metadata |> 
  select(c(1,10,11,17))

# Data pre processing

metadata_modified<- metadata_subset |> 
  rename(tissue=characteristics_ch1, metastasis=characteristics_ch1.1) |> 
  mutate(tissue= gsub("tissue: ", "", tissue)) |> 
  mutate(metastasis= gsub("metastasis: ", "", metastasis))
  

# Reshaping data
counts_data_long<-Counts_data |> 
  rename(gene=X) |> 
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm")

# Joining data 

counts_final_data<-counts_data_long |> 
  left_join(metadata_modified, by=c("sample"="description"))
  
# Export data
write.csv(counts_final_data,"GSE183947_counts.csv", row.names = F)








