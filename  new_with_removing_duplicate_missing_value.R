#Load required packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)


# Import dataframe 

data<-read.csv("feature_wide.csv")
gg_miss_var(data)

data<-na.omit(data)

res<-getGEO("GSE10245", GSEMatrix = TRUE)
class(res)
metadata<- pData(phenoData(res[[1]]))


# select and rename variables of metadata

metadata_modified<-metadata |> 
  rename(Sample=geo_accession) |> 
  select(2,10,34)



head(metadata_modified)

#Clean data and Removing duplicate value

clean_data<-data |> 
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") |> 
  distinct(Gene_Symbol, .keep_all = TRUE)

# Reshape clean data

clean_data_long<-clean_data |> 
  pivot_longer(cols = c(-ID,-Gene_Symbol), names_to = "Sample", values_to = "fpkm")


# Join data

data_join<-clean_data_long |> 
  left_join(metadata_modified,by="Sample") |> 
  select(-1) |> 
  rename(dissease_state=)

# Export data
write.csv(data_join,"data_join.csv", row.names = FALSE)
write.csv(metadata,"metadata.csv",row.names = FALSE)

 