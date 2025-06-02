# Load packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)

# Import data
express<-read.csv("express.csv")
feature<-read.csv("feature.csv")
metadata<-read.csv("metadata.csv")

# Data pre processing 

express_subsetting<-express |> 
  select(c(1,4, 63, 122))
gene_symbol<-feature |> 
  select(1,12)

# Joining the data

counts_final_data<-gene_symbol |> 
  left_join(express_subsetting,by="X")

head(counts_final_data)

# Creating a new column
total_counts<-counts_final_data |> 
  mutate(fpkm=rowMeans(across(c(GSM803617,GSM803676,GSM803735)), na.rm = T))
  
# Creating final data frame

modified_counts<-total_counts |> 
  filter(fpkm >=13) |> 
  select(c(2,6))

# Export the data
write.csv(modified_counts,"modified_counts.csv", row.names = F)













