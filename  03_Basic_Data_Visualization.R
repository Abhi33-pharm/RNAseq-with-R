#Load required packages 
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(GEOquery)
library(readxl)
library(naniar)

# Import data 

data<-read.csv("GSE183947_counts.csv")

# Visualization structure
ggplot(data=, aes(x=, y=))+ geom_type()

#  1. Bar plot

data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=sample, y=fpkm, fill = tissue))+
  geom_col()

# 2.Box plot


data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=metastasis, y=fpkm, fill = tissue))+
  geom_boxplot()

# 3. violin plot

data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=metastasis, y=fpkm, fill = tissue))+
  geom_violin()

# 4. Histogram

data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=fpkm, fill = tissue))+
  geom_histogram()

# Split 

data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=fpkm, fill = tissue))+
  geom_histogram()+
  facet_wrap(~tissue)

#5.Density plot 

data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=fpkm, fill = tissue))+
  geom_density()

# Split figure


data |> 
  filter(gene=="BRCA1") |>
  ggplot(aes(x=fpkm, fill = tissue))+
  geom_density()+
  facet_wrap(~tissue)


# 6. Scatter plot

data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |>
  spread(key = gene, value = fpkm) |> 
  ggplot(aes(x=BRCA1, y=BRCA2, colour =tissue))+
  geom_point()
  

# Add stats


data |> 
  filter(gene=="BRCA1" | gene=="BRCA2") |>
  spread(key = gene, value = fpkm) |> 
  ggplot(aes(x=BRCA1, y=BRCA2, colour =tissue))+
  geom_point()+
  geom_smooth(method = "lm", se= FALSE)


# Heatmap

gene_of_interest<- c("BRCA1", "BRCA2", "MYCN", "TP53")

data |> 
  filter(gene %in% gene_of_interest) |> 
  ggplot(aes(x=sample, y=gene, fill = fpkm))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")



















  
  
  
  
  
  
  
  