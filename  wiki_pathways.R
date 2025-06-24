library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# WIKI pathway

gene <- "2261"

WIKI<-enrichWP(gene = gene,
               organism = "Homo sapiens")
dotplot(WIKI)+
  scale_fill_gradientn(colours = c("red", "yellow"))

barplot(WIKI)+
  scale_fill_gradient(low = "red",high = "yellow")
 