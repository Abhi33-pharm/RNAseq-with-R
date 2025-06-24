# Load Packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


# Reactome 

gene<-"2261"
reactome<-enrichPathway(gene = gene,
                        organism = "human",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
head(reactome)
dotplot(reactome)

options(enrichplot.colours = c("red", "blue"))


