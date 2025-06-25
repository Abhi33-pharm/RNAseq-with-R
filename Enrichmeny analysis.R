# Installation of Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")


# Load Packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


# Upload Data
data <- read_table("Data/table_similar(1).txt") |>
  dplyr::select(GeneSymbol)

# Data Formatting
data_entrez <- bitr(data$GeneSymbol,
     fromType = "SYMBOL",
     toType = "ENTREZID",
     OrgDb = "org.Hs.eg.db")

# Gene Ontology

GO_BP <- enrichGO(data_entrez$ENTREZID,
         OrgDb = "org.Hs.eg.db",
         ont = "MF",
         pvalueCutoff = 0.05,
         pAdjustMethod = "BH")
dotplot(GO_BP)


# KEGG Enrichment

kegg <- enrichKEGG(gene = data_entrez$ENTREZID,
           organism = "hsa",
           pvalueCutoff = 0.05)

dotplot(kegg)


# Wikipathway

WIKI <- enrichWP(gene = data_entrez$ENTREZID,
         organism = "Homo sapiens")
dotplot(WIKI)


# Reactome

Reactome <- enrichPathway(data_entrez$ENTREZID,
              organism = "human",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH")
dotplot(Reactome)

options(enrichplot.colours = c("red", "blue"))
