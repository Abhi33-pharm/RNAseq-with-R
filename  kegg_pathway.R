# Load packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
# kegg pathway

gene <- "2261"  # Entrez ID for FGFR3

kegg_result <- enrichKEGG(gene = gene,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)
head(kegg_result)

barplot(kegg_result, showCategory = 10)
dotplot(kegg_result)


