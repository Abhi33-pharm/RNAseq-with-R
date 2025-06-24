library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# Go biological process

# FGFR3 Entrez ID

gene <- "2261"

# GO Biological Process (BP)
ego_bp <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

head(ego_bp)
barplot(ego_bp, showCategory = 10, title = "GO Biological Process",
        color = "pvalue", "p.adjust")+
  scale_fill_gradient(low = "yellow",high = "blue")

dotplot(ego_bp, showCategory = 10, title = "GO Biological Process",
        color = "pvalue", "p.adjust")+
  scale_fill_gradient(low = "yellow",high = "blue")



# GO Molecular Function (MF)
ego_mf <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pvalueCutoff = 0.05)

barplot(ego_mf, showCategory = 10, title = "GO Molecular Function")

# GO Cellular Component (CC)

ego_cc <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "cc",
                   pvalueCutoff = 0.05)

head(ego_cc)

barplot(ego_cc, showCategory = 10, title = "GO Molecular Function",
        color = "pvalue", "p.adjust")+
  scale_fill_gradient(low = "blue",high = "red")
