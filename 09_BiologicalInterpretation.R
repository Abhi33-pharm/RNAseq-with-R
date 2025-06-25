# Load libraries
library(recount)
library(DESeq2)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(clusterProfiler)
library(tidyverse)
library(ggpubr)

# Problem: Effects of EWS-FLI1 KD in Ewing sarcoma
# RQ: What is the effect of knocking down EWS-FLI1 in Ewing sarcoma? How does it change cellular behavior?

# The dataset we use here is downloaded from recount2
# Find a project of interest (note: it's hit or miss, the database is not limitless)
project_info <- abstract_search(query = "Ewing sarcoma")

# select study
download_study("SRP015989")

# load the data
load("SRP015989/rse_gene.Rdata")

# examine the read counts
counts_data <- assay(rse_gene)

# examine the coldata
col_data <- as.data.frame(colData(rse_gene))

# fix the colData to give a column with the appropriate groups
col_data$condition <- c(rep("shCTR", 3), rep("ShEF1", 4))


# check colnames of count_data and rownames of col_data
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))

# are they in same order
all(colnames(counts_data) == rownames(col_data))

# perform DESeq2 analysis using DESeqDataSetFromMatrix (not suitable for this condition)
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# normalize the data
rld <- rlog(dds)
plotPCA(rld)

# perform DEG analysis
dds <- DESeq(dds)

# save as results
res <- results(dds)

# -- plotMA
plotMA(res)

# LFC shrink
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)


# -- plotMA with resNorm
plotMA(resNorm)

# Make a DF
resdf <- as.data.frame(resNorm)

# -- convert ENSG to gene symbol
str(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)

ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                      keys = keys(EnsDb.Hsapiens.v86), 
                      columns = c("SYMBOL"))

# -- wrangle the resdf and join with ens2sym map
anno_results <- resdf |> 
  rownames_to_column() |> 
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) |> 
  dplyr::select(-rowname) |> 
  inner_join(y = ens2sym, by = "GENEID") 

# -- volcano plot
EnhancedVolcano(
  anno_results, 
  lab = anno_results$SYMBOL,
  x = "log2FoldChange", 
  y = "padj"
)

# -- Get the over-expressed genes
anno_results |> 
  dplyr::filter(padj < 0.01 & log2FoldChange > 2) |> 
  write_csv(file = "tables/over_expressed_genes.csv")

# -- Get the under-expressed genes
anno_results |> 
  dplyr::filter(padj < 0.01 & log2FoldChange < -2) |> 
  write_csv(file = "tables/under_expressed_genes.csv")



