# packages for data manipulation
library(tidyverse)

# packages for DEG
library(airway)
library(DESeq2)

# packages for visualization 
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)

# package for poison distance calculation
library(PoiClaClu)

# package for gene filter 
library(genefilter)

# packages for gene annotations 
library(org.Hs.eg.db)
library(AnnotationDbi)

# know our data 
?airway

# load data from "airway" dataset 
data("airway")
class(airway)

# get the raw counts_data 
counts_data <- assay(airway)

# get the col_data / metadata  
col_data <- as.data.frame(colData(airway))

# making sure the row names in col_data matches to column names in counts_data
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))

# are they in the same order?
all(colnames(counts_data) == rownames(col_data))

# Construct a DESeqDataSetFromMatrix data object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data, 
  colData = col_data, 
  design = ~dex
)

# Plot heatmap of poison distances between samples 
# Poison distance for raw (non-normalized) count data
# Use Euclidean distance for data normalized by regularized-logarithm transformation (rlog) or variance stablization transfromation (vst)
counts(dds)
t(counts(dds))

poisd <- PoissonDistance(t(counts(dds)))
sample_poison_matrix <- as.matrix(poisd$dd)
rownames(sample_poison_matrix) <- paste(dds$dex, sep = "-")
colnames(sample_poison_matrix) <- NULL 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sample_poison_matrix, 
  clustering_distance_rows = poisd$dd, 
  clustering_distance_cols = poisd$dd, 
  color = colors
)

# perform differential gene expression analysis 
dds <- DESeq(dds)

# obtain results 
res <- results(dds)
names(res)
res_df <- as.data.frame(res)

# plot average expression versus log2 fold change - points are colored blue if Padj < 0.1
plotMA(res)
plotMA(res, ylim = c(-4,4))

# plot histogram of P-values
hist(res$pvalue, breaks = 20, col = "grey50", border = "white")

# plot histogram of padj values
hist(res$padj, breaks = 20, col = "grey50", border = "white")

# plot histogram of P-values - improved version by filtering out genes with very low expression levels
hist(res$pvalue[res$baseMean > 1], breaks = 20, col = "grey50", border = "white")

# Add gene annotation to results 
rownames(res)
colnames(res)

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              rownames(res), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")
# check colnames 
names(anno)

# convert results into data frame 
res <- cbind(ENSEMBL = rownames(res), res)
res_df <- as.data.frame(res)
names(res_df)

# join annotation tables to results 
anno_results <- left_join(res_df, anno, by = "ENSEMBL")

# Volcano Plot
# Default cutoffs are log2FC > |2| and adjusted P-value < 0.05
EnhancedVolcano(
  anno_results, 
  lab = anno_results$SYMBOL, 
  x = "log2FoldChange", 
  y = "padj"
)

# Add custom log2FC and adjusted P-value cutoffs and size of points and labels 
EnhancedVolcano(
  anno_results, 
  lab = anno_results$SYMBOL, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  title = "Untreated vs. Treated"
  
)

 
# Adjust axis limits 
EnhancedVolcano(
  anno_results, 
  lab = anno_results$SYMBOL, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  xlim = c(-5, 5), 
  ylim = c(0, -log10(10e-10)), 
  title = "Untreated vs. Treated"
)

# Modify border and remove grid lines 
EnhancedVolcano(
  anno_results, 
  lab = anno_results$SYMBOL, 
  x = "log2FoldChange", 
  y = "padj", 
  pCutoff = 0.001, 
  FCcutoff = 2, 
  pointSize = 1.5, 
  labSize = 3.0, 
  xlim = c(-5, 5), 
  ylim = c(0, -log10(10e-10)), 
  border = "full", 
  borderWidth = 1.5, 
  borderColour = "black", 
  gridlines.major = FALSE, 
  title = "Untreated vs. Treated"
)

# perform regularized-logarithm transformation (rlog) on the data
rld <- rlog(dds)

# plot Principal Component Analysis
plotPCA(rld, intgroup = "dex")

# subset genes 
results_sig <- anno_results[which(anno_results$padj < 0.01 & abs(anno_results$log2FoldChange) >= 1 & anno_results$baseMean >= 20), ]

# DE genes with strongest downregulation (head)
downregulation <- head(results_sig[order(results_sig$log2FoldChange ), ])

# DE genes with strongest upregulation(tail)
upregulation <- tail(results_sig[order(results_sig$log2FoldChange ), ])

# plot expression of individual genes 
# gene with largest positive log2FC 
plotCounts(dds, gene = which.max(anno_results$log2FoldChange), intgroup = "dex")
plotCounts(dds, gene = which.min(anno_results$log2FoldChange), intgroup = "dex")

# specific gene of interest 
plotCounts(dds, gene = "ENSG00000127954", intgroup = "dex")