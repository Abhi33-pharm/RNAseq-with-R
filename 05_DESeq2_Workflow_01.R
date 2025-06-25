# install required bioconductor packages
BiocManager::install(c("airway", "DESeq2", "EnhancedVolcano"))

# install required packages
install.packages(c("tidyverse", "ggthemes"))

# load required packages
library(tidyverse)
library(airway)
library(DESeq2)

# get info about the data
?airway

# get the data
data(airway)

# counts table
counts_data <- assay(airway)

# metadata (coldata)
col_data <- as.data.frame(colData(airway))

# making sure the row names in `col_data` (metadata) matches to columns names in `counts_data`
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))

# are they in same order?
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) == rownames(col_data))

# construct a DESeqDataSetFromMatrix data object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~dex # condition
)

# pre-filtering: removing rows with low gene counts
# keep rows that have at least 10 reads total
rowSums(counts(dds))
rowSums(counts(dds)) >= 10

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

# reference category ~ set the factor level
dds$dex <- relevel(dds$dex, ref = "untrt")

# perform differential gene expression analysis
dds <- DESeq(dds)

# save as results
res <- results(dds)
res

# exploring the results
summary(res)

# working with alpha (significance level)
# if alpha = 0.1 or 10% ~ 90% CI (default)
# if alpha = 0.01 or 1% ~ 99% CI
res_0.01 <- results(dds, alpha = 0.01)
summary(res_0.01)

# if alpha = 0.05 or 5% ~ 95% CI
res_0.05 <- results(dds, alpha = 0.05)
summary(res_0.05)

# results name
resultsNames(dds)

# contrast
contrast_res <- results(dds, contrast = c("dex", "trt", "untrt"))
contrast_res
summary(contrast_res)










