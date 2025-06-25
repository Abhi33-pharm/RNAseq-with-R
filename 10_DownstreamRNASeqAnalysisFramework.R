# List of CRAN packages to be installed
packages_to_install <- c(
  "tidyverse",   # Collection of packages for data science
  "ggpubr",      # Simplifying plotting with ggplot2
  "msigdbr",     # Molecular signatures database
  "pheatmap",    # Pretty heatmaps
  "PoiClaClu",    # Poisson Clustering Model
  "RColorBrewer"  # Colors
)

# Install the packages using install.packages()
install.packages(packages_to_install)

# Ensure the BiocManager package is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages to be installed
bioconductor_packages <- c(
  "DESeq2",           # Differential gene expression analysis
  "EnhancedVolcano",  # Creating Volcano plots
  "recount",          # Accessing recount project data
  "apeglm",           # Accurate estimation of GLM coefficients
  "EnsDb.Hsapiens.v86", # Ensemble database for Homo sapiens
  "clusterProfiler",  # Statistical analysis and visualization of functional profiles
  "org.Hs.eg.db",     # Genome wide annotation for Human
  "AnnotationDbi",    # Annotation database interface
  "genefilter"        # Methods for filtering genes
)

# Install the Bioconductor packages
BiocManager::install(bioconductor_packages)

# Load libraries for data manipulation and visualization
library(tidyverse)       # Collection of data science packages (ggplot2, dplyr, etc.)
library(ggpubr)          # Convenient interface for ggplot2
library(RColorBrewer)    # Color palettes for graphics

# Load libraries for high-throughput genomic data analysis
library(recount)         # Access and analyze recount project data
library(DESeq2)          # Differential gene expression analysis
library(EnhancedVolcano) # Enhanced volcano plots for visualizing differential expression

# Load libraries for genomic data annotation and management
library(EnsDb.Hsapiens.v86) # Ensemble-based annotations for Homo sapiens
library(org.Hs.eg.db)       # Genome wide annotation for Human, from OrgDb
library(AnnotationDbi)      # Interface to multiple annotation databases

# Load libraries for pathway and gene set analysis
library(msigdbr)        # Molecular signatures database
library(clusterProfiler) # Statistical analysis and visualization of functional profiles

# Load libraries for data visualization and filtering in bioinformatics
library(pheatmap)       # Pretty heatmaps
library(genefilter)     # Methods for filtering genes based on expression levels

# Load library for specific statistical methods
library(PoiClaClu)      # Poisson Clustering Model


# Search for projects related to "Ewing sarcoma" using a hypothetical abstract search function
project_info <- abstract_search(query = "Ewing sarcoma")

# View the search results in a viewer window to select a relevant study
View(project_info)

# Manually select a study identifier based on review of project information
selected_study <- "SRP015989"

# Download the data for the selected study using a hypothetical download function
download_study(selected_study)

# Load the downloaded gene-level RNA-seq expression data
load("SRP015989/rse_gene.Rdata")

# Extract the read count matrix from the loaded SummarizedExperiment object
counts_data <- assay(rse_gene)

# View the count matrix to inspect its structure and data
View(counts_data)

# Convert and view the sample metadata (colData) to DataFrame to prepare for DESeq2 analysis
col_data <- as.data.frame(colData(rse_gene))

# View the column data to check sample information and experimental design
View(col_data)

# Manually define experimental conditions in the dataset for differential expression analysis
# shCTR and shEF1 represent different experimental groups
rse_gene$condition <- c(rep("shCTR", 3), rep("shEF1", 4))

# Convert the 'condition' column to a factor to ensure proper handling in DESeq2
rse_gene$condition <- as.factor(rse_gene$condition)

# Create a DESeqDataSet from the SummarizedExperiment object
# This function prepares the data for differential expression analysis using DESeq2
# 'design = ~condition' specifies that the analysis should consider the 'condition' variable
dds <- DESeqDataSet(rse_gene, design = ~condition)


# Calculate Poisson distances between samples using raw count data
# Poisson distance is suitable for raw count data as it accounts for the nature of count variability
poisd <- PoissonDistance(t(counts(dds)))

# Convert the distance object to a matrix format for visualization
samplePoisDistMatrix <- as.matrix(poisd$dd)

# Label the rows of the matrix with a combination of 'dex' and 'cell' sample attributes
# This labeling provides clear identification of samples in the heatmap
rownames(samplePoisDistMatrix) <- paste(dds$condition, sep=" - ")

# Remove column names for cleaner appearance
colnames(samplePoisDistMatrix) <- NULL

# Define a color palette using RColorBrewer to enhance visual appeal
# 'Blues' color palette is reversed to have dark colors represent smaller distances
colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot the heatmap using pheatmap function
# The distances used for clustering rows and columns are set directly from the Poisson distance calculation
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


# Run the DESeq function to perform differential expression analysis
# This function will use the information in 'dds' and the design formula to estimate size factors,
# dispersions, and perform negative binomial GLM fitting and Wald statistics
dds <- DESeq(dds)

# Extract results from the DESeq2 object
# The 'contrast' argument specifies which conditions to compare; here comparing "shEF1" vs "shCTR"
res <- results(dds, contrast = c("condition", "shEF1", "shCTR"))

# Ensure that the 'res' object contains the results from DESeq2 analysis
# Example DESeq2 results fetching (make sure to replace with your actual analysis if different):
# res <- results(dds, contrast = c("condition", "shEF1", "shCTR"))
# Plot an MA plot using the DESeq2 results
# 'ylim' sets the y-axis limits to range between -4 and 4, which helps in focusing on significant changes
plotMA(res, ylim = c(-4, 4))

# Ensure that 'res' (or another variable) contains the results from DESeq2 analysis
# Example DESeq2 results fetching (replace or confirm as needed):
# res <- results(dds, contrast = c("condition", "shEF1", "shCTR"))

# Plot a histogram of p-values from the DESeq2 results
# 'breaks' controls the number of bins in the histogram
# 'col' sets the color of the bars, and 'border' defines the color of the bin borders
hist(res$pvalue, breaks=20, col="grey50", border="white")


# Plot a histogram of p-values, excluding genes with a base mean expression level below 1
# This helps to remove genes that are lowly expressed and less likely to be biologically significant
# 'breaks' defines the number of bins in the histogram, improving the granularity of the data representation
# 'col' sets the fill color of the bars, and 'border' defines the color of the bin borders for better visual distinction
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")

# Print structure of the EnsDb.Hsapiens.v86 database
str(EnsDb.Hsapiens.v86)

# Print available columns in the EnsDb.Hsapiens.v86 database
columns(EnsDb.Hsapiens.v86)

# Retrieve and preview available keys in the EnsDb.Hsapiens.v86 database
keys_preview <- head(keys(EnsDb.Hsapiens.v86))
print(keys_preview)

# Select SYMBOL annotations for a subset of keys to ensure the operation completes in a reasonable time
sample_keys <- head(keys(EnsDb.Hsapiens.v86), 100)  # Adjust number as needed
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = sample_keys,
                                 columns = c("SYMBOL"), keytype = "GENEID")


# Convert rownames to a column, clean up the gene IDs, and join with ens2sym map
resdf <- resdf %>%
  # Convert rownames to a column named 'rowname'
  rownames_to_column(var = "rowname") %>%
  # Create 'GENEID' by removing version numbers from 'rowname'
  mutate(GENEID = gsub("\\..+", "", rowname)) %>%
  # Remove the old 'rowname' column as it's no longer needed
  select(-rowname) %>%
  # Join with the ens2sym dataframe on 'GENEID'
  inner_join(ens2sym, by = "GENEID")


# Perform a regularized log transformation of count data
# This transformation stabilizes variance across the range of mean values
rld <- rlog(dds)

# Plot Principal Component Analysis (PCA) to examine data clustering and batch effects
# PCA plot helps to visualize the overall effect of experimental conditions on gene expression
plotPCA(rld)


# Plot MA plot to visualize the differences between conditions
# MA plot shows log ratios (M) against average expression values (A), highlighting significant changes
plotMA(res)



