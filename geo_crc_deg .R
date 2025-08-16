# Install Packages 
BiocManager::install("hugene10sttranscriptcluster.db", force = TRUE)

# Load libraries
library(GEOquery)    
library(limma)       
library(ggplot2)     
library(pheatmap)    
library(dplyr)       
library(tidyr)
library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)

#Load data

gse <- getGEO("GSE44076", GSEMatrix = TRUE, getGPL = TRUE)
expr_data <- exprs(gse[[1]])  # Expression matrix
pdata <- pData(gse[[1]])      # Phenotype data

head(pdata)


# Check the platform used (e.g., GPL6244 for Hugene 1.0 ST array)
gpl <- annotation(gse[[1]])
print(gpl)  # Verify the platform ID

gpl <- getGEO("GPL13667")
mapping <- Table(gpl)[, c("ID", "Gene Symbol")]

# Map probe IDs to gene symbols

head(rownames(expr_data))
 


# Check sample groups

table(pdata$characteristics_ch1)

# Subsetting samples of interest

healthy_samples<-which(pdata$characteristics_ch1=="sample type: Mucosa")
tumor_samples<-which(pdata$characteristics_ch1=="sample type: Tumor")

# Creating a subset expression matrix
expr_subset <- expr_data[, c(healthy_samples, tumor_samples)]
group_labels <- factor(c(rep("Healthy", length(healthy_samples)), 
                         rep("Tumor", length(tumor_samples))),
                       levels = c("Healthy", "Tumor"))


#Differential expression analysis

# Design matrix
design <- model.matrix(~ 0 + group_labels)
colnames(design) <- levels(group_labels)

# Fit linear model
fit <- lmFit(expr_subset, design)
contrasts <- makeContrasts(Tumor_vs_Healthy = Tumor - Healthy, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


# Extract results
de_genes <- topTable(fit2, number = Inf, adjust.method = "BH", p.value = 0.05)
head(de_genes)
de_genes$GeneSymbol <- mapping[match(rownames(de_genes), mapping$ID), "Gene Symbol"]


