
BiocManager::install(c("maftools","SummarizedExperiment","BSgenome", "Biostrings"),force = TRUE)


#Load and install packages

library(TCGAbiolinks)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(DESeq2)
library(maftools)
library(SummarizedExperiment)
library(ggplot2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Get a list of projects 

GDC_projects<-getGDCprojects()
getProjectSummary("TCGA-COAD")


# Building a query

query_tcga<-GDCquery(project ="TCGA-COAD",
         data.category = "Transcriptome Profiling")

output_querY_tcga<-getResults(query_tcga)

#Build a query to retrieve the gene expression data 

unlink("GDCdata", recursive = TRUE)

query_tcga<-GDCquery(project ="TCGA-COAD",
                     data.category = "Transcriptome Profiling",
                     experimental.strategy = "RNA-Seq",
                     workflow.type ="STAR - Counts",
                     access = "open",
                     barcode = c("TCGA-AZ-6601-11A-01R-1774-07","TCGA-AA-3511-11A-01R-1839-07",
                                 "TCGA-A6-2680-11A-01R-A32Z-07","TCGA-AA-3712-11A-01R-1723-07",
                                 "TCGA-A6-2686-11A-01R-A32Z-07","TCGA-AZ-6605-11A-01R-1839-07",
                                 "TCGA-AA-3516-11A-01R-A32Z-07","TCGA-AA-3525-11A-01R-A32Z-07",
                                 "TCGA-A6-2682-11A-01R-A32Z-07","TCGA-AA-3662-11A-01R-1723-07",
                                 "TCGA-AA-3688-01A-01R-0905-07","TCGA-G4-6298-01A-11R-1723-07",
                                 "TCGA-AA-3672-01A-01R-0905-07","TCGA-G4-6314-01A-11R-1723-07",
                                 "TCGA-A6-2682-01A-01R-1410-07","TCGA-AA-3562-01A-02R-0821-07",
                                 "TCGA-AA-3979-01A-01R-1022-07","TCGA-AA-3524-01A-02R-0821-07",
                                 "TCGA-DM-A285-01A-11R-A16W-07","TCGA-A6-4105-01A-02R-1774-07",
                                 "TCGA-AA-3976-01A-01R-1022-07","TCGA-AD-6964-01A-11R-1928-07"))

output_query<-getResults(query_tcga)



#Download data 

GDCdownload(query_tcga,directory = "GDCdata")

data<-GDCprepare(query_tcga,directory = "GDCdata")

# Extract count matrix and clinical data

count_data <- assay(data)
clinical_data <- colData(data)
clinical_data <- as.data.frame(colData(data))
colnames(clinical_data)

# Add gene annotation and gene symbols
# First, get the rowData which contains gene information
gene_info <- rowData(data)

# Convert to data frame
gene_info_df <- as.data.frame(gene_info)
# If ENSEMBL IDs are in the format ENSG00000133703.13, we can remove the version number
gene_info_df$ensembl_id <- sub("\\..*", "", gene_info_df$gene_id)

# First get all possible ENSEMBL IDs from the database
all_ensembl <- keys(org.Hs.eg.db, keytype="ENSEMBL")

# See how many of your IDs match
sum(gene_info_df$ensembl_id %in% all_ensembl)

# For any remaining NAs, you could try:
gene_info_df$symbol <- ifelse(is.na(gene_info_df$symbol),
                              gene_info_df$gene_id,  # fallback to original ID
                              gene_info_df$symbol)

# Ensure count data has gene_id column
count_df <- as.data.frame(count_data)
count_df$gene_id <- rownames(count_df)

# Merge
merged_data <- merge(gene_info_df, count_df, by = "gene_id")



# First, create a design matrix based on sample type (tumor vs normal)
# Extract sample type from barcode (01A = tumor, 11A = normal)
clinical_data$sample_type <- ifelse(substr(clinical_data$barcode, 14, 15) == "01", 
                                    "tumor", "normal")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = clinical_data,
                              design = ~ sample_type)

# Pre-filtering: remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("sample_type", "tumor", "normal"))
res <- as.data.frame(res)

# Add gene symbols to results
res$symbol <- gene_info_df$symbol[match(rownames(res), gene_info_df$gene_id)]

# Define my genes of interest

genes_of_interest<-c("PIK3CA","STAT3","MMP9","PARP1","GSK3B","SRC","ESR1","EGFR","AKT1","ERBB2")


volcano_plot <- EnhancedVolcano(res,
                                lab = res$symbol,
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                title = 'Tumor vs Normal in TCGA-COAD',
                                pCutoff = 0.05,
                                FCcutoff = 1,
                                pointSize = 3.0,
                                labSize = 4.0,
                                colAlpha = 0.7,
                                legendPosition = 'right',
                                drawConnectors = TRUE,
                                widthConnectors = 0.5,
                                selectLab = genes_of_interest,
                                xlim = c(-10, 10),
                                ylim = c(0, max(-log10(res$pvalue), na.rm = TRUE) + 2),
                                caption = paste('Total genes =', nrow(res),
                                                '| Significant genes =', 
                                                sum(res$pvalue < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)))

print(volcano_plot)
