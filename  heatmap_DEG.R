# Load packages

library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(RColorBrewer)


#load geo data

gset<-getGEO("GSE110223",GSEMatrix =TRUE, AnnotGPL=TRUE)

# Extract phenodata 
pheno<-pData(gset[[1]])

if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]



# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "01010101010101010101010101"
sml <- strsplit(gsms, split="")[[1]]


# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Assign samples to groups 

gs <- factor(sml)
groups <- make.names(c("Normal","CRC"))
levels(gs) <- groups
gset$group <- gs

# Set up design matrix

design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B",number = 50000)

colnames(tT)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Platform_SPOTID","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# heatmap

#Genes of interest 

genes_of_interest<-c("PIK3CA","STAT3","MMP9","PARP1", "GSK3B", "SRC","ESR1","EGFR","AKT1","ERBB2")


# Extract expression data
expr_data <- exprs(gset)

# Check if we need to map probe IDs to gene symbols
if(!any(rownames(expr_data) %in% genes_of_interest)) {
  # Map probes to gene symbols (adjust based on  platform)
  gene_symbols <- fData(gset)$Gene.symbol
  rownames(expr_data) <- gene_symbols
}

# Subset expression data for your genes
heatmap_data <- expr_data[rownames(expr_data) %in% genes_of_interest, ]

# Remove duplicates by keeping highest expressed probe
heatmap_data <- heatmap_data[!duplicated(rownames(heatmap_data)), ]

# Check which genes were found
missing_genes <- setdiff(genes_of_interest, rownames(heatmap_data))
if(length(missing_genes) > 0) {
  message("These genes were not found: ", paste(missing_genes, collapse=", "))
}

# Scale data by row (z-score normalization)
scaled_data <- t(scale(t(heatmap_data)))

# Prepare sample annotations
annotation_df <- data.frame(
  Tissue_Type = gset$group
)
rownames(annotation_df) <- colnames(scaled_data)

# Define colors for annotations
ann_colors <- list(
  Tissue_Type = c(Normal = "blue", CRC = "red")
)

# Create the heatmap
pheatmap(
  scaled_data,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,  # Often cleaner with many samples
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  main = "Expression of 10 Key Genes in Colorectal Cancer",
  border_color = NA,
  gaps_col = cumsum(table(gset$group)))# Adds space between groups



