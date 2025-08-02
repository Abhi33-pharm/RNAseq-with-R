# Load Packages 

library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(EnhancedVolcano)
library(tidyverse)
library(dplyr)

#Load geo series 

gse_crc <- getGEO("GSE35896", GSEMatrix = TRUE)[[1]]  # CRC dataset
gse_ctrl <- getGEO("GSE9254", GSEMatrix = TRUE)[[1]] # Normal dataset

# Checking platform compatibility
platform_crc <- annotation(gse_crc)
platform_ctrl <- annotation(gse_ctrl)

if (platform_crc != platform_ctrl) {
  stop("Datasets are from different platforms and not directly comparable.")
}

# Combining the expression sets
expr_crc <- exprs(gse_crc)
expr_ctrl <- exprs(gse_ctrl)

# Combining and labelling
combined_expr <- cbind(expr_crc, expr_ctrl)
group <- factor(c(rep("CRC", ncol(expr_crc)), rep("Control", ncol(expr_ctrl))))



# Get platform ID from your dataset
platform_id <- annotation(gse_crc)  #  "GPL570"

# Downloading the platform data
gpl <- getGEO(platform_id, AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Inspect the column names to find gene symbol column
head(colnames(gpl_table))  # Look for something like "Gene Symbol" or "GENE_SYMBOL"



# Differential Expression 

# Design matrix
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "CRCvsControl")

# Fit model
fit <- lmFit(combined_expr, design)
fit <- eBayes(fit)

# Get top table
deg_results <- topTable(fit, coef = "CRCvsControl", number = Inf)

#Merging DEG Results with Gene Symbols

# Giving rownames of deg_results are probe IDs

deg_results$probe_id <- rownames(deg_results)

# Adjusting names  (e.g., "Gene Symbol" to "GENE_SYMBOL")
symbol_map <- gpl_table[, c("ID", "Gene symbol")]
colnames(symbol_map) <- c("probe_id", "Gene_symbol")

# Merge with DEG table
deg_annotated <- merge(deg_results, symbol_map, by = "probe_id", all.x = TRUE)

# Cleaning up multiple symbols (if present as "TP53 /// XYZ")
deg_annotated$Gene_symbol <- sapply(strsplit(deg_annotated$Gene_symbol, " /// "), `[`, 1)

# Mark genes of interest
genes_of_interest <- c(" PIK3CA", "STAT3", "MMP9", "PARP1", "CTNNB1", "SRC", "ESR1", "EGFR", "AKT1","ERBB2")
deg_annotated$Highlight <- ifelse(deg_annotated$Gene_symbol %in% genes_of_interest, "TargetGene",
                                  ifelse(deg_annotated$adj.P.Val < 0.05 & abs(deg_annotated$logFC) > 1, "Significant", "Not Significant"))


# Filter duplicates - keep most significant (lowest adj.P.Val)
deg_annotated_unique <- deg_annotated |> 
  arrange(adj.P.Val) |>   # Sort by significance
  distinct(Gene_symbol, .keep_all = TRUE)  # Keep first occurrence per gene

# Createing Volcano Plot using EnhancedVolcano
volcano_plot <- EnhancedVolcano(
  deg_annotated_unique,
  lab = deg_annotated_unique$Gene_symbol,
  x = 'logFC',
  y = 'adj.P.Val',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 3.0,
  title = 'CRC vs Control - Differential Expression',
  subtitle = 'Volcano plot of differentially expressed genes',
  caption = 'FC cutoff = 1; p-value cutoff = 0.05',
  colAlpha = 0.7,
  legendLabels = c('Not sig.', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
  legendPosition = 'right',
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  selectLab = genes_of_interest,
  boxedLabels = FALSE,
  max.overlaps = Inf
)

# Save the plot
ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

# Display the plot
print(volcano_plot)


# Box Plots for Genes of Interest

# Preparing data for box plots
expr_data <- as.data.frame(t(combined_expr))
expr_data$Group <- group

# Getting probe IDs for genes of interest
gene_probes <- deg_annotated[deg_annotated$Gene_symbol %in% genes_of_interest, c("probe_id", "Gene_symbol")]

# Creating a list to store plots
plot_list <- list()

# Custom color palette
group_colors <- c("CRC" = "#E64B35", "Control" = "#4DBBD5")

# Create box plots for each gene
for(i in 1:nrow(gene_probes)) {
  probe <- gene_probes$probe_id[i]
  gene <- gene_probes$Gene_symbol[i]
  
  if(!probe %in% colnames(expr_data)) next
  
  # Calculate y-axis limits
  y_values <- expr_data[[probe]]
  y_range <- range(y_values)
  y_limit <- c(y_range[1] - 0.1*diff(y_range), y_range[2] + 0.2*diff(y_range))
  
  p <- ggplot(expr_data, aes(x = Group, y = .data[[probe]], fill = Group)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = 21,
      outlier.size = 2,
      outlier.fill = "white",
      outlier.color = "blue"
    ) +
    geom_jitter(
      width = 0.15,
      size = 1,
      alpha = 0.6,
      color = "black"
    ) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = gene,
      y = "Normalized Expression",
      x = ""
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(size = 11),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      aspect.ratio = 1.2
    ) +
    coord_cartesian(ylim = y_limit) +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      label.x = 1.4,
      label.y = y_limit[2] * 0.95,
      size = 4,
      bracket.size = 0.5
    )
  
  # Add to plot list
  plot_list[[gene]] <- p
  
  # Save individual plots
  ggsave(
    paste0("boxplot_", gene, ".png"),
    plot = p,
    width = 3.5,
    height = 4,
    dpi = 600,
    bg = "white"
  )
}

# Combine plots into a multi-panel figure (if desired)
combined_plots <- ggarrange(plotlist = plot_list, ncol = 3, nrow = ceiling(length(plot_list)/3))
ggsave(
  "combined_boxplots.png",
  plot = combined_plots,
  width = 10,
  height = 3 * ceiling(length(plot_list)/3),
  dpi = 600,
  bg = "white"
)


#Heatmap
# Get unique genes of interest removing double symbols
deg_annotated_unique <- deg_annotated %>% 
  filter(Gene_symbol %in% genes_of_interest) %>% 
  group_by(Gene_symbol) |> 
  arrange(adj.P.Val) |>  # Keeping most significant probe
  dplyr::slice(1) |>     # Explicitly use dplyr's slice
  ungroup() |> 
  as.data.frame()        # Convert to regular data frame

# Extracting expression data
heatmap_data <- combined_expr[deg_annotated_unique$probe_id, ]
row_annot <- deg_annotated_unique[, c("Gene_symbol", "logFC", "adj.P.Val")]
rownames(heatmap_data) <- row_annot$Gene_symbol

# Creating sample annotation
sample_annot <- data.frame(Group = group)
rownames(sample_annot) <- colnames(heatmap_data)

# Creating gene annotation
gene_annot <- data.frame(
  LogFC = row_annot$logFC,
  FDR = -log10(row_annot$adj.P.Val)
)
rownames(gene_annot) <- row_annot$Gene_symbol

# Set colors
group_colors <- list(
  Group = c("CRC" = "#E64B35", "Control" = "#4DBBD5")
)

# Create heatmap
heatmap_goi <- pheatmap(
  heatmap_data,
  scale = "row",  # Z-score normalization by row
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = sample_annot,
  annotation_row = gene_annot,
  annotation_colors = group_colors,
  fontsize_row = 9,
  fontsize_col = 8,
  cellwidth = NA,
  cellheight = 12,
  main = "Expression Heatmap of Target Genes (Z-scores)",
  gaps_col = cumsum(table(group)),  # Adding gap between groups
  treeheight_row = 20,
  treeheight_col = 20,
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_row = TRUE,
  annotation_names_col = TRUE
)

# Save heatmap
ggsave(
  "heatmap_genes_of_interest.png",
  plot = heatmap_goi,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)


# Correlation Analysis Heatmap of Genes of Interest

# Extraction of expression data for genes of interest
corr_genes <- deg_annotated |> 
  filter(Gene_symbol %in% genes_of_interest) |> 
  group_by(Gene_symbol) |> 
  arrange(adj.P.Val) |>   # Keep most significant probe per gene
  dplyr::slice(1) |>  
  ungroup() |> 
  as.data.frame()

corr_data <- combined_expr[corr_genes$probe_id, ]
rownames(corr_data) <- corr_genes$Gene_symbol

# Calculation of correlation matrix
corr_matrix <- cor(t(corr_data), method = "pearson")  # Pearson correlation between genes

# Creating annotation for genes (using logFC and FDR)
gene_annot_corr <- data.frame(
  LogFC = corr_genes$logFC,
  FDR = -log10(corr_genes$adj.P.Val)
)
rownames(gene_annot_corr) <- corr_genes$Gene_symbol

# Creating correlation heatmap
corr_heatmap <- pheatmap(
  corr_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = gene_annot_corr,
  annotation_col = gene_annot_corr,
  annotation_colors = list(
    LogFC = colorRampPalette(c("blue", "white", "red"))(100),
    FDR = colorRampPalette(c("white", "darkgreen"))(100)
  ),
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 8,
  main = "Correlation Analysis of Target Genes",
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 30,
  cellheight = 30
)

# Saving correlation heatmap
ggsave(
  "correlation_heatmap.png",
  plot = corr_heatmap,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)


