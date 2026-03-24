# Regenerate Figure 4A with samples grouped by condition
# Control samples together, Abx samples together (no clustering)

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE")

# Load data
sig_TE <- read.csv("sig_TE.csv", row.names = 1)
count_table <- read.delim("TEtranscripts_out.cntTable", header = TRUE, row.names = 1)

# Create metadata
metadata <- data.frame(
  Sample = c("Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5", "Sample_6"),
  Type = factor(c("Abx", "Control", "Control", "Abx", "Control", "Abx"),
                levels = c("Control", "Abx"))
)
rownames(metadata) <- colnames(count_table)

cat("\n=== Metadata ===\n")
print(metadata)

# ============================================================================
# FIGURE 4A: Heatmap of significant TEs ordered by fold change
# ============================================================================

# Order by fold change (lowest to highest for corrected data)
sig_TE_ordered <- sig_TE[order(sig_TE$log2FoldChange, decreasing = FALSE), ]

# Get normalized counts for these TEs
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = metadata,
                              design = ~ Type)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Extract normalized counts for significant TEs
sig_TE_counts <- normalized_counts[sig_TE_ordered$GeneID, ]

# Use simple gene names for better visualization
rownames(sig_TE_counts) <- sig_TE_ordered$simple_geneName

# Reorder columns to group by condition: Control, Control, Control, Abx, Abx, Abx
control_samples <- rownames(metadata)[metadata$Type == "Control"]
abx_samples <- rownames(metadata)[metadata$Type == "Abx"]
sample_order <- c(control_samples, abx_samples)

sig_TE_counts_ordered <- sig_TE_counts[, sample_order]

cat("\n=== Sample order in heatmap ===\n")
cat("Control samples:", control_samples, "\n")
cat("Abx samples:", abx_samples, "\n")

# Log2 transform (add pseudocount to avoid log(0))
log2_counts <- log2(sig_TE_counts_ordered + 1)

# Z-score normalization for better heatmap visualization
z_score_counts <- t(scale(t(log2_counts)))

# Create annotation for samples (in the new order)
annotation_col <- data.frame(
  Group = metadata[sample_order, "Type"],
  row.names = sample_order
)

# Color palette
ann_colors <- list(
  Group = c("Control" = "#4575b4", "Abx" = "#d73027")
)

# Generate heatmap with NO column clustering
pdf("Figure_4A_heatmap_significant_TEs.pdf", width = 8, height = 14)
pheatmap(z_score_counts,
         cluster_rows = FALSE,  # Keep order by fold change
         cluster_cols = FALSE,  # NO clustering - keep Control|Abx grouping
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Significant TEs (ordered by log2FC)\nControl vs Abx")
dev.off()

cat("\n=== Figure 4A regenerated! ===\n")
cat("Samples are now grouped: Control | Control | Control | Abx | Abx | Abx\n")
cat("TEs are ordered from most downregulated (-1.79) to most upregulated (+1.01)\n")
