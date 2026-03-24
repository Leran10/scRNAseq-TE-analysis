# Figure 4B: TE Category Distribution - ALL TEs detected in scRNA-seq
# Shows ALL TEs (not just significant ones) separately for Control and Abx groups
# Bar chart with Control and Abx side by side

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE")

# Read count table (contains ALL TEs and genes)
count_table <- read.delim("TEtranscripts_out.cntTable", header = TRUE, row.names = 1)

# Create metadata
metadata <- data.frame(
  Sample = colnames(count_table),
  Type = factor(c("Abx", "Control", "Control", "Abx", "Control", "Abx"),
                levels = c("Control", "Abx"))
)
rownames(metadata) <- colnames(count_table)

cat("\n=== Data loaded ===\n")
cat("Total features (genes + TEs):", nrow(count_table), "\n")

# ============================================================================
# Parse TE information
# ============================================================================

# TEs have format "Name:Family:Class", genes have format "ENSMUSG..."
# Filter for TEs only (features with colons)
te_features <- rownames(count_table)[grepl(":", rownames(count_table))]
te_count_table <- count_table[te_features, ]

cat("Total TEs in count table:", nrow(te_count_table), "\n")

# Parse class and family
parse_TE_info <- function(gene_id) {
  parts <- strsplit(as.character(gene_id), ":")[[1]]
  return(data.frame(
    TE_name = ifelse(length(parts) >= 1, parts[1], "Unknown"),
    TE_family = ifelse(length(parts) >= 2, parts[2], "Unknown"),
    TE_class = ifelse(length(parts) >= 3, parts[3], "Unknown")
  ))
}

te_info <- do.call(rbind, lapply(te_features, parse_TE_info))
te_info$GeneID <- te_features

# ============================================================================
# Categorize TEs
# ============================================================================

categorize_TE <- function(family, class) {
  # ERVs (LTR retrotransposons with specific families)
  if (family == "ERVK") return("ERVK")
  if (family == "ERV1") return("ERV1")
  if (grepl("ERVL", family)) return("ERVL")
  if (family == "ERVB") return("ERVB")

  # Other major TE types
  if (class == "LINE") return("LINE")
  if (class == "SINE") return("SINE")
  if (class == "DNA") return("DNA Transposon")
  if (class == "LTR" && !grepl("ERV", family)) return("LTR (other)")

  # Everything else
  return("Other")
}

te_info$Category <- mapply(categorize_TE, te_info$TE_family, te_info$TE_class)

cat("\n=== TE Categories ===\n")
category_summary <- table(te_info$Category)
print(category_summary)

# ============================================================================
# Normalize counts and calculate group means
# ============================================================================

# Run DESeq2 normalization to get normalized counts
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = metadata,
                              design = ~ Type)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Get normalized counts for TEs
te_normalized <- normalized_counts[te_features, ]

# Calculate mean expression per group
te_data <- data.frame(
  GeneID = te_features,
  Category = te_info$Category,
  Control_mean = rowMeans(te_normalized[, metadata$Type == "Control"]),
  Abx_mean = rowMeans(te_normalized[, metadata$Type == "Abx"])
)

# ============================================================================
# Detect TEs in each group (expression > threshold)
# ============================================================================

threshold <- 10  # Minimum normalized count to consider "detected"

te_data$Detected_Control <- te_data$Control_mean > threshold
te_data$Detected_Abx <- te_data$Abx_mean > threshold

cat("\n=== Detection Summary ===\n")
cat("TEs detected in Control (mean count >", threshold, "):", sum(te_data$Detected_Control), "\n")
cat("TEs detected in Abx (mean count >", threshold, "):", sum(te_data$Detected_Abx), "\n")

# ============================================================================
# Count TEs by category for each group
# ============================================================================

control_counts <- te_data %>%
  filter(Detected_Control) %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Group = "Control")

abx_counts <- te_data %>%
  filter(Detected_Abx) %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Group = "Abx")

# Combine
group_counts <- rbind(control_counts, abx_counts)

# Reorder categories for better visualization
category_order <- c("ERVK", "ERV1", "ERVL", "ERVB", "LTR (other)",
                   "LINE", "SINE", "DNA Transposon", "Other")
group_counts$Category <- factor(group_counts$Category,
                               levels = category_order)

# Make sure Group is a factor with correct order
group_counts$Group <- factor(group_counts$Group, levels = c("Control", "Abx"))

cat("\n=== TE Category Distribution by Group ===\n")
print(group_counts)

# ============================================================================
# Create grouped bar chart (Control and Abx side by side)
# ============================================================================

# Define colors for categories
category_colors <- c(
  "ERVK" = "#e41a1c",      # Red
  "ERV1" = "#377eb8",      # Blue
  "ERVL" = "#4daf4a",      # Green
  "ERVB" = "#984ea3",      # Purple
  "LTR (other)" = "#ff7f00", # Orange
  "LINE" = "#a65628",      # Brown
  "SINE" = "#f781bf",      # Pink
  "DNA Transposon" = "#999999", # Gray
  "Other" = "#666666"      # Dark gray
)

pdf("Figure_4B_TE_category_distribution_all.pdf", width = 10, height = 7)
ggplot(group_counts, aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Control" = "#66c2a5", "Abx" = "#fc8d62"),
                    labels = c("Control", "Antibiotic-treated")) +
  labs(title = "TE Category Distribution Across Sample Groups",
       subtitle = paste0("All TEs detected in scRNA-seq (mean normalized count > ", threshold, ")"),
       x = "TE Category",
       y = "Number of TEs Detected",
       fill = "Group") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top"
  )
dev.off()

# ============================================================================
# Save data
# ============================================================================

write.csv(group_counts, "Figure_4B_category_counts_by_group.csv", row.names = FALSE)
write.csv(te_data, "TE_classification_with_expression_ALL.csv", row.names = FALSE)

# Summary statistics
cat("\n=== Summary by Category ===\n")
summary_wide <- group_counts %>%
  pivot_wider(names_from = Group, values_from = Count, values_fill = 0)
print(summary_wide)

cat("\n=== Files Generated ===\n")
cat("1. Figure_4B_TE_category_distribution_all.pdf\n")
cat("2. Figure_4B_category_counts_by_group.csv\n")
cat("3. TE_classification_with_expression_ALL.csv\n")
cat("\n=== DONE! ===\n")
