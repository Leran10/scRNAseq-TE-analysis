# Figure 4C: ERV Class Distribution - ALL ERVs detected in scRNA-seq
# Shows ALL ERVs (not just significant ones) separately for Control and Abx groups

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
# Classify ERVs into Class I/II/III
# ============================================================================

classify_ERV <- function(family, name, class) {
  # Only classify if it's an LTR retrotransposon
  if (class != "LTR") {
    return("Not_ERV")
  }

  # Class I: Gammaretrovirus-like (ERV1, MuLV-related)
  if (grepl("ERV1", family)) {
    return("Class_I_Gamma")
  }

  # Class II: Betaretrovirus-like (ERVK, MMTV, IAP-related)
  if (grepl("ERVK", family)) {
    return("Class_II_Beta")
  }

  # Class III: Spumaretrovirus-like (ERVL, MaLR)
  if (grepl("ERVL|MaLR", family)) {
    return("Class_III_Spuma")
  }

  # Gypsy is not an ERV but a different LTR retrotransposon
  if (grepl("Gypsy", family)) {
    return("Not_ERV")
  }

  # Other LTR elements that are not clearly classified ERVs
  return("LTR_Other")
}

# Apply classification
te_info$ERV_Class <- mapply(classify_ERV, te_info$TE_family, te_info$TE_name, te_info$TE_class)

# Filter to keep only ERVs (exclude Not_ERV and LTR_Other)
erv_only <- te_info[te_info$ERV_Class %in% c("Class_I_Gamma", "Class_II_Beta", "Class_III_Spuma"), ]

cat("\n=== ERV Classification ===\n")
cat("Total ERVs (Class I/II/III) detected:", nrow(erv_only), "\n")
cat("  Class I (Gamma):", sum(erv_only$ERV_Class == "Class_I_Gamma"), "\n")
cat("  Class II (Beta):", sum(erv_only$ERV_Class == "Class_II_Beta"), "\n")
cat("  Class III (Spuma):", sum(erv_only$ERV_Class == "Class_III_Spuma"), "\n")

# ============================================================================
# Normalize counts and calculate group means
# ============================================================================

# Subset count table to ERVs only
erv_counts <- te_count_table[erv_only$GeneID, ]

# Run DESeq2 normalization to get normalized counts
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = metadata,
                              design = ~ Type)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Get normalized counts for ERVs
erv_normalized <- normalized_counts[erv_only$GeneID, ]

# Calculate mean expression per group
erv_data <- data.frame(
  GeneID = erv_only$GeneID,
  ERV_Class = erv_only$ERV_Class,
  Control_mean = rowMeans(erv_normalized[, metadata$Type == "Control"]),
  Abx_mean = rowMeans(erv_normalized[, metadata$Type == "Abx"])
)

# ============================================================================
# Detect ERVs in each group (expression > threshold)
# ============================================================================

threshold <- 10  # Minimum normalized count to consider "detected"

erv_data$Detected_Control <- erv_data$Control_mean > threshold
erv_data$Detected_Abx <- erv_data$Abx_mean > threshold

cat("\n=== Detection Summary ===\n")
cat("ERVs detected in Control (mean count >", threshold, "):", sum(erv_data$Detected_Control), "\n")
cat("ERVs detected in Abx (mean count >", threshold, "):", sum(erv_data$Detected_Abx), "\n")

# ============================================================================
# Count ERVs by class for each group
# ============================================================================

control_counts <- erv_data %>%
  filter(Detected_Control) %>%
  group_by(ERV_Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Group = "Control")

abx_counts <- erv_data %>%
  filter(Detected_Abx) %>%
  group_by(ERV_Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Group = "Abx")

# Combine
group_counts <- rbind(control_counts, abx_counts)

# Clean up class names
group_counts$ERV_Class <- factor(group_counts$ERV_Class,
                                levels = c("Class_I_Gamma", "Class_II_Beta", "Class_III_Spuma"),
                                labels = c("Class I (Gamma-like)",
                                          "Class II (Beta-like)",
                                          "Class III (Spuma-like)"))

# Calculate percentages within each group
group_counts <- group_counts %>%
  group_by(Group) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1)) %>%
  ungroup()

cat("\n=== ERV Class Distribution by Group ===\n")
print(group_counts)

# ============================================================================
# Create pie charts
# ============================================================================

# Overall pie chart (all ERVs)
overall_counts <- erv_data %>%
  group_by(ERV_Class) %>%
  summarise(Count = n(), .groups = "drop")

overall_counts$ERV_Class <- factor(overall_counts$ERV_Class,
                                   levels = c("Class_I_Gamma", "Class_II_Beta", "Class_III_Spuma"),
                                   labels = c("Class I (Gamma-like)",
                                             "Class II (Beta-like)",
                                             "Class III (Spuma-like)"))

overall_counts$Percentage <- round(overall_counts$Count / sum(overall_counts$Count) * 100, 1)

pdf("Figure_4C_ERV_class_distribution_overall.pdf", width = 8, height = 8)
ggplot(overall_counts, aes(x = "", y = Count, fill = ERV_Class)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Class I (Gamma-like)" = "#66c2a5",
                               "Class II (Beta-like)" = "#fc8d62",
                               "Class III (Spuma-like)" = "#8da0cb")) +
  labs(title = "ERV Class Distribution - Overall",
       subtitle = "All ERVs detected in scRNA-seq (regardless of significance)",
       fill = "ERV Class") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom")
dev.off()

# Control pie chart
control_data <- group_counts[group_counts$Group == "Control", ]

pdf("Figure_4C_ERV_class_distribution_Control.pdf", width = 8, height = 8)
ggplot(control_data, aes(x = "", y = Count, fill = ERV_Class)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Class I (Gamma-like)" = "#66c2a5",
                               "Class II (Beta-like)" = "#fc8d62",
                               "Class III (Spuma-like)" = "#8da0cb")) +
  labs(title = "ERV Class Distribution - Control",
       subtitle = paste0("ERVs with mean normalized count > ", threshold),
       fill = "ERV Class") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom")
dev.off()

# Abx pie chart
abx_data <- group_counts[group_counts$Group == "Abx", ]

pdf("Figure_4C_ERV_class_distribution_Abx.pdf", width = 8, height = 8)
ggplot(abx_data, aes(x = "", y = Count, fill = ERV_Class)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Class I (Gamma-like)" = "#66c2a5",
                               "Class II (Beta-like)" = "#fc8d62",
                               "Class III (Spuma-like)" = "#8da0cb")) +
  labs(title = "ERV Class Distribution - Antibiotic-treated",
       subtitle = paste0("ERVs with mean normalized count > ", threshold),
       fill = "ERV Class") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom")
dev.off()

# ============================================================================
# Save data
# ============================================================================

write.csv(overall_counts, "Figure_4C_ERV_class_counts_overall.csv", row.names = FALSE)
write.csv(group_counts, "Figure_4C_ERV_class_counts_by_group.csv", row.names = FALSE)
write.csv(erv_data, "ERV_classification_with_expression_ALL.csv", row.names = FALSE)

cat("\n=== Files Generated ===\n")
cat("1. Figure_4C_ERV_class_distribution_overall.pdf\n")
cat("2. Figure_4C_ERV_class_distribution_Control.pdf\n")
cat("3. Figure_4C_ERV_class_distribution_Abx.pdf\n")
cat("4. Figure_4C_ERV_class_counts_overall.csv\n")
cat("5. Figure_4C_ERV_class_counts_by_group.csv\n")
cat("6. ERV_classification_with_expression_ALL.csv\n")
cat("\n=== DONE! ===\n")
