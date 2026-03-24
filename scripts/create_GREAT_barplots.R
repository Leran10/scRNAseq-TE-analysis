# Create GREAT Analysis Bar Plots
# Shows top enriched GO terms (BP, MF, CC) for upregulated and downregulated TEs

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/leranwang/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE")

cat("\n=== Creating GREAT Analysis Bar Plots ===\n\n")

# Function to process GREAT results and get top terms
process_great_results <- function(bp_file, mf_file, cc_file, top_n = 10) {

  # Read the three GO category files
  bp <- read.csv(bp_file, stringsAsFactors = FALSE)
  mf <- read.csv(mf_file, stringsAsFactors = FALSE)
  cc <- read.csv(cc_file, stringsAsFactors = FALSE)

  # Add category label
  bp$Category <- "Biological Process"
  mf$Category <- "Molecular Function"
  cc$Category <- "Cellular Component"

  # Combine all
  all_go <- rbind(bp, mf, cc)

  # Calculate -log10(p-value) for visualization
  all_go$NegLog10P <- -log10(all_go$Hyper_Adjp_BH + 1e-300)  # Add small value to avoid log(0)

  # Get top N from each category by adjusted p-value
  top_terms <- all_go %>%
    group_by(Category) %>%
    arrange(Hyper_Adjp_BH) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    arrange(Category, Hyper_Adjp_BH)

  # Shorten term names if too long
  top_terms$name_short <- ifelse(nchar(top_terms$name) > 50,
                                  paste0(substr(top_terms$name, 1, 47), "..."),
                                  top_terms$name)

  return(top_terms)
}

# Process upregulated TEs
cat("Processing upregulated TEs results...\n")
up_results <- process_great_results(
  "GREAT_upregulated_TEs_GO_GO Biological Process.csv",
  "GREAT_upregulated_TEs_GO_GO Molecular Function.csv",
  "GREAT_upregulated_TEs_GO_GO Cellular Component.csv",
  top_n = 10
)

# Process downregulated TEs
cat("Processing downregulated TEs results...\n")
down_results <- process_great_results(
  "GREAT_downregulated_TEs_GO_GO Biological Process.csv",
  "GREAT_downregulated_TEs_GO_GO Molecular Function.csv",
  "GREAT_downregulated_TEs_GO_GO Cellular Component.csv",
  top_n = 10
)

# Create function to make bar plot
create_great_barplot <- function(data, title) {

  # Reorder factor levels for plotting
  data$Category <- factor(data$Category,
                         levels = c("Biological Process",
                                   "Molecular Function",
                                   "Cellular Component"))

  # Create unique identifier for each term for ordering
  data$term_id <- paste(data$Category, data$name, sep = "_")

  # Order by category and p-value within category
  data <- data %>%
    arrange(Category, Hyper_Adjp_BH) %>%
    mutate(term_id = factor(term_id, levels = unique(term_id)))

  # Create plot
  p <- ggplot(data, aes(x = NegLog10P, y = reorder(name_short, NegLog10P), fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Biological Process" = "#E41A1C",
      "Molecular Function" = "#377EB8",
      "Cellular Component" = "#4DAF4A"
    )) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = title,
      x = "-log10(Adjusted P-value)",
      y = "",
      fill = "GO Category"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "lightgray"),
      legend.position = "none",
      panel.grid.major.y = element_blank()
    )

  return(p)
}

# Create plots
cat("\nCreating bar plots...\n")

# Upregulated TEs plot
p_up <- create_great_barplot(
  up_results,
  "GREAT Analysis: Upregulated TEs in Antibiotic-Treated Samples"
)

ggsave("Figure_4H_GREAT_upregulated_TEs_barplot.pdf",
       plot = p_up,
       width = 10,
       height = 12)

# Downregulated TEs plot
p_down <- create_great_barplot(
  down_results,
  "GREAT Analysis: Downregulated TEs in Antibiotic-Treated Samples"
)

ggsave("Figure_4H_GREAT_downregulated_TEs_barplot.pdf",
       plot = p_down,
       width = 10,
       height = 12)

# Save summary tables
write.csv(up_results, "Figure_4H_GREAT_upregulated_summary.csv", row.names = FALSE)
write.csv(down_results, "Figure_4H_GREAT_downregulated_summary.csv", row.names = FALSE)

cat("\n=== Summary ===\n")
cat("\nUpregulated TEs - Top enriched terms by category:\n")
cat("\nBiological Process:\n")
print(up_results %>% filter(Category == "Biological Process") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\nMolecular Function:\n")
print(up_results %>% filter(Category == "Molecular Function") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\nCellular Component:\n")
print(up_results %>% filter(Category == "Cellular Component") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\n\nDownregulated TEs - Top enriched terms by category:\n")
cat("\nBiological Process:\n")
print(down_results %>% filter(Category == "Biological Process") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\nMolecular Function:\n")
print(down_results %>% filter(Category == "Molecular Function") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\nCellular Component:\n")
print(down_results %>% filter(Category == "Cellular Component") %>%
        select(name, Hyper_Fold_Enrichment, Hyper_Adjp_BH) %>% head(5))

cat("\n=== Files Generated ===\n")
cat("1. Figure_4H_GREAT_upregulated_TEs_barplot.pdf\n")
cat("2. Figure_4H_GREAT_downregulated_TEs_barplot.pdf\n")
cat("3. Figure_4H_GREAT_upregulated_summary.csv\n")
cat("4. Figure_4H_GREAT_downregulated_summary.csv\n")

cat("\n=== DONE! ===\n")
