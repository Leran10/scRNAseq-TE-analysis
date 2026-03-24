# scRNA-seq Transposable Element (TE) Analysis

This repository contains the analysis scripts, data, and outputs for Figure 4 and Supplemental Figure 4 of the manuscript, focusing on transposable element expression changes in single-cell RNA-seq data from antibiotic-treated versus control mouse bone marrow samples.

## Repository Structure

```
.
├── scripts/          # R scripts for TE analysis
├── data/             # Input data and significant TE lists
├── outputs/          # Generated figures (PDFs)
├── methods/          # Methods section documentation
└── README.md         # This file
```

## Scripts

### Main Analysis Scripts

1. **TE_analysis.R**
   - Initial TE quantification using TEtranscripts output
   - DESeq2 differential expression analysis
   - Identification of 42 significant TEs (36 downregulated, 6 upregulated)
   - Generates volcano plots

2. **regenerate_figure_4A_grouped.R**
   - Creates heatmap of significant TEs
   - Samples grouped by condition (Control vs Antibiotic-treated)
   - Z-score normalization and visualization
   - **Output:** `Figure_4A_heatmap_significant_TEs.pdf`

3. **create_GREAT_barplots.R**
   - Processes GREAT enrichment analysis results
   - Creates bar plots for GO term enrichment
   - Separate analysis for upregulated and downregulated TEs
   - **Outputs:**
     - `Figure_4H_GREAT_upregulated_TEs_barplot.pdf`
     - `Figure_4H_GREAT_downregulated_TEs_barplot.pdf`

4. **regenerate_figure_4B_all_TEs.R**
   - TE category distribution across treatment groups
   - Categorizes TEs into ERVK, ERV1, ERVL, LINE, SINE, DNA transposons, etc.
   - Grouped bar chart comparing Control vs Antibiotic-treated
   - **Output:** `Figure_4B_TE_category_distribution_all.pdf`

5. **regenerate_figure_4C_all_ERVs.R**
   - ERV class distribution (Class I, II, III)
   - Pie charts for overall, control, and antibiotic-treated samples
   - **Outputs:**
     - `Figure_4C_ERV_class_distribution_overall.pdf`
     - `Figure_4C_ERV_class_distribution_Control.pdf`
     - `Figure_4C_ERV_class_distribution_Abx.pdf`

## Data Files

- **sig_TE.csv** - List of 42 significant TEs with fold changes and annotations
- **TEtranscripts_out_gene_TE_analysis.txt** - Full DESeq2 results for all TEs and genes

## Software Requirements

- **R** (version 4.0+)
- **TEtranscripts** v2.2.3
- **R packages:**
  - DESeq2 v1.48.1
  - tidyverse
  - pheatmap
  - ggplot2
  - RColorBrewer

## Analysis Pipeline

1. **TE Quantification**
   - Run TEtranscripts on STAR-aligned BAM files
   - Input: scRNA-seq BAM files, gene GTF, TE GTF
   - Output: Count table with genes and TEs

2. **Differential Expression**
   - Run `TE_analysis.R` for DESeq2 analysis
   - Identify significant TEs (padj < 0.05)

3. **Figure Generation**
   - `regenerate_figure_4A_grouped.R` → Figure 4A (heatmap)
   - `create_GREAT_barplots.R` → Figure 4H (GREAT enrichment)
   - `regenerate_figure_4B_all_TEs.R` → Supplemental Figure 4A-B (TE categories)
   - `regenerate_figure_4C_all_ERVs.R` → Supplemental Figure 4C-D (ERV classes)

## Key Findings

- **42 significantly differentially expressed TEs** identified
  - 36 downregulated in antibiotic-treated samples
  - 6 upregulated in antibiotic-treated samples
- **ERVK family** most affected (17 significant TEs)
- **ERV families** (ERVK, ERV1, ERVL) account for majority of dysregulated TEs
- GREAT analysis reveals enrichment in immune-related GO terms

## Methods

Detailed methods documentation is available in the `methods/` directory:
- `Methods_TE_Analysis.txt` - Clean version for manuscript
- `Methods_TE_Analysis_with_script_references.txt` - Detailed version with script references
- Word document versions (.docx) also available

## Contact

For questions or issues, please contact:
- Leran Wang (leran.wang@wustl.edu)
