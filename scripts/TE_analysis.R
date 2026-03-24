if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("data.table")
install.packages("tidyverse")
.rs.restartR() 



library("data.table")
library("tidyverse")
library("DESeq2")

TE_id <- read.table("~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_gene_id_check.csv",stringsAsFactors = TRUE,sep = ";") %>%
         separate(into = c("gene_id","transcript_id","family_id","class_id"),col = V1,sep = ";")




data <- read.table("~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TEtranscripts_out.cntTable",header=T,row.names=1)


rownames(data)[!rownames(data) %like% "ENSMUSG"]

unique(rownames(data)[!rownames(data) %like% "ENSMUSG"])
#78316     6


# sample1: HMFJGDSXY_AATGCCATGA-TACGTAATGC_L003_R1.fastq.gz    VNAM
# sample2: HMFJGDSXY_GCACTGAGAA-TATGCGTGAA_L002_R1.fastq.gz    non-VNAM
# sample3: HMFJGDSXY_GCCCGATGGA-AATCGTCTAG_L002_R1.fastq.gz    non-VNAM
# sample4: HMFJGDSXY_GCCTTCGGTA-CCAACGATTT_L003_R1.fastq.gz    VNAM
# sample5: HMFJGDSXY_TATTGAGGCA-CAGGTAAGTG_L002_R1.fastq.gz    non-VNAM
# sample6: HMFJGDSXY_TCGGCTCTAC-CCGATGGTCT_L003_R1.fastq.gz    VNAM


metadata <- read.delim("~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/data/Samplemap.txt") %>%
            filter(File.Name %like% "_R1.") %>%
            select(File.Name,Sample.Name,Type) %>%
            mutate(rownames = ifelse(File.Name %like% "AATGCCATGA","sample1_Aligned.out.bam.T",
                                     ifelse(File.Name %like% "GCACTGAGAA","sample2_Aligned.out.bam.C",
                                            ifelse(File.Name %like% "GCCCGATGGA","sample3_Aligned.out.bam.C",
                                                   ifelse(File.Name %like% "GCCTTCGGTA","sample4_Aligned.out.bam.T",
                                                          ifelse(File.Name %like% "TATTGAGGCA","sample5_Aligned.out.bam.C",
                                                                 "sample6_Aligned.out.bam.T")))))) %>%
            column_to_rownames("rownames")



metadata <- metadata[match(colnames(data), rownames(metadata)), ]

metadata$Type <- factor(metadata$Type,levels = c("non-VNAM","VNAM"))

# DESEQ2
ds <- DESeqDataSetFromMatrix(data, metadata, ~Type)
dds <- DESeq(ds)
res <- results(dds) %>%
       data.frame(.) %>%
       rownames_to_column("GeneID") %>%
       mutate(Gene_type = ifelse(GeneID %like% "ENSMUSG","Genes","TEs"))

View(res %>% filter(Gene_type == "TEs"))


resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TEtranscripts_out_gene_TE_analysis.txt", sep="\t",quote=F,row.names=F)

sig.TE <- unique(resSig %>% filter(Gene_type == "TEs") %>% pull(GeneID))


dim(resSig %>% filter(Gene_type == "TEs") %>% filter(log2FoldChange > 0)) # 6
resSig %>% filter(Gene_type == "TEs") %>% filter(log2FoldChange > 0) %>% pull(GeneID)
# [1] "ERVB4_2-I_MM-int:ERVK:LTR" "ERVB5_1-LTR_MM:ERVK:LTR"   "ERVB7_3-LTR_MM:ERVK:LTR"  
#[4] "MamGypLTR1b:Gypsy:LTR"     "RLTR27:ERVK:LTR"           "RLTR44C:ERVK:LTR"   





dim(resSig %>% filter(Gene_type == "TEs") %>% filter(log2FoldChange < 0)) # 36
resSig %>% filter(Gene_type == "TEs") %>% filter(log2FoldChange < 0) %>% pull(GeneID)
# [1] "Arthur1A:hAT-Tip100:DNA"    "ERVB2_1-I_MM-int:ERVK:LTR"  "IAPLTR3:ERVK:LTR"          
# [4] "Kanga2_a:TcMar-Tc2:DNA"     "L1M4a2:L1:LINE"             "L1ME2z:L1:LINE"            
# [7] "LTR16D:ERVL:LTR"            "LTR16E1:ERVL:LTR"           "LTR52-int:ERVL:LTR"        
# [10] "LTR80A:ERVL:LTR"            "MER117:hAT-Charlie:DNA"     "MER129:LTR:LTR?"           
# [13] "MER4CL34:ERV1:LTR"          "MER53:hAT:DNA"              "MER67B:ERV1:LTR"           
# [16] "MLT1H2:ERVL-MaLR:LTR"       "MLT1J2:ERVL-MaLR:LTR"       "MLTR25C:ERVK:LTR"          
# [19] "MLTR31C_MM:ERVK:LTR"        "MMTV-int:ERVK:LTR"          "MamRep1879:hAT-Tip100:DNA" 
# [22] "MuLV-int:ERV1:LTR"          "ORSL-2a:hAT-Tip100:DNA"     "RLTR17D_Mm:ERVK:LTR"       
# [25] "RLTR1F_Mm:ERVK:LTR"         "RLTR20C2_MM:ERVK:LTR"       "RLTR26C_MM:ERVK:LTR"       
# [28] "RLTR3_Mm:ERVK:LTR"          "RLTR44-int:ERVK:LTR"        "RLTR48B:ERV1:LTR"          
# [31] "RLTR4_MM-int:ERV1:LTR"      "RLTR4_Mm:ERV1:LTR"          "RLTR9C:ERVK:LTR"           
# [34] "Tigger11a:TcMar-Tigger:DNA" "UCON27:UCON27:Unknown"      "X5A_LINE:CR1:LINE"  


resSig.mod <- resSig %>% filter(Gene_type == "TEs") %>% 
       mutate(expression = ifelse(log2FoldChange < 0,"highly_expressed_in_nonVNAM","highly_expressed_in_VNAM")) %>%
       mutate(simple_geneName = sapply(str_split(GeneID,":"),"[",1))

write.csv(resSig.mod,"~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/sig_TE.csv")

View(resSig.mod)
resSig.mod$GeneID

library(plotly)
library(ggrepel)

res <- res %>%
       mutate(annotation = ifelse(GeneID %in% sig.TE,"yes","no")) %>%
       mutate(simple_geneName = sapply(str_split(GeneID,":"),"[",1))



subset(res,annotation == "yes") %>% pull(GeneID)



WT.Naive <- ggplot(data = res,
                                      aes(x = log2FoldChange,
                                          y = -log10(pvalue),
                                          color = "grey",
                                          label = simple_geneName)) +
  geom_point(data = subset(res,res$padj > 0.05 & res$Gene_type == "Genes"),
             color = "grey",show.legend = FALSE) +
  geom_point(data = subset(res,res$padj > 0.05 & res$Gene_type == "TEs"),
             color = "black",size = 1) +
  geom_point(data = subset(res,res$padj < 0.05 & res$Gene_type == "Genes"),
             color = "orange",size = 1) +
  geom_point(data = subset(res,res$padj < 0.05 & res$Gene_type == "TEs"),
             color = "red",size = 1) +
  geom_vline(xintercept = 0,lty = 2)+
  geom_hline(yintercept = -log10(0.05)) +
 # geom_text_repel(data = subset(res,annotation == "yes"),label=subset(res,annotation == "yes") %>% pull(GeneID),show.legend = FALSE,color = "black") +
  labs(title = paste0("VNAM vs non-VNAM")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        legend.position = "none")

ggplotly(WT.Naive)

library(htmlwidgets)

# Save the plot as an HTML file
saveWidget(ggplotly(WT.Naive), file = "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_volcano_plot.html")




# TSNE plot for the significant genes



##### Per manuscript, below code is to refine the plots per Arushna's request #####

# 1.Can you please annotate the following ERVs that I have highlighted in the table below on the volcano plot.
# I have validated these by qPCR that they are indeed lower in VNAM.FYI you also generated this table that I have attached.
# 
# 2.Additionally, can we add a threshold line for significance. I see a dark line I am guessing that might be the one. For clarity, can you make it a dashed line.
# 
# 3.Can you please zoom in to the plot. We can show the y-axis up to the 50 mark and x-axis up to the + 5- and -5-fold change.


res <- res %>%
       mutate(GeneID = ifelse(GeneID == "MER129:LTR:LTR?","MER129:LTR:LTR",GeneID))



WT.Naive.refine <- ggplot(data = res,
                   aes(x = log2FoldChange,
                       y = -log10(pvalue),
                       color = "grey",
                       label = simple_geneName)) +
  geom_point(data = subset(res,res$padj > 0.05 & res$Gene_type == "Genes"),
             color = "grey",show.legend = FALSE) +
  geom_point(data = subset(res,res$padj > 0.05 & res$Gene_type == "TEs"),
             color = "#A6A1A1",size = 1) +
  geom_point(data = subset(res,res$padj < 0.05 & res$Gene_type == "Genes"),
             color = "#A6A1A1",size = 1) +
  geom_point(data = subset(res,res$padj < 0.05 & res$Gene_type == "TEs"),
             color = "red",size = 1) +
  geom_point(
    data = subset(res, GeneID %in% c("IAPLTR3:ERVK:LTR", "L1M4a2:L1:LINE", "L1ME2z:L1:LINE","MER129:LTR:LTR","MER4CL34:ERV1:LTR","MLT1J2:ERVL-MaLR:LTR","MMTV-int:ERVK:LTR","MuLV-int:ERV1:LTR")),
    aes(x = log2FoldChange, y = -log10(pvalue)),   # use your plot's x/y (e.g., y = neglog10p)
    shape = 21, fill = "red", color = "black",  # black border, hollow fill
    stroke = 0.8,                        # border thickness
    size = 1.8,                          # match your main point size
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(res, GeneID %in% c("IAPLTR3:ERVK:LTR", "L1M4a2:L1:LINE", "L1ME2z:L1:LINE","MER129:LTR:LTR","MER4CL34:ERV1:LTR","MLT1J2:ERVL-MaLR:LTR","MMTV-int:ERVK:LTR","MuLV-int:ERV1:LTR")),
    aes(x = log2FoldChange, y = -log10(pvalue), label = simple_geneName),  # <-- match to your x/y columns
    color = "black",           # <-- text color
    segment.color = "black",   # (optional) label line color
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 123
  ) +
  geom_vline(xintercept = 0,lty = 2)+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  # geom_text_repel(data = subset(res,annotation == "yes"),label=subset(res,annotation == "yes") %>% pull(GeneID),show.legend = FALSE,color = "black") +
  labs(title = paste0("Control vs Abx")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        legend.position = "bottom")

ggsave(WT.Naive.refine,file = "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_volcano_plot_refine.pdf",dpi = 600, width = 8,height = 6)
ggsave(WT.Naive.refine,file = "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_volcano_plot_refine_simpleGeneID.pdf",dpi = 600, width = 8,height = 6)



library(htmlwidgets)

# Save the plot as an HTML file
saveWidget(ggplotly(WT.Naive.refine), file = "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_volcano_plot_refine.html")



res.subset <- res %>%
              filter(pvalue >= 1E-50 & log2FoldChange)


WT.Naive.refine.subset <- ggplot(data = res.subset,
                          aes(x = log2FoldChange,
                              y = -log10(pvalue),
                              color = "grey",
                              label = GeneID)) +
  geom_point(data = subset(res.subset,res.subset$padj > 0.05 & res.subset$Gene_type == "Genes"),
             color = "grey",show.legend = FALSE) +
  geom_point(data = subset(res.subset,res.subset$padj > 0.05 & res.subset$Gene_type == "TEs"),
             color = "#A6A1A1",size = 1) +
  geom_point(data = subset(res.subset,res.subset$padj < 0.05 & res.subset$Gene_type == "Genes"),
             color = "#A6A1A1",size = 1) +
  geom_point(data = subset(res.subset,res.subset$padj < 0.05 & res.subset$Gene_type == "TEs"),
             color = "red",size = 1) +
  geom_point(
    data = subset(res.subset, GeneID %in% c("IAPLTR3:ERVK:LTR", "L1M4a2:L1:LINE", "L1ME2z:L1:LINE","MER129:LTR:LTR","MER4CL34:ERV1:LTR","MLT1J2:ERVL-MaLR:LTR","MMTV-int:ERVK:LTR","MuLV-int:ERV1:LTR")),
    aes(x = log2FoldChange, y = -log10(pvalue)),   # use your plot's x/y (e.g., y = neglog10p)
    shape = 21, fill = "red", color = "black",  # black border, hollow fill
    stroke = 0.8,                        # border thickness
    size = 1.8,                          # match your main point size
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(res.subset, GeneID %in% c("IAPLTR3:ERVK:LTR", "L1M4a2:L1:LINE", "L1ME2z:L1:LINE","MER129:LTR:LTR","MER4CL34:ERV1:LTR","MLT1J2:ERVL-MaLR:LTR","MMTV-int:ERVK:LTR","MuLV-int:ERV1:LTR")),
    aes(x = log2FoldChange, y = -log10(pvalue), label = GeneID),  # <-- match to your x/y columns
    color = "black",           # <-- text color
    segment.color = "black",   # (optional) label line color
    min.segment.length = 0,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.3,     # space around points
    max.overlaps = Inf,
    force = 2,               # stronger repulsion
    seed = 123
  ) +
  geom_vline(xintercept = 0,lty = 2)+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
  # geom_text_repel(data = subset(res,annotation == "yes"),label=subset(res,annotation == "yes") %>% pull(GeneID),show.legend = FALSE,color = "black") +
  labs(title = paste0("VNAM vs non-VNAM")) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        legend.position = "bottom")

ggsave(WT.Naive.refine.subset,file = "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Forest/TE/TE_volcano_plot_refine_subset.pdf",dpi = 600, width = 8,height = 5)

