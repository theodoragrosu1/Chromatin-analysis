library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(ggpubr)

#this script uses output from 3_reannotation_of_peaks.R script 
#output: correlation analysis with gene expression

annotated_peaks_wt_me3_saf <- annotated_peaks_WT_me3 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_a7_me3_saf <- annotated_peaks_A7_me3 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_me3_saf <- annotated_peaks_C9_me3 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()


WT_counts_me3 <- featureCounts(".\\merged_bams\\WT_me3_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_wt_me3_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "WT_counts_me3" = "value")) %>%
  na.omit()


A7_counts_me3 <- featureCounts(".\\merged_bams\\A7_me3_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_a7_me3_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "A7_counts_me3" = "value")) %>%
  na.omit()

C9_counts_me3 <- featureCounts(".\\merged_bams\\C9_me3_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_c9_me3_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "C9_counts_me3" = "value")) %>%
  na.omit()


C9_counts_me3$GeneID <- as.character(C9_counts_me3$GeneID)
A7_counts_me3$GeneID <- as.character(A7_counts_me3$GeneID)
WT_counts_me3$GeneID <- as.character(WT_counts_me3$GeneID)


gene_lengths <- annotated_peaks_WT_me3 %>%
  rbind(annotated_peaks_A7_me3) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts_me3, A7_counts_me3, by = "GeneID") %>%
  full_join(C9_counts_me3, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts_me3 <- as.numeric(counts_matrix$WT_counts_me3)
counts_matrix$C9_counts_me3 <- as.numeric(counts_matrix$C9_counts_me3)
counts_matrix$A7_counts_me3 <- as.numeric(counts_matrix$A7_counts_me3)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts_me3", "A7_counts_me3", "C9_counts_me3")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)

RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts_me3 - RPKM$WT_counts_me3
RPKM$logFC_C9_WT <- RPKM$C9_counts_me3 - RPKM$WT_counts_me3

RPKM$geneID <- rownames(RPKM)


expression_wt_vs_a7 <- read.csv(".\\normalised_expression_jurkat_wt_vs_a7.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                                 LogFC>0 ~ "Upregulated"))
expression_wt_vs_c9 <- read.csv(".\\normalised_expression_jurkat_wt_vs_c9.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))

RPKM_A7 <- RPKM[,-5]

me3_and_rna <- left_join(RPKM_A7, expression_wt_vs_a7, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(me3_and_rna$logFC_A7_WT, me3_and_rna$LogFC, method = "pearson")


a7_me3_gene_expr <- ggplot(me3_and_rna, aes(x = logFC_A7_WT, y = LogFC)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  geom_point(size = 1.5, alpha = 0.3, color = "#8B008B") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "H3K27me3 A7 - WT", y = "mRNA A7 - WT") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R = ", round(cor_results$estimate, 2),
                          "\np = ", format.pval(cor_results$p.value, digits = 2)),
           size = 4.5, color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

a7_me3_vs_expr <- left_join(RPKM_A7, expression_wt_vs_a7, by = "geneID") %>%
  drop_na(LogFC) %>%
  #dplyr::filter(geneID %in% unique(c(annotated_peaks_wt_me3_saf$GeneID, annotated_peaks_c5_me3_saf$GeneID))) %>%
  ggplot(aes(x = logFC_A7_WT, y = LogFC)) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 1, size = 5) + 
  labs(x = "H3K27me3 LFC\nA7 - WT", y = "mRNA LFC\nA7 - WT") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null")


ggsave(a7_me3_gene_expr, file = ".\\scatter_density_a7_me3.pdf", width = 7, height = 6, dpi = 300)

RPKM_C9 <- RPKM[,-4]

me3_and_rna <- left_join(RPKM_C9, expression_wt_vs_c9, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(me3_and_rna$logFC_C9_WT, me3_and_rna$LogFC, method = "pearson")


c9_me3_gene_expr <- ggplot(me3_and_rna, aes(x = logFC_C9_WT, y = LogFC)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  geom_point(size = 1.5, alpha = 0.3, color = "#FFA500") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "H3K27me3 C9 - WT", y = "mRNA C9 - WT") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R = ", round(cor_results$estimate, 2),
                          "\np = ", format.pval(cor_results$p.value, digits = 2)),
           size = 4.5, color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")
ggsave(c9_me3_gene_expr, file = ".\\scatter_density_c9_me3.pdf", width = 7, height = 6, dpi = 300)

