library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from reannotation_of_peaks.R script 
#ooutput: correlation analysis with gene expression

annotated_peaks_wt_ub_saf <- annotated_peaks_WT_ub %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_a7_ub_saf <- annotated_peaks_A7_ub %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_ub_saf <- annotated_peaks_C9_ub %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()


WT_counts_ub <- featureCounts(".\\merged_bams\\WT_ub_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_wt_ub_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "WT_counts_ub" = "value")) %>%
  na.omit()


A7_counts_ub <- featureCounts(".\\merged_bams\\A7_ub_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_a7_ub_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "A7_counts_ub" = "value")) %>%
  na.omit()

C9_counts_ub <- featureCounts(".\\merged_bams\\C9_ub_deduped.bam", # count reads
                               annot.inbuilt = "hg38",
                               annot.ext = annotated_peaks_c9_ub_saf,
                               isGTFAnnotationFile = F,
                               isPairedEnd = T,
                               nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "C9_counts_ub" = "value")) %>%
  na.omit()


C9_counts_ub$GeneID <- as.character(C9_counts_ub$GeneID)
A7_counts_ub$GeneID <- as.character(A7_counts_ub$GeneID)
WT_counts_ub$GeneID <- as.character(WT_counts_ub$GeneID)


gene_lengths <- annotated_peaks_WT_ub %>%
  rbind(annotated_peaks_A7_ub) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts_ub, A7_counts_ub, by = "GeneID") %>%
  full_join(C9_counts_ub, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts_ub <- as.numeric(counts_matrix$WT_counts_ub)
counts_matrix$C9_counts_ub <- as.numeric(counts_matrix$C9_counts_ub)
counts_matrix$A7_counts_ub <- as.numeric(counts_matrix$A7_counts_ub)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts_ub", "A7_counts_ub", "C9_counts_ub")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)

RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts_ub - RPKM$WT_counts_ub
RPKM$logFC_C9_WT <- RPKM$C9_counts_ub - RPKM$WT_counts_ub

RPKM$geneID <- rownames(RPKM)


expression_wt_vs_a7 <- read.csv(".\\RNAseq\\normalised_expression_jurkat_wt_vs_a7.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))
expression_wt_vs_c9 <- read.csv(".\\RNAseq\\normalised_expression_jurkat_wt_vs_c9.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))

RPKM_A7 <- RPKM[,-5]

ub_and_rna <- left_join(RPKM_A7, expression_wt_vs_a7, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(ub_and_rna$logFC_A7_WT, ub_and_rna$LogFC, method = "pearson")


a7_ub_gene_expr <- ggplot(ub_and_rna, aes(x = logFC_A7_WT, y = LogFC)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  geom_point(size = 1.5, alpha = 0.3, color = "#8B008B") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "H2AK119ub A7 - WT", y = "mRNA A7 - WT") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R = ", round(cor_results$estimate, 2),
                          "\np = ", format.pval(cor_results$p.value, digits = 2)),
           size = 4.5, color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

ggsave(a7_ub_gene_expr, file = ".\\scatter_density_a7_ub.pdf", width = 7, height = 6, dpi = 300)

RPKM_C9 <- RPKM[,-4]

ub_and_rna <- left_join(RPKM_C9, expression_wt_vs_c9, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(ub_and_rna$logFC_C9_WT, ub_and_rna$LogFC, method = "pearson")


c9_ub_gene_expr <- ggplot(ub_and_rna, aes(x = logFC_C9_WT, y = LogFC)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  geom_point(size = 1.5, alpha = 0.3, color = "#FFA500") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "H2AK119ub C9 - WT", y = "mRNA C9 - WT") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("R = ", round(cor_results$estimate, 2),
                          "\np = ", format.pval(cor_results$p.value, digits = 2)),
           size = 4.5, color = "black") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")
ggsave(c9_ub_gene_expr, file = ".\\scatter_density_c9_ub.pdf", width = 7, height = 6, dpi = 300)

