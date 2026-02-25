library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)
library(dplyr)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Cairo)
library(colorspace)

#Function to annotate peak files
annotate_peaks <- function(peak_file_name, condition) {
  annotated_peaks <- annotatePeak(peak_file_name[[condition]], tssRegion = c(-3000, 3000),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
  annotated_peaks_df <- as.data.frame(annotated_peaks)
  annotated_peaks_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_df$annotation) #This merges all Intron/Exon as one annotation, independent of intron/exon number
  annotated_peaks_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_df$annotation)
  return(annotated_peaks_df)
}


#load annotations
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "path\\ATAC-peaks\\"

#load peak files
peak_files <- list(WT_ATAC_peaks = paste0(dir, "WT_ATAC_peaks.bed"),
                   A7_ATAC_peaks = paste0(dir, "A7_ATAC_peaks.bed"),
                   C9_ATAC_peaks = paste0(dir, "C9_ATAC_peaks.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_ATAC <- annotate_peaks(peak_files, "WT_ATAC_peaks") 
annotated_peaks_A7_ATAC <- annotate_peaks(peak_files, "A7_ATAC_peaks")
annotated_peaks_C9_ATAC <- annotate_peaks(peak_files, "C9_ATAC_peaks")

annotated_peaks_wt_atac_saf <- annotated_peaks_WT_ATAC %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_a7_atac_saf <- annotated_peaks_A7_ATAC %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_atac_saf <- annotated_peaks_C9_ATAC %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()


WT_counts_ATAC <- featureCounts("path\\ATAC-peaks\\ATACseq_bam_files\\WT.mononucleosomal.sorted.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_atac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("WT_counts_ATAC" = "value") %>%
  na.omit()

A7_counts_ATAC <- featureCounts("path\\ATAC-peaks\\ATACseq_bam_files\\A7.mononucleosomal.sorted.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_atac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))  %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("A7_counts_ATAC" = "value") %>%
  na.omit()

C9_counts_ATAC <- featureCounts("path\\ATAC-peaks\\ATACseq_bam_files\\C9.mononucleosomal.sorted.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_atac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))%>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("C9_counts_ATAC" = "value") %>%
  na.omit()


C9_counts_ATAC$GeneID <- as.character(C9_counts_ATAC$GeneID)
A7_counts_ATAC$GeneID <- as.character(A7_counts_ATAC$GeneID)
WT_counts_ATAC$GeneID <- as.character(WT_counts_ATAC$GeneID)


gene_lengths <- annotated_peaks_WT_ATAC %>%
  rbind(annotated_peaks_A7_ATAC) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  dplyr::rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts_ATAC, A7_counts_ATAC, by = "GeneID") %>%
  full_join(C9_counts_ATAC, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts_ATAC <- as.numeric(counts_matrix$WT_counts_ATAC)
counts_matrix$C9_counts_ATAC <- as.numeric(counts_matrix$C9_counts_ATAC)
counts_matrix$A7_counts_ATAC <- as.numeric(counts_matrix$A7_counts_ATAC)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts_ATAC", "A7_counts_ATAC", "C9_counts_ATAC")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)


RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts_ATAC - RPKM$WT_counts_ATAC
RPKM$logFC_C9_WT <- RPKM$C9_counts_ATAC - RPKM$WT_counts_ATAC
RPKM$geneID <- rownames(RPKM)


write.csv(RPKM, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\ATAC_rpkm.csv", quote = FALSE, row.names = FALSE)

open_clones <- RPKM %>%
  dplyr::filter (logFC_A7_WT > 0) %>%
  dplyr::filter (logFC_C9_WT > 0)
open_clones$geneID <- rownames(open_clones)

open_wt <- RPKM %>%
  dplyr::filter (logFC_A7_WT < 0) %>%
  dplyr::filter (logFC_C9_WT< 0)
open_wt$geneID <- rownames(open_wt)

open_clones_avg <- open_clones %>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)

open_wt_avg <- open_wt %>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)

RPKM_avg <- RPKM %>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)


##import RNAseq
normalised_expression_wt_vs_clones <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

normalised_expression_wt_vs_clones$geneID <- as.character(normalised_expression_wt_vs_clones$geneID)
open_clones_avg$geneID <- as.character(open_clones_avg$geneID)
open_wt_avg$geneID <- as.character(open_wt_avg$geneID)

genes_with_open_atac_clones <- normalised_expression_wt_vs_clones %>%
  inner_join(open_clones_avg, by = "geneID")
median_log2FC <- median(genes_with_open_atac_clones$logFC)

plot_data <- genes_with_open_atac_clones %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Expression by RNAseq", avg_logFC = "Opening in clones"))


library(ggplot2)
ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_boxplot(aes(fill = type), alpha = 0.5, width = 0.6) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs ATAC Opening") +
  theme_classic() +
  theme(legend.position = "none")


#####
genes_with_open_atac_wt <- normalised_expression_wt_vs_clones %>%
  inner_join(open_wt_avg, by = "geneID")
median_log2FC <- median(genes_with_open_atac_wt$logFC)

plot_data <- genes_with_open_atac_wt %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Expression by RNAseq", avg_logFC = "Closing in clones"))


library(ggplot2)
ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_boxplot(aes(fill = type), alpha = 0.5, width = 0.6) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs ATAC Opening") +
  theme_classic() +
  theme(legend.position = "none")


####better question: at which median chromatin logFC do we see upregulation in our clones?
upregulated_genes_clones <- normalised_expression_wt_vs_clones %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val <= 0.1)

upregulated_genes_wt <- normalised_expression_wt_vs_clones %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val <= 0.1)

upreg_genes_with_atac_clones <- upregulated_genes_clones %>%
  inner_join(RPKM_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_atac_clones$avg_logFC)

upreg_no_ATAC_change_clones <- upreg_genes_with_atac_clones %>%
  dplyr::filter(avg_logFC < median_log2FC) 

write.csv(upreg_no_ATAC_change_clones, "path\\upreg_clones_no_ATAC_change.csv")

library(WebGestaltR)
upreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                 "geneontology_Molecular_Function_noRedundant"),
                              interestGeneType = "genesymbol",
                              referenceGeneType = "genesymbol",
                              interestGene = upreg_no_ATAC_change_clones$geneID,
                              referenceGene = normalised_expression_wt_vs_clones$geneID,
                              sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones\n no ATAC change") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 




####better question: at which median chromatin logFC do we see downregulation in our clones?
upregulated_genes_wt <- normalised_expression_wt_vs_clones %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val <= 0.1)

upreg_genes_with_atac_wt <- upregulated_genes_wt %>%
  inner_join(RPKM_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_atac_wt$avg_logFC)

upreg_no_ATAC_change_wt <- upreg_genes_with_atac_wt %>%
  dplyr::filter(avg_logFC > median_log2FC) 

downreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                 "geneontology_Molecular_Function_noRedundant"),
                              interestGeneType = "genesymbol",
                              referenceGeneType = "genesymbol",
                              interestGene = upreg_no_ATAC_change_wt$geneID,
                              referenceGene = normalised_expression_wt_vs_clones$geneID,
                              sigMethod = "fdr", fdrThr = 0.1)

write.csv(upreg_no_ATAC_change_wt, "path\\upreg_wt_no_ATAC_change.csv")
