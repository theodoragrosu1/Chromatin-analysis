library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)
library(dplyr)

#this script uses output from reannotation_of_peaks.R script

annotated_peaks_wt_SUZ12_saf <- annotated_peaks_WT_SUZ12 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_a7_SUZ12_saf <- annotated_peaks_A7_SUZ12 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_SUZ12_saf <- annotated_peaks_C9_SUZ12 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

WT_counts_SUZ12 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\WT_SUZ12\\WT_SUZ12_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_SUZ12_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("WT_counts_SUZ12" = "value") %>%
  na.omit()


A7_counts_SUZ12 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\A7_SUZ12\\A7_SUZ12_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_SUZ12_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))  %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("A7_counts_SUZ12" = "value") %>%
  na.omit()

C9_counts_SUZ12 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\C9_SUZ12\\C9_SUZ12_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_SUZ12_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))%>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("C9_counts_SUZ12" = "value") %>%
  na.omit()


C9_counts_SUZ12$GeneID <- as.character(C9_counts_SUZ12$GeneID)
A7_counts_SUZ12$GeneID <- as.character(A7_counts_SUZ12$GeneID)
WT_counts_SUZ12$GeneID <- as.character(WT_counts_SUZ12$GeneID)


gene_lengths <- annotated_peaks_WT_SUZ12 %>%
  rbind(annotated_peaks_A7_SUZ12) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  dplyr::rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts_SUZ12, A7_counts_SUZ12, by = "GeneID") %>%
  full_join(C9_counts_SUZ12, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts_SUZ12 <- as.numeric(counts_matrix$WT_counts_SUZ12)
counts_matrix$C9_counts_SUZ12 <- as.numeric(counts_matrix$C9_counts_SUZ12)
counts_matrix$A7_counts_SUZ12 <- as.numeric(counts_matrix$A7_counts_SUZ12)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts_SUZ12", "A7_counts_SUZ12", "C9_counts_SUZ12")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)


RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts_SUZ12 - RPKM$WT_counts_SUZ12
RPKM$logFC_C9_WT <- RPKM$C9_counts_SUZ12 - RPKM$WT_counts_SUZ12

RPKM$geneID <- rownames(RPKM)

write.csv(RPKM, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
####################questions
SUZ12_RPKM <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_rpkm.csv")
SUZ12_and_rna <- left_join(SUZ12_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

SUZ12_up_a7 <- SUZ12_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

SUZ12_up_c9 <- SUZ12_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

SUZ12_down_a7 <- SUZ12_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

SUZ12_down_c9 <- SUZ12_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

wt_SUZ12_lost <- c(SUZ12_down_a7, SUZ12_down_c9) #here: wt_SUZ12_lost means the SUZ12 sites that were lost in the clones
clones_SUZ12_gained <- c(SUZ12_up_a7, SUZ12_up_c9)

venn.diagram(list("A7 SUZ12 \u2191" = SUZ12_up_a7,
                 "C9 SUZ12 \u2191" = SUZ12_up_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 2, cat.cex = 1.6, print.mode = "raw", margin = 0.2,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#ffaf42", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\rna_SUZ12_up_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 SUZ12 \u2193" = SUZ12_down_a7, 
                  "C9 SUZ12 \u2193" = SUZ12_down_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\rna_SUZ12_down_genes_up_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 SUZ12 \u2193" = SUZ12_down_a7, 
                  "C9 SUZ12 \u2193" = SUZ12_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\rna_SUZ12_down_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 SUZ12 \u2191" = SUZ12_up_a7,
                  "C9 SUZ12 \u2191" = SUZ12_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#ffaf42", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\rna_SUZ12_up_genes_up_overlaps.png", disable.logging = TRUE)


SUZ12_down_genes_up_common <- intersect(intersect(SUZ12_down_a7, SUZ12_down_c9), rna_up_genes) 
SUZ12_up_genes_down_common <- intersect(intersect(SUZ12_up_a7,SUZ12_up_c9), rna_down_genes) 

SUZ12_down_genes_down_common <- intersect(intersect(SUZ12_down_a7, SUZ12_down_c9), rna_down_genes) 
SUZ12_up_genes_up_common <- intersect(intersect(SUZ12_up_a7,SUZ12_up_c9), rna_up_genes)

WT_ATAC_WT_SUZ12 <- intersect(WT_ATAC_peaks, wt_SUZ12_lost) #half of SUZ12 peaks are at WT ATAC regions

#get H3K27me3
me3_RPKM <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv")
me3_and_rna <- left_join(me3_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()


me3_up_a7 <- me3_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

me3_up_c9 <- me3_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

me3_down_a7 <- me3_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

me3_down_c9 <- me3_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

wt_me3_lost <- c(me3_down_a7, me3_down_c9) #here: wt_me3_lost means the me3 sites that were lost in the clones
clones_me3_gained <- c(me3_up_a7, me3_up_c9)

wt_suz12_wt_me3 <- intersect(wt_SUZ12_lost, wt_me3_lost)

wt_atac_wt_suz12_wt_me3 <- intersect(WT_ATAC_peaks, wt_suz12_wt_me3) #most of the peaks that have me3 and SUZ12 also have open chromatin


venn.diagram(list("WT SUZ12 \u2191" = wt_SUZ12_lost, 
                  "WT H3K27me3 \u2191" = wt_me3_lost,
                  "WT ATAC \u2191" = WT_ATAC_peaks),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\WT_ATAC_SUZ12_WT_me3.png", disable.logging = TRUE)


venn.diagram(list("WT SUZ12 \u2191" = wt_SUZ12_lost, 
                  "WT H3K27me3 \u2191" = wt_me3_lost,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\SUZ12_me3_WT_genes_up_clones_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7&C9 SUZ12 \u2191" = clones_SUZ12_gained,
                  "A7&C9 H3K27me3 \u2191" = clones_me3_gained,
                  "A7&C9 RNA \u2193" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.15,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#ffaf42", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\SUZ12_me3_clones_genes_down_clones_overlaps.png", disable.logging = TRUE)


####get intersect for heatmaps
WT_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\WT_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()

A7_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\A7_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()
C9_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\C9_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()

clones_ATAC_peaks <- c(A7_ATAC_peaks, C9_ATAC_peaks)


#get H3K27ac
ac_RPKM <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\macs2_peaks\\H3K27ac\\contrasts\\ac_rpkm.csv")
ac_and_rna <- left_join(ac_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()


ac_up_a7 <- ac_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

ac_up_c9 <- ac_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

ac_down_a7 <- ac_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

ac_down_c9 <- ac_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

wt_ac_lost <- c(ac_down_a7, ac_down_c9) #here: wt_ac_lost means the ac sites that were lost in the clones
clones_ac_gained <- c(ac_up_a7, ac_up_c9)

venn.diagram(list("WT SUZ12 \u2191" = wt_SUZ12_lost,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "A7 ATAC peaks \u2191" = A7_ATAC_peaks,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#fe6e00", "#42A5F5", "#e27602", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\SUZ12_in_WT_H3K27ac_ATAC_gained_A7_overlaps_RNA.png", disable.logging = TRUE)


venn.diagram(list("WT SUZ12 \u2191" = wt_SUZ12_lost,
                  "C9 H3K27ac \u2191" = ac_up_c9,
                  "C9 ATAC peaks \u2191" = C9_ATAC_peaks,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#fe6e00", "#4FC3F7", "#f1b04c", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\SUZ12_lost_WT_H3K27ac_ATAC_gained_C9_overlaps_RNA.png", disable.logging = TRUE)

WT_unique_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\WT_unique_atac.csv") %>%
  drop_na()
venn.diagram(list("A7 SUZ12 \u2191" = SUZ12_up_a7,
                  "C9 SUZ12 \u2191" = SUZ12_up_c9,
                  "WT ATAC peaks \u2191" = WT_unique_ATAC_peaks$SYMBOL),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#1138be", "#ffae42"),
             "C:\\Users\\tgrosu\\Downloads\\SUZ12_up_clones_ATAC_peaks_WT.png", disable.logging = TRUE)

