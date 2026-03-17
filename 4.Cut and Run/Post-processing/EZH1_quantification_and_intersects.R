library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from reannotation_of_peaks.R script

annotated_peaks_wt_EZH1_saf <- annotated_peaks_WT_EZH1 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_a7_EZH1_saf <- annotated_peaks_A7_EZH1 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_EZH1_saf <- annotated_peaks_C9_EZH1 %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

WT_counts_EZH1 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\WT_EZH1\\WT_EZH1_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_EZH1_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("WT_counts_EZH1" = "value") %>%
  na.omit()


A7_counts_EZH1 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\A7_EZH1\\A7_EZH1_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_EZH1_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))  %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("A7_counts_EZH1" = "value") %>%
  na.omit()

C9_counts_EZH1 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\deduped\\C9_EZH1\\C9_EZH1_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_EZH1_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))%>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("C9_counts_EZH1" = "value") %>%
  na.omit()


C9_counts_EZH1$GeneID <- as.character(C9_counts_EZH1$GeneID)
A7_counts_EZH1$GeneID <- as.character(A7_counts_EZH1$GeneID)
WT_counts_EZH1$GeneID <- as.character(WT_counts_EZH1$GeneID)


gene_lengths <- annotated_peaks_WT_EZH1 %>%
  rbind(annotated_peaks_A7_EZH1) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  dplyr::rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts_EZH1, A7_counts_EZH1, by = "GeneID") %>%
  full_join(C9_counts_EZH1, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts_EZH1 <- as.numeric(counts_matrix$WT_counts_EZH1)
counts_matrix$C9_counts_EZH1 <- as.numeric(counts_matrix$C9_counts_EZH1)
counts_matrix$A7_counts_EZH1 <- as.numeric(counts_matrix$A7_counts_EZH1)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts_EZH1", "A7_counts_EZH1", "C9_counts_EZH1")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)


RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts_EZH1 - RPKM$WT_counts_EZH1
RPKM$logFC_C9_WT <- RPKM$C9_counts_EZH1 - RPKM$WT_counts_EZH1

RPKM$geneID <- rownames(RPKM)

write.csv(RPKM, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\EZH1\\EZH1_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
####################questions
EZH1_RPKM <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\EZH1\\EZH1_rpkm.csv")
EZH1_and_rna <- left_join(EZH1_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

EZH1_up_a7 <- EZH1_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

EZH1_up_c9 <- EZH1_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

EZH1_down_a7 <- EZH1_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

EZH1_down_c9 <- EZH1_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

wt_EZH1_peaks <- c(EZH1_down_a7, EZH1_down_c9) #here: wt_EZH1_lost means the EZH1 sites that were lost in the clones
clones_EZH1_gained <- c(EZH1_up_a7, EZH1_up_c9)

venn.diagram(list("A7 EZH1 \u2191" = EZH1_up_a7,
                 "C9 EZH1 \u2191" = EZH1_up_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 2, cat.cex = 1.6, print.mode = "raw", margin = 0.2,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#ffaf42", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\rna_EZH1_up_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 EZH1 \u2193" = EZH1_down_a7, 
                  "C9 EZH1 \u2193" = EZH1_down_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\rna_EZH1_down_genes_up_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 EZH1 \u2193" = EZH1_down_a7, 
                  "C9 EZH1 \u2193" = EZH1_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffe86d", "#ffcf2f", "#6E1372"),
             "C:\\Users\\tgrosu\\Downloads\\rna_EZH1_down_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 EZH1 \u2191" = EZH1_up_a7,
                  "C9 EZH1 \u2191" = EZH1_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#ffaf42", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\rna_EZH1_up_genes_up_overlaps.png", disable.logging = TRUE)


EZH1_down_genes_up_common <- intersect(intersect(EZH1_down_a7, EZH1_down_c9), rna_up_genes) 
EZH1_up_genes_down_common <- intersect(intersect(EZH1_up_a7,EZH1_up_c9), rna_down_genes) 

EZH1_down_genes_down_common <- intersect(intersect(EZH1_down_a7, EZH1_down_c9), rna_down_genes) 
EZH1_up_genes_up_common <- intersect(intersect(EZH1_up_a7,EZH1_up_c9), rna_up_genes)

#non-canonical EZH1 is shown to be able to recruit H3K27ac to chromatin and PolII for active transcription

#quantification of H3K27ac script has to also be run for this part
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
ac_gained_clones <- c(ac_up_a7, ac_up_c9)

intersect_EZH1_wt_ac_clones <- intersect(wt_EZH1_peaks, ac_gained_clones)
intersect_EZH1_clones_ac_clones <- intersect(clones_EZH1_gained, ac_gained_clones)

venn.diagram(list("WT EZH1 \u2191" = wt_EZH1_peaks,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#42A5F5", "#4FC3F7"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_in_WT_H3K27ac_gained_clones_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7&C9 EZH1 \u2191" = clones_EZH1_gained,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffcf2f", "#42A5F5", "#4FC3F7"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_in_clones_H3K27ac_gained_clones_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7&C9 EZH1 \u2191" = clones_EZH1_gained,
                  "A7 H3K27ac \u2193" = ac_down_a7,
                  "C9 H3K27ac \u2193" = ac_down_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffcf2f", "#42A5F5", "#4FC3F7"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_in_clones_H3K27ac_lost_clones_overlaps.png", disable.logging = TRUE)


####get intersect for heatmaps
ATAC_RPKM <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\ATAC_rpkm.csv")
atac_and_rna <- left_join(ATAC_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

atac_up_a7 <- atac_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

atac_up_c9 <- atac_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

atac_down_a7 <- atac_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

atac_down_c9 <- atac_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)
atac_gained_clones <- c(atac_up_a7, atac_up_c9)

venn.diagram(list("WT EZH1 \u2191" = wt_EZH1_peaks,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "A7 ATAC peaks \u2191" = atac_up_a7,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#fe6e00", "#42A5F5", "#e27602", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_in_WT_H3K27ac_ATAC_gained_A7_overlaps_RNA.png", disable.logging = TRUE)


venn.diagram(list("WT EZH1 \u2191" = wt_EZH1_peaks,
                  "C9 H3K27ac \u2191" = ac_up_c9,
                  "C9 ATAC peaks \u2191" = atac_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#fe6e00", "#4FC3F7", "#f1b04c", "#E16036"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_lost_WT_H3K27ac_ATAC_gained_C9_overlaps_RNA.png", disable.logging = TRUE)


venn.diagram(list("A7&C9 EZH1 \u2191" = clones_EZH1_gained,
                  "A7 ATAC \u2193" = atac_down_a7,
                  "C9 ATAC \u2193"= atac_down_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#fe6e00", "#4FC3F7", "#f1b04c"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_lost_WT_ATAC.png", disable.logging = TRUE)



#hypothesis: in WT we might be predominantly seeing non-canonical EZH1 functions, where it corresponds a lot 
#to gene activation and RNA PolII recruitment
#we should hopefully see the A7 and C9 EZH1 peaks relating to WT ATAC regions (that are lost in A7 and C9)
WT_unique_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\WT_unique_atac.csv") %>%
  drop_na()
venn.diagram(list("A7 EZH1 \u2191" = EZH1_up_a7,
                  "A7 ATAC \u2193" = atac_down_a7),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#173518", "#ffae42"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_up_A7_ATAC_peaks_A7_down.png", disable.logging = TRUE)

venn.diagram(list("C9 EZH1 \u2191" = EZH1_up_c9,
                  "C9 ATAC \u2193"= atac_down_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#1138be", "#fe6e00"),
             "C:\\Users\\tgrosu\\Downloads\\EZH1_up_C9_ATAC_peaks_C9_down.png", disable.logging = TRUE)
