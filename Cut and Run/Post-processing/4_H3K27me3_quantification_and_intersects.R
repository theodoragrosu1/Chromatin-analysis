library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from 3_reannotation_of_peaks.R script, because the 2_annotation_contrasts_peaks script only treats unique peaks
#input files: reannotated peaks (WT_unique, but then A7 and C9 that also contain the common peaks between them)
#input files: also the deduped bam files for each condition
#input files: ATAC peaks
##########this script also has to be run for other marks if intersect is of interest
#output files: quantification and RPKM normalization at every peak present (unique fot WT, and A7 and C9) based on read depth
#output files: GSEA analysis
 
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


WT_counts_me3 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\merged_bams\\WT_me3_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_me3_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "WT_counts_me3" = "value")) %>%
  na.omit()


A7_counts_me3 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\merged_bams\\A7_me3_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_me3_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "A7_counts_me3" = "value")) %>%
  na.omit()

C9_counts_me3 <- featureCounts("C:\\Users\\tgrosu\\OneDrive\\Desktop\\merged_bams\\C9_me3_deduped.bam", # count reads
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

write.csv(RPKM, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv", quote = FALSE, row.names = FALSE)

#load analysed RNAseq (DEGs)
wt_vs_clones_diff_expr <- read.csv(".\\Jurkat_wt_v_clones_all_genes (1).csv") %>%  #this is the list of DEGs (common between A7 and C9) that was published in Lefeivre et al (2025, Blood Advances); GSE240152
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))

####################Analysis using the RPKM normalization 
me3_RPKM <- read.csv(".\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv")
me3_and_rna <- left_join(me3_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()


rna_up_genes <- me3_and_rna %>%
  dplyr::filter(logFC >= 1) %>%
  #dplyr::filter(adj.P.Val < 0.1) %>%
  pull(geneID)

me3_up_a7 <- me3_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

me3_up_c9 <- me3_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

rna_down_genes <- me3_and_rna %>%
  dplyr::filter(logFC <= -1) %>%
  #dplyr::filter(adj.P.Val < 0.1) %>%
  pull(geneID)

me3_down_a7 <- me3_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

me3_down_c9 <- me3_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

wt_me3_lost <- c(me3_down_a7, me3_down_c9) #here: wt_me3_lost means the me3 sites that were lost in the clones

venn.diagram(list("A7 H3K27me3 \u2191" = me3_up_a7,
                  "C9 H3K27me3 \u2191" = me3_up_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 2, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#6E1372"),
             ".\\rna_me3_up_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2193" = me3_down_a7, 
                  "C9 H3K27me3 \u2193" = me3_down_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#E16036"),
             ".\\rna_me3_down_genes_up_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2193" = me3_down_a7, 
                  "C9 H3K27me3 \u2193" = me3_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#6E1372"),
             ".\\rna_me3_down_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2191" = me3_up_a7,
                  "C9 H3K27me3 \u2191" = me3_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#E16036"),
             ".\\rna_me3_up_genes_up_overlaps.png", disable.logging = TRUE)


me3_down_genes_up_common <- intersect(intersect(me3_down_a7, me3_down_c9), rna_up_genes) 
me3_up_genes_down_common <- intersect(intersect(me3_up_a7,me3_up_c9), rna_down_genes) 

me3_down_genes_down_common <- intersect(intersect(me3_down_a7, me3_down_c9), rna_down_genes) 
me3_up_genes_up_common <- intersect(intersect(me3_up_a7,me3_up_c9), rna_up_genes)

#perform GSEA 

library(WebGestaltR)
upreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                "geneontology_Molecular_Function_noRedundant"),
                             interestGeneType = "genesymbol",
                             referenceGeneType = "genesymbol",
                             interestGene = me3_down_genes_up_common,
                             referenceGene = wt_vs_clones_diff_expr$geneID,
                             sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:20,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones\n with decreased H3K27me3") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

####overlap between datasets to get genes of interest: 
#quantification of H3K27ac script has to also be run for this part

ac_gained_clones <- c(ac_up_a7, ac_up_c9)

intersect_me3_wt_ac_clones <- intersect(wt_me3_lost, ac_gained_clones) %>%
  as.matrix()


venn.diagram(list("WT H3K27me3 \u2193" = wt_me3_lost,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#42A5F5", "#4FC3F7"),
             ".\\H3K27me3_lost_WT_H3K27ac_gained_clones_overlaps.png", disable.logging = TRUE)

####get intersect for heatmaps with the ATAC peaks
WT_ATAC_peaks <- read.csv(".\\ATAC-peaks\\WT_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()

A7_ATAC_peaks <- read.csv(".\\ATAC-peaks\\A7_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()
C9_ATAC_peaks <- read.csv(".\\ATAC-peaks\\C9_unique_atac.csv") %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit()

clones_ATAC_peaks <- c(A7_ATAC_peaks, C9_ATAC_peaks)

venn.diagram(list("WT H3K27me3 \u2193" = wt_me3_lost,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "A7 ATAC peaks \u2191" = A7_ATAC_peaks),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#42A5F5", "#e27602"),
             ".\\H3K27me3_lost_WT_H3K27ac_ATAC_gained_A7_overlaps.png", disable.logging = TRUE)

venn.diagram(list("WT H3K27me3 \u2193" = wt_me3_lost,
                  "C9 H3K27ac \u2191" = ac_up_c9,
                  "C9 ATAC peaks \u2191" = C9_ATAC_peaks),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#4FC3F7", "#f1b04c"),
             ".\\H3K27me3_lost_WT_H3K27ac_ATAC_gained_C9_overlaps.png", disable.logging = TRUE)


genes_to_plot <- intersect(intersect(intersect(intersect(ac_up_a7, ac_up_c9), wt_me3_lost), A7_ATAC_peaks), C9_ATAC_peaks)

normalised_expression_jurkat <- read.csv(".\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv(".\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]


normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) %>%
  dplyr::filter(logFC > 3) #makes it more manageable, otherwise there are way too many genes for this to be visualized nicely in a heatmap

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))


normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")


#heatmap of selected genes
library(pheatmap)
library(RColorBrewer)
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))

set.seed(1234)
upreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                 "geneontology_Molecular_Function_noRedundant"),
                              interestGeneType = "genesymbol",
                              referenceGeneType = "genesymbol",
                              interestGene = genes_to_plot,
                              referenceGene = wt_vs_clones_diff_expr$geneID,
                              sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "PRC2 targets") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank())


venn.diagram(list("A7 H3K27me3 \u2191" = me3_up_a7,
                  "C9 H3K27me3 \u2191" = me3_up_c9,
                  "A7&C9 ATAC peaks \u2193" = WT_ATAC_peaks,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 1, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#c95b0c", "#6E1372"),
             ".\\H3K27me3_gained_clones_ATAC_WT_lost_overlaps.png", disable.logging = TRUE)


venn.diagram(list("A7&C9 H3K27me3 \u2193" = wt_me3_lost, 
                  "A7&C9 ATAC peaks \u2191" = clones_ATAC_peaks,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.5, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#c95b0c", "#E16036"),
             ".\\H3K27me3_lost_clones_ATAC_WT_gained_overlaps.png", disable.logging = TRUE)

ub_gained_clones <- c(ub_up_a7,ub_up_c9)
wt_ub_lost <- c(ub_down_a7, ub_down_c9) #here: wt_ub_lost means the ub sites that are lost in the clones

intersect_me3_wt_ub_wt <- intersect(wt_me3_lost, wt_ub_lost) %>%
  as.matrix()


venn.diagram(list("A7 H2AK119ub \u2191" = ub_up_a7,
                  "C9 H2AK119ub \u2191" = ub_up_c9,
                  "A7&C9 ATAC peaks \u2193" = WT_ATAC_peaks,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 1, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#c95b0c", "#6E1372"),
             ".\\H2AK119ub_gained_clones_ATAC_WT_lost_overlaps.png", disable.logging = TRUE)

