library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from reannotation_of_peaks.R script

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


WT_counts_me3 <- featureCounts("path\\merged_bams\\WT_me3_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_me3_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("WT_counts_me3" = "value") %>%
  na.omit()


A7_counts_me3 <- featureCounts("path\\merged_bams\\A7_me3_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_me3_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))  %>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("A7_counts_me3" = "value") %>%
  na.omit()

C9_counts_me3 <- featureCounts("path\\merged_bams\\C9_me3_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_me3_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value"))%>%
  dplyr::rename("GeneID" = "Var1") %>%
  dplyr::rename("C9_counts_me3" = "value") %>%
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
  dplyr::rename("GeneID" = "SYMBOL")

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

write.csv(RPKM, file = "path\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
####################questions
me3_RPKM <- read.csv("path\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv")
me3_and_rna <- left_join(me3_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

me3_up_a7 <- me3_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

me3_up_c9 <- me3_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
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
             lwd = 0, cex = 2, cat.cex = 1.6, print.mode = "raw", margin = 0.2,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#6E1372"),
             "path\\rna_me3_up_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2193" = me3_down_a7, 
                  "C9 H3K27me3 \u2193" = me3_down_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#E16036"),
             "path\\rna_me3_down_genes_up_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2193" = me3_down_a7, 
                  "C9 H3K27me3 \u2193" = me3_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#6E1372"),
             "path\\rna_me3_down_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H3K27me3 \u2191" = me3_up_a7,
                  "C9 H3K27me3 \u2191" = me3_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#E16036"),
             "path\\rna_me3_up_genes_up_overlaps.png", disable.logging = TRUE)


me3_down_genes_up_common <- intersect(intersect(me3_down_a7, me3_down_c9), rna_up_genes) 
me3_up_genes_down_common <- intersect(intersect(me3_up_a7,me3_up_c9), rna_down_genes) 

me3_down_genes_down_common <- intersect(intersect(me3_down_a7, me3_down_c9), rna_down_genes) 
me3_up_genes_up_common <- intersect(intersect(me3_up_a7,me3_up_c9), rna_up_genes)


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
  labs(y = "GO term", title = "Genes with absent H3K27me3\n and increased H3K27ac in clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

####overlap between datasets to get genes of interest: 

ac_RPKM <- read.csv("path\\macs2_peaks\\H3K27ac\\contrasts\\ac_rpkm.csv")
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

intersect_me3_wt_ac_clones <- intersect(wt_me3_lost, ac_gained_clones)


venn.diagram(list("WT H3K27me3 \u2193" = wt_me3_lost,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9, 
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5, 0.5), fill = c("#173518", "#42A5F5", "#4FC3F7", "#E16036"),
             "path\\H3K27me3_lost_WT_H3K27ac_gained_clones_overlaps_RNA.png", disable.logging = TRUE)


venn.diagram(list("WT H3K27me3 \u2193" = wt_me3_lost,
                  "A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.2,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#42A5F5", "#4FC3F7"),
             "path\\H3K27me3_lost_WT_H3K27ac_gained_clones_overlaps.png", disable.logging = TRUE)



#BMP-17 signature Z-scores
genes_to_plot <- c("S100A4", "LGALS1", "FAM30A", "IGLL1", "MEF2C", "CTSW", "PRDX1",
                   "HCST", "HSH2D", "KLF2", "VAMP8", "HOPX", "CYBA", "MT-ND3",
                   "ENO1", "PEBP1", "CD44")

normalised_expression_jurkat <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]


normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) 

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


#annotate the increase in /h3K27ac peaks, where H3K27me3 used to be in WT

dir <- "path\\peaks\\"

#load peak files
peak_files <- list(A7_ac_WT_me3 = paste0(dir, "H3K27ac_replaces_H3K27me3_A7.bed"),
                   C9_ac_WT_me3 = paste0(dir, "H3K27ac_replaces_H3K27me3_C9.bed"))

#Annotate peaks for each condition and save files
annotated_peaks_A7_ac_WT_me3 <- annotate_peaks(peak_files, "A7_ac_WT_me3") 
annotated_peaks_C9_ac_WT_me3 <- annotate_peaks(peak_files, "C9_ac_WT_me3")

#Count number of peaks in each type of annotation for peaks
A7_ac_WT_me3_summary <- as.data.frame(table(annotated_peaks_A7_ac_WT_me3$annotation))
colnames(A7_ac_WT_me3_summary) <- c("Annotation", "A7_ac_WT_me3")

C9_ac_WT_me3_summary <- as.data.frame(table(annotated_peaks_C9_ac_WT_me3$annotation))
colnames(C9_ac_WT_me3_summary) <- c("Annotation", "C9_ac_WT_me3")

summary_annotation_me3 <- left_join(A7_ac_WT_me3_summary, C9_ac_WT_me3_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_me3$Annotation <- factor(summary_annotation_me3$Annotation, 
                                            levels = c("3' UTR", "Exon", 
                                                       "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                       "Distal Intergenic", "Intron"))

summary_annotation_me3$Condition <- factor(summary_annotation_me3$Condition, 
                                           levels = c("A7_ac_WT_me3", "C9_ac_WT_me3"))

#plot
annotated_me3_peaks_plot <- summary_annotation_me3 %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 3000)) +
  labs(x = "Sample") +
  ylab("No. H3K27ac peaks replacing H3K27me3")
