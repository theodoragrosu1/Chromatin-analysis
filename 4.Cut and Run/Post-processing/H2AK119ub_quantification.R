library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from 5_0_annotate_peaks_peak_width.R

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


WT_counts <- featureCounts("path\\merged_bams\\WT_ub_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_ub_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename(c("GeneID" = "Var1", "WT_counts" = "value")) %>%
  na.omit()


A7_counts <- featureCounts("path\\merged_bams\\A7_ub_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_ub_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename(c("GeneID" = "Var1", "A7_counts" = "value")) %>%
  na.omit()

C9_counts <- featureCounts("path\\merged_bams\\C9_ub_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_ub_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename(c("GeneID" = "Var1", "C9_counts" = "value")) %>%
  na.omit()


C9_counts$GeneID <- as.character(C9_counts$GeneID)
A7_counts$GeneID <- as.character(A7_counts$GeneID)
WT_counts$GeneID <- as.character(WT_counts$GeneID)


gene_lengths <- annotated_peaks_WT_ub %>%
  rbind(annotated_peaks_A7_ub) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  dplyr::rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts, A7_counts, by = "GeneID") %>%
  full_join(C9_counts, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts <- as.numeric(counts_matrix$WT_counts)
counts_matrix$C9_counts <- as.numeric(counts_matrix$C9_counts)
counts_matrix$A7_counts <- as.numeric(counts_matrix$A7_counts)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts", "A7_counts", "C9_counts")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)

RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_A7_WT <- RPKM$A7_counts - RPKM$WT_counts
RPKM$logFC_C9_WT <- RPKM$C9_counts - RPKM$WT_counts

RPKM$geneID <- rownames(RPKM)

write.csv(RPKM, file = "path\\macs2_peaks\\H2AK119ub\\contrasts\\ub_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
####################questions
RPKM_ub <- read.csv("path\\macs2_peaks\\H2AK119ub\\contrasts\\ub_rpkm.csv")
ub_and_rna <- left_join(RPKM_ub, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

rna_down_genes <- ub_and_rna %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

rna_up_genes <- ub_and_rna %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

ub_up_a7 <- ub_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

ub_up_c9 <- ub_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

ub_down_a7 <- ub_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

ub_down_c9 <- ub_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)




library(WebGestaltR)

ub_down_genes_up_common <- intersect(intersect(ub_down_a7, ub_down_c9), rna_up_genes) 
ub_up_genes_down_common <- intersect(intersect(ub_up_a7,ub_up_c9), rna_down_genes) 

ub_down_genes_down_common <- intersect(intersect(ub_down_a7, ub_down_c9), rna_down_genes) 
ub_up_genes_up_common <- intersect(intersect(ub_up_a7,ub_up_c9), rna_up_genes)  

upreg_genes_ub_down_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                "geneontology_Molecular_Function_noRedundant"),
                             interestGeneType = "genesymbol",
                             referenceGeneType = "genesymbol",
                             interestGene = ub_up_genes_up_common,
                             referenceGene = wt_vs_clones_diff_expr$geneID,
                             sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_ub_down_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:30,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones\n with decreased H2AK119ub") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 


venn.diagram(list("A7 H2AK119ub \u2191" = ub_up_a7,
                  "C9 H2AK119ub \u2191" = ub_up_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 2, cat.cex = 1.3, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#6E1372"),
             "path\\rna_ub_up_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H2AK119ub \u2193" = ub_down_a7, 
                  "C9 H2AK119ub \u2193" = ub_down_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#E16036"),
             "path\\rna_ub_down_genes_up_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H2AK119ub \u2193" = ub_down_a7, 
                  "C9 H2AK119ub \u2193" = ub_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#6E1372"),
             "path\\rna_ub_down_genes_down_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 H2AK119ub \u2191" = ub_up_a7,
                  "C9 H2AK119ub \u2191" = ub_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89", "#E16036"),
             "path\\rna_ub_up_genes_up_overlaps.png", disable.logging = TRUE)


