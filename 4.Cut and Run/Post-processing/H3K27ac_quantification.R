library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from reannotation_of_peaks.R script

annotated_peaks_wt_ac_saf <- annotated_peaks_WT_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

####for A7 and C9, we used the "_ac_specific_peaks.bed" files that include the clone-specific peaks, but also the peaks from the overlap between the clones
annotated_peaks_a7_ac_saf <- annotated_peaks_A7_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_ac_saf <- annotated_peaks_C9_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()


WT_counts <- featureCounts("path\\merged_bams\\WT_ac_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_ac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename(c("GeneID" = "Var1", "WT_counts" = "value")) %>%
  na.omit()


A7_counts <- featureCounts("path\\merged_bams\\A7_ac_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_a7_ac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  dplyr::rename(c("GeneID" = "Var1", "A7_counts" = "value")) %>%
  na.omit()

C9_counts <- featureCounts("path\\merged_bams\\C9_ac_deduped_filtered.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_ac_saf,
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


gene_lengths <- annotated_peaks_WT_ac %>%
  rbind(annotated_peaks_A7_ac) %>%
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

write.csv(RPKM, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\ac_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
####################questions
ac_RPKM <- read.csv("path\\macs2_peaks\\H3K27ac\\contrasts\\ac_rpkm.csv")
ac_and_rna <- left_join(ac_RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  #dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

ac_up_a7 <- ac_and_rna %>%
  dplyr::filter(logFC_A7_WT >= 0.5) %>%
  pull(geneID)%>%
  as.matrix()

ac_up_c9 <- ac_and_rna %>%
  dplyr::filter(logFC_C9_WT >= 0.5) %>%
  pull(geneID) 

rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  #dplyr::filter(adj.P.Val < 0.2) %>%
  pull(geneID)

ac_down_a7 <- ac_and_rna %>%
  dplyr::filter(logFC_A7_WT <= -0.5) %>%
  pull(geneID)

ac_down_c9 <- ac_and_rna %>%
  dplyr::filter(logFC_C9_WT <= -0.5) %>%
  pull(geneID)

venn.diagram(list("A7 H3K27ac \u2191" = ac_up_a7,
                  "C9 H3K27ac \u2191" = ac_up_c9,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#ff5a89","#E16036"),
             "path\\rna_ac_up_genes_up_overlaps.png", disable.logging = TRUE)



venn.diagram(list("A7 H3K27ac \u2193" = ac_down_a7, 
                  "C9 H3K27ac \u2193" = ac_down_c9,
                  "A7&C9 RNA \u2193" = rna_down_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#9ECAE1", "#2171B6", "#6E1372"),
             "path\\rna_ac_down_genes_down_overlaps.png", disable.logging = TRUE)


ac_down_genes_down_common <- intersect(intersect(ac_down_a7, ac_down_c9), rna_down_genes) 
ac_up_genes_up_common <- intersect(intersect(ac_up_a7,ac_up_c9), rna_up_genes) 


library(WebGestaltR)


downreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                 "geneontology_Molecular_Function_noRedundant"),
                              interestGeneType = "genesymbol",
                              referenceGeneType = "genesymbol",
                              interestGene = ac_down_genes_down_common,
                              referenceGene = wt_vs_clones_diff_expr$geneID,
                              sigMethod = "fdr", fdrThr = 0.1)
 

upreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                "geneontology_Molecular_Function_noRedundant"),
                             interestGeneType = "genesymbol",
                             referenceGeneType = "genesymbol",
                             interestGene = ac_up_genes_up_common,
                             referenceGene = wt_vs_clones_diff_expr$geneID,
                             sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = downreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:20,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones\n with increased H3K27ac") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 


GSEA.chosen = upreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones\n with increased H3K27ac") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 
