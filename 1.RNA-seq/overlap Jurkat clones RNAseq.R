#input file: normalized tpm from RNAseq for each sample
#input file: pairwise comparison of each EZH2 clone individually versus WT
#Output: overlap between clones of downregulated and upregulated genes 
#Output: Venn diagrams
#Output: GO enrichment for the overlapped genes

library(tidyverse)
library(VennDiagram)
library(WebGestaltR)
library(dplyr)
library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)

upregulated_a7 <- read.csv("#path to the differential expression analysis of clone A7 vs WT") %>%
  dplyr::filter(LogFC >= 1) %>%
  dplyr::pull("geneID")


upregulated_c9 <- read.csv("#path to the differential expression analysis of clone C9 vs WT") %>%
  dplyr::filter(LogFC >= 1) %>%
  dplyr::pull("geneID")

downregulated_a7 <- read.csv("#path to the differential expression analysis of clone A7 vs WT") %>%
  dplyr::filter(LogFC <= -1) %>%
  dplyr::pull("geneID")

downregulated_c9 <- read.csv("#path to the differential expression analysis of clone C9 vs WT") %>%
  dplyr::filter(LogFC <= -1) %>%
  dplyr::pull("geneID")

genes.down <- list(Jurkat_A7 = downregulated_a7, Jurkat_C9 = downregulated_c9)
genes.up <- list(Jurkat_A7 = upregulated_a7, Jurkat_C9 = upregulated_c9)

venn.diagram(genes.down, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#a9d6e5", "#014f86"), 
             "genes.down.venn.a7.c9.tiff")
venn.diagram(genes.up, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#ffccd5", "#c9184a"), 
             "genes.up.venn.a7.c9.tiff")


write.table(intersect(downregulated_a7, downregulated_c9), row.names = FALSE, quote = FALSE, "#path and name")
write.table(intersect(upregulated_a7, upregulated_c9), row.names = FALSE, quote = FALSE, col.names = FALSE, "#path and name")


normalised_expression_jurkat_vs_clones <- read.csv("#path to normalized count file")

upreg_genes_overlap_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                         "geneontology_Molecular_Function_noRedundant"),
                                      interestGeneType = "genesymbol",
                                      referenceGeneType = "genesymbol",
                                      interestGene = intersect(upregulated_a7, upregulated_c9),
                                      referenceGene = normalised_expression_jurkat_vs_clones$geneID,
                                      sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_overlap_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:10,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Overlap upregulated genes in clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

down_genes_overlap_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                        "geneontology_Molecular_Function_noRedundant"),
                                     interestGeneType = "genesymbol",
                                     referenceGeneType = "genesymbol",
                                     interestGene = intersect(downregulated_a7, downregulated_c9),
                                     referenceGene = normalised_expression_jurkat_vs_clones$geneID,
                                     sigMethod = "fdr", fdrThr = 0.5)
GSEA.chosen = down_genes_overlap_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:10,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Overlap downregulated genes in clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 
