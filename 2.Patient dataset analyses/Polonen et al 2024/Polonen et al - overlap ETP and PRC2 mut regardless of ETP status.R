#input file: genes up- and downregulated in all the conditions of interest
#output file: Venn diagram with overlaps
#output file: GO enrichment of overlapped genes

library(tidyverse)
library(VennDiagram)
library(WebGestaltR)
library(dplyr)

all_etp <- read.csv("#path")
all_prc2_mut <- read.csv("#path")

upregulated_patients_etp <- read.csv("#path") %>%
  dplyr::filter(condition == "Genes_upregulated_in_ETP") %>%
  dplyr::pull("geneID")

upregulated_prc2_mut <- read.csv("#path") %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::pull("SYMBOL")

downregulated_patients_etp <- read.csv("#path") %>%
  dplyr::filter(condition == "Genes_upregulated_in_T-ALL") %>%
  dplyr::pull("geneID")

downregulated_prc2_mut <- read.csv("#path") %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::pull("SYMBOL")

down_ezh2_ko <- read.csv("#path to Jurkat DEGs") %>%
  dplyr::filter(logFC <= -1) %>%
  dplyr::filter(adj.P.Val >= 0.05) %>%
  dplyr::pull("geneID")

up_ezh2_ko <- read.csv("#path to Jurkat DEGs") %>%
  dplyr::filter(logFC >= 1) %>%
  dplyr::filter(adj.P.Val >= 0.05) %>%
  dplyr::pull("geneID")

genes.down <- list(ETP = downregulated_patients_etp, PRC2_mut = downregulated_prc2_mut)
genes.up <- list(ETP = upregulated_patients_etp, PRC2_mut = upregulated_prc2_mut)

venn.diagram(genes.down, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#a9d6e5", "#014f86"), 
             "genes.down.venn.etp.tiff")
venn.diagram(genes.up, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#ffccd5", "#c9184a"), 
             "genes.up.venn.etp.tiff")


write.table(intersect(downregulated_patients_etp, downregulated_prc2_mut), row.names = FALSE, quote = FALSE, "downregulated_genes_etp_intersect_prc2_mut.csv")
write.table(intersect(upregulated_patients_etp, upregulated_prc2_mut), row.names = FALSE, quote = FALSE, col.names = FALSE, "upregulated_genes_etp_intersect_jurkat.csv")


upreg_genes_overlap_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                         "geneontology_Molecular_Function_noRedundant"),
                                      interestGeneType = "genesymbol",
                                      referenceGeneType = "genesymbol",
                                      interestGene = intersect(upregulated_patients_etp, upregulated_prc2_mut),
                                      referenceGene = intersect(all_etp$geneID, all_prc2_mut$SYMBOL),
                                      sigMethod = "fdr", fdrThr = 0.999)
GSEA.chosen = upreg_genes_overlap_go[1:20,]  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Overlap ETP patients and PRC2 mut patients") + #title to be changed depending on which collection is run
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
                                     interestGene = intersect(downregulated_patients_etp, downregulated_prc2_mut),
                                     referenceGene = intersect(all_etp$geneID, all_prc2_mut$SYMBOL),
                                     sigMethod = "fdr", fdrThr = 0.99)
GSEA.chosen = down_genes_overlap_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Overlap ETP patients and EZH2 KO cells") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 



#Venn DIagram with overlap for patient samples and cell lines
genes.down <- list(ETP = downregulated_patients_etp, PRC2_mut = downregulated_prc2_mut, EZH2_KO = down_ezh2_ko)
genes.up <- list(ETP = upregulated_patients_etp, PRC2_mut = upregulated_prc2_mut, EZH2_KO = up_ezh2_ko)
venn.diagram(genes.down, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5, 0.5), fill = c("#a9d6e5", "#014f86", "#08519c"), 
             "#path\\genes.down.venn.etp.prc2mut.ezh2ko.tiff")
venn.diagram(genes.up, lwd = 0, cex = 1, cat.cex = 0.5, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#c9184a", "#ff69b4"), 
             "#path\\genes.up.venn.etp.prc2mut.ezh2ko.tiff")
