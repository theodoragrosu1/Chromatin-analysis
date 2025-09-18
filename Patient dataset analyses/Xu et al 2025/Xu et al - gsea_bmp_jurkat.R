#input files: list of genes expressed in BMP versus T-specified from paper
#input files: list of normalized expression from cell line RNAseq
#output: GSEA with BMP# signature as backgorund

###load packages ----
library(tidyverse)
library(matrixStats)
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot)
library(GSEABase)
library(limma)
library(ggplot2)
library(readxl)


normalised_expression_jurkat_vs_clones <- read.csv("#path to normalized expression from cell line RNAseq", header = TRUE)

data.df <- normalised_expression_jurkat_vs_clones 

signature_bmp <- read_xlsx("#path to csv from paper", sheet = 8) %>%
  dplyr::select(comparison, geneID, avg_log2FC) %>%
  dplyr::mutate(comparison = case_when(avg_log2FC < -0.1 ~ 'upregulated_in_BMP',
                                       avg_log2FC > 0.1 ~ 'upregulated_in_T_spec')) 

write.csv(signature_bmp, file = "#save BMP signature file")

signature_bmp_anno <- signature_bmp %>%
  dplyr::select(comparison, geneID)

#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(data.df, geneID, logFC) 
data.gsea <- data.df.sub$logFC 
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)

#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=signature_bmp_anno, verbose=FALSE, maxGSSize = 1000, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.999, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)

#plot for individual signature
gseaplot2(GSEA.res, 
          geneSetID = c(1), #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = TRUE)#, 
title = GSEA.res$Description[1])

#plots significant GSEA
GSEA.chosen = GSEA.df# [1:30, ]  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=NES, y= fct_reorder(ID, NES), fill = NES) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "C2 - WT vs clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

