#input file: normalized expression for cell line RNAseq
#input file: stemness-related genes
#output: signature 


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

normalised_expression_jurkat_vs_clones <- read.csv("#path to normalized cell lines RNAseq", header = TRUE)

data.df <- normalised_expression_jurkat_vs_clones 


signature_stemness <- read.csv("#path to paper csv", header = TRUE) %>%
  dplyr::select(comparison, geneID)

#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(data.df, geneID, logFC) 
data.gsea <- data.df.sub$logFC 
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)

#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=signature_stemness, verbose=FALSE, maxGSSize = 1000, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.999, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)

#plot for individual signature
gseaplot2(GSEA.res, 
          geneSetID = c(1), #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = TRUE)#, 
title = GSEA.res$Description[1])

