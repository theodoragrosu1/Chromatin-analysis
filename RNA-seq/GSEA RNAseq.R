#This script performs differential gene expression analysis and GSEA on Jurkat EZH2-WT vs Jurkat-EZH2-KO
#Input: kallisto abundance for each sample - available on GEO
#Input: study design file available with the code
#Output: volcano plot, GSEA plots
#Output: differentially expressed genes comparing WT with both A7 and C9 EZH2-KO clones


library(tidyverse)
library(matrixStats)
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot)
library(GSEABase)
library(limma)
library(ggplot2)


normalised_expression_jurkat_vs_clones <- read.csv("#path to csv containing normalised counts", header = TRUE)

data.df <- normalised_expression_jurkat_vs_clones %>%
  dplyr::filter(abs(logFC) > 0.5)



# GSEA for C2----
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      collection = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 



#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(data.df, geneID, logFC)
data.gsea <- data.df.sub$logFC
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)
print(data.gsea)
#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.05, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)

write.csv(GSEA.df, file = "#path and name of file to be saved") #change name depending on which collection is run

#plot for individual signature
gseaplot2(GSEA.res, 
          geneSetID = c(18), #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = FALSE,
          title = "PRC2 mutant versus PRC2 wild-type: cell lines")

#plots significant GSEA
GSEA.chosen = GSEA.df[1:30,] #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
  
  ggplot(GSEA.chosen) + aes(x=NES, y= fct_reorder(ID, NES), fill = NES) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "PRC2 mutant - PRC2 wild-type cell lines") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 
