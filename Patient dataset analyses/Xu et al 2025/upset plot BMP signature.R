#input file: genes that are differentially expressed in the different conditions 
#output: upset plots for both up- and downregulated genes across conditions

library(UpSetR)
library(dplyr)

set_genes_up_when_ezh2_lost <- read.csv("#DEG RNAseq cell lines (both clones together") %>%
  dplyr::filter(logFC >= 1) %>%
  dplyr::filter(adj.P.Val >= 0.05) %>%
  dplyr::select("geneID")

set_genes_up_when_prc2_lost_polonen <- read.csv("#PRC2 altered DEGs in patient dataset from Polonen et al") %>%
  dplyr::filter(adj.P.Val < 0.1) %>% 
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::select("SYMBOL") %>%
  mutate(gs_name = "Genes_up_upon_PRC2_alt_Polonen") %>%
  dplyr::select(c("gs_name", "SYMBOL"))

upregulated_patients_etp <- read.csv("#DEGs in ETPs from Polonen et al") %>%
  dplyr::filter(condition == "Genes_upregulated_in_ETP") %>%
  dplyr::select("geneID")

set_genes_up_when_bmp <- read.csv("#genes in BMP Xu et al set")

UpSetR::upset(fromList(list("EZH2 KO up" = set_genes_up_when_ezh2_lost$geneID, "PRC2 altered up" = set_genes_up_when_prc2_lost_polonen$SYMBOL, 
                            "BMP up" = set_genes_up_when_bmp$geneID)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("EZH2 KO up", "PRC2 altered up", "BMP up"), 
                                  color = "#E68D74", active = T)))
overlap <- intersect(intersect(intersect(set_genes_up_when_ezh2_lost$geneID, set_genes_up_when_prc2_lost_polonen$SYMBOL), upregulated_patients_etp$geneID), set_genes_up_when_bmp$geneID)
overlap

overlap <- intersect(intersect(set_genes_up_when_ezh2_lost$geneID, set_genes_up_when_prc2_lost_polonen$SYMBOL), set_genes_up_when_bmp$geneID)

UpSetR::upset(fromList(list("EZH2 KO up" = set_genes_up_when_ezh2_lost$geneID, "PRC2 altered up" = set_genes_up_when_prc2_lost_polonen$SYMBOL, 
                            "ETP up" = upregulated_patients_etp$geneID)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("EZH2 KO up", "ETP up", "PRC2 altered up"), 
                                  color = "#E68D74", active = T),
                             list(query = intersects, 
                                  params = list("EZH2 KO up", "PRC2 altered up", "ETP up"), 
                                  color = "#E68D74", active = T)))

#same analysis as above but for downregulated genes
set_genes_down_when_ezh2_lost <- read.csv("#DEG RNAseq cell lines (both clones together") %>%
  dplyr::filter(logFC <= -1) %>%
  dplyr::filter(adj.P.Val >= 0.05) %>%
  dplyr::select("geneID")

set_genes_down_when_prc2_lost_polonen <- read.csv("#PRC2 altered DEGs in patient dataset from Polonen et al") %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select("SYMBOL")

set_genes_down_when_etp <-  read.csv("#DEGs in ETPs from Polonen et al") %>%
  dplyr::filter(condition == "Genes_upregulated_in_T-ALL") %>%
  dplyr::select("geneID")

set_genes_down_when_bmp <- read.csv("#genes in T-specified set in Xu et al") 

UpSetR::upset(fromList(list("EZH2 KO down" = set_genes_down_when_ezh2_lost$geneID, "PRC2 altered down" = set_genes_down_when_prc2_lost_polonen$SYMBOL, 
                            "ETP down" = set_genes_down_when_etp$geneID)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("EZH2 KO down", "ETP down", "PRC2 altered down"), 
                                  color = "#266CA9", active = T),
                             list(query = intersects, 
                                  params = list("EZH2 KO down", "PRC2 altered down", "ETP down"), 
                                  color = "#266CA9", active = T)))
overlap_down <- intersect(intersect(intersect(set_genes_down_when_ezh2_lost$geneID, set_genes_down_when_prc2_lost_liu$geneID), set_genes_down_when_etp_liu_et_al$geneID), set_genes_down_when_bmp$geneID)
overlap_down


UpSetR::upset(fromList(list("EZH2 KO" = set_genes_down_when_ezh2_lost$geneID, "PRC2 altered" = set_genes_down_when_prc2_lost_polonen$SYMBOL, 
                            "BMP" = set_genes_down_when_bmp$geneID)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("EZH2 KO", "PRC2 altered", "BMP"), 
                                  color = "#266CA9", active = T)))

overlap_down <- intersect(intersect(set_genes_down_when_ezh2_lost$geneID, set_genes_down_when_prc2_lost_polonen$SYMBOL), set_genes_down_when_bmp$geneID)
