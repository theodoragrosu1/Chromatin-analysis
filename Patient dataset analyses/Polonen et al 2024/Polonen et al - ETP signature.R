#input file: raw counts RNAseq (Polonen et al, 2024)
#input file: clinical data from sequenced patients (Polonen et al, 2024)
#output file: differentially expressed genes in ETP-ALL versus T-ALL patients
#output file: signature of genes typical for ETP-ALL
#output file: GSEA on the ETP-ALL vs T-ALL genes

library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase)
library(Biobase)
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(gplots)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(tximport)
library(ggrepel)

#read patient samples RNAseq raw counts
IDs <- read.csv("#path to csv file", sep = '\t') %>%
  dplyr::rename("geneID" = "X")

IDs$geneID <- gsub("\\..*","",IDs$geneID)


annots <- select(org.Hs.eg.db, keys=IDs$geneID, 
                 columns="SYMBOL", keytype="ENSEMBL")

annots %>%
  drop_na()

annots <- annots %>%
  dplyr::rename("geneID" = "ENSEMBL")

rna_seq_counts <- left_join(IDs, annots, by = "geneID")

rna_seq_counts<- rna_seq_counts[,-1]


#R crashes when grouping and summarizing across 6000 genes, so simplify by only summarizing genes that appear twice
summarized_data <- rna_seq_counts %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

summarized_data %>%
  drop_na()
single_row_groups <- rna_seq_counts %>%
  group_by(SYMBOL) %>%
  dplyr::filter(n() == 1) %>%
  ungroup()

rna_seq_counts_filtered <- bind_rows(summarized_data, single_row_groups) %>%
  dplyr::filter(duplicated(SYMBOL))
rownames(rna_seq_counts_filtered) <- NULL
rownames(rna_seq_counts_filtered) <- rna_seq_counts_filtered$'SYMBOL'
rna_seq_counts_filtered <- rna_seq_counts_filtered[,-1] %>%
  as.data.frame()

#ETP status must be added to study design 

#read patient ID xlsx
status <- read_xlsx("#path to clinical data from xlsx", sheet =2) %>%
  dplyr::select(c("USI", "ETP.STATUS")) %>%
  dplyr::rename("ID"="USI") %>%
  drop_na()

#for the analysis, select only ETP and Non-ETP
status <- status %>%
  dplyr::filter(ETP.STATUS!="Unknown") %>%
  dplyr::filter(ETP.STATUS!="Near-ETP")


#remove patients from the RNA seq counts file which do not have their ETP status annotated

rna_seq_counts_filtered_status <- rna_seq_counts_filtered %>%
  dplyr::select(status$ID)

#perform differential expression analysis
DGEList <- DGEList(rna_seq_counts_filtered_status)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=110 #user defined - depends on studydesign, eliminates lowly expressed genes --> here, 110 patients are ETP
DGEList.filtered <- DGEList[keepers,]
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>%
  cpm(log=TRUE) %>%
  as_tibble(rownames = "geneID")

conditions <- status$ETP.STATUS %>%
  factor()
design <- model.matrix(~0 + conditions)
colnames(design) <- c("ETP", "Non_ETP")

v.DEGList.patients.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.DEGList.patients.filtered.norm, design)
contrast.matrix <- makeContrasts(differences = Non_ETP - ETP, levels=design)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust = "BH", coef=1, number=40000, sort.by="p")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 5, ymin = -log10(0.1), ymax = 155, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -5, ymin = -log10(0.1), ymax = 155, alpha=.2, fill="#2C467A") +
  #annotate(geom = "text", x = -8, y = 2, label = "649", size = 7, colour = "#2C467A") + #labels are manually added after running summary(results.t.all downstream)
  #annotate(geom = "text", x = 8, y = 2, label = "582", size = 7, colour = "#BE684D") +
  labs(title="Volcano plot",
       subtitle = "T ALL vs ETP-ALL", #change
       caption="adj.p.value >0.05\n|logFC|≥1" ) + #change
  theme_bw(base_size = 20) +
  geom_label_repel(data=dplyr::filter(myTopHits.df, abs(-log10(adj.P.Val))>0.01 & abs(logFC) > 1), aes(label=geneID), max.overlaps = 30) +
  geom_label_repel(data=dplyr::filter(myTopHits.df, geneID %in% c("MEF2C", "CCND2", "LYL1", "SMARCC2")), aes(label=geneID), max.overlaps = 30)

vplot

signature_etp <- myTopHits.df %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(logFC >= 1 | logFC <= -1) %>%
  mutate(condition = case_when(logFC >= 1 ~ "Genes_upregulated_in_T-ALL",
                               logFC <= -1 ~ "Genes_upregulated_in_ETP")) %>%
  dplyr::select(condition, geneID, logFC)

data.df <- read.csv("/Users/theo/Desktop/Jurkat_wt_v_clones_all_genes (2).csv", header = TRUE)
data.df.sub <- dplyr::select(data.df, geneID, logFC)
data.gsea <- data.df.sub$logFC
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=signature_etp, verbose=FALSE, maxGSSize = 1000, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 1, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)

gseaplot2(GSEA.res, geneSetID = c(1), pvalue_table = TRUE)

write.csv(signature_etp, file = "#path to where the csv should be saved", quote = FALSE, row.names = FALSE)


#run GSEA on ETP signature
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      collection = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol)

data.df <- signature_etp

data.df.sub <- dplyr::select(data.df, geneID, logFC)
data.gsea <- data.df.sub$logFC
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)
print(data.gsea)
#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.05, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)


#plots significant GSEA
GSEA.chosen = GSEA.df[1:20,] #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

ggplot(GSEA.chosen) + aes(x=NES, y= fct_reorder(ID, NES), fill = NES) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "Non-ETP - ETP paediatric T-ALL") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 


gseaplot2(GSEA.res, 
          geneSetID = c(2, 33), #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = FALSE,
          title = "PRC2 mutant versus PRC2 wild-type paediatric T-ALL")
