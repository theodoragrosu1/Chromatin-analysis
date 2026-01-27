#input file: raw counts RNAseq (Polonen et al, 2024)
#input file: clinical data from sequenced patients (Polonen et al, 2024)
#output file: differentially expressed genes in ETP-ALL versus T-ALL patients and volcano plot
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

#do the same with PRC2 altered cases

#mutations sheet14
prc2_mutations <- read_xlsx("#path to clinical data", sheet =15) %>%
  dplyr::filter(gene == "EZH2" | gene == "SUZ12" | gene == "EED") %>%
  dplyr::select(sample, gene, mutation_class, ExonicFunc) %>%
  dplyr::rename("ID" = "sample")

ID_all <- read_xlsx("#path to clinical data", sheet =2) %>%
  dplyr::select("USI") %>%
  dplyr::rename("ID"="USI") %>%
  drop_na()

prc2_mutations_status <- left_join(prc2_mutations, ID_all, by = "ID") %>%
  drop_na()

#cna sheet 19
prc2_cna <- read_xlsx("#path to clinical data", sheet = 20) %>%
  dplyr::filter(chrom == "chr7" | chrom == "chr11" | chrom == "chr17")


#grch37
#ezh2 chr7:148504464-148581441
#eed chr11:85955806-85989785
#suz12 chr17:30264044-30328057

prc2_cna <- prc2_cna %>%
  dplyr::mutate(prc2_status = case_when(chrom == "chr7" & loc.start<148504464 & loc.end>148581441 & LogRatio < 0.7 ~ "ezh2_loss",
                                        chrom == "chr11" & loc.start<85955806 & loc.end>85989785 & LogRatio < 0.7 ~ "eed_loss",
                                        chrom == "chr17" & loc.start<30264044 & loc.end>30328057 & LogRatio < 0.7 ~ "suz12_loss")) %>%
  drop_na()


prc2_lost <- prc2_cna %>%
  dplyr::select(c(sample, prc2_status)) %>%
  dplyr::filter(prc2_status != "NA") 

prc2_lost <- prc2_lost %>%
  dplyr::rename("ID" = "sample")

prc2_lost_status <- left_join(prc2_lost, ID_all, by = "ID") %>%
  drop_na()

#read patient samples RNAseq raw counts
IDs <- read.csv("#path to raw counts", sep = '\t') %>%
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

#attribute PRC2 status to patients in the RNAseq data
rna_seq_counts_filtered_status_prc2_lost <- rna_seq_counts_filtered %>%
  dplyr::select(prc2_lost_status$ID)

rna_seq_counts_filtered_status_prc2_mut <- rna_seq_counts_filtered %>%
  dplyr::select(prc2_mutations_status$ID)

#make study design
study_design <- ID_all %>%
  dplyr::mutate(prc2_status = case_when(ID %in% unique(c(prc2_mutations_status$ID, prc2_cna$sample)) ~ "prc2_mut", TRUE ~ "prc2_wt"))


study_design <- study_design %>%
  dplyr:: filter(ID %in% colnames(rna_seq_counts_filtered))

study_design <- study_design %>%
  dplyr::filter(ID %in% colnames(rna_seq_counts_filtered))

rna_seq_counts_filtered_1 <- rna_seq_counts_filtered %>%
  dplyr::select(study_design$ID) %>%
  as.matrix

#perform differential expression analysis
DGEList <- DGEList(rna_seq_counts_filtered_1)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=217 #user defined - depends on studydesign, eliminates lowly expressed genes --> here, 216 patients have PRC2 alterations
DGEList.filtered <- DGEList[keepers,]
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>%
  cpm(log=TRUE) %>%
  as_tibble(rownames = "SYMBOL")


group.t.all <- study_design$prc2_status
group.t.all <- factor(group.t.all)
design.t.all <- model.matrix(~0 + group.t.all)
colnames(design.t.all) <- levels(group.t.all)

#v.DEGList.t.all.filtered.norm <- voom(DGEList.t.all.filtered.norm, design.t.all, plot = TRUE)
fit.t.all <- lmFit(log2.cpm.filtered.norm, design.t.all)
contrast.matrix.t.all <- makeContrasts(differences = prc2_mut - prc2_wt,
                                       levels=design.t.all)

fits.t.all <- contrasts.fit(fit.t.all, contrast.matrix.t.all)
ebFit.t.all <- eBayes(fits.t.all)
myTopHits.t.all <- topTable(ebFit.t.all, adjust ="BH", coef=1, number=40000, sort.by="p")

myTopHits.t.all.df <- myTopHits.t.all
rownames(myTopHits.t.all.df) <- NULL



hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      collection = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol)

data.df <- myTopHits.t.all %>%
  dplyr::filter(abs(logFC) > 0.05) 

data.df.sub <- dplyr::select(data.df, SYMBOL, logFC)
data.gsea <- data.df.sub$logFC
names(data.gsea) <- as.character(data.df.sub$SYMBOL)
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
  labs(y = "ID", title = "PRC2 mutant - PRC2 wild-type paediatric T-ALL") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 


gseaplot2(GSEA.res, 
          geneSetID = c(4, 12), #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = FALSE,
          title = "PRC2 mutant versus PRC2 wild-type paediatric T-ALL")




#volcano plot DEGs
upregulated_in_mut <- myTopHits.t.all.df %>%
  dplyr::filter(logFC > 1 & adj.P.Val<= 0.05)

upregulated_in_wt <- myTopHits.t.all.df %>%
  dplyr::filter(logFC < -1 & adj.P.Val<= 0.05)

vplot <- ggplot(myTopHits.t.all.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", SYMBOL)) +
  geom_point(size=2) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 0.5, xmax = 5, ymin = -log10(10), ymax = 55, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -0.5, xmax = -5, ymin = -log10(10), ymax = 55, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -4, y = 5, label = "1245 in PRC2 wt", size = 5, colour = "#2C467A") + #labels are manually added after running summary(results.t.all downstream)
  annotate(geom = "text", x = 4, y = 5, label = "1026 in PRC2 mut", size = 5, colour = "#BE684D") +
  labs(title="Volcano plot",
       subtitle = "T ALL prc2 wt & mut",
       caption="adj.P.Val > 0.1\n|logFC|>1" ) +
  theme_bw(base_size = 20) +
  geom_label_repel(data=dplyr::filter(myTopHits.t.all.df, abs(-log10(adj.P.Val))>0.1 & abs(logFC) > 2.2), aes(label=SYMBOL), max.overlaps = 45) +
  geom_label_repel(data=dplyr::filter(myTopHits.t.all.df, SYMBOL %in% c("MEF2C", "CCND2", "LYL1", "CD1A")), aes(label=SYMBOL), max.overlaps = 45)

vplot

results.t.all <- decideTests(ebFit.t.all, method="global", adjust.method="BH", p.value=0.1, lfc=0.5)
summary(results.t.all) #results are added manually to vplot upstream


#check expression of individual genes
EZH1_t_all_rna_seq <- rna_seq_counts %>%
  dplyr::filter(SYMBOL == "EZH1") %>%
  dplyr::select(-SYMBOL) %>%
  pivot_longer(cols = everything(), values_to = "RNA_expression", names_to = "patient_id") %>%
  mutate(group = case_when(patient_id %in% pull(study_design[study_design[, "prc2_status"] == "prc2_mut", "ID"], "ID") ~"prc2_mut",
                           patient_id %in% pull(study_design[study_design[, "prc2_status"] == "prc2_wt", "ID"], "ID") ~ "prc2_wt")) %>%
  drop_na()

EZH1_t_all_rna_seq %>%
  ggplot(aes(x = group, y = RNA_expression)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  theme_bw(base_size = 20) +
  labs(title = "EZH1 expression")

study_design %>%
  dplyr::filter(prc2_status == "prc2_mut") %>%
  nrow()
