library(tidyverse)
library(readxl)
library(Cairo)
library(GSEABase) 
library(clusterProfiler)
library(enrichplot)
library(limma)
library(stats)

clinical <- read_xlsx(".\\41591_2022_2112_MOESM3_ESM.xlsx", sheet = 1, skip = 2) #data from supplemental of Lee et al
drug_screen <- read_xlsx(".\\41591_2022_2112_MOESM3_ESM.xlsx", sheet = 2, skip = 2)

drug_screen_t_all <- read_xlsx(".\\41591_2022_2112_MOESM3_ESM.xlsx", sheet = 1, skip = 2) %>%
  dplyr::filter(Immunophenotype == "T") %>%
  drop_na(Asparaginase_normalized) %>%
  dplyr::select(c("Patient ID", "Molecular subtype", "Asparaginase_normalized")) %>%
  left_join(drug_screen) %>%
  dplyr::select(c("Patient ID", "Molecular subtype", "Asparaginase_normalized")) %>%
  drop_na()

#downloaded from https://permalinks.stjude.cloud/permalinks/all-pharmacotype
Lee_et_al_fpkm <- read.csv(".\\pharmacotyping_ped_rnaseq_fpkm_ALLids_0823.csv") %>%
  dplyr::select(-GeneID)

#keep only samples present in the drug screen data
Lee_et_al_fpkm <- Lee_et_al_fpkm %>%
  dplyr::select(any_of(c("GeneName", drug_screen_t_all$`Patient ID`)))

#fitler drug screen data using rna-seq samples
drug_screen_t_all <- drug_screen_t_all %>%
  dplyr::filter(`Patient ID` %in% colnames(Lee_et_al_fpkm)) %>%
  as.data.frame()

ggplot(drug_screen_t_all, aes(x = Asparaginase_normalized)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey80")+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_minimal(base_size = 20)

histogram_asparaginase_treatment <- ggplot(drug_screen_t_all, aes(x = Asparaginase_normalized)) + 
  geom_histogram(colour="black", fill="grey80")+ 
  #geom_vline(xintercept = 0.25, linetype = "dashed", color = "#FF6666", size = 1) +
  annotate("rect", xmin = 0.75, xmax = 1.02, ymin = 0, ymax = 10, alpha=.2, fill="#FF6666") +
  annotate("rect", xmin = -0.00, xmax = 0.30, ymin = 0, ymax = 10, alpha=.2, fill="#007EA7") +
  annotate(geom = "text", x = 0.85, y = 9, label = "resistant", size = 5, colour = "#FF6666") +
  annotate(geom = "text", x = 0.15, y = 9, label = "sensitive", size = 5, colour = "#007EA7") +
  #geom_vline(xintercept = 0.75, linetype = "dashed", color = "#FF6666", size = 1) +
  theme_minimal(base_size = 20) + 
  labs(x = "Normalized Asparaginase sensitivity", y = "No. samples", title = "Lee et al. Asparaginase treatment")

ggsave("./asp_res_vs_sens_gsea_signatures_from_cell_lines/histogram_asp_treatment.svg", 
       plot = histogram_asp_treatment, width = 16, height = 10, 
       dpi = 1000, units = "cm") 


mean_asp <- mean(drug_screen_t_all$Asparaginase_normalized)

#stratify
asp_resistant <- drug_screen_t_all[drug_screen_t_all$Asparaginase_normalized > 0.75, "Patient ID"]
asp_sensitive <- drug_screen_t_all[drug_screen_t_all$Asparaginase_normalized < 0.25, "Patient ID"]


#build study design matrix
study_design <- drug_screen_t_all %>% 
  dplyr::mutate(Asparaginase_resistance = case_when(`Patient ID` %in% asp_resistant ~ "Resistant",
                                           `Patient ID` %in% asp_sensitive ~ "Sensitive")) %>%
  dplyr::select(c("Patient ID", "Asparaginase_resistance")) %>%
  na.omit()


#build groups resistant vs sensitive
group <- study_design$Asparaginase_resistance
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#filter lowly expressed genes and log2 normalise

rownames(Lee_et_al_fpkm) <- Lee_et_al_fpkm$GeneName
Lee_et_al_fpkm <- Lee_et_al_fpkm %>%
  dplyr::select(-GeneName)

#select only the ones present in study design
Lee_et_al_fpkm <- Lee_et_al_fpkm %>%
  dplyr::select(study_design$`Patient ID`)

keepers <- rowSums(Lee_et_al_fpkm>0.2)>= 4 # user defined - depends on studydesign #27 chosen here because it's the smallest group
Lee_et_al_fpkm_filtered <- Lee_et_al_fpkm[keepers, ]
log_Lee_et_al_fpkm <- log2(Lee_et_al_fpkm_filtered+0.1)

#fit model
fit <- lmFit(log_Lee_et_al_fpkm, design)
contrast.matrix <- makeContrasts(differences = Resistant - Sensitive,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits, trend = TRUE)
top_hits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="p")

top_hits$gene <- rownames(top_hits)
top_hits_df <- top_hits %>%
  as_tibble()

Asparaginase_signature <- top_hits_df %>%
  dplyr::mutate(comparison = case_when(logFC < -0.5 ~ 'Asparaginase_sensitive',
                                       logFC > 0.5 ~ 'Asparaginase_resistant')) %>%
  dplyr::select("comparison", "gene")
  
#create gene sets from existing RNA-Seq on PRC2 WT vs PRC2-depleted
#A. from cell lines
#B. from patient data

set_genes_up_when_prc2_lost <- read.csv(".\\Jurkat_wt_v_clones_all_genes (2).csv") %>%
  dplyr::filter(adj.P.Val < 0.1) %>% 
  dplyr::filter(logFC >= 0.5) %>%
  mutate(gs_name = "Genes_up_upon_EZH2_KO_Jurkat") %>%
  dplyr::select(c("gs_name", "geneID", "logFC"))

set_genes_down_when_prc2_lost <- read.csv(".\\Jurkat_wt_v_clones_all_genes (2).csv") %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::filter(logFC <= -0.5) %>%
  mutate(gs_name = "Genes_down_upon_EZH2_KO_Jurkat") %>%
  dplyr::select(c("gs_name", "geneID", "logFC"))

set_genes_from_cell_lines <- rbind(set_genes_up_when_prc2_lost, set_genes_down_when_prc2_lost)
# the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis

mydata.t.all.gsea <- set_genes_from_cell_lines$logFC
names(mydata.t.all.gsea) <- as.character(set_genes_from_cell_lines$geneID)
mydata.t.all.gsea <- sort(mydata.t.all.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
###Run on cell line gene sets
set.seed(1234)
myGSEA.t.all.res <- GSEA(mydata.t.all.gsea, TERM2GENE=Asparaginase_signature, verbose=FALSE, seed = TRUE, nPermSimple = 1000, eps = 0, pvalueCutoff = 0.99)
myGSEA.t.all.df <- as_tibble(myGSEA.t.all.res@result)

dir.create("asp_res_vs_sens_gsea_signatures_from_cell_lines")
write.csv(myGSEA.t.all.df, file = ".\\GSEA_df_Asparaginase_Jurkat_signatures.csv", quote = FALSE, row.names = FALSE)

#save figures
#change svg and .svg with pdf and .pdf, respectively to change format of output figure

pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[1], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 1, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = FALSE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[1]))



pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[2], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 2, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = FALSE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[2]))
dev.off()

gseaplot2(myGSEA.t.all.res, 
          geneSetID = c(1,2), #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases"))

#ETP signature from Polonen et al
set_genes_up_when_etp <- read.csv(".\\Polonen et al RNAseq analysis\\rscripts\\t_all_etp_v_non-etp_signature.csv") %>%
  dplyr::filter(logFC >= 1.5) %>% 
  dplyr::select(geneID) %>%
  mutate(gs_name = "Genes_up_in_ETP") %>%
  dplyr::select(c("gs_name", "geneID"))

set_genes_down_when_etp <- read.csv(".\\Polonen et al RNAseq analysis\\rscripts\\t_all_etp_v_non-etp_signature.csv") %>%
  dplyr::filter(logFC <= -1.5) %>%
  dplyr::select(geneID) %>%
  mutate(gs_name = "Genes_down_in_ETP") %>%
  dplyr::select(c("gs_name", "geneID"))


set_genes_in_etp <- rbind(set_genes_up_when_etp, set_genes_down_when_etp)
colnames(set_genes_in_etp) <- c("gs_name", "gene")
  
# the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.t.all.gsea <- top_hits$logFC
names(mydata.t.all.gsea) <- as.character(rownames(top_hits))
mydata.t.all.gsea <- sort(mydata.t.all.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
myGSEA.t.all.res <- GSEA(mydata.t.all.gsea, TERM2GENE=set_genes_in_etp, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 1)
myGSEA.t.all.df <- as_tibble(myGSEA.t.all.res@result)

write.csv(myGSEA.t.all.df, file = "./asp_res_vs_sens_gsea_signatures_from_cell_lines/GSEA_df_Liu_et_al_signatures.csv", quote = FALSE, row.names = FALSE)

#save figures
#change svg and .svg with pdf and .pdf, respectively to change format of output figure
pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[1], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 1, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[1]))
dev.off()


pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[2], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 2, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[2]))
dev.off()

gseaplot2(myGSEA.t.all.res, 
          geneSetID = c(1,2), #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases"))

#PRc2mut signature from Polonen et al
set_genes_up_when_prc2mut <- read.csv(".\\Polonen et al RNAseq analysis\\rscripts\\t_all_prc2_wt_v_mut.csv") %>%
  dplyr::filter(adj.P.Val < 0.1) %>% 
  dplyr::filter(logFC >= 1) %>%
  dplyr::select(SYMBOL) %>%
  mutate(gs_name = "Genes_up_upon_PRC2_mut") %>%
  dplyr::select(c("gs_name", "SYMBOL"))

set_genes_down_when_prc2mut <- read.csv(".\\Polonen et al RNAseq analysis\\rscripts\\t_all_prc2_wt_v_mut.csv") %>%
  dplyr::filter(adj.P.Val < 0.1) %>% 
  dplyr::filter(logFC <= -1) %>%
  dplyr::select(SYMBOL) %>%
  mutate(gs_name = "Genes_down_upon_PRC2_mut") %>%
  dplyr::select(c("gs_name", "SYMBOL"))


set_genes_in_prc2_mut <- rbind(set_genes_up_when_prc2mut, set_genes_down_when_prc2mut)
colnames(set_genes_in_prc2_mut) <- c("gs_name", "gene")

# the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.t.all.gsea <- top_hits$logFC
names(mydata.t.all.gsea) <- as.character(rownames(top_hits))
mydata.t.all.gsea <- sort(mydata.t.all.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
myGSEA.t.all.res <- GSEA(mydata.t.all.gsea, TERM2GENE=set_genes_in_prc2_mut, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 1)
myGSEA.t.all.df <- as_tibble(myGSEA.t.all.res@result)

write.csv(myGSEA.t.all.df, file = "./asp_res_vs_sens_gsea_signatures_from_cell_lines/GSEA_df_Liu_et_al_signatures.csv", quote = FALSE, row.names = FALSE)

#save figures
#change svg and .svg with pdf and .pdf, respectively to change format of output figure
pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[1], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 1, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[1]))
dev.off()


pdf(paste0("./asp_res_vs_sens_gsea_signatures_from_cell_lines/", myGSEA.t.all.res$Description[2], ".pdf"), width = 8, height = 5)
gseaplot2(myGSEA.t.all.res, 
          geneSetID = 2, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = TRUE, 
          title = paste("Asparaginase resistant vs sensitive cases -", myGSEA.t.all.res$Description[2]))
dev.off()
