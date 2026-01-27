library(edgeR)
library(readr)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(glmnet)
library(ggplot2)
library(tibble)
library(rsample)
library(matrixStats)
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot)
library(GSEABase)


counts <- read_tsv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\TALL_X01_counts.tsv")
counts <- as.data.frame(counts)

rownames(counts) <- counts[[1]]
counts <- counts[, -1]

clinical_data <- readxl::read_excel("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\41586_2024_7807_MOESM4_ESM.xlsx", sheet = "ST1_Clinical_Data")
metadata <- clinical_data %>%
  dplyr::select(1, `ETP.STATUS`) %>%
  dplyr::rename(sampleID = 1, condition = `ETP.STATUS`) %>%
  dplyr::filter(condition %in% c("ETP", "Non-ETP"))

metadata$condition <- factor(metadata$condition, levels = c("Non-ETP", "ETP"))
counts <- counts[, metadata$sampleID]

geneIDs <- gsub("\\..*", "", rownames(counts))
annots <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = geneIDs,
                                columns = "SYMBOL",
                                keytype = "ENSEMBL")

annots <- annots %>%
  tidyr::drop_na() %>%
  dplyr::rename(geneID = ENSEMBL)

glimpse(annots)

counts$geneID <- geneIDs
counts <- counts %>%
  left_join(annots, by = "geneID") %>%
  filter(!is.na(SYMBOL))


counts <- counts %>%
  tibble::as_tibble() %>%
  dplyr::select(-geneID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), sum))

counts <- as.data.frame(counts)
rownames(counts) <- counts$SYMBOL
counts <- counts[, -1]


#perform differential expression analysis
DGEList <- DGEList(counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=110 #user defined - depends on studydesign, eliminates lowly expressed genes --> here, 110 patients are ETP
DGEList.filtered <- DGEList[keepers,]
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>%
  cpm(log=TRUE) %>%
  as.data.frame(rownames = "geneID")

X <- t(log2.cpm.filtered.norm)        # samples x genes

mrd <- clinical_data %>%
  filter(`ETP.STATUS` %in% c("ETP", "Non-ETP")) %>%
  dplyr::select(USI, `Day.29.MRD`)

Y <- mrd$`Day.29.MRD`[match(metadata$sampleID, mrd$USI)]

MRD_binary <- as.integer(Y > 0.01)  

table(MRD_binary)

valid_idx <- !is.na(MRD_binary)
X_filtered <- X[valid_idx, ]
Y_filtered <- MRD_binary[valid_idx]

dim(X_filtered)
length(Y_filtered)

set.seed(42) # saw this as a standard value to choose in data analysis for reproducibility

data_combined <- as.data.frame(X_filtered)
data_combined$MRD_binary <- Y_filtered

split_obj <- initial_split(data_combined, prop = 0.8, strata = "MRD_binary")
train_data <- training(split_obj)
test_data <- testing(split_obj)

X_train <- as.matrix(train_data[, -ncol(train_data)])
Y_train <- train_data$MRD_binary

X_test <- as.matrix(test_data[, -ncol(test_data)])
Y_test <- test_data$MRD_binary

cvfit_logistic <- cv.glmnet(X_train,Y_train,family = "binomial",alpha = 1)

plot(cvfit_logistic)
title("Logistic LASSO Curve\nPredicting MRD > 0.01", line = 2.5)

lasso_coef <- coef(cvfit_logistic, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_values <- lasso_coef[which(lasso_coef != 0)]

print(lasso_genes)

selected_genes_df <- data.frame(
  Gene = lasso_genes,
  Coefficient = as.numeric(lasso_values)
)

write_csv(selected_genes_df, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\LASSO_Selected_Genes_MRD_NatureThreshold.csv")

# Predict on unseen test data
pred_prob <- predict(cvfit_logistic, newx = X_test, s = "lambda.min", type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

print(table(Predicted = pred_class, Actual = Y_test))

# AUROC on test data
library(pROC)

roc_obj <- pROC::roc(Y_test, as.numeric(pred_prob))
print(roc_obj$auc)

plot(roc_obj, print.auc = TRUE, col = "blue",
     main = "AUROC: Logistic LASSO Predicting MRD > 0.01 (Test Data)")

ci.auc(roc_obj)

#run GSEAs on selected genes
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      collection = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol)

#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(selected_genes_df, Gene, Coefficient)
data.gsea <- data.df.sub$Coefficient
names(data.gsea) <- as.character(data.df.sub$Gene)
data.gsea <- sort(data.gsea, decreasing = TRUE)
print(data.gsea)
#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.99, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)
GSEA.chosen = GSEA.df[] #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

library(tidyverse)
ggplot(GSEA.chosen) + aes(x=NES, y= fct_reorder(ID, NES), fill = NES) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "Selected genes by Lasso regression") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

#perform GO enrichment analysis
reference <- as.data.frame(log2.cpm.filtered.norm) %>%
  dplyr::filter(keep == "TRUE")

library(WebGestaltR)
go_selected_genes <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                "geneontology_Molecular_Function_noRedundant"),
                             interestGeneType = "genesymbol",
                             referenceGeneType = "genesymbol",
                             interestGene = selected_genes_df$Gene,
                             referenceGene = rownames(reference),
                             sigMethod = "fdr", fdrThr = 0.999)


#do the same analysis for PRC2 mutated patients
#mutations sheet14
ezh2_mutations <- read_xlsx("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\41586_2024_7807_MOESM4_ESM.xlsx", sheet =15) %>%
  dplyr::filter(gene == "EZH2") %>%
  dplyr::select(sample, gene, mutation_class, ExonicFunc) %>%
  dplyr::rename("ID" = "sample")

ID_all <- read_xlsx("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\41586_2024_7807_MOESM4_ESM.xlsx", sheet =2) %>%
  dplyr::select("USI") %>%
  dplyr::rename("ID"="USI") %>%
  drop_na()

ezh2_mutations_status <- left_join(ezh2_mutations, ID_all, by = "ID") %>%
  drop_na()

#cna sheet 19
ezh2_cna <- read_xlsx("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\41586_2024_7807_MOESM4_ESM.xlsx", sheet = 20) %>%
  dplyr::filter(chrom == "chr7")


#grch37
#ezh2 chr7:148504464-148581441
#eed chr11:85955806-85989785
#suz12 chr17:30264044-30328057

ezh2_cna <- ezh2_cna %>%
  dplyr::mutate(ezh2_status = case_when(chrom == "chr7" & loc.start<148504464 & loc.end>148581441 & LogRatio < 0.7 ~ "ezh2_loss")) %>%
  drop_na()


ezh2_lost <- ezh2_cna %>%
  dplyr::select(c(sample, ezh2_status)) %>%
  dplyr::filter(ezh2_status != "NA") 

ezh2_lost <- ezh2_lost %>%
  dplyr::rename("ID" = "sample")

ezh2_lost_status <- left_join(ezh2_lost, ID_all, by = "ID") %>%
  drop_na()


study_design <- ID_all %>%
  dplyr::mutate(ezh2_status = case_when(ID %in% unique(c(ezh2_mutations_status$ID, ezh2_cna$sample)) ~ "ezh2_mut", TRUE ~ "ezh2_wt"))

group.t.all <- study_design$ezh2_status
group.t.all <- factor(group.t.all)
design.t.all <- model.matrix(~0 + group.t.all)
colnames(design.t.all) <- levels(group.t.all)

counts <- read_tsv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\Polonen et al RNAseq analysis\\TALL_X01_counts.tsv")
counts <- as.data.frame(counts)

rownames(counts) <- counts[[1]]
counts <- counts[, -1]

metadata <- study_design %>%
  dplyr::rename(condition = 'ezh2_status')

metadata$condition <- factor(metadata$condition, levels = c("ezh2_wt", "ezh2_mut"))
counts <- counts[, metadata$ID]

geneIDs <- gsub("\\..*", "", rownames(counts))
annots <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = geneIDs,
                                columns = "SYMBOL",
                                keytype = "ENSEMBL")

annots <- annots %>%
  tidyr::drop_na() %>%
  dplyr::rename(geneID = ENSEMBL)

glimpse(annots)

counts$geneID <- geneIDs
counts <- counts %>%
  left_join(annots, by = "geneID") %>%
  dplyr::filter(!is.na(SYMBOL))


counts <- counts %>%
  tibble::as_tibble() %>%
  dplyr::select(-geneID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), sum))

counts <- as.data.frame(counts)
rownames(counts) <- counts$SYMBOL
counts <- counts[, -1]


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = as.data.frame(metadata),
                              design = ~ 1)

keep <- rowSums(counts(dds) >= 10) >= 101 #101 PRC2-altered patients
dds <- dds[keep, ]

vsd <- vst(dds)
vst_matrix <- assay(vsd)  # genes x samples
X <- t(vst_matrix)        # samples x genes

mrd <- clinical_data %>%
  dplyr::select(USI, `Day.29.MRD`)

Y <- mrd$`Day.29.MRD`[match(metadata$ID, mrd$USI)]

MRD_binary <- as.integer(Y > 0.01)  

table(MRD_binary)

valid_idx <- !is.na(MRD_binary)
X_filtered <- X[valid_idx, ]
Y_filtered <- MRD_binary[valid_idx]

dim(X_filtered)
length(Y_filtered)

set.seed(42) # saw this as a standard value to choose in data analysis for reproducibility

data_combined <- as.data.frame(X_filtered)
data_combined$MRD_binary <- Y_filtered

split_obj <- initial_split(data_combined, prop = 0.8, strata = "MRD_binary")
train_data <- training(split_obj)
test_data <- testing(split_obj)

X_train <- as.matrix(train_data[, -ncol(train_data)])
Y_train <- train_data$MRD_binary

X_test <- as.matrix(test_data[, -ncol(test_data)])
Y_test <- test_data$MRD_binary

cvfit_logistic <- cv.glmnet(
  X_train,
  Y_train,
  family = "binomial",
  alpha = 1  
)

plot(cvfit_logistic)
title("Logistic LASSO Curve\nPredicting MRD > 0.01", line = 2.5)

lasso_coef <- coef(cvfit_logistic, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_values <- lasso_coef[which(lasso_coef != 0)]

print(lasso_genes)

selected_genes_df <- data.frame(
  Gene = lasso_genes,
  Coefficient = as.numeric(lasso_values)
)

write_csv(selected_genes_df, "LASSO_Selected_Genes_MRD_NatureThreshold.csv")

# Predict on unseen test data
pred_prob <- predict(cvfit_logistic, newx = X_test, s = "lambda.min", type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

print(table(Predicted = pred_class, Actual = Y_test))

# AUROC on test data
library(pROC)

roc_obj <- pROC::roc(Y_test, as.numeric(pred_prob))
print(roc_obj$auc)

plot(roc_obj, print.auc = TRUE, col = "blue",
     main = "AUROC: Logistic LASSO Predicting MRD > 0.01 (Test Data)")

ci.auc(roc_obj)


#run GSEAs on selected genes
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      collection = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol)

#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(selected_genes_df, Gene, Coefficient)
data.gsea <- data.df.sub$Coefficient
names(data.gsea) <- as.character(data.df.sub$Gene)
data.gsea <- sort(data.gsea, decreasing = TRUE)
print(data.gsea)
#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.99, pAdjustMethod = "BH")
GSEA.df <- as_tibble(GSEA.res@result)
GSEA.chosen = GSEA.df[] #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

library(tidyverse)
ggplot(GSEA.chosen) + aes(x=NES, y= fct_reorder(ID, NES), fill = NES) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "Selected genes by Lasso regression") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

#perform GO enrichment analysis
reference <- as.data.frame(keep) %>%
  dplyr::filter(keep == "TRUE")

library(WebGestaltR)
go_selected_genes <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                    "geneontology_Molecular_Function_noRedundant"),
                                 interestGeneType = "genesymbol",
                                 referenceGeneType = "genesymbol",
                                 interestGene = selected_genes_df$Gene,
                                 referenceGene = rownames(reference),
                                 sigMethod = "fdr", fdrThr = 0.999)
GSEA.chosen = go_selected_genes  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Downregulated genes in clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 
