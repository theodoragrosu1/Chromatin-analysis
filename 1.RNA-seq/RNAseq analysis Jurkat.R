library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase)
library(gprofiler2) 
library(gplots)
library(enrichplot)
library(ggrepel)

###Jurkat wt vs clones
#This script performs differential gene expression analysis and GSEA on Jurkat EZH2-WT vs Jurkat-EZH2-KO
#Input: kallisto abundance for each sample - available on GEO
#Input: study design file available with the code
#Output: volcano plot, GSEA plots
#Output: differentially expressed genes comparing WT with both A7 and C9 EZH2-KO clones
#Output: files for downstream PROGENy analysis

###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(limma)
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot)


### read data----
design <- read_tsv("study_design_jurkat.txt") 
path <- file.path("#folder for raw data", design$sample, "abundance.tsv") # set file paths to mapped data

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) %>% #gene symbols 
  as_tibble() %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select("target_id", "gene_name")

Tx_gene <- tximport(path, #imports the data for all samples
                    type = "kallisto",
                    tx2gene = Tx,
                    txOut = FALSE, #data represented at gene level rather than transcript
                    countsFromAbundance = "lengthScaledTPM", #transcripts per million
                    ignoreTxVersion = TRUE) 


### preprocessing----
sample_labels <- design$sample
DGEList <- DGEList(Tx_gene$counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.df) <- c("geneID", sample_labels)

###tidy data
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                  names_to = "samples", # name of new column
                                  values_to = "expression") # name of new column storing all the data



###filter data
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=3 #user defined - depends on studydesign, eliminates lowly expressed genes
DGEList.filtered <- DGEList[keepers,]

###normalize data
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>% 
  cpm(log=TRUE) %>% 
  as_tibble(rownames = "geneID")

colnames(log2.cpm.filtered.norm) <- c("geneID", sample_labels)


### multivariate analysis ----
data.df <- mutate(log2.cpm.filtered.norm,
                  jurkat.wt.AVG = (A13 + A11 + A12)/3,
                  jurkat.clones.AVG = (A10 + A14 + A15 + A16 + A17 + A18)/6,
                  LogFC = (jurkat.clones.AVG - jurkat.wt.AVG)) %>% 
  mutate_if(is.numeric, round, 2)


group <- design$phenotype
group <- factor(group)


write.csv(data.df[ ,c(1:10)], file = "normalised_expression_jurkat_wt_vs_clones.csv", row.names = FALSE)


### Differential gene expression analysis ----
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE) #models the mean-variance trend
fit <- lmFit(v.DEGList.filtered.norm, design) #fits linear model
contrast.matrix <- makeContrasts(differences = Jurkat_c - Jurkat_wt,
                                 levels=design) #contrast matrix

fits <- contrasts.fit(fit, contrast.matrix) #contrasts between the 2 groups
ebFit <- eBayes(fits)
TopHits <- topTable(ebFit, adjust ="BH", coef=1, number=20000, sort.by="logFC")

TopHits.df <- TopHits %>%
  as_tibble(rownames = "geneID")

write.csv(TopHits.df, file = "#path and name to save differential expression analysis", quote = FALSE, row.names = FALSE)


###diff genes - volcano plot ----
results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(results) #results are added manually to vplot upstream

vplot <- ggplot(TopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  annotate("rect", xmin = 1, xmax = 11, ymin = -log10(0.05), ymax = 2.1, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -11, ymin = -log10(0.05), ymax = 2.1, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -10, y = 1.4, label = "10", size = 7, colour = "#2C467A") + #labels are manually added after running summary(results)
  annotate(geom = "text", x = 10, y = 1.4, label = "12", size = 7, colour = "#BE684D") +
  labs(title="Jurkat wt & clones",
       subtitle = "Volcano plot") +
  theme_bw(base_size = 15)
vplot



diff_expr_jurkat_vs_clones <- read.csv("#path to the file containing differential analysis")

upregulated_genes <- diff_expr_jurkat_vs_clones %>%
  dplyr::filter(logFC >= 1) %>%
  dplyr::filter(adj.P.Val > 0.05) %>%
  dplyr::pull(geneID)
downregulated_genes <- diff_expr_jurkat_vs_clones %>%
  dplyr::filter(logFC <= -1) %>%
  dplyr::filter(adj.P.Val > 0.05) %>%
  dplyr::pull(geneID)

upreg_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                         "geneontology_Molecular_Function_noRedundant"),
                                      interestGeneType = "genesymbol",
                                      referenceGeneType = "genesymbol",
                                      interestGene = upregulated_genes,
                                      referenceGene = diff_expr_jurkat_vs_clones$geneID,
                                      sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = upreg_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen[1:30,]) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "Upregulated genes in clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

down_genes_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                        "geneontology_Molecular_Function_noRedundant"),
                                     interestGeneType = "genesymbol",
                                     referenceGeneType = "genesymbol",
                                     interestGene = downregulated_genes,
                                     referenceGene = normalised_expression_jurkat_vs_clones$geneID,
                                     sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = down_genes_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
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
