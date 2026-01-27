library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggraph)
library(igraph)
library(WebGestaltR)
library(tidyverse)
library(viper)
library(dorothea)

#input files: normalized RNAseq expression and TRRUST data from the reference database of human transcriptional regulatory interactions https://www.nature.com/articles/srep11432

Jurkat_DEGs <- read.csv(".\\Jurkat_wt_v_clones_all_genes (1).csv")
trrust_df <- read.csv(".\\trrust_rawdata.human.tsv", sep="\t") 
trrust_df <- trrust_df[,-6]
trrust_df <- trrust_df[,-5]

filtered_trrust <- trrust_df %>%
  dplyr::filter(Target %in% Jurkat_DEGs$geneID)

#count how many DEGs are regulated by each TF - this gives you TFs that regulate the most DEGs in the list
tf_enrichment <- filtered_trrust %>%
  group_by(TF) %>%
  summarize(DEG_targets = list(Target), overlap_count = n()) %>%
  arrange(desc(overlap_count))


tf_enrichment_top <- tf_enrichment %>%
  top_n(30)

go_analysis <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                         "geneontology_Molecular_Function_noRedundant"),
                                      interestGeneType = "genesymbol",
                                      referenceGeneType = "genesymbol",
                                      interestGene = tf_enrichment$TF,
                                      referenceGene = Jurkat_DEGs$geneID,
                                      sigMethod = "fdr", fdrThr = 0.01)
GSEA.chosen = go_analysis[1:30,]  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "TFs enriched in DEGs clones - WT") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 

#overlap with the TFs known to appear in ATACseq
atac_tf_a7 <- read.table(".\\bindetect_results_a7.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
atac_tf_c9 <- read.table(".\\bindetect_results_c9.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
atac_tf_overlap <- merge(atac_tf_a7, atac_tf_c9, by = "motif_id") 

atac_tf_overlap <- atac_tf_overlap[,-c(2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19)]

intersect <- intersect(tf_enrichment$TF, atac_tf_overlap$motif_id)
intersect[1:10]


go_analysis <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                              "geneontology_Molecular_Function_noRedundant"),
                           interestGeneType = "genesymbol",
                           referenceGeneType = "genesymbol",
                           interestGene = intersect,
                           referenceGene = atac_tf_overlap$motif_id,
                           sigMethod = "fdr", fdrThr = 0.05)
GSEA.chosen = go_analysis[1:30,]  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
ggplot(GSEA.chosen) + aes(x=-log10(FDR), y= fct_reorder(description, -FDR), fill = FDR) +
  geom_bar(stat = "identity") +
  labs(y = "GO term", title = "TFs enriched in DEGs clones - WT overlapped with DE TFs in ATACseq") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 


#regulon activity using VIPER and dorothea
expression_matrix <-  read.csv(".\\normalised_expression_jurkat_wt_vs_clones.csv", row.names = 1)

data.df <- mutate(expression_matrix,
                  jurkat.wt.AVG = (A13 + A11 + A12)/3,
                  jurkat.A7.AVG = (A10 + A14 + A15)/3,
                  jurkat.C9.AVG = (A16 + A17 + A18)/3) %>% 
  mutate_if(is.numeric, round, 2)
expression_matrix <- data.df[, c(10,11,12)]

study_design <- read_tsv(".\\study_Design.txt")

all(colnames(expression_matrix) %in% study_design$sample_id)

# Load regulons (human)
data(dorothea_hs, package = "dorothea")

# Use high-confidence TF-target interactions (A, B, C)
regulons_df <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

# Make sure gene names are uppercase to match regulon format
rownames(expression_matrix) <- toupper(rownames(expression_matrix))

# Create list of regulons manually
regulons <- lapply(split(regulons_df, regulons_df$tf), function(tf_df) {
  tf_targets <- tf_df %>%
    select(target, mor) %>%
    distinct() %>%
    tibble::column_to_rownames("target")
  
  # Ensure numeric
  tfmode <- as.numeric(tf_targets$mor)
  names(tfmode) <- rownames(tf_targets)
  
  # Assign likelihood = 1 for all targets
  likelihood <- rep(1, length(tfmode))
  names(likelihood) <- names(tfmode)
  
  list(tfmode = tfmode, likelihood = likelihood)
})


tf_activity <- viper(expression_matrix, regulons, method = "scale")

# Merge TF activity with metadata
tf_activity_df <- as.data.frame(tf_activity)
tf_activity_df$TF <- rownames(tf_activity_df)

tf_long <- tf_activity %>%
  as.data.frame() %>%
  rownames_to_column(var = "TF") %>%
  pivot_longer(-TF, names_to = "sample_id", values_to = "activity") %>%
  left_join(study_design, by = "sample_id")


anova_results <- tf_long %>%
  group_by(TF) %>%
  summarise(
    p_value = summary(aov(activity ~ condition))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  )

tukey_results <- tf_long %>%
  group_by(TF) %>%
  group_modify(~ {
    model <- aov(activity ~ condition, data = .x)
    tukey <- TukeyHSD(model)
    res <- as.data.frame(tukey$condition)
    res$comparison <- rownames(res)
    res$TF <- unique(.x$TF)
    res
  }) %>%
  ungroup()

A7.TF <- tukey_results %>%
  dplyr::filter(comparison == "WT-A7") %>%
  dplyr:: select(TF, diff) %>%
  dplyr::filter(abs(diff) > 1)


top_tfs_a7 <- A7.TF %>% 
  mutate(abs_mean_diff = abs(diff)) %>%
  arrange(desc(abs_mean_diff)) %>%
  slice_head(n = 25)


tf_order <- top_tfs_a7 %>%
  mutate(
    mean_group = ifelse(diff >= 0, "Positive mean", "Negative mean")
  ) %>%
  arrange(mean_group, desc(diff))

# Order TFs within groups descending by overall mean
tf_order_pos <- tf_order %>%
  filter(mean_group == "Positive mean") %>%
  arrange(desc(diff)) %>%
  mutate(TF_ordered = factor(TF, levels = TF))

tf_order_neg <- tf_order %>%
  filter(mean_group == "Negative mean") %>%
  arrange(desc(diff)) %>%
  mutate(TF_ordered = factor(TF, levels = TF))

# Combine back
tf_order <- bind_rows(tf_order_pos, tf_order_neg)


plot_data <- tf_order %>%
  mutate(
    TF = factor(TF, levels = levels(tf_order$TF_ordered)),
    mean_group = factor(mean_group, levels = c("Positive mean", "Negative mean"))
  )

ggplot(plot_data, aes(x = TF, y = diff, fill = diff)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_gradient(low = "blue", high = "red")+
  labs(
    title = "Top 25 Differentially Active TFs",
    x = "WT - A7",
    y = "Mean TF Activity"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



C9.TF <- tukey_results %>%
  dplyr::filter(comparison == "WT-C9") %>%
  dplyr:: select(TF, diff) %>%
  dplyr::filter(abs(diff) > 1)


top_tfs_c9 <- C9.TF %>% 
  mutate(abs_mean_diff = abs(diff)) %>%
  arrange(desc(abs_mean_diff)) %>%
  slice_head(n = 25)


tf_order.c9 <- top_tfs_c9 %>%
  mutate(
    mean_group = ifelse(diff >= 0, "Positive mean", "Negative mean")
  ) %>%
  arrange(mean_group, desc(diff))

# Order TFs within groups descending by overall mean
tf_order_pos.c9 <- tf_order.c9 %>%
  filter(mean_group == "Positive mean") %>%
  arrange(desc(diff)) %>%
  mutate(TF_ordered = factor(TF, levels = TF))

tf_order_neg.c9 <- tf_order.c9 %>%
  filter(mean_group == "Negative mean") %>%
  arrange(desc(diff)) %>%
  mutate(TF_ordered = factor(TF, levels = TF))

# Combine back
tf_order.c9 <- bind_rows(tf_order_pos.c9, tf_order_neg.c9)


plot_data.c9 <- tf_order.c9 %>%
  mutate(
    TF = factor(TF, levels = levels(tf_order.c9$TF_ordered)),
    mean_group = factor(mean_group, levels = c("Positive mean", "Negative mean"))
  )

ggplot(plot_data.c9, aes(x = TF, y = diff, fill = diff)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_gradient(low = "blue", high = "red")+
  labs(
    title = "Top 25 Differentially Active TFs",
    x = "WT - C9",
    y = "Mean TF Activity"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


