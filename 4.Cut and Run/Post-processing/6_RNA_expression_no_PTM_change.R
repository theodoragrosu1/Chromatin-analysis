library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)
library(dplyr)
library(Cairo)
library(ggplot2)
library(colorspace)

#this code imports RPKM normalized Cut and Run and sees which DEGs cannot be explained by changes in PTMs
RPKM_me3 <- read.csv("path\\me3_rpkm.csv") %>%
  drop_na()
RPKM_ac <- read.csv("path\\ac_rpkm.csv") %>%
  drop_na()
RPKM_ub <- read.csv("path\\ub_rpkm.csv") %>%
  drop_na()

RPKM_me3_avg <- RPKM_me3%>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)

RPKM_ac_avg <- RPKM_ac %>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)

RPKM_ub_avg <- RPKM_ub %>%
  dplyr::mutate(avg_logFC = (logFC_A7_WT + logFC_C9_WT)/2)

me3_lost_clones <- RPKM_me3_avg %>%
  dplyr::filter(avg_logFC < 0)
me3_gained_clones <- RPKM_me3_avg %>%
  dplyr::filter(avg_logFC > 0)
ac_lost_clones <- RPKM_ac_avg %>%
  dplyr::filter(avg_logFC < 0)
ac_gained_clones <- RPKM_ac_avg %>%
  dplyr::filter(avg_logFC >0)
ub_lost_clones <- RPKM_ub_avg %>%
  dplyr::filter(avg_logFC < 0)
ub_gained_clones <- RPKM_ub_avg %>%
  dplyr::filter(avg_logFC >0)

normalised_expression_wt_vs_clones <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

normalised_expression_wt_vs_clones$geneID <- as.character(normalised_expression_wt_vs_clones$geneID)
me3_lost_clones$geneID <- as.character(me3_lost_clones$geneID)
me3_gained_clones$geneID <- as.character(me3_gained_clones$geneID)
ac_lost_clones$geneID <- as.character(ac_lost_clones$geneID)
ac_gained_clones$geneID <- as.character(ac_gained_clones$geneID)
ub_lost_clones$geneID <- as.character(ub_lost_clones$geneID)
ub_gained_clones$geneID <- as.character(ub_gained_clones$geneID)

upregulated_genes_clones <- normalised_expression_wt_vs_clones %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val <= 0.1)
write.csv(upregulated_genes_clones, file = "path\\upregulated_genes_clones.csv")

downregulated_genes_clones <- normalised_expression_wt_vs_clones %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val <= 0.1)
write.csv(downregulated_genes_clones, file = "path\\downregulated_genes_clones.csv")

upreg_genes_with_me3_wt <- upregulated_genes_clones %>%
  inner_join(RPKM_me3_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_me3_wt$avg_logFC)

plot_data <- upreg_genes_with_me3_wt %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Upregulated in clones", avg_logFC = "H3K27me3"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27me3") +
  theme_classic() +
  theme(legend.position = "none")

upreg_no_me3_change_wt <- upreg_genes_with_me3_wt %>%
  dplyr::filter(avg_logFC > median_log2FC)  #there are 23 genes where H3K27me3 does not change but the genes are upregulated in the clones


downreg_genes_with_me3_clones <- downregulated_genes_clones %>%
  inner_join(RPKM_me3_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(downreg_genes_with_me3_clones$avg_logFC)

plot_data <- downreg_genes_with_me3_clones %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Downregulated in clones", avg_logFC = "H3K27me3"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27me3") +
  theme_classic() +
  theme(legend.position = "none")


downreg_no_me3_change_clones <- downreg_genes_with_me3_clones %>%
  dplyr::filter(avg_logFC > median_log2FC)  #there are 5 genes where H3K27me3 doesn't change but the genes are downregulated in the clones


####at which H3K27ac does gene expression changes
upreg_genes_with_ac_wt <- upregulated_genes_clones %>%
  inner_join(RPKM_ac_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_ac_wt$avg_logFC)

plot_data <- upreg_genes_with_ac_wt %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Upregulated in clones", avg_logFC = "H3K27ac"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27ac") +
  theme_classic() +
  theme(legend.position = "none")

upreg_no_ac_change_clones <- upreg_genes_with_ac_wt %>%
  dplyr::filter(avg_logFC < median_log2FC)  #there are 31 genes where H3K27ac does not change but the genes are upregulated in the clones

downreg_genes_with_ac_wt <- downregulated_genes_clones %>%
  inner_join(RPKM_me3_avg, by = "geneID") %>%
  drop_na()
median_log2FC <- median(downreg_genes_with_ac_wt$avg_logFC)

plot_data <- downreg_genes_with_ac_wt %>%
  dplyr::select(geneID, logFC, avg_logFC) %>%
  pivot_longer(
    cols = c(logFC, avg_logFC),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Downregulated in clones", avg_logFC = "H3K27ac"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27ac") +
  theme_classic() +
  theme(legend.position = "none")

downreg_no_ac_change_wt <- downreg_genes_with_ac_wt %>%
  dplyr::filter(avg_logFC > median_log2FC)  #there are 5 genes where H3K27ac doesn't change but the genes are downregulated in the clones

