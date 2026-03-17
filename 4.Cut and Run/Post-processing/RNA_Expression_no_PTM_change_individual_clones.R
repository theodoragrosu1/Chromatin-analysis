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
RPKM_me3 <- read.csv("path\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv") %>%
  drop_na()
RPKM_ac <- read.csv("path\\macs2_peaks\\H3K27ac\\contrasts\\ac_rpkm.csv") %>%
  drop_na()
RPKM_ub <- read.csv("path\\macs2_peaks\\H2AK119ub\\contrasts\\ub_rpkm.csv") %>%
  drop_na()

RPKM_me3_a7 <- RPKM_me3[,c(-3,-5)]
RPKM_me3_c9 <- RPKM_me3[,c(-2,-4)]

RPKM_ac_a7 <- RPKM_ac[,c(-3,-5)]
RPKM_ac_c9 <- RPKM_ac[,c(-2,-4)]

me3_lost_a7 <- RPKM_me3_a7 %>%
  dplyr::filter(logFC_A7_WT < 0)
me3_lost_c9 <- RPKM_me3_c9 %>%
  dplyr::filter(logFC_C9_WT < 0)
me3_gained_a7 <- RPKM_me3_a7 %>%
  dplyr::filter(logFC_A7_WT > 0)
me3_gained_c9 <- RPKM_me3_c9 %>%
  dplyr::filter(logFC_C9_WT > 0)

ac_lost_a7 <- RPKM_ac_a7 %>%
  dplyr::filter(logFC_A7_WT < 0)
ac_lost_c9 <- RPKM_ac_c9 %>%
  dplyr::filter(logFC_C9_WT < 0)
ac_gained_a7 <- RPKM_ac_a7 %>%
  dplyr::filter(logFC_A7_WT > 0)
ac_gained_c9 <- RPKM_ac_c9 %>%
  dplyr::filter(logFC_C9_WT > 0)

normalised_expression_wt_vs_a7 <- read.csv("path\\normalised_expression_jurkat_wt_vs_a7.csv")

normalised_expression_wt_vs_c9 <- read.csv("path\\normalised_expression_jurkat_wt_vs_c9.csv")

normalised_expression_wt_vs_a7$geneID <- as.character(normalised_expression_wt_vs_a7$geneID)
normalised_expression_wt_vs_c9$geneID <- as.character(normalised_expression_wt_vs_c9$geneID)

me3_lost_a7$geneID <- as.character(me3_lost_a7$geneID)
me3_lost_c9$geneID <- as.character(me3_lost_c9$geneID)
me3_gained_a7$geneID <- as.character(me3_gained_a7$geneID)
me3_gained_c9$geneID <- as.character(me3_gained_c9$geneID)

ac_lost_a7$geneID <- as.character(ac_lost_a7$geneID)
ac_lost_c9$geneID <- as.character(ac_lost_c9$geneID)
ac_gained_a7$geneID <- as.character(ac_gained_a7$geneID)
ac_gained_c9$geneID <- as.character(ac_gained_c9$geneID)

###################################################################
upregulated_genes_a7 <- normalised_expression_wt_vs_a7 %>%
  dplyr::filter(LogFC >= 0.5)

upregulated_genes_c9 <- normalised_expression_wt_vs_c9 %>%
  dplyr::filter(LogFC >= 0.5)

downregulated_genes_a7 <- normalised_expression_wt_vs_a7 %>%
  dplyr::filter(LogFC <= -0.5)

downregulated_genes_c9 <- normalised_expression_wt_vs_c9 %>%
  dplyr::filter(LogFC <= -0.5)

upreg_genes_with_me3_a7 <- upregulated_genes_a7 %>%
  inner_join(RPKM_me3_a7, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_me3_a7$logFC_A7_WT)

plot_data <- upreg_genes_with_me3_a7 %>%
  dplyr::select(geneID, LogFC, logFC_A7_WT) %>%
  pivot_longer(
    cols = c(LogFC, logFC_A7_WT),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Upregulated in A7", logFC_A7_WT = "H3K27me3"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27me3") +
  theme_classic() +
  theme(legend.position = "none")

upreg_no_me3_change_a7 <- upreg_genes_with_me3_a7 %>%
  dplyr::filter(logFC_A7_WT > median_log2FC)  #there are 324 genes where H3K27me3 does not change but the genes are upregulated in the clones


upreg_genes_with_me3_c9 <- upregulated_genes_c9 %>%
  inner_join(RPKM_me3_c9, by = "geneID") %>%
  drop_na()
median_log2FC <- median(upreg_genes_with_me3_c9$logFC_C9_WT)

plot_data <- upreg_genes_with_me3_c9 %>%
  dplyr::select(geneID, LogFC, logFC_C9_WT) %>%
  pivot_longer(
    cols = c(LogFC, logFC_C9_WT),
    names_to = "type",
    values_to = "logFC"
  ) %>%
  dplyr::mutate(type = recode(type, logFC = "Upregulated in C9", logFC_C9_WT = "H3K27me3"))

ggplot(plot_data, aes(x = type, y = logFC, group = geneID)) +
  geom_line(alpha = 0.3, color = "gray") +  # lines connecting same gene
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "black", fill = "pink") +
  geom_hline(yintercept = median_log2FC, linetype = "dashed", color = "red", size = 1) +
  ylab("logFC") +
  xlab("") +
  ggtitle("Gene-wise Expression vs H3K27me3") +
  theme_classic() +
  theme(legend.position = "none")

upreg_no_me3_change_c9 <- upreg_genes_with_me3_c9 %>%
  dplyr::filter(logFC_C9_WT > median_log2FC)  #there are 533 genes where H3K27me3 does not change but the genes are upregulated in the clones


#######################
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

