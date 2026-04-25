library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# 1. Load manifest — skip first 7 rows (Illumina header junk)
manifest <- fread("path\\MethylationEPIC v2.0 Files\\EPIC-8v2-0_A2.csv", skip = 7)

# Keep probe ID and hg38 coordinates only
probe_coords <- manifest[, .(IlmnID, CHR, MAPINFO, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, UCSC_RefGene_Name)]
probe_coords <- probe_coords[!is.na(MAPINFO) & CHR != ""]

mean_beta_values <- read.csv("path\\Mean_Beta_values_Sesame_QCDPB_Jurkat_cell_line_930659CpGs(Mean_Beta_values_Sesame_QCDPB_J).csv")  # columns: cg_id, beta
head(mean_beta_values)
dim(mean_beta_values)
colnames(mean_beta_values)

mean_beta_values <- mean_beta_values %>%
  mutate(
    delta_A7 = Jurkat.A7 - Jurkat.WT,
    delta_C9 = Jurkat.C9 - Jurkat.WT
  ) %>%
  na.omit()
###hypermethylated regions
hypermethylated_A7_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_A7 > 0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.A7", "delta_A7"))
row.names(hypermethylated_A7_vs_WT) <- hypermethylated_A7_vs_WT$Name

hypermethylated_C9_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_C9 > 0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.C9", "delta_C9"))
row.names(hypermethylated_C9_vs_WT) <- hypermethylated_C9_vs_WT$Name

hypermethylated_clones_common <- left_join(hypermethylated_C9_vs_WT, hypermethylated_A7_vs_WT, by = "Name") %>%
  na.omit()

head(hypermethylated_clones_common)

##hypomethylated regions
hypomethylated_A7_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_A7 < -0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.A7", "delta_A7"))
row.names(hypomethylated_A7_vs_WT) <- hypomethylated_A7_vs_WT$Name

hypomethylated_C9_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_C9 < -0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.C9", "delta_C9"))
row.names(hypomethylated_C9_vs_WT) <- hypomethylated_C9_vs_WT$Name

hypomethylated_clones_common <- left_join(hypomethylated_C9_vs_WT, hypomethylated_A7_vs_WT, by = "Name") %>%
  na.omit()
head(hypomethylated_clones_common)


###annotate regions based on proximity to genes
probe_coords_clean <- probe_coords %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, UCSC_RefGene_Group, .keep_all = TRUE) 

annotation_genes_clean <- probe_coords_clean %>%
  mutate(
    region_priority = case_when(
      grepl("TSS200", UCSC_RefGene_Group) ~ "TSS200",
      grepl("TSS1500", UCSC_RefGene_Group) ~ "TSS1500",
      grepl("5UTR", UCSC_RefGene_Group) ~ "5UTR",
      grepl("3UTR", UCSC_RefGene_Group) ~ "3UTR",
      grepl("exon", UCSC_RefGene_Group) ~ "Exon",
      TRUE ~ "intron/intergenic"
    )
  ) %>%
  select(Name, region_priority)


mean_beta_values_annotated <- mean_beta_values %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

###Analysis for separate
hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##Anaylsis for common sites
hypermethylated_clones_common_annotated <- hypermethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_clones_common_annotated <- hypomethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##############is there any correlation between TSS methylation and gene expression
####annotate and extract the genes that have changed beta value at TSS200 only!!!
hypermethylated_clones_common_annotated_tss <- left_join(hypermethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS200", UCSC_RefGene_Group)) 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  mutate(UCSC_RefGene_Name = str_remove(UCSC_RefGene_Name, ";.*"))

normalised_expression_jurkat <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]
genes_to_plot <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name

normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) 

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

#normalised_expression_jurkat <- normalised_expression_jurkat %>%
  #dplyr::filter(abs(logFC) > 2.5)

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))


normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")

library(pheatmap)
library(RColorBrewer)
library(ggrepel)
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = T, cluster_rows = T, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))

hypermethylated_clones_common_annotated_tss$geneID <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name


expression_wt_vs_a7 <- read.csv("path\\normalised_expression_jurkat_wt_vs_a7.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))
expression_wt_vs_c9 <- read.csv("path\\normalised_expression_jurkat_wt_vs_c9.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))
me_and_rna_a7 <- left_join(hypermethylated_clones_common_annotated_tss, expression_wt_vs_a7, by = "geneID") %>%
  na.omit()
me_and_rna_c9 <- left_join(hypermethylated_clones_common_annotated_tss, expression_wt_vs_c9, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(me_and_rna_a7$delta_A7, me_and_rna_a7$LogFC, method = "pearson")

a7_me3_vs_expr <- me_and_rna_a7 %>%
  drop_na(LogFC) %>%
  ggplot(aes(x = delta_A7, y = LogFC)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 0.05, size = 3) + 
  labs(x = "\u0394 \u03B2 A7 - WT\n at TSS200", y = "mRNA LFC\nA7 - WT") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null") + 
  labs(caption = "Correlation between hypermethylation and gene expression",
)+
  geom_label_repel(data=filter(me_and_rna_a7, LogFC >= 1),
                   aes(label=geneID), max.overlaps = 25) +
  geom_label_repel(data=filter(me_and_rna_a7, LogFC <= -1),
                   aes(label=geneID), max.overlaps = 25)


cor_results <- cor.test(me_and_rna_c9$delta_C9, me_and_rna_c9$LogFC, method = "pearson")

c9_me3_vs_expr <- me_and_rna_c9 %>%
  drop_na(LogFC) %>%
  ggplot(aes(x = delta_C9, y = LogFC)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 0.05, size = 3) + 
  labs(x = "\u0394 \u03B2 C9 - WT\nat TSS200", y = "mRNA LFC\nC9 - WT") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null") + 
  labs(caption = "Correlation between hypermethylation and gene expression",
  )+
  geom_label_repel(data=filter(me_and_rna_c9, LogFC >= 1),
                   aes(label=geneID), max.overlaps = 25) +
  geom_label_repel(data=filter(me_and_rna_c9, LogFC <= -1),
                   aes(label=geneID), max.overlaps = 25)




####annotate and extract the genes that have changed beta value at TSS200 only!!!
hypomethylated_clones_common_annotated_tss <- left_join(hypomethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS200", UCSC_RefGene_Group)) 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  mutate(UCSC_RefGene_Name = str_remove(UCSC_RefGene_Name, ";.*"))

normalised_expression_jurkat <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]
genes_to_plot <- hypomethylated_clones_common_annotated_tss$UCSC_RefGene_Name

normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) 

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

#normalised_expression_jurkat <- normalised_expression_jurkat %>%
#dplyr::filter(abs(logFC) > 2.5)

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))


normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")

library(pheatmap)
library(RColorBrewer)
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = T, cluster_rows = T, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))



####which hypermethylated and which hypomethylated regions have H3K27me3 in clones?
me3_regions <- read.csv("path\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv")
genes <- manifest[, .(IlmnID, UCSC_RefGene_Name)]
genes_clean <- genes %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, UCSC_RefGene_Name, .keep_all = TRUE) 
mean_beta_values_annotated <- mean_beta_values %>%   
  left_join(genes_clean, by = "Name")
mean_beta_values_annotated_tss <- left_join(mean_beta_values_annotated, probe_coords_clean, by = "Name") 
mean_beta_values_annotated_tss <- mean_beta_values_annotated_tss %>%
  dplyr::filter(grepl("TSS", UCSC_RefGene_Group)) 
mean_beta_values_annotated_tss$geneID <- mean_beta_values_annotated_tss$UCSC_RefGene_Name.x
me3_hyperme <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  na.omit() 
library(tidyverse)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(ggplot2)

cor_results <- cor.test(me3_hyperme$logFC_A7_WT, me3_hyperme$delta_A7, method = "pearson")


a7_me3_vs_expr <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  drop_na(logFC_A7_WT) %>%
  ggplot(aes(x = logFC_A7_WT, y = delta_A7)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 1, size = 5) + 
  labs(x = "H3K27me3 LFC\nA7 - WT", y = "\u0394 \u03B2 A7 - WT\nat promoter regions") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null")

cor_results <- cor.test(me3_hyperme$logFC_C9_WT, me3_hyperme$delta_C9, method = "pearson")


c9_me3_vs_expr <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  drop_na(logFC_C9_WT) %>%
  ggplot(aes(x = logFC_C9_WT, y = delta_C9)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 1, size = 5) + 
  labs(x = "H3K27me3 LFC\nC9 - WT", y = "\u0394 \u03B2 C9 - WT\nat promoter regions") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null")


########################are there overlaps between hypermethylated regions and me3 lost regions?
###try venn diagram for visualization
hypermethylated_clones_common_annotated_tss <- left_join(hypermethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS", UCSC_RefGene_Group)) 

hypomethylated_clones_common_annotated_tss <- left_join(hypomethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS", UCSC_RefGene_Group)) 
wt_me3_lost <- me3_regions %>%
  dplyr::filter(logFC_A7_WT < 0 & logFC_C9_WT < 0)
clones_me3_gained <- me3_regions %>%
  dplyr::filter(logFC_A7_WT > 0 & logFC_C9_WT > 0)
hypermethylated_clones_common_annotated_tss$geneID <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name
hypomethylated_clones_common_annotated_tss$geneID <- hypomethylated_clones_common_annotated_tss$UCSC_RefGene_Name

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  pull(geneID)
rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  pull(geneID)
######visualize overlaps
library(VennDiagram)
venn.diagram(list("Hypermethylated promoters\nin both clones" = hypermethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 lost regions\nin both clones" = wt_me3_lost$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypermethylated_genes_me3_lost.png", disable.logging = TRUE)

venn.diagram(list("Hypomethylated promoters\nin both clones" = hypomethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 gained regions\nin both clones" = clones_me3_gained$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypomethylated_genes_me3_gained.png", disable.logging = TRUE)

###there cannot be a proper gain without EZH2?
##################look at ubiquitination: there should be positive correlation?
#####H2AK119ub recruits DNMTs (i think)










###annotate regions based on CGIs
probe_coords_clean <- probe_coords %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, Relation_to_UCSC_CpG_Island, .keep_all = TRUE) 

annotation_genes_clean <- probe_coords_clean %>%
  mutate(
    region_priority = case_when(
      grepl("Island", Relation_to_UCSC_CpG_Island) ~ "Island",
      grepl("S_Shelf", Relation_to_UCSC_CpG_Island) ~ "S_Shelf",
      grepl("S_Shore", Relation_to_UCSC_CpG_Island) ~ "S_Shore",
      grepl("N_Shelf", Relation_to_UCSC_CpG_Island) ~ "N_Shelf",
      grepl("N_Shore", Relation_to_UCSC_CpG_Island) ~ "N_Shore",
      TRUE ~ "Open sea"
    )
  ) %>%
  dplyr::select(Name, region_priority)


mean_beta_values_annotated <- mean_beta_values %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

###Analysis for separate
hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 


##Anaylsis for common sites
hypermethylated_clones_common_annotated <- hypermethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_clones_common_annotated <- hypomethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##############is there any correlation between CpG methylation and gene expression
####annotate and extract the genes that have changed beta value at CGI only!!!
hypermethylated_clones_common_annotated_tss <- left_join(hypermethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("Island", region_priority)) 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  mutate(UCSC_RefGene_Name = str_remove(UCSC_RefGene_Name, ";.*"))

normalised_expression_jurkat <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]
genes_to_plot <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name

normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) 

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

#normalised_expression_jurkat <- normalised_expression_jurkat %>%
#dplyr::filter(abs(logFC) > 2.5)

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))


normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")

library(pheatmap)
library(RColorBrewer)
library(ggrepel)
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))

hypermethylated_clones_common_annotated_tss$geneID <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name


expression_wt_vs_a7 <- read.csv("path\\normalised_expression_jurkat_wt_vs_a7.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))
expression_wt_vs_c9 <- read.csv("path\\normalised_expression_jurkat_wt_vs_c9.csv") %>%
  dplyr::mutate(regulation = case_when(LogFC<0 ~ "Downregulated", 
                                       LogFC>0 ~ "Upregulated"))
me_and_rna_a7 <- left_join(hypermethylated_clones_common_annotated_tss, expression_wt_vs_a7, by = "geneID") %>%
  na.omit()
me_and_rna_c9 <- left_join(hypermethylated_clones_common_annotated_tss, expression_wt_vs_c9, by = "geneID") %>%
  na.omit()

cor_results <- cor.test(me_and_rna_a7$delta_A7, me_and_rna_a7$LogFC, method = "pearson")

a7_me3_vs_expr <- me_and_rna_a7 %>%
  drop_na(LogFC) %>%
  ggplot(aes(x = delta_A7, y = LogFC)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 0.05, size = 3) + 
  labs(x = "\u0394 \u03B2 A7 - WT", y = "mRNA LFC\nA7 - WT") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null") + 
  labs(caption = "Correlation between hypermethylation and gene expression\nA7 EZH2 KO, Island methylation",
  )+
  geom_label_repel(data=filter(me_and_rna_a7, LogFC >= 2.5),
                   aes(label=geneID), max.overlaps = 25) +
  geom_label_repel(data=filter(me_and_rna_a7, LogFC <= -2.5),
                   aes(label=geneID), max.overlaps = 25)


cor_results <- cor.test(me_and_rna_c9$delta_C9, me_and_rna_c9$LogFC, method = "pearson")

c9_me3_vs_expr <- me_and_rna_c9 %>%
  drop_na(LogFC) %>%
  ggplot(aes(x = delta_C9, y = LogFC)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 0.05, size = 3) + 
  labs(x = "\u0394 \u03B2 C9 - WT", y = "mRNA LFC\nC9 - WT") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null") + 
  labs(caption = "Correlation between hypermethylation and gene expression\nC9 EZH2 KO, Island methylation",
  )+
  geom_label_repel(data=filter(me_and_rna_c9, LogFC >= 5),
                   aes(label=geneID), max.overlaps = 25) +
  geom_label_repel(data=filter(me_and_rna_c9, LogFC <= -2.5),
                   aes(label=geneID), max.overlaps = 25)



####annotate and extract the genes that have changed beta value at CpG island only!!!
hypomethylated_clones_common_annotated_tss <- left_join(hypomethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("Island", region_priority)) 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  mutate(UCSC_RefGene_Name = str_remove(UCSC_RefGene_Name, ";.*"))

normalised_expression_jurkat <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]
genes_to_plot <- hypomethylated_clones_common_annotated_tss$UCSC_RefGene_Name

normalised_expression_jurkat <- all%>%
  dplyr::filter(geneID %in% genes_to_plot) 

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

#normalised_expression_jurkat <- normalised_expression_jurkat %>%
#dplyr::filter(abs(logFC) > 2.5)

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))


normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")

library(pheatmap)
library(RColorBrewer)
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))

########################are there overlaps between hypermethylated regions and me3 lost regions?
###try venn diagram for visualization
hypermethylated_clones_common_annotated_tss <- left_join(hypermethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("Island", Relation_to_UCSC_CpG_Island)) 

hypomethylated_clones_common_annotated_tss <- left_join(hypomethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("Island", Relation_to_UCSC_CpG_Island)) 
wt_me3_lost <- me3_regions %>%
  dplyr::filter(logFC_A7_WT < 0 & logFC_C9_WT < 0)
clones_me3_gained <- me3_regions %>%
  dplyr::filter(logFC_A7_WT > 0 & logFC_C9_WT > 0)
hypermethylated_clones_common_annotated_tss$geneID <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name
hypomethylated_clones_common_annotated_tss$geneID <- hypomethylated_clones_common_annotated_tss$UCSC_RefGene_Name

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  pull(geneID)
rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  pull(geneID)
######visualize overlaps
library(VennDiagram)
venn.diagram(list("Hypermethylated CpG\nin both clones" = hypermethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 lost regions\nin both clones" = wt_me3_lost$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypermethylated_cpg_genes_me3_lost.png", disable.logging = TRUE)

venn.diagram(list("Hypomethylated CpG\nin both clones" = hypomethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 gained regions\nin both clones" = clones_me3_gained$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypomethylated_cpg_genes_me3_gained.png", disable.logging = TRUE)


####which hypermethylated and which hypomethylated regions have H3K27me3 in clones?
me3_regions <- read.csv("path\\macs2_peaks\\H3K27me3\\contrasts\\me3_rpkm.csv")
genes <- manifest[, .(IlmnID, UCSC_RefGene_Name)]
genes_clean <- genes %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, UCSC_RefGene_Name, .keep_all = TRUE) 
mean_beta_values_annotated <- mean_beta_values %>%   
  left_join(genes_clean, by = "Name")
mean_beta_values_annotated_tss <- left_join(mean_beta_values_annotated, probe_coords_clean, by = "Name") 
mean_beta_values_annotated_tss <- mean_beta_values_annotated_tss %>%
  dplyr::filter(grepl("Island", Relation_to_UCSC_CpG_Island)) 
mean_beta_values_annotated_tss$geneID <- mean_beta_values_annotated_tss$UCSC_RefGene_Name.x
me3_hyperme <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  na.omit() 
library(tidyverse)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(ggplot2)

cor_results <- cor.test(me3_hyperme$logFC_A7_WT, me3_hyperme$delta_A7, method = "pearson")


a7_me3_vs_expr <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  drop_na(logFC_A7_WT) %>%
  ggplot(aes(x = logFC_A7_WT, y = delta_A7)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 1, size = 5) + 
  labs(x = "H3K27me3 LFC\nA7 - WT", y = "\u0394 \u03B2 A7 - WT\nat CpG islands") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null")

cor_results <- cor.test(me3_hyperme$logFC_C9_WT, me3_hyperme$delta_C9, method = "pearson")


c9_me3_vs_expr <- left_join(mean_beta_values_annotated_tss, me3_regions, by = "geneID") %>%
  drop_na(logFC_C9_WT) %>%
  ggplot(aes(x = logFC_C9_WT, y = delta_C9)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointdensity(size = 2) +
  scale_color_viridis() + 
  stat_cor(aes(label = paste("R==", ..r.., "*paste(\";\")~~~p ==",
                             10^(log10(..p..) %% 1), "%*% 10^",
                             floor(log10(..p..)))), label.x = 1, size = 5) + 
  labs(x = "H3K27me3 LFC\nC9 - WT", y = "\u0394 \u03B2 C9 - WT\nat CpG islands") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_classic(base_size = 24) + 
  theme(legend.position = "null")


########################are there overlaps between hypermethylated regions and me3 lost regions?
###try venn diagram for visualization
hypermethylated_clones_common_annotated_tss <- left_join(hypermethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypermethylated_clones_common_annotated_tss <- hypermethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS", UCSC_RefGene_Group)) 

hypomethylated_clones_common_annotated_tss <- left_join(hypomethylated_clones_common_annotated, probe_coords_clean, by = "Name") 
hypomethylated_clones_common_annotated_tss <- hypomethylated_clones_common_annotated_tss %>%
  dplyr::filter(grepl("TSS", UCSC_RefGene_Group)) 
wt_me3_lost <- me3_regions %>%
  dplyr::filter(logFC_A7_WT < 0 & logFC_C9_WT < 0)
clones_me3_gained <- me3_regions %>%
  dplyr::filter(logFC_A7_WT > 0 & logFC_C9_WT > 0)
hypermethylated_clones_common_annotated_tss$geneID <- hypermethylated_clones_common_annotated_tss$UCSC_RefGene_Name
hypomethylated_clones_common_annotated_tss$geneID <- hypomethylated_clones_common_annotated_tss$UCSC_RefGene_Name

wt_vs_clones_diff_expr <- read.csv("path\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))

rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  pull(geneID)
rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  pull(geneID)
######visualize overlaps
library(VennDiagram)
venn.diagram(list("Hypermethylated promoters\nin both clones" = hypermethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 lost regions\nin both clones" = wt_me3_lost$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypermethylated_genes_me3_lost.png", disable.logging = TRUE)

venn.diagram(list("Hypomethylated promoters\nin both clones" = hypomethylated_clones_common_annotated_tss$geneID,
                  "H3K27me3 gained regions\nin both clones" = clones_me3_gained$geneID),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#173518"),
             "path\\overlap_hypomethylated_genes_me3_gained.png", disable.logging = TRUE)




