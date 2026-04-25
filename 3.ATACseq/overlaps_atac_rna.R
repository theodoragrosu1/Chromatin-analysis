library(VennDiagram)
library(WebGestaltR)
library(UpSetR)
library(tidyverse)
library(Cairo)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(colorspace)
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\"

constrasts <- c("A7_peaks_only", "C9_peaks_only", "WT_peaks_only", 
                "WT_A7_C9_peaks_overlap")

contrasts_peak_files <- list(A7_peaks_only = paste0(dir, "A7_ATAC_peaks.bed"), 
                             C9_peaks_only = paste0(dir, "C9_ATAC_peaks.bed"), 
                             WT_peaks_only = paste0(dir, "WT_ATAC_peaks.bed"),
                             WT_A7_C9_peaks_overlap = paste0(dir, "A7_C9_WT_merged.bed"))

annotate_peaks_func <- function(peak_files_list, chosen_df)  {
  annotated_peaks <- annotatePeak(peak_files_list[[chosen_df]], tssRegion = c(-3000, 3000),
                                  TxDb = edb, annoDb = "org.Hs.eg.db") %>%
    as.data.frame()
  
  annotated_peaks$annotation <- gsub("^Intron.*", "Intron", annotated_peaks$annotation)
  annotated_peaks$annotation <- gsub("^Exon.*", "Exon", annotated_peaks$annotation)
  return(annotated_peaks)
}

A7_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "A7_peaks_only")
C9_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "C9_peaks_only")
WT_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "WT_peaks_only")
WT_A7_C9_peaks_overlap_annotated <- annotate_peaks_func(contrasts_peak_files, "WT_A7_C9_peaks_overlap")

#write.csv(WT_peaks_only_annotated, file = paste0(dir, "\\WT_unique_atac.csv"), quote = FALSE, row.names = FALSE)
#write.csv(A7_peaks_only_annotated, file = paste0(dir, "\\A7_unique_atac.csv"), quote = FALSE, row.names = FALSE)
#write.csv(C9_peaks_only_annotated, file = paste0(dir, "\\C9_unique_atac.csv"), quote = FALSE, row.names = FALSE)
#write.csv(WT_A7_C9_peaks_overlap_annotated, file = paste0(dir, "\\WT_A7_C9_overlap_atac.csv"), quote = FALSE, row.names = FALSE)


A7_peaks_only_annotated_summary <- as.data.frame(table(A7_peaks_only_annotated$annotation))
colnames(A7_peaks_only_annotated_summary) <- c("Annotation", "A7_peaks_only")

C9_peaks_only_annotated_summary <- as.data.frame(table(C9_peaks_only_annotated$annotation))
colnames(C9_peaks_only_annotated_summary) <- c("Annotation", "C9_peaks_only")

WT_peaks_only_annotated_summary <- as.data.frame(table(WT_peaks_only_annotated$annotation))
colnames(WT_peaks_only_annotated_summary) <- c("Annotation", "WT_peaks_only")

WT_A7_C9_peaks_overlap_annotated_summary <- as.data.frame(table(WT_A7_C9_peaks_overlap_annotated$annotation))
colnames(WT_A7_C9_peaks_overlap_annotated_summary) <- c("Annotation", "WT_A7_C9_peaks_overlap")


A7_peaks_only_annotated_summary_perc <- A7_peaks_only_annotated_summary %>%
  mutate(A7_percentage = A7_peaks_only/sum(A7_peaks_only_annotated_summary$A7_peaks_only)*100) %>%
  dplyr::select(-A7_peaks_only)

C9_peaks_only_annotated_summary_perc <- C9_peaks_only_annotated_summary %>%
  mutate(C9_percentage = C9_peaks_only/sum(C9_peaks_only_annotated_summary$C9_peaks_only)*100) %>%
  dplyr::select(-C9_peaks_only)

WT_peaks_only_annotated_summary_perc <- WT_peaks_only_annotated_summary %>%
  mutate(WT_percentage = WT_peaks_only/sum(WT_peaks_only_annotated_summary$WT_peaks_only)*100) %>%
  dplyr::select(-WT_peaks_only)

WT_A7_C9_peaks_overlap_annotated_summary_perc <- WT_A7_C9_peaks_overlap_annotated_summary %>%
  mutate(WT_A7_C9_percentage = WT_A7_C9_peaks_overlap/sum(WT_A7_C9_peaks_overlap_annotated_summary$WT_A7_C9_peaks_overlap)*100) %>%
  dplyr::select(-WT_A7_C9_peaks_overlap)


summary_annotation_perc <- left_join(A7_peaks_only_annotated_summary_perc, C9_peaks_only_annotated_summary_perc, by = "Annotation") %>%
  left_join(WT_peaks_only_annotated_summary_perc, by = "Annotation") %>%
  left_join(WT_A7_C9_peaks_overlap_annotated_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))
summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition, 
                                            levels = c("WT_A7_C9_percentage", 
                                                       "C9_percentage", "A7_percentage", "WT_percentage"))

annotated_peaks_plot <- summary_annotation_perc %>%
  ggplot(aes(x = Condition, y = percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               "#FEB95F", "#CC79A7",
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5)),
                    guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  labs(x = "Sample") +
  ylab("Accessible regions (%)")

ggsave(filename = "./ATAC_Seq/plots/annotated_constrast_atac_peaks_percentage.png", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)

####Overlap between peaks and RNA expression#########
wt_vs_a7_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_a7.csv") %>%
  dplyr::filter(abs(LogFC) >= 0.5) 
a7_down <- wt_vs_a7_diff_expr %>%
  dplyr::filter(LogFC < -0.5) %>%
  dplyr::pull(geneID)
a7_up <- wt_vs_a7_diff_expr %>%
  dplyr::filter(LogFC > 0.5) %>%
  dplyr::pull(geneID)

wt_vs_c9_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_c9.csv") %>%
  dplyr::filter(abs(LogFC) > 0.5) 
c9_down <- wt_vs_c9_diff_expr %>%
  dplyr::filter(LogFC < -0.5) %>%
  dplyr::pull(geneID)
c9_up <- wt_vs_c9_diff_expr %>%
  dplyr::filter(LogFC > 0.5) %>%
  dplyr::pull(geneID)


a7_open_promoters <- A7_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

a7_open <- A7_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()


c9_open_promoters <- C9_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c9_open <- C9_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

WT_minus_A7_annotated <- anti_join(WT_peaks_only_annotated, A7_peaks_only_annotated, by = "geneId")

a7_closed_promoters <- WT_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)" | annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()
a7_closed <- WT_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

#add to change thresholds # | annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)" 
WT_minus_C9_annotated <- anti_join(WT_peaks_only_annotated, C9_peaks_only_annotated, by = "geneId")
c9_closed_promoters <- WT_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)"| annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c9_closed <- WT_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

svg(file="C:\\Users\\tgrosu\\OneDrive\\Desktop\\ATAC-peaks\\Upset_a7_RNA_ATAC.svg", width = 15, height = 6)
UpSetR::upset(fromList(list("A7 RNA \u2191" = a7_up, "A7 RNA \u2193" = a7_down, 
                    "A7 promoter accessibilty \u2191" = a7_open_promoters, "A7 promoter accessibilty \u2193" = a7_closed_promoters)),
      order.by = "freq", text.scale = 2, point.size = 3,
      queries = list(list(query = intersects, 
                          params = list("A7 RNA \u2191", "A7 promoter accessibilty \u2191"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA \u2193", "A7 promoter accessibilty \u2193"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA \u2193"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter accessibilty \u2193"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA \u2191"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter accessibilty \u2191"), 
                          color = "#6679AB", active = T)),
      sets.bar.color = c("#E68D74", "#6679AB", "#E68D74", "#6679AB"))


dev.off()


svg(file="C:\\Users\\tgrosu\\OneDrive\\Desktop\\ATAC-peaks\\Upset_c9_RNA_ATAC.svg", width = 15, height = 6)

UpSetR::upset(fromList(list("C9 RNA \u2191" = c9_up, "C9 RNA \u2193" = c9_down, 
                            "C9 promoter accessibilty \u2191" = c9_open_promoters, "C9 promoter accessibilty \u2193" = c9_closed_promoters)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("C9 RNA \u2191", "C9 promoter accessibilty \u2191"), 
                                  color = "#E68D74", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2193", "C9 promoter accessibilty \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 promoter accessibilty \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2191"), 
                                  color = "#E68D74", active = T),
                             list(query = intersects, 
                                  params = list("C9 promoter accessibilty \u2191"), 
                                  color = "#E68D74", active = T)),
              sets.bar.color = c("#E68D74", "#6679AB", "#E68D74", "#6679AB"))

dev.off()

setwd("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ATAC-peaks\\")
venn.diagram(list("A7 RNA \u2191\nA7 promoter accessibilty \u2191" = intersect(a7_up, a7_open_promoters), 
                  "C9 RNA \u2191\nC9 promoter accessibilty \u2191" = intersect(c9_up, c9_open_promoters)), 
             lwd = 0, cex = 1.5, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), cat.pos = c(340, 0), fill = c("#ffccd5", "#c9184a"),
             "atac_open_rna_up_a7_c9_overlaps.png", disable.logging = TRUE)


genes_to_plot <- intersect(intersect(a7_up, a7_open_promoters), intersect(c9_up, c9_open_promoters))


fisher.test(matrix(c(length(union(a7_up, a7_open_promoters))-length(union(intersect(a7_up, a7_open_promoters),intersect(c9_up, c9_open_promoters))), 
                     length(setdiff(intersect(a7_up, a7_open_promoters), intersect(c9_up, c9_open_promoters))), 
                     length(setdiff(intersect(c9_up, c9_open_promoters), intersect(a7_up, a7_open_promoters))), 
                     length(intersect(intersect(a7_up, a7_open_promoters), intersect(c9_up, c9_open_promoters)))), nrow = 2), alternative = "greater")


fisher.test(matrix(c(length(union(a7_down, a7_closed_promoters))-length(union(intersect(a7_down, a7_closed_promoters),intersect(c9_down, c9_closed_promoters))), 
                     length(setdiff(intersect(a7_down, a7_closed_promoters), intersect(c9_down, c9_closed_promoters))), 
                     length(setdiff(intersect(c9_down, c9_closed_promoters), intersect(a7_down, a7_closed_promoters))), 
                     length(intersect(intersect(a7_down, a7_closed_promoters), intersect(c9_down, c9_closed_promoters)))), nrow = 2), alternative = "greater")

library(pheatmap)
library(RColorBrewer)

normalised_expression_jurkat <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]

normalised_expression_jurkat <- all %>%
  dplyr::filter(geneID %in% genes_to_plot) %>%
  dplyr::filter(logFC > 5.5) #this filters out the many many entries to a more manageable 34

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))

normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")

#heatmap of selected genes
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))


venn.diagram(list("A7 RNA \u2193\nA7 promoter accessibilty \u2193" = intersect(a7_down, a7_closed_promoters), 
                  "C9 RNA \u2193\nC9 promoter accessibilty \u2193" = intersect(c9_down, c9_closed_promoters)), 
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), cat.pos = c(340, 20), fill = c("#a9d6e5", "#014f86"),
             "atac_closed_rna_down_a7_c9_overlaps.png", disable.logging = TRUE)


#-------from here:---------
genes_to_plot <- intersect(intersect(a7_down, a7_closed_promoters), intersect(c9_down, c9_closed_promoters))

normalised_expression_jurkat <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_clones.csv")
fold_changes <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_clones_all_genes (1).csv")

all <- left_join(normalised_expression_jurkat, fold_changes, by = "geneID")
all <- all[,-c(12:16)]

filter_negative_sums <- function(all) {
  # Calculate sums for specified columns
  sum1 <- rowMeans(all[, c("A10", "A14", "A15")])
  sum2 <- rowMeans(all[, c("A16", "A17", "A18")])
  
  # Filter rows where both sums are less than 0
  filtered_df <- all[sum1 < 0 & sum2 < 0, ]
  
  return(filtered_df)
}
filtered_df <- filter_negative_sums(all)

normalised_expression_jurkat <- filtered_df %>%
  dplyr::filter(geneID %in% genes_to_plot) 

rownames(normalised_expression_jurkat) <- normalised_expression_jurkat$geneID

normalised_expression_jurkat <- normalised_expression_jurkat %>%
  dplyr::select(-c("geneID", "logFC"))

normalised_expression_jurkat <- normalised_expression_jurkat[, c(2, 3, 4, 1, 5, 6, 7, 8, 9)]
display_labels <- c("WT", "WT", "WT", "A7", "A7", "A7", "C9", "C9", "C9")


#heatmap of selected genes
pheatmap(normalised_expression_jurkat, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20, labels_col = display_labels,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))



list_sets <- c("A7 RNA up", "A7 RNA down", 
                  "A7 promoter open", "A7 promoter closed",
                  "C9 RNA up", "C9 RNA down", 
                  "C9 promoter open", "C9 promoter closed")

list_sets <- factor(list_sets, levels=list_sets)


upset(fromList(list("A7 RNA up" = a7_up, "A7 RNA down" = a7_down, 
                    "A7 promoter open" = a7_open_promoters, "A7 promoter closed" = a7_closed_promoters,
                    "C9 RNA up" = c9_up, "C9 RNA down" = c9_down, 
                    "C9 promoter open" = c9_open_promoters, "C9 promoter closed" = c9_closed_promoters)),
      order.by = "freq", text.scale = 1.5, nsets = 8, nintersects = 55, sets = list_sets,
      queries = list(list(query = intersects, 
                          params = list("A7 RNA up", "A7 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA down", "A7 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA up", "C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA down", "C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed", "C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "C9 promoter open", "A7 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open", "A7 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "C9 promoter open", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed", "C9 promoter closed", "A7 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed", "C9 promoter closed", "A7 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "C9 promoter open", "A7 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter closed", "A7 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter closed", "A7 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open", "A7 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "A7 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("A7 promoter open", "C9 RNA up"), 
                          color = "#E68D74", active = T)),
      sets.bar.color = c("#E68D74", "#E68D74", "#6679AB", "#6679AB",
                         "#E68D74", "#6679AB", "#6679AB", "#E68D74"))

venn.diagram(list("A7 promoter\nopen" = a7_open_promoters, "C9 promoter\nopen" = c9_open_promoters), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "open_atac_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 promoter\nclosed" = a7_closed_promoters, "C9 promoter\nclosed" = c9_closed_promoters), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(200, 170), fill = c("#a9d6e5", "#014f86"),
             "closed_atac_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 promoter\nopen" = a7_open_promoters, "A7 RNA up" = a7_up), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "atac_open_rna_up_a7_overlaps.png", disable.logging = TRUE)

venn.diagram(list("A7 promoter\nclosed" = a7_closed_promoters, "A7 RNA down" = a7_down), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(195, 175), fill = c("#a9d6e5", "#014f86"),
             "atac_closed_rna_down_a7_overlaps.png", disable.logging = TRUE)


venn.diagram(list("C9 promoter\nopen" = c9_open_promoters, "C9 RNA up" = c9_up), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "atac_open_rna_up_c9_overlaps.png", disable.logging = TRUE)

venn.diagram(list("C9 promoter\nclosed" = c9_closed_promoters, "C9 RNA down" = c9_down), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(195, 175), fill = c("#a9d6e5", "#014f86"),
             "atac_closed_rna_down_c9_overlaps.png", disable.logging = TRUE)


####Overlap between peaks and RNA expression######### take the unchanged RNA expression for the DNAme list
wt_vs_a7_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_a7.csv")  

wt_vs_c9_diff_expr <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_c9.csv") 


a7_unchanged <- wt_vs_a7_diff_expr %>%
  dplyr::filter(abs(LogFC) < 0.5) %>%
  dplyr::filter(abs(LogFC) > 0) %>%
  dplyr::pull(geneID)

c9_unchanged <- wt_vs_c9_diff_expr %>%
  dplyr::filter(abs(LogFC) < 0.5) %>%
  dplyr::filter(abs(LogFC) > 0) %>%
  dplyr::pull(geneID)

wt_open_promoters <- WT_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)"| annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

a7_open_promoters <- A7_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)"| annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

a7_open <- setdiff(a7_open_promoters, wt_open_promoters)

c9_open_promoters <- C9_peaks_only_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)"| annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c9_open <- setdiff(c9_open_promoters, wt_open_promoters)

WT_minus_A7_annotated <- anti_join(WT_peaks_only_annotated, A7_peaks_only_annotated, by = "geneId")

a7_closed_promoters <- WT_minus_A7_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)" | annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()
a7_closed <- WT_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

#add to change thresholds # | annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)" 
WT_minus_C9_annotated <- anti_join(WT_peaks_only_annotated, C9_peaks_only_annotated, by = "geneId")
c9_closed_promoters <- WT_minus_C9_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)"| annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c9_closed <- WT_peaks_only_annotated %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()


dir <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\"
a7_open_promoter_a7_unchanged <- intersect(a7_open, a7_unchanged) %>%
  write.csv(file = paste0(dir, "\\a7_open_a7_unchanged.csv"), quote = FALSE, row.names = FALSE)
c9_open_promoter_c9_unchanged <- intersect(c9_open, c9_unchanged) %>%
  write.csv(file = paste0(dir, "\\c9_open_c9_unchanged.csv"), quote = FALSE, row.names = FALSE)
clone_open_clones_unchanged <- intersect(a7_open_promoter_a7_unchanged, c9_open_promoter_c9_unchanged)

a7_closed_promoter_a7_unchanged <- intersect(a7_closed_promoters, a7_unchanged)%>%
  write.csv(file = paste0(dir, "\\a7_closed_a7_unchanged.csv"), quote = FALSE, row.names = FALSE)
c9_closed_promoter_c9_unchanged <- intersect(c9_closed_promoters, c9_unchanged)%>%
  write.csv(file = paste0(dir, "\\c9_closed_c9_unchanged.csv"), quote = FALSE, row.names = FALSE)
clone_closed_clones_unchanged <- intersect(a7_closed_promoter_a7_unchanged, c9_closed_promoter_c9_unchanged)

#might have to quantify ATACseq for this analysis

