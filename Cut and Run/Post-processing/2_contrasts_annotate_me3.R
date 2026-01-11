#input files: peak files from upset_overlaps_me3.R since we have both unique peaks but also the overlapped peaks
#output files: .csv file with genomic annotation of the peaks in every condition for downstream analyses
#output files: Venn diagrams for intersect between Cut and run peaks, and also ATAC peaks

library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(Cairo)
library(colorspace)

#Function to annotate peak files
annotate_peaks <- function(peak_file_name, condition) {
  annotated_peaks <- annotatePeak(peak_file_name[[condition]], tssRegion = c(-3000, 3000),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
  annotated_peaks_df <- as.data.frame(annotated_peaks)
  annotated_peaks_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_df$annotation) #This merges all Intron/Exon as one annotation, independent of intron/exon number
  annotated_peaks_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_df$annotation)
  return(annotated_peaks_df)
}


#load annotations
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- ".\\macs2_peaks\\H3K27me3\\contrasts\\"

#load peak files
peak_files <- list(A7_C9_common_me3 = paste0(dir, "A7_c9_me3_overlap_peaks.bed"),
                   WT_unique_me3 = paste0(dir, "WT_me3_unique_peaks.bed"),
                   A7_unique_me3 = paste0(dir, "A7_me3_unique_peaks.bed"),
                   C9_unique_me3 = paste0(dir, "C9_me3_unique_peaks.bed"), 
                   WT_A7_C9_overlap_me3 = paste0(dir, "wt_A7_c9_me3_overlap_peaks.bed"),
                   all_me3_peaks = paste0(dir, "all_me3_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_me3 <- annotate_peaks(peak_files, "WT_unique_me3") 
annotated_peaks_A7_unique_me3 <- annotate_peaks(peak_files, "A7_unique_me3")
annotated_peaks_C9_unique_me3 <- annotate_peaks(peak_files, "C9_unique_me3")
annotated_peaks_A7_C9_common_me3 <- annotate_peaks(peak_files, "A7_C9_common_me3") 
annotated_peaks_WT_A7_C9_overlap_me3 <- annotate_peaks(peak_files, "WT_A7_C9_overlap_me3")
annotated_all_me3_peaks <- annotate_peaks(peak_files, "all_me3_peaks")

dir.create(".\\macs2_peaks\\H3K27me3\\contrasts\\annotated_peaks")

write.csv(annotated_peaks_WT_unique_me3, file = paste0(dir, ".\\annotated_peaks\\WT_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_unique_me3, file = paste0(dir, ".\\annotated_peaks\\A7_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_me3, file = paste0(dir, ".\\annotated_peaks\\C9_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_C9_common_me3, file = paste0(dir, ".\\annotated_peaks\\A7_C9_common_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_A7_C9_overlap_me3, file = paste0(dir, ".\\annotated_peaks\\WT_A7_C9_overlap_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_all_me3_peaks, file = paste0(dir, ".\\annotated_peaks\\annotated_all_me3.csv"), quote = FALSE, row.names = FALSE)


#Count number of peaks in each type of annotation for H3K27me3 peaks
WT_unique_me3_summary <- as.data.frame(table(annotated_peaks_WT_unique_me3$annotation))
colnames(WT_unique_me3_summary) <- c("Annotation", "WT_unique_me3")

A7_unique_me3_summary <- as.data.frame(table(annotated_peaks_A7_unique_me3$annotation))
colnames(A7_unique_me3_summary) <- c("Annotation", "A7_unique_me3")

C9_unique_me3_summary <- as.data.frame(table(annotated_peaks_C9_unique_me3$annotation))
colnames(C9_unique_me3_summary) <- c("Annotation", "C9_unique_me3")

A7_C9_common_me3_summary <- as.data.frame(table(annotated_peaks_A7_C9_common_me3$annotation))
colnames(A7_C9_common_me3_summary) <- c("Annotation", "A7_C9_common_me3")

WT_A7_C9_overlap_me3_summary <- as.data.frame(table(annotated_peaks_WT_A7_C9_overlap_me3$annotation))
colnames(WT_A7_C9_overlap_me3_summary) <- c("Annotation", "WT_A7_C9_overlap_me3")

summary_annotation_me3 <- left_join(WT_unique_me3_summary, A7_unique_me3_summary, by = "Annotation") %>%
  left_join(C9_unique_me3_summary, by = "Annotation") %>%
  #left_join(A7_C9_common_me3_summary, by = "Annotation") %>%
  #left_join(WT_A7_C9_overlap_me3_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_me3$Annotation <- factor(summary_annotation_me3$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_me3$Condition <- factor(summary_annotation_me3$Condition, 
                                          levels = c(#"WT_A7_C9_overlap_me3", "A7_C9_common_me3", 
                                                     "C9_unique_me3", "A7_unique_me3", "WT_unique_me3"))

#plot
annotated_me3_peaks_plot <- summary_annotation_me3 %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 60000)) +  #change the limits based on the larger category for the x-axis
  labs(x = "Sample") +
  ylab("No. H3K27me3 peaks")


ggsave(dir, filename = ".\\annotated_peaks\\annotated_contrast_peaks_me3_unique.pdf", plot = annotated_me3_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


WT_unique_me3_summary_perc <- WT_unique_me3_summary %>%
  mutate(WT_unique_me3_percentage = WT_unique_me3/sum(WT_unique_me3_summary$WT_unique_me3)*100) %>%
  dplyr::select(-WT_unique_me3) 

A7_unique_me3_summary_perc <- A7_unique_me3_summary %>%
  mutate(A7_unique_me3_percentage = A7_unique_me3/sum(A7_unique_me3_summary$A7_unique_me3)*100) %>%
  dplyr::select(-A7_unique_me3) 

C9_unique_me3_summary_perc <- C9_unique_me3_summary %>%
  mutate(C9_unique_me3_percentage = C9_unique_me3/sum(C9_unique_me3_summary$C9_unique_me3)*100) %>%
  dplyr::select(-C9_unique_me3) 

A7_C9_common_me3_summary_perc <- A7_C9_common_me3_summary %>%
  mutate(A7_C9_common_me3_percentage = A7_C9_common_me3/sum(A7_C9_common_me3_summary$A7_C9_common_me3)*100) %>%
  dplyr::select(-A7_C9_common_me3) 

WT_A7_C9_overlap_me3_summary_perc <- WT_A7_C9_overlap_me3_summary %>%
  mutate(WT_A7_C9_overlap_me3_percentage = WT_A7_C9_overlap_me3/sum(WT_A7_C9_overlap_me3_summary$WT_A7_C9_overlap_me3)*100) %>%
  dplyr::select(-WT_A7_C9_overlap_me3) 

summary_annotation_perc <- left_join(WT_unique_me3_summary_perc, A7_unique_me3_summary_perc, by = "Annotation") %>%
  left_join(C9_unique_me3_summary_perc, by = "Annotation") %>%
  left_join(A7_C9_common_me3_summary_perc, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_me3_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -c("Annotation"), names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))

summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition,
                                            levels = c("WT_A7_C9_overlap_me3_percentage",
                                                       "A7_C9_common_me3_percentage",
                                                       "C9_unique_me3_percentage",
                                                       "A7_unique_me3_percentage",
                                                       "WT_unique_me3_percentage"))

annotated_peaks_plot_perc <- summary_annotation_perc %>%
  ggplot(aes(x = Condition, y = percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00",  
                               "#FEB95F", "#CC79A7", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5)), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  #scale_x_discrete(labels=c("WT_percentage" = "WT", "C9_percentage" = "C9")) +
  labs(x = "Sample") +
  ylab("H3K27me3 peaks (%)") + 
  theme(panel.grid = element_blank()) 

ggsave(dir, filename = ".\\annotated_peaks\\annotated_constrast_me3_peaks_percentage.pdf", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)

#annotate all peaks found in the 3 cell lines to see if there are any peaks in promoter regions that overlap and to get Venn diagrams
dir1 <- ".\\macs2_peaks\\H3K27me3\\"
dir2 <- ".\\macs2_peaks\\H2AK119ub\\"
dir3 <- ".\\macs2_peaks\\H3K27ac\\"

#load peak files
peak_files <- list(WT_all_me3 = paste0(dir1, "WT_me3_peaks.bed"),
                   A7_all_me3 = paste0(dir1, "A7_me3_peaks.bed"),
                   C9_all_me3 = paste0(dir1, "C9_me3_peaks.bed"),
                   WT_all_ub = paste0(dir2, "WT_ub_peaks.bed"),
                   A7_all_ub = paste0(dir2, "A7_ub_peaks.bed"),
                   C9_all_ub = paste0(dir2, "C9_ub_peaks.bed"),
                   WT_all_ac = paste0(dir3, "WT_ac_peaks.bed"),
                   A7_all_ac = paste0(dir3, "A7_ac_peaks.bed"),
                   C9_all_ac = paste0(dir3, "C9_ac_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_all_me3 <- annotate_peaks(peak_files, "WT_all_me3") 
annotated_peaks_A7_all_me3 <- annotate_peaks(peak_files, "A7_all_me3")
annotated_peaks_C9_all_me3 <- annotate_peaks(peak_files, "C9_all_me3")

annotated_peaks_WT_all_ub <- annotate_peaks(peak_files, "WT_all_ub") 
annotated_peaks_A7_all_ub <- annotate_peaks(peak_files, "A7_all_ub")
annotated_peaks_C9_all_ub <- annotate_peaks(peak_files, "C9_all_ub")

annotated_peaks_WT_all_ac <- annotate_peaks(peak_files, "WT_all_ac") 
annotated_peaks_A7_all_ac <- annotate_peaks(peak_files, "A7_all_ac")
annotated_peaks_C9_all_ac <- annotate_peaks(peak_files, "C9_all_ac")


WT_me3_promoters <- annotated_peaks_WT_all_me3 %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit

A7_me3_promoters <- annotated_peaks_A7_all_me3 %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit


C9_me3_promoters <- annotated_peaks_C9_all_me3 %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit

WT_ub_promoters <- annotated_peaks_WT_all_ub %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit

A7_ub_promoters <- annotated_peaks_A7_all_ub %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit


C9_ub_promoters <- annotated_peaks_C9_all_ub %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit

WT_ac_promoters <- annotated_peaks_WT_all_ac %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit

A7_ac_promoters <- annotated_peaks_A7_all_ac %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  na.omit


C9_ac_promoters <- annotated_peaks_C9_all_ac %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit



library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

venn.diagram(list("WT H3K27me3\n promoters" = WT_me3_promoters,
                  "A7 H3K27me3\n promoters" = A7_me3_promoters,
                  "C9 H3K27me3\n promoters" = C9_me3_promoters),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#173518", "#2e6930", "#449e48"),
             ".\\H3K27me3_overlaps.png", disable.logging = TRUE)

venn.diagram(list("WT H2AK119ub\n promoters" = WT_ub_promoters,
                  "A7 H2AK119ub\n promoters" = A7_ub_promoters,
                  "C9 H2AK119ub\n promoters" = C9_ub_promoters),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#4a2574", "#924ed1", "#bc83bd"),
             ".\\H2AK119ub_promoters_overlaps.png", disable.logging = TRUE)


venn.diagram(list("WT H3K27me3\n promoters" = WT_me3_promoters,
                  "WT H2AK119ub\n promoters" = WT_ub_promoters),
             lwd = 0, cex = 2, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#173518", "#4a2574"),
             ".\\WT_H3K27me3_H2AK119ub_promoters_overlaps.png", disable.logging = TRUE)


venn.diagram(list("WT H3K27ac\n promoters" = WT_ac_promoters,
                  "A7 H3K27ac\n promoters" = A7_ac_promoters,
                  "C9 H3K27ac\n promoters" = C9_ac_promoters),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#0D47A1", "#42A5F5", "#4FC3F7"),
             ".\H3K27ac_promoters_overlaps.png", disable.logging = TRUE)



WT_ATAC_peaks <- read.csv(".\\ATAC-peaks\\WT_unique_atac.csv") %>%
  na.omit
A7_ATAC_peaks <- read.csv(".\\ATAC-peaks\\A7_unique_atac.csv") %>%
  na.omit
C9_ATAC_peaks <- read.csv(".\\ATAC-peaks\\C9_unique_atac.csv") %>%
  na.omit

WT_ac_ATAC_peaks <- intersect(WT_unique_ac_promoters, WT_ATAC_peaks$SYMBOL)
A7_ac_ATAC_peaks <- intersect(A7_unique_ac_promoters, A7_ATAC_peaks$SYMBOL)
C9_ac_ATAC_peaks <- intersect(C9_unique_ac_promoters, C9_ATAC_peaks$SYMBOL)

wt_vs_clones_diff_expr <- read.csv(".\\Jurkat_wt_v_clones_all_genes (1).csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))
rna_up_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  pull(geneID)

rna_down_genes <- wt_vs_clones_diff_expr %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  pull(geneID)

venn.diagram(list("A7 ATAC peaks \u2191" = A7_ATAC_peaks$SYMBOL,
                  "A7 H3K27ac \u2191" = A7_unique_ac_promoters,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#e27602", "#42A5F5", "#E16036"),
             ".\\rna_up_A7_ac_A7_ATAC_overlaps.png", disable.logging = TRUE)

venn.diagram(list("C9 ATAC peaks \u2191" = C9_ATAC_peaks$SYMBOL,
                  "C9 H3K27ac \u2191" = C9_unique_ac_promoters,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#f1b04c", "#4FC3F7", "#E16036"),
             ".\\rna_up_C9_ac_C9_ATAC_overlaps.png", disable.logging = TRUE)

##now intersect H3K27me3 promoter peaks with ATAC peaks
WT_unique_me3_promoters <- read.csv(".\\macs2_peaks\\H3K27me3\\contrasts\\annotated_peaks\\WT_unique_me3.csv") %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit

A7_unique_me3_promoters <- read.csv(".\\macs2_peaks\\H3K27me3\\contrasts\\annotated_peaks\\A7_unique_me3.csv") %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit

C9_unique_me3_promoters <- read.csv(".\\macs2_peaks\\H3K27me3\\contrasts\\annotated_peaks\\C9_unique_me3.csv") %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  dplyr::pull(SYMBOL) %>%
  unique()%>%
  na.omit

venn.diagram(list("A7 ATAC peaks \u2191" = A7_ATAC_peaks$SYMBOL,
                  "WT H3K27me3 \u2193" = WT_unique_me3_promoters,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#e27602", "#173518", "#E16036"),
             ".\\rna_up_WT_me3_unique_A7_ATAC_overlaps.png", disable.logging = TRUE)

venn.diagram(list("C9 ATAC peaks \u2191" = C9_ATAC_peaks$SYMBOL,
                  "WT H3K27me3 \u2193" = WT_unique_me3_promoters,
                  "A7&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#f1b04c", "#173518", "#E16036"),
             ".\\rna_up_WT_me3_unique_C9_ATAC_overlaps.png", disable.logging = TRUE)

a <- intersect(intersect(A7_ATAC_peaks$SYMBOL,WT_unique_me3_promoters), rna_up_genes)
b <- intersect(intersect(C9_ATAC_peaks$SYMBOL,WT_unique_me3_promoters), rna_up_genes)


