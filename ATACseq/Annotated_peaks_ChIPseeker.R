#input files: bed files from peaks for each condition
#output: annotated genomic regions for each bed file
#output: annotated genomic files for common peaks too

library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(colorspace)

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "./ATAC_Seq/"

#peak .bed files are created with HMMRATAC in WSL where we also used the bedtools package to compute overlapped peaks
constrasts <- c("Jurkat_A7_peaks_only", "Jurkat_C9_peaks_only", "Jurkat_WT_peaks_only", 
                "Jurkat_WT_A7_C9_peaks_overlap", "Jurkat_A7_C9_overlap_peaks")

contrasts_peak_files <- list(Jurkat_A7_peaks_only = paste0(dir, "Jurkat_A7_peaks_only.bed"), 
                             Jurkat_C9_peaks_only = paste0(dir, "Jurkat_C9_peaks_only.bed"), 
                             Jurkat_WT_peaks_only = paste0(dir, "Jurkat_WT_peaks_only.bed"),
                             Jurkat_WT_A7_C9_peaks_overlap = paste0(dir, "Jurkat_WT_A7_C9_peaks_overlap.bed"),
                             Jurkat_A7_C9_overlap_peaks = paste0(dir, "Jurkat_A7_C9_overlap_peaks.bed"))

#]function annotates all bed files indepenedntly without having to run them one by one (can be rewritten to run one by one)
annotate_peaks_func <- function(peak_files_list, chosen_df)  {
  annotated_peaks <- annotatePeak(peak_files_list[[chosen_df]], tssRegion = c(-3000, 3000),
                                  TxDb = edb, annoDb = "org.Hs.eg.db") %>%
    as.data.frame()
  
  annotated_peaks$annotation <- gsub("^Intron.*", "Intron", annotated_peaks$annotation)
  annotated_peaks$annotation <- gsub("^Exon.*", "Exon", annotated_peaks$annotation)
  return(annotated_peaks)
}

Jurkat_A7_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "Jurkat_A7_peaks_only")
Jurkat_C9_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "Jurkat_C9_peaks_only")
Jurkat_WT_peaks_only_annotated <- annotate_peaks_func(contrasts_peak_files, "Jurkat_WT_peaks_only")
Jurkat_WT_A7_C9_peaks_overlap_annotated <- annotate_peaks_func(contrasts_peak_files, "Jurkat_WT_A7_C9_peaks_overlap")
Jurkat_A7_C9_overlap_peaks_annotated <- annotate_peaks_func(contrasts_peak_files, "Jurkat_A7_C9_overlap_peaks")

#save annotated peaks
create.dir("./ATAC_Seq/annotated_peaks")
write.csv(Jurkat_WT_peaks_only_annotated, file = paste0(dir, "./annotated_peaks/WT_unique_atac.csv"), quote = FALSE, row.names = FALSE)
write.csv(Jurkat_A7_peaks_only_annotated, file = paste0(dir, "./annotated_peaks/A7_unique_atac.csv"), quote = FALSE, row.names = FALSE)
write.csv(Jurkat_C9_peaks_only_annotated, file = paste0(dir, "./annotated_peaks/C9_unique_atac.csv"), quote = FALSE, row.names = FALSE)
write.csv(Jurkat_A7_C9_overlap_peaks_annotated, file = paste0(dir, "./annotated_peaks/A7_C9_common_atac.csv"), quote = FALSE, row.names = FALSE)
write.csv(Jurkat_WT_A7_C9_peaks_overlap_annotated, file = paste0(dir, "./annotated_peaks/WT_A7_C9_overlap_atac.csv"), quote = FALSE, row.names = FALSE)


Jurkat_A7_peaks_only_annotated_summary <- as.data.frame(table(Jurkat_A7_peaks_only_annotated$annotation))
colnames(Jurkat_A7_peaks_only_annotated_summary) <- c("Annotation", "Jurkat_A7_peaks_only")

Jurkat_C9_peaks_only_annotated_summary <- as.data.frame(table(Jurkat_C9_peaks_only_annotated$annotation))
colnames(Jurkat_C9_peaks_only_annotated_summary) <- c("Annotation", "Jurkat_C9_peaks_only")

Jurkat_WT_peaks_only_annotated_summary <- as.data.frame(table(Jurkat_WT_peaks_only_annotated$annotation))
colnames(Jurkat_WT_peaks_only_annotated_summary) <- c("Annotation", "Jurkat_WT_peaks_only")

Jurkat_WT_A7_C9_peaks_overlap_annotated_summary <- as.data.frame(table(Jurkat_WT_A7_C9_peaks_overlap_annotated$annotation))
colnames(Jurkat_WT_A7_C9_peaks_overlap_annotated_summary) <- c("Annotation", "Jurkat_WT_A7_C9_peaks_overlap")

Jurkat_A7_C9_overlap_peaks_annotated_summary <- as.data.frame(table(Jurkat_A7_C9_overlap_peaks_annotated$annotation))
colnames(Jurkat_A7_C9_overlap_peaks_annotated_summary) <- c("Annotation", "Jurkat_A7_C9_overlap_peaks")


Jurkat_A7_peaks_only_annotated_summary_perc <- Jurkat_A7_peaks_only_annotated_summary %>%
  mutate(A7_percentage = Jurkat_A7_peaks_only/sum(Jurkat_A7_peaks_only_annotated_summary$Jurkat_A7_peaks_only)*100) %>%
  select(-Jurkat_A7_peaks_only)

Jurkat_C9_peaks_only_annotated_summary_perc <- Jurkat_C9_peaks_only_annotated_summary %>%
  mutate(C9_percentage = Jurkat_C9_peaks_only/sum(Jurkat_C9_peaks_only_annotated_summary$Jurkat_C9_peaks_only)*100) %>%
  select(-Jurkat_C9_peaks_only)

Jurkat_WT_peaks_only_annotated_summary_perc <- Jurkat_WT_peaks_only_annotated_summary %>%
  mutate(WT_percentage = Jurkat_WT_peaks_only/sum(Jurkat_WT_peaks_only_annotated_summary$Jurkat_WT_peaks_only)*100) %>%
  select(-Jurkat_WT_peaks_only)

Jurkat_WT_A7_C9_peaks_overlap_annotated_summary_perc <- Jurkat_WT_A7_C9_peaks_overlap_annotated_summary %>%
  mutate(WT_A7_C9_percentage = Jurkat_WT_A7_C9_peaks_overlap/sum(Jurkat_WT_A7_C9_peaks_overlap_annotated_summary$Jurkat_WT_A7_C9_peaks_overlap)*100) %>%
  select(-Jurkat_WT_A7_C9_peaks_overlap)

Jurkat_A7_C9_overlap_peaks_annotated_summary_perc <- Jurkat_A7_C9_overlap_peaks_annotated_summary %>%
  mutate(A7_C9_percentage = Jurkat_A7_C9_overlap_peaks/sum(Jurkat_A7_C9_overlap_peaks_annotated_summary$Jurkat_A7_C9_overlap_peaks)*100) %>%
  select(-Jurkat_A7_C9_overlap_peaks)

summary_annotation_perc <- left_join(Jurkat_A7_peaks_only_annotated_summary_perc, Jurkat_C9_peaks_only_annotated_summary_perc, by = "Annotation") %>%
  left_join(Jurkat_WT_peaks_only_annotated_summary_perc, by = "Annotation") %>%
  left_join(Jurkat_WT_A7_C9_peaks_overlap_annotated_summary_perc, by = "Annotation") %>%
  left_join(Jurkat_A7_C9_overlap_peaks_annotated_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))
summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition, 
                                            levels = c("WT_A7_C9_percentage", "A7_C9_percentage",
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

ggsave(filename = "./ATAC_Seq/plots/annotated_constrast_atac_peaks_percentage.pdf", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)