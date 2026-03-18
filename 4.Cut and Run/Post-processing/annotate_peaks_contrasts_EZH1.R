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

dir <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\EZH1\\EZH1_filtered\\contrasts\\"

#load peak files
peak_files <- list(A7_C9_common_EZH1 = paste0(dir, "a7_c9_EZH1_overlap_peaks.bed"),
                   WT_unique_EZH1 = paste0(dir, "WT_EZH1_unique_peaks.bed"),
                   A7_unique_EZH1 = paste0(dir, "A7_EZH1_unique_peaks.bed"),
                   C9_unique_EZH1 = paste0(dir, "C9_EZH1_unique_peaks.bed"), 
                   WT_A7_C9_overlap_EZH1 = paste0(dir, "wt_a7_c9_EZH1_overlap_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_EZH1 <- annotate_peaks(peak_files, "WT_unique_EZH1") 
annotated_peaks_A7_unique_EZH1 <- annotate_peaks(peak_files, "A7_unique_EZH1")
annotated_peaks_C9_unique_EZH1 <- annotate_peaks(peak_files, "C9_unique_EZH1")
annotated_peaks_A7_C9_common_EZH1 <- annotate_peaks(peak_files, "A7_C9_common_EZH1") 
annotated_peaks_WT_A7_C9_overlap_EZH1 <- annotate_peaks(peak_files, "WT_A7_C9_overlap_EZH1")

dir.create("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\EZH1\\EZH1_filtered\\contrasts\\annotated_peaks")

write.csv(annotated_peaks_WT_unique_EZH1, file = paste0(dir, ".\\annotated_peaks\\WT_unique_EZH1.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_unique_EZH1, file = paste0(dir, ".\\annotated_peaks\\A7_unique_EZH1.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_EZH1, file = paste0(dir, ".\\annotated_peaks\\C9_unique_EZH1.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_C9_common_EZH1, file = paste0(dir, ".\\annotated_peaks\\A7_C9_common_EZH1.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_A7_C9_overlap_EZH1, file = paste0(dir, ".\\annotated_peaks\\WT_A7_C9_overlap_EZH1.csv"), quote = FALSE, row.names = FALSE)


#Count number of peaks in each type of annotation for EZH1 peaks
WT_unique_EZH1_summary <- as.data.frame(table(annotated_peaks_WT_unique_EZH1$annotation))
colnames(WT_unique_EZH1_summary) <- c("Annotation", "WT_unique_EZH1")

A7_unique_EZH1_summary <- as.data.frame(table(annotated_peaks_A7_unique_EZH1$annotation))
colnames(A7_unique_EZH1_summary) <- c("Annotation", "A7_unique_EZH1")

C9_unique_EZH1_summary <- as.data.frame(table(annotated_peaks_C9_unique_EZH1$annotation))
colnames(C9_unique_EZH1_summary) <- c("Annotation", "C9_unique_EZH1")

A7_C9_common_EZH1_summary <- as.data.frame(table(annotated_peaks_A7_C9_common_EZH1$annotation))
colnames(A7_C9_common_EZH1_summary) <- c("Annotation", "A7_C9_common_EZH1")

WT_A7_C9_overlap_EZH1_summary <- as.data.frame(table(annotated_peaks_WT_A7_C9_overlap_EZH1$annotation))
colnames(WT_A7_C9_overlap_EZH1_summary) <- c("Annotation", "WT_A7_C9_overlap_EZH1")


summary_annotation_EZH1 <- left_join(WT_unique_EZH1_summary, A7_unique_EZH1_summary, by = "Annotation") %>%
  left_join(C9_unique_EZH1_summary, by = "Annotation") %>%
  left_join(A7_C9_common_EZH1_summary, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_EZH1_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_EZH1$Annotation <- factor(summary_annotation_EZH1$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_EZH1$Condition <- factor(summary_annotation_EZH1$Condition, 
                                          levels = c("WT_A7_C9_overlap_EZH1", "A7_C9_common_EZH1", 
                                                     "C9_unique_EZH1", "A7_unique_EZH1", "WT_unique_EZH1"))

#plot
annotated_EZH1_peaks_plot <- summary_annotation_EZH1 %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 8000)) +
  labs(x = "Sample") +
  ylab("No. EH1 peaks")


ggsave(dir, filename = ".\\annotated_peaks\\annotated_contrast_peaks_EZH1.pdf", plot = annotated_EZH1_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


