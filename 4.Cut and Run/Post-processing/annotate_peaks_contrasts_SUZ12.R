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

dir <- "path\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\"

#load peak files
peak_files <- list(A7_C9_common_SUZ12 = paste0(dir, "a7_c9_SUZ12_overlap_peaks.bed"),
                   WT_unique_SUZ12 = paste0(dir, "WT_SUZ12_unique_peaks.bed"),
                   A7_unique_SUZ12 = paste0(dir, "A7_SUZ12_unique_peaks.bed"),
                   C9_unique_SUZ12 = paste0(dir, "C9_SUZ12_unique_peaks.bed"), 
                   WT_A7_C9_overlap_SUZ12 = paste0(dir, "wt_a7_c9_SUZ12_overlap_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_SUZ12 <- annotate_peaks(peak_files, "WT_unique_SUZ12") 
annotated_peaks_A7_unique_SUZ12 <- annotate_peaks(peak_files, "A7_unique_SUZ12")
annotated_peaks_C9_unique_SUZ12 <- annotate_peaks(peak_files, "C9_unique_SUZ12")
annotated_peaks_A7_C9_common_SUZ12 <- annotate_peaks(peak_files, "A7_C9_common_SUZ12") 
annotated_peaks_WT_A7_C9_overlap_SUZ12 <- annotate_peaks(peak_files, "WT_A7_C9_overlap_SUZ12")

dir.create("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\annotated_peaks")

write.csv(annotated_peaks_WT_unique_SUZ12, file = paste0(dir, ".\\annotated_peaks\\WT_unique_SUZ12.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_unique_SUZ12, file = paste0(dir, ".\\annotated_peaks\\A7_unique_SUZ12.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_SUZ12, file = paste0(dir, ".\\annotated_peaks\\C9_unique_SUZ12.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_C9_common_SUZ12, file = paste0(dir, ".\\annotated_peaks\\A7_C9_common_SUZ12.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_A7_C9_overlap_SUZ12, file = paste0(dir, ".\\annotated_peaks\\WT_A7_C9_overlap_SUZ12.csv"), quote = FALSE, row.names = FALSE)


#Count number of peaks in each type of annotation for SUZ12 peaks
WT_unique_SUZ12_summary <- as.data.frame(table(annotated_peaks_WT_unique_SUZ12$annotation))
colnames(WT_unique_SUZ12_summary) <- c("Annotation", "WT_unique_SUZ12")

A7_unique_SUZ12_summary <- as.data.frame(table(annotated_peaks_A7_unique_SUZ12$annotation))
colnames(A7_unique_SUZ12_summary) <- c("Annotation", "A7_unique_SUZ12")

C9_unique_SUZ12_summary <- as.data.frame(table(annotated_peaks_C9_unique_SUZ12$annotation))
colnames(C9_unique_SUZ12_summary) <- c("Annotation", "C9_unique_SUZ12")

A7_C9_common_SUZ12_summary <- as.data.frame(table(annotated_peaks_A7_C9_common_SUZ12$annotation))
colnames(A7_C9_common_SUZ12_summary) <- c("Annotation", "A7_C9_common_SUZ12")

WT_A7_C9_overlap_SUZ12_summary <- as.data.frame(table(annotated_peaks_WT_A7_C9_overlap_SUZ12$annotation))
colnames(WT_A7_C9_overlap_SUZ12_summary) <- c("Annotation", "WT_A7_C9_overlap_SUZ12")


summary_annotation_SUZ12 <- left_join(WT_unique_SUZ12_summary, A7_unique_SUZ12_summary, by = "Annotation") %>%
  left_join(C9_unique_SUZ12_summary, by = "Annotation") %>%
  left_join(A7_C9_common_SUZ12_summary, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_SUZ12_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_SUZ12$Annotation <- factor(summary_annotation_SUZ12$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_SUZ12$Condition <- factor(summary_annotation_SUZ12$Condition, 
                                          levels = c("WT_A7_C9_overlap_SUZ12", "A7_C9_common_SUZ12", 
                                                     "C9_unique_SUZ12", "A7_unique_SUZ12", "WT_unique_SUZ12"))

#plot
annotated_SUZ12_peaks_plot <- summary_annotation_SUZ12 %>%
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
  ylab("No. SUZ12 peaks")


ggsave(dir, filename = ".\\annotated_peaks\\annotated_contrast_peaks_SUZ12.pdf", plot = annotated_SUZ12_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


