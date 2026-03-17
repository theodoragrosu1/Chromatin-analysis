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

dir <- "path\\macs2_peaks\\H2AK119ub\\contrasts\\"

#load peak files
peak_files <- list(A7_C9_common_ub = paste0(dir, "A7_c9_ub_overlap_peaks.bed"),
                   WT_unique_ub = paste0(dir, "WT_ub_unique_peaks.bed"),
                   A7_unique_ub = paste0(dir, "A7_ub_unique_peaks.bed"),
                   C9_unique_ub = paste0(dir, "C9_ub_unique_peaks.bed"), 
                   WT_A7_C9_overlap_ub = paste0(dir, "wt_A7_c9_ub_overlap_peaks.bed"),
                   all_ub_peaks = paste0(dir, "all_ub_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_ub <- annotate_peaks(peak_files, "WT_unique_ub") 
annotated_peaks_A7_unique_ub <- annotate_peaks(peak_files, "A7_unique_ub")
annotated_peaks_C9_unique_ub <- annotate_peaks(peak_files, "C9_unique_ub")
annotated_peaks_A7_C9_common_ub <- annotate_peaks(peak_files, "A7_C9_common_ub") 
annotated_peaks_WT_A7_C9_overlap_ub <- annotate_peaks(peak_files, "WT_A7_C9_overlap_ub")
annotated_peaks_all_ub <- annotate_peaks(peak_files, "all_ub_peaks")

dir.create("path\\macs2_peaks\\H2AK119ub\\contrasts\\annotated_peaks")

write.csv(annotated_peaks_WT_unique_ub, file = paste0(dir, ".\\annotated_peaks\\WT_unique_ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_unique_ub, file = paste0(dir, ".\\annotated_peaks\\A7_unique_ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_ub, file = paste0(dir, ".\\annotated_peaks\\C9_unique_ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_C9_common_ub, file = paste0(dir, ".\\annotated_peaks\\A7_C9_common_ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_A7_C9_overlap_ub, file = paste0(dir, ".\\annotated_peaks\\WT_A7_C9_overlap_ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_all_ub, file = paste0(dir, ".\\annotated_peaks\\all_ub_peaks.csv"), quote = FALSE, row.names = FALSE)

#Count number of peaks in each type of annotation for H2AK119ub peaks
WT_unique_ub_summary <- as.data.frame(table(annotated_peaks_WT_unique_ub$annotation))
colnames(WT_unique_ub_summary) <- c("Annotation", "WT_unique_ub")

A7_unique_ub_summary <- as.data.frame(table(annotated_peaks_A7_unique_ub$annotation))
colnames(A7_unique_ub_summary) <- c("Annotation", "A7_unique_ub")

C9_unique_ub_summary <- as.data.frame(table(annotated_peaks_C9_unique_ub$annotation))
colnames(C9_unique_ub_summary) <- c("Annotation", "C9_unique_ub")

A7_C9_common_ub_summary <- as.data.frame(table(annotated_peaks_A7_C9_common_ub$annotation))
colnames(A7_C9_common_ub_summary) <- c("Annotation", "A7_C9_common_ub")

WT_A7_C9_overlap_ub_summary <- as.data.frame(table(annotated_peaks_WT_A7_C9_overlap_ub$annotation))
colnames(WT_A7_C9_overlap_ub_summary) <- c("Annotation", "WT_A7_C9_overlap_ub")

summary_annotation_ub <- left_join(WT_unique_ub_summary, A7_unique_ub_summary, by = "Annotation") %>%
  left_join(C9_unique_ub_summary, by = "Annotation") %>%
  #left_join(A7_C9_common_ub_summary, by = "Annotation") %>%
  #left_join(WT_A7_C9_overlap_ub_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_ub$Annotation <- factor(summary_annotation_ub$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_ub$Condition <- factor(summary_annotation_ub$Condition, 
                                          levels = c(#"WT_A7_C9_overlap_ub", "A7_C9_common_ub", 
                                                     "C9_unique_ub", "A7_unique_ub", "WT_unique_ub"))

#plot
annotated_ub_peaks_plot <- summary_annotation_ub %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 12000)) +
  labs(x = "Sample") +
  ylab("No. H2AK119ub peaks")


ggsave(dir, filename = ".\\annotated_contrast_peaks_ub.pdf", plot = annotated_ub_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


WT_unique_ub_summary_perc <- WT_unique_ub_summary %>%
  mutate(WT_unique_ub_percentage = WT_unique_ub/sum(WT_unique_ub_summary$WT_unique_ub)*100) %>%
  dplyr::select(-WT_unique_ub) 

A7_unique_ub_summary_perc <- A7_unique_ub_summary %>%
  mutate(A7_unique_ub_percentage = A7_unique_ub/sum(A7_unique_ub_summary$A7_unique_ub)*100) %>%
  dplyr::select(-A7_unique_ub) 

C9_unique_ub_summary_perc <- C9_unique_ub_summary %>%
  mutate(C9_unique_ub_percentage = C9_unique_ub/sum(C9_unique_ub_summary$C9_unique_ub)*100) %>%
  dplyr::select(-C9_unique_ub) 

A7_C9_common_ub_summary_perc <- A7_C9_common_ub_summary %>%
  mutate(A7_C9_common_ub_percentage = A7_C9_common_ub/sum(A7_C9_common_ub_summary$A7_C9_common_ub)*100) %>%
  dplyr::select(-A7_C9_common_ub) 

WT_A7_C9_overlap_ub_summary_perc <- WT_A7_C9_overlap_ub_summary %>%
  mutate(WT_A7_C9_overlap_ub_percentage = WT_A7_C9_overlap_ub/sum(WT_A7_C9_overlap_ub_summary$WT_A7_C9_overlap_ub)*100) %>%
  dplyr::select(-WT_A7_C9_overlap_ub) 

summary_annotation_perc <- left_join(WT_unique_ub_summary_perc, A7_unique_ub_summary_perc, by = "Annotation") %>%
  left_join(C9_unique_ub_summary_perc, by = "Annotation") %>%
  left_join(A7_C9_common_ub_summary_perc, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_ub_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -c("Annotation"), names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))

summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition,
                                            levels = c("WT_A7_C9_overlap_ub_percentage",
                                                       "A7_C9_common_ub_percentage",
                                                       "C9_unique_ub_percentage",
                                                       "A7_unique_ub_percentage",
                                                       "WT_unique_ub_percentage"))

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
  ylab("H2AK119ub peaks (%)") + 
  theme(panel.grid = element_blank()) 

ggsave(dir, filename = ".\\annotated_constrast_ub_peaks_percentage.pdf", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)
