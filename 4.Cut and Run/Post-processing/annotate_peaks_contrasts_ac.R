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

dir <- "path\\macs2_peaks\\H3K27ac\\contrasts\\"

#load peak files
peak_files <- list(A7_C9_common_ac = paste0(dir, "A7_c9_ac_overlap_peaks.bed"),
                   WT_unique_ac = paste0(dir, "WT_ac_unique_peaks.bed"),
                   A7_unique_ac = paste0(dir, "A7_ac_unique_peaks.bed"),
                   C9_unique_ac = paste0(dir, "C9_ac_unique_peaks.bed"), 
                   WT_A7_C9_overlap_ac = paste0(dir, "wt_A7_c9_ac_overlap_peaks.bed"),
                   all_ac_peaks = paste0(dir, "all_ac_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_ac <- annotate_peaks(peak_files, "WT_unique_ac") 
annotated_peaks_A7_unique_ac <- annotate_peaks(peak_files, "A7_unique_ac")
annotated_peaks_C9_unique_ac <- annotate_peaks(peak_files, "C9_unique_ac")
annotated_peaks_A7_C9_common_ac <- annotate_peaks(peak_files, "A7_C9_common_ac") 
annotated_peaks_WT_A7_C9_overlap_ac <- annotate_peaks(peak_files, "WT_A7_C9_overlap_ac")
annotated_all_ac_peaks <- annotate_peaks(peak_files, "all_ac_peaks")

dir.create("path\\macs2_peaks\\H3K27ac\\contrasts\\annotated_peaks")

write.csv(annotated_peaks_WT_unique_ac, file = paste0(dir, ".\\annotated_peaks\\WT_unique_ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_unique_ac, file = paste0(dir, ".\\annotated_peaks\\A7_unique_ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_ac, file = paste0(dir, ".\\annotated_peaks\\C9_unique_ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_A7_C9_common_ac, file = paste0(dir, ".\\annotated_peaks\\A7_C9_common_ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_A7_C9_overlap_ac, file = paste0(dir, ".\\annotated_peaks\\WT_A7_C9_overlap_ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_all_ac_peaks, file = paste0(dir, ".\\annotated_peaks\\annotated_all_ac.csv"), quote = FALSE, row.names = FALSE)


#Count number of peaks in each type of annotation for H3K27ac peaks
WT_unique_ac_summary <- as.data.frame(table(annotated_peaks_WT_unique_ac$annotation))
colnames(WT_unique_ac_summary) <- c("Annotation", "WT_unique_ac")

A7_unique_ac_summary <- as.data.frame(table(annotated_peaks_A7_unique_ac$annotation))
colnames(A7_unique_ac_summary) <- c("Annotation", "A7_unique_ac")

C9_unique_ac_summary <- as.data.frame(table(annotated_peaks_C9_unique_ac$annotation))
colnames(C9_unique_ac_summary) <- c("Annotation", "C9_unique_ac")

A7_C9_common_ac_summary <- as.data.frame(table(annotated_peaks_A7_C9_common_ac$annotation))
colnames(A7_C9_common_ac_summary) <- c("Annotation", "A7_C9_common_ac")

WT_A7_C9_overlap_ac_summary <- as.data.frame(table(annotated_peaks_WT_A7_C9_overlap_ac$annotation))
colnames(WT_A7_C9_overlap_ac_summary) <- c("Annotation", "WT_A7_C9_overlap_ac")



summary_annotation_ac <- left_join(WT_unique_ac_summary, A7_unique_ac_summary, by = "Annotation") %>%
  left_join(C9_unique_ac_summary, by = "Annotation") %>%
  left_join(A7_C9_common_ac_summary, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_ac_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_ac$Annotation <- factor(summary_annotation_ac$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_ac$Condition <- factor(summary_annotation_ac$Condition, 
                                          levels = c("WT_A7_C9_overlap_ac", "A7_C9_common_ac", 
                                                     "C9_unique_ac", "A7_unique_ac", "WT_unique_ac"))

#plot
annotated_ac_peaks_plot <- summary_annotation_ac %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 30000)) +
  labs(x = "Sample") +
  ylab("No. H3K27ac peaks")


ggsave(dir, filename = ".\\annotated_peaks\\annotated_contrast_peaks_ac.pdf", plot = annotated_ac_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


WT_unique_ac_summary_perc <- WT_unique_ac_summary %>%
  mutate(WT_unique_ac_percentage = WT_unique_ac/sum(WT_unique_ac_summary$WT_unique_ac)*100) %>%
  dplyr::select(-WT_unique_ac) 

A7_unique_ac_summary_perc <- A7_unique_ac_summary %>%
  mutate(A7_unique_ac_percentage = A7_unique_ac/sum(A7_unique_ac_summary$A7_unique_ac)*100) %>%
  dplyr::select(-A7_unique_ac) 

C9_unique_ac_summary_perc <- C9_unique_ac_summary %>%
  mutate(C9_unique_ac_percentage = C9_unique_ac/sum(C9_unique_ac_summary$C9_unique_ac)*100) %>%
  dplyr::select(-C9_unique_ac) 

A7_C9_common_ac_summary_perc <- A7_C9_common_ac_summary %>%
  mutate(A7_C9_common_ac_percentage = A7_C9_common_ac/sum(A7_C9_common_ac_summary$A7_C9_common_ac)*100) %>%
  dplyr::select(-A7_C9_common_ac) 

WT_A7_C9_overlap_ac_summary_perc <- WT_A7_C9_overlap_ac_summary %>%
  mutate(WT_A7_C9_overlap_ac_percentage = WT_A7_C9_overlap_ac/sum(WT_A7_C9_overlap_ac_summary$WT_A7_C9_overlap_ac)*100) %>%
  dplyr::select(-WT_A7_C9_overlap_ac) 

summary_annotation_perc <- left_join(WT_unique_ac_summary_perc, A7_unique_ac_summary_perc, by = "Annotation") %>%
  left_join(C9_unique_ac_summary_perc, by = "Annotation") %>%
  left_join(A7_C9_common_ac_summary_perc, by = "Annotation") %>%
  left_join(WT_A7_C9_overlap_ac_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -c("Annotation"), names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))

summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition,
                                            levels = c("WT_A7_C9_overlap_ac_percentage",
                                                       "A7_C9_common_ac_percentage",
                                                       "C9_unique_ac_percentage",
                                                       "A7_unique_ac_percentage",
                                                       "WT_unique_ac_percentage"))

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
  ylab("H3K27ac peaks (%)") + 
  theme(panel.grid = element_blank()) 

ggsave(dir, filename = ".\\annotated_peaks\\annotated_constrast_ac_peaks_percentage.pdf", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)
