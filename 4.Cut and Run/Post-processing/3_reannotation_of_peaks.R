#input files: unique WT peaks as given by 1_upset_overlaps_me3.R and then A7 and C9 peaks that also include the overlap between themselves, to not miss out on any shared peaks
#output files: annotated peaks for all conditions

library(dplyr)
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
peak_files <- list(WT_me3_peaks = paste0(dir, "WT_me3_unique_peaks.bed"),
                   A7_me3_peaks = paste0(dir, "a7_me3_specific_peaks.bed"),
                   C9_me3_peaks = paste0(dir, "a7_me3_specific_peaks.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_me3 <- annotate_peaks(peak_files, "WT_me3_peaks") 
annotated_peaks_A7_me3 <- annotate_peaks(peak_files, "A7_me3_peaks")
annotated_peaks_C9_me3 <- annotate_peaks(peak_files, "C9_me3_peaks")


dir_ac <- ".\\macs2_peaks\\H3K27ac\\contrasts\\"

#load peak files
peak_files <- list(WT_ac_peaks = paste0(dir_ac, "WT_ac_unique_peaks.bed"),
                   A7_ac_peaks = paste0(dir_ac, "a7_ac_specific_peaks.bed"),
                   C9_ac_peaks = paste0(dir_ac, "c9_ac_specific_peaks.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_ac <- annotate_peaks(peak_files, "WT_ac_peaks") 
annotated_peaks_A7_ac <- annotate_peaks(peak_files, "A7_ac_peaks")
annotated_peaks_C9_ac <- annotate_peaks(peak_files, "C9_ac_peaks")


dir_ub <- ".\\macs2_peaks\\H2AK119ub\\contrasts\\"
#load peak files
peak_files <- list(WT_ub_peaks = paste0(dir_ub, "WT_ub_unique_peaks.bed"),
                   A7_ub_peaks = paste0(dir_ub, "a7_ub_specfic_peaks.bed"),
                   C9_ub_peaks = paste0(dir_ub, "c9_ub_specfic_peaks.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_ub <- annotate_peaks(peak_files, "WT_ub_peaks") 
annotated_peaks_A7_ub <- annotate_peaks(peak_files, "A7_ub_peaks")
annotated_peaks_C9_ub <- annotate_peaks(peak_files, "C9_ub_peaks")










