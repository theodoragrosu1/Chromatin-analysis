library(dplyr)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(Cairo)
library(colorspace)
library(dplyr)

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

dir <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\EZH1\\"

#load peak files
peak_files <- list(WT_EZH1_peaks = paste0(dir, "WT_EZH1_peaks_filtered.bed"),
                   A7_EZH1_peaks = paste0(dir, "A7_EZH1_peaks_filtered.bed"),
                   C9_EZH1_peaks = paste0(dir, "C9_EZH1_peaks_filtered.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_EZH1 <- annotate_peaks(peak_files, "WT_EZH1_peaks") 
annotated_peaks_A7_EZH1 <- annotate_peaks(peak_files, "A7_EZH1_peaks")
annotated_peaks_C9_EZH1 <- annotate_peaks(peak_files, "C9_EZH1_peaks")

######here: change for SUZ12
dir_ac <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\macs2_peaks\\H3K27ac\\contrasts\\"

#load peak files
peak_files <- list(WT_ac_peaks = paste0(dir_ac, "WT_ac_unique_peaks.bed"),
                   A7_ac_peaks = paste0(dir_ac, "a7_ac_specific_peaks.bed"),
                   C9_ac_peaks = paste0(dir_ac, "c9_ac_specific_peaks.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_ac <- annotate_peaks(peak_files, "WT_ac_peaks") 
annotated_peaks_A7_ac <- annotate_peaks(peak_files, "A7_ac_peaks")
annotated_peaks_C9_ac <- annotate_peaks(peak_files, "C9_ac_peaks")


###load ATAC peaks
WT_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\WT_unique_atac.csv") %>%
  na.omit
A7_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\A7_unique_atac.csv") %>%
  na.omit
C9_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\C9_unique_atac.csv") %>%
  na.omit








