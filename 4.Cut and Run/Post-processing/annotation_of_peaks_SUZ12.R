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

dir <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\"

#load peak files
peak_files <- list(WT_SUZ12_peaks = paste0(dir, "WT_SUZ12_peaks_filtered.bed"),
                   A7_SUZ12_peaks = paste0(dir, "A7_SUZ12_peaks_filtered.bed"),
                   C9_SUZ12_peaks = paste0(dir, "C9_SUZ12_peaks_filtered.bed"))
#Annotate peaks for each condition and save files
annotated_peaks_WT_SUZ12 <- annotate_peaks(peak_files, "WT_SUZ12_peaks") 
annotated_peaks_A7_SUZ12 <- annotate_peaks(peak_files, "A7_SUZ12_peaks")
annotated_peaks_C9_SUZ12 <- annotate_peaks(peak_files, "C9_SUZ12_peaks")


###load ATAC peaks
WT_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\WT_unique_atac.csv") %>%
  na.omit
A7_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\A7_unique_atac.csv") %>%
  na.omit
C9_ATAC_peaks <- read.csv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\ETP-ALL_chapter\\ATAC-peaks\\C9_unique_atac.csv") %>%
  na.omit








