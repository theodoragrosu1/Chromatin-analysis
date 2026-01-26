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

dir_me3 <- ".\\macs2_peaks\\H3K27me3\\"
dir_ub <- ".\\macs2_peaks\\H2AK119ub\\"

#load peak files and reannotate all WT me3 and ub peaks to get all PcG loci
peak_files <- list(WT_me3_peaks_all = paste0(dir_me3, "WT_me3_peaks.bed"))
peak_files_ub <- list(WT_ub_peaks_all = paste0(dir_ub, "WT_ub_peaks.bed"))

#Annotate peaks for each condition and save files
annotated_peaks_WT_me3_all <- annotate_peaks(peak_files, "WT_me3_peaks_all") 
annotated_peaks_WT_ub_all <- annotate_peaks(peak_files_ub, "WT_ub_peaks_all") 

WT_me3_all_promoters <- annotated_peaks_WT_me3_all %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  na.omit

write.table(WT_me3_all_promoters,".\\macs2_peaks\\H3K27me3\\contrasts\\annotated_peaks\\WT_me3_all_promoters.bed",row.names = F,col.names = F,quote = FALSE,sep="\t")

WT_ub_all_promoters <- annotated_peaks_WT_ub_all %>%
  dplyr::filter(annotation %in% c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)")) %>%
  na.omit

write.table(WT_ub_all_promoters,".\\macs2_peaks\\H2AK119ub\\contrasts\\annotated_peaks\\WT_ub_all_promoters.bed",row.names = F,col.names = F,quote = FALSE,sep="\t")

library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

venn.diagram(list("WT H3K27me3\n promoters" = WT_me3_all_promoters$SYMBOL,
                  "WT H2AK119ub\n promoters" = WT_ub_all_promoters$SYMBOL),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#173518", "#4a2574"),
             ".\\WT_PcGloci.png", disable.logging = TRUE)

