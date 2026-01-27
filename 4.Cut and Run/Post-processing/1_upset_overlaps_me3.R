#input files: MACS2 peak files that have only chr name, start, and end (first 3 columns, MACS2 outputs more)
#output files: unique peaks in every condition (WT, A7, and C9) and also overlaps between the conditions

library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(tidyverse)
library(ChIPseeker)

fs <- dir(".\\macs2_peaks\\H3K27me3\\", 
          pattern = "*_peaks.bed", full.names = TRUE)



venn <- vennCount(fs, maxgap = 1000, FUN = min)


upset_themes_fix <- lapply(ComplexUpset::upset_themes, function(.ele){
  lapply(.ele, function(.e){
    do.call(theme, .e[names(.e) %in% names(formals(theme))])
  })
})

upsetPlot(venn,
          themes = list(default=theme_bw()))

combinations <- venn$combinations
expInput <- venn$counts

plotdata <- combinations[rep(rownames(combinations), expInput), ] %>%
  as.data.frame()
colnames(plotdata)
names(plotdata) <- c("Jurkat A7 H3K27me3", "Jurkat C9 H3K27me3", "Jurkat WT H3K27me3")

plotdata <- plotdata %>%
  dplyr::select("Jurkat A7 H3K27me3", "Jurkat C9 H3K27me3", "Jurkat WT H3K27me3")


sets_order <- c("Jurkat C9 H3K27me3", "Jurkat A7 H3K27me3", "Jurkat WT H3K27me3")


library(UpSetR)

upset(
  plotdata,
  sets = sets_order,   # specify which columns are sets
  keep.order = TRUE,
  order.by = "freq",                    # order intersections by frequency
  main.bar.color = c("black", "#FFA500", "darkgrey", "darkgrey", "darkgrey", "#8B008B",  "darkgrey") ,        # intersection bars color
  sets.bar.color = c("#FFA500", "#8B008B", "black"),  # A7, C9, WT,
  text.scale = 1.5
)


wt_me3_unique_peaks <- venn$overlapList[["001"]][["WT_me3_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c9_me3_unique_peaks <- venn$overlapList[["010"]][["C9_me3_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_me3_unique_peaks <- venn$overlapList[["100"]][["A7_me3_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

wt_a7_c9_me3_overlap_peaks <- venn$overlapList[["111"]][["WT_me3_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_c9_me3_overlap_peaks <- venn$overlapList[["110"]][["A7_me3_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

dir.create(".\\macs2_peaks\\H3K27me3\\contrasts\\")
write.table(wt_me3_unique_peaks, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\WT_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_me3_unique_peaks, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\A7_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c9_me3_unique_peaks, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\C9_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wt_a7_c9_me3_overlap_peaks, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\wt_a7_c9_me3_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_c9_me3_overlap_peaks, file = ".\\macs2_peaks\\H3K27me3\\contrasts\\a7_c9_me3_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_me3_unique_peaks, c9_me3_unique_peaks) %>%
  rbind(a7_c9_me3_overlap_peaks) %>%
  write.table(file = ".\macs2_peaks\\H3K27me3\\contrasts\\a7_c9_me3_union_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_me3_unique_peaks, a7_c9_me3_overlap_peaks) %>%
  write.table(file = ".\\macs2_peaks\\H3K27me3\\contrasts\\a7_me3_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(c9_me3_unique_peaks, a7_c9_me3_overlap_peaks) %>%
  write.table(file = ".\\macs2_peaks\\H3K27me3\\contrasts\\c9_me3_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


rbind(wt_me3_unique_peaks, a7_me3_unique_peaks) %>%
  rbind(c9_me3_unique_peaks) %>%
  rbind(a7_c9_me3_overlap_peaks) %>%
  rbind(wt_a7_c9_me3_overlap_peaks) %>%
  write.table(file = ".\\macs2_peaks\\H3K27me3\\contrasts\\all_me3_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
