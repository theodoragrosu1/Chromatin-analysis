library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(tidyverse)
library(ChIPseeker)

fs <- dir("path\\macs2_peaks\\H3K27ac\\", 
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
names(plotdata) <- c("Jurkat A7 H3K27ac", "Jurkat C9 H3K27ac", "Jurkat WT H3K27ac")

plotdata <- plotdata %>%
  dplyr::select("Jurkat A7 H3K27ac", "Jurkat C9 H3K27ac", "Jurkat WT H3K27ac")


sets_order <- c("Jurkat C9 H3K27ac", "Jurkat A7 H3K27ac", "Jurkat WT H3K27ac")


library(UpSetR)

upset(
  plotdata,
  sets = sets_order,   # specify which columns are sets
  keep.order = TRUE,
  order.by = "freq",                    # order intersections by frequency
  main.bar.color = c("#8B008B", "black","darkgrey", "darkgrey", "#FFA500", "darkgrey", "darkgrey") ,        # intersection bars color
  sets.bar.color = c("#FFA500", "#8B008B", "black"),  # A7, C9, WT,
  text.scale = 1.5
)


wt_ac_unique_peaks <- venn$overlapList[["001"]][["WT_ac_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c9_ac_unique_peaks <- venn$overlapList[["010"]][["C9_ac_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_ac_unique_peaks <- venn$overlapList[["100"]][["A7_ac_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

wt_a7_c9_ac_overlap_peaks <- venn$overlapList[["111"]][["WT_ac_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_c9_ac_overlap_peaks <- venn$overlapList[["110"]][["A7_ac_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

dir.create("path\\macs2_peaks\\H3K27ac\\contrasts\\")
write.table(wt_ac_unique_peaks, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\WT_ac_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_ac_unique_peaks, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\A7_ac_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c9_ac_unique_peaks, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\C9_ac_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wt_a7_c9_ac_overlap_peaks, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\wt_a7_c9_ac_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_c9_ac_overlap_peaks, file = "path\\macs2_peaks\\H3K27ac\\contrasts\\a7_c9_ac_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_ac_unique_peaks, c9_ac_unique_peaks) %>%
  rbind(a7_c9_ac_overlap_peaks) %>%
  write.table(file = "path\\macs2_peaks\\H3K27ac\\contrasts\\a7_c9_ac_union_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_ac_unique_peaks, a7_c9_ac_overlap_peaks) %>%
  write.table(file = "path\\macs2_peaks\\H3K27ac\\contrasts\\a7_ac_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(c9_ac_unique_peaks, a7_c9_ac_overlap_peaks) %>%
  write.table(file = "path\\macs2_peaks\\H3K27ac\\contrasts\\c9_ac_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


rbind(wt_ac_unique_peaks, a7_ac_unique_peaks) %>%
  rbind(c9_ac_unique_peaks) %>%
  rbind(a7_c9_ac_overlap_peaks) %>%
  rbind(wt_a7_c9_ac_overlap_peaks) %>%
  write.table(file = "path\\macs2_peaks\\H3K27ac\\contrasts\\all_ac_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

