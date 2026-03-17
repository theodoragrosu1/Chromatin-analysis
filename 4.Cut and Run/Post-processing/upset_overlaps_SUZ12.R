library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(tidyverse)
library(ChIPseeker)

fs <- dir("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\", 
          pattern = "*_filtered.bed", full.names = TRUE)


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
names(plotdata) <- c("Jurkat A7 SUZ12", "Jurkat C9 SUZ12", "Jurkat WT SUZ12")

plotdata <- plotdata %>%
  dplyr::select("Jurkat A7 SUZ12", "Jurkat C9 SUZ12", "Jurkat WT SUZ12")


sets_order <- c("Jurkat C9 SUZ12", "Jurkat A7 SUZ12", "Jurkat WT SUZ12")


library(UpSetR)

upset(
  plotdata,
  sets = sets_order,   # specify which columns are sets
  keep.order = TRUE,
  order.by = "freq",                    # order intersections by frequency
  main.bar.color = c("#8B008B","black", "#FFA500", "darkgrey", "darkgrey", "darkgrey", "darkgrey") ,        # intersection bars color
  sets.bar.color = c("#FFA500", "#8B008B", "black"),  # A7, C9, WT,
  text.scale = 1.5
)


wt_SUZ12_unique_peaks <- venn$overlapList[["001"]][["WT_SUZ12_peaks_filtered.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c9_SUZ12_unique_peaks <- venn$overlapList[["010"]][["C9_SUZ12_peaks_filtered.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_SUZ12_unique_peaks <- venn$overlapList[["100"]][["A7_SUZ12_peaks_filtered.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

wt_a7_c9_SUZ12_overlap_peaks <- venn$overlapList[["111"]][["WT_SUZ12_peaks_filtered.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

a7_c9_SUZ12_overlap_peaks <- venn$overlapList[["110"]][["A7_SUZ12_peaks_filtered.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

dir.create("C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\")
write.table(wt_SUZ12_unique_peaks, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\WT_SUZ12_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_SUZ12_unique_peaks, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\A7_SUZ12_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c9_SUZ12_unique_peaks, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\C9_SUZ12_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wt_a7_c9_SUZ12_overlap_peaks, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\wt_a7_c9_SUZ12_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(a7_c9_SUZ12_overlap_peaks, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\a7_c9_SUZ12_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_SUZ12_unique_peaks, c9_SUZ12_unique_peaks) %>%
  rbind(a7_c9_SUZ12_overlap_peaks) %>%
  write.table(file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\a7_c9_SUZ12_union_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(a7_SUZ12_unique_peaks, a7_c9_SUZ12_overlap_peaks) %>%
  write.table(file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\a7_SUZ12_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(c9_SUZ12_unique_peaks, a7_c9_SUZ12_overlap_peaks) %>%
  write.table(file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\c9_SUZ12_specific_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


rbind(wt_SUZ12_unique_peaks, a7_SUZ12_unique_peaks) %>%
  rbind(c9_SUZ12_unique_peaks) %>%
  rbind(a7_c9_SUZ12_overlap_peaks) %>%
  rbind(wt_a7_c9_SUZ12_overlap_peaks) %>%
  write.table(file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\cutrun_peaks\\SUZ12\\SUZ12_filtered\\contrasts\\all_SUZ12_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
