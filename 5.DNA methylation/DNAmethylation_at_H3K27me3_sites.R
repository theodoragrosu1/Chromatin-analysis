library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# 1. Load manifest — skip first 7 rows (Illumina header junk)
manifest <- fread("path\\EPIC-8v2-0_A2.csv", skip = 7)

# Keep probe ID and hg38 coordinates only
probe_coords <- manifest[, .(IlmnID, CHR, MAPINFO)]
probe_coords <- probe_coords[!is.na(MAPINFO) & CHR != ""]


# Load your mean beta values (adjust filename/format)
betas <- read.csv("path\\Beta_values_SeSAMe_QCDPB_Jurkat_KO_cell_lines_937690probes(Beta_values_SeSAMe_QCDPB_Jurkat).csv")  # columns: cg_id, beta
head(betas)
dim(betas)
colnames(betas)

# What does your manifest look like after loading?
head(manifest[, 1:5])
colnames(manifest)[1:10]

merged <- merge(betas, manifest[, c("IlmnID", "CHR", "MAPINFO")], 
                by.x = "Probe", 
                by.y = "IlmnID")
head(merged)

# Strip suffix to get base CpG ID
merged$cg_base <- sub("_.*", "", merged$Probe)

# Keep one row per CpG (average if duplicated)
merged_dedup <- merged %>%
  group_by(cg_base, CHR, MAPINFO) %>%
  summarise(
    beta_WT = mean(Jurkat.WT, na.rm = TRUE),
    beta_A7 = mean(Jurkat.A7, na.rm = TRUE),
    beta_C9 = mean(Jurkat.C9, na.rm = TRUE),
    .groups = "drop"
  )

dim(merged_dedup)  # should be ~800k rows
head(merged_dedup)

library(GenomicRanges)

probe_gr <- makeGRangesFromDataFrame(merged_dedup,
                                     seqnames.field = "CHR",
                                     start.field    = "MAPINFO",
                                     end.field      = "MAPINFO",
                                     keep.extra.columns = TRUE)

seqlevels(probe_gr) <- paste0("chr", seqlevels(probe_gr))
probe_gr$category <- NA

peaks_atac_only <- import("path\\peak_files\\WT_ATAC_peaks.bed")
hits_atac_only  <- findOverlaps(probe_gr, peaks_atac_only)
probe_gr$category[queryHits(hits_atac_only)] <- "ATAC_only"

peaks_atac_me3 <- import("path\\peaks\\WT_ATAC_H3K27me3.bed")
hits_atac_me3  <- findOverlaps(probe_gr, peaks_atac_me3)
probe_gr$category[queryHits(hits_atac_me3)] <- "ATAC_H3K27me3"

peaks_atac_suz12_me3 <- import("path\\peaks\\WT_ATAC_SUZ12_H3K27me3.bed")
hits_atac_suz12_me3  <- findOverlaps(probe_gr, peaks_atac_suz12_me3)
probe_gr$category[queryHits(hits_atac_suz12_me3)] <- "ATAC_SUZ12_H3K27me3"

peaks_atac_suz12 <- import("path\\peaks\\WT_ATAC_SUZ12.bed")
hits_atac_suz12  <- findOverlaps(probe_gr, peaks_atac_suz12_me3)
probe_gr$category[queryHits(hits_atac_suz12)] <- "ATAC_SUZ12"

peaks_suz12_me3 <- import("path\\peaks\\WT_SUZ12_H3K27me3.bed")
hits_suz12_me3  <- findOverlaps(probe_gr, peaks_atac_suz12_me3)
probe_gr$category[queryHits(hits_suz12_me3)] <- "SUZ12_H3K27me3"

table(probe_gr$category, useNA = "always")

plot_df <- as.data.frame(probe_gr) %>%
  filter(!is.na(category)) %>%
  tidyr::pivot_longer(cols = c(beta_WT, beta_A7, beta_C9),
                      names_to = "sample", values_to = "beta") %>%
  filter(!is.na(beta)) %>%
  filter(sample == "beta_WT") %>%   # keep only WT
  mutate(category = factor(category, levels = c("ATAC_only", "ATAC_SUZ12",
                                                "ATAC_H3K27me3", "ATAC_SUZ12_H3K27me3",
                                                "SUZ12_H3K27me3")))

ggplot(plot_df, aes(x = category, y = beta, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.08, outlier.size = 0.3, fill = "white") +
  scale_fill_manual(values = c("#6baed6","#74c476","#fd8d3c","#8B0000","#9e9ac8")) +
  labs(y = "Mean CpG beta value", x = NULL,
       title = "CpG methylation by chromatin state in Jurkat WT") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylim(0, 1)

ggsave("path\\methylation_WT_only.pdf", width = 8, height = 5)



# Add a mean KO beta column
merged_dedup$beta_KO_mean <- rowMeans(
  cbind(merged_dedup$beta_A7, merged_dedup$beta_C9), 
  na.rm = TRUE)

# Rebuild probe_gr with the new column included
probe_gr <- makeGRangesFromDataFrame(merged_dedup,
                                     seqnames.field = "CHR",
                                     start.field    = "MAPINFO",
                                     end.field      = "MAPINFO",
                                     keep.extra.columns = TRUE)
seqlevels(probe_gr) <- paste0("chr", seqlevels(probe_gr))

probe_gr$category_ko <- NA
peaks_lost     <- import("path\\peaks\\H3K27me3_lost_both_clones.bed")
hits <- findOverlaps(probe_gr, peaks_lost)
probe_gr$category_ko[queryHits(hits)] <- "Lost in KO"

peaks_retained <- import("path\\peaks\\H3K27me3_retained_both_clones.bed")
hits <- findOverlaps(probe_gr, peaks_retained)
probe_gr$category_ko[queryHits(hits)] <- "Retained in KO"

peaks_gained   <- import("path\\peaks\\H3K27me3_gained_both_clones.bed")
hits <- findOverlaps(probe_gr, peaks_gained)
probe_gr$category_ko[queryHits(hits)] <- "Gained in KO"

table(probe_gr$category_ko, useNA = "always")

dodge <- position_dodge(width = 0.8)

# --- A7 ---
probe_gr$category_A7 <- NA

suppressWarnings({
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_lost_A7.bed"))
  probe_gr$category_A7[queryHits(hits)] <- "Lost"
  
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_retained_A7.bed"))
  probe_gr$category_A7[queryHits(hits)] <- "Retained"
  
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_gained_A7.bed"))
  probe_gr$category_A7[queryHits(hits)] <- "Gained"
})

# --- C9 ---
probe_gr$category_C9 <- NA

suppressWarnings({
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_lost_C9.bed"))
  probe_gr$category_C9[queryHits(hits)] <- "Lost"
  
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_retained_C9.bed"))
  probe_gr$category_C9[queryHits(hits)] <- "Retained"
  
  hits <- findOverlaps(probe_gr, import("path\\peaks\\H3K27me3_gained_C9.bed"))
  probe_gr$category_C9[queryHits(hits)] <- "Gained"
})

table(probe_gr$category_A7, useNA = "always")
table(probe_gr$category_C9, useNA = "always")


df <- as.data.frame(probe_gr) %>%
  filter(beta_WT >= 0 & beta_WT <= 1,
         beta_A7 >= 0 & beta_A7 <= 1,
         beta_C9 >= 0 & beta_C9 <= 1)

# A7 dataframe
df_A7 <- df %>%
  filter(!is.na(category_A7)) %>%
  tidyr::pivot_longer(cols = c(beta_WT, beta_A7),
                      names_to = "sample", values_to = "beta") %>%
  mutate(clone    = "A7",
         category = factor(category_A7, levels = c("Lost", "Retained", "Gained")),
         sample   = recode(sample, "beta_WT" = "WT", "beta_A7" = "EZH2 KO"))

# C9 dataframe
df_C9 <- df %>%
  filter(!is.na(category_C9)) %>%
  tidyr::pivot_longer(cols = c(beta_WT, beta_C9),
                      names_to = "sample", values_to = "beta") %>%
  mutate(clone    = "C9",
         category = factor(category_C9, levels = c("Lost", "Retained", "Gained")),
         sample   = recode(sample, "beta_WT" = "WT", "beta_C9" = "EZH2 KO"))

plot_df <- bind_rows(df_A7, df_C9)


# Focus on just Lost regions for the key story
plot_lost <- bind_rows(
  df_A7 %>% filter(category == "Lost") %>% mutate(clone = "A7"),
  df_C9 %>% filter(category == "Lost") %>% mutate(clone = "C9")
)



# Lost
p_lost <- ggplot(bind_rows(df_A7, df_C9) %>% filter(category == "Lost"),
                 aes(x = beta, fill = sample, colour = sample)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~clone) +
  scale_fill_manual(values   = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  scale_colour_manual(values = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  labs(x = "Mean CpG beta value", y = "Density", fill = NULL, colour = NULL,
       title = "Lost H3K27me3 regions") +
  theme_classic() + xlim(0, 1)

# Retained
p_retained <- ggplot(bind_rows(df_A7, df_C9) %>% filter(category == "Retained"),
                     aes(x = beta, fill = sample, colour = sample)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~clone) +
  scale_fill_manual(values   = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  scale_colour_manual(values = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  labs(x = "Mean CpG beta value", y = "Density", fill = NULL, colour = NULL,
       title = "Retained H3K27me3 regions") +
  theme_classic() + xlim(0, 1)

# Gained
p_gained <- ggplot(bind_rows(df_A7, df_C9) %>% filter(category == "Gained"),
                   aes(x = beta, fill = sample, colour = sample)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~clone) +
  scale_fill_manual(values   = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  scale_colour_manual(values = c("WT" = "#2166ac", "EZH2 KO" = "#d73027")) +
  labs(x = "Mean CpG beta value", y = "Density", fill = NULL, colour = NULL,
       title = "Gained H3K27me3 regions") +
  theme_classic() + xlim(0, 1)

# Save separately
ggsave("path\\density_lost.pdf",     p_lost,     width = 8, height = 4)
ggsave("path\\density_retained.pdf", p_retained, width = 8, height = 4)
ggsave("path\\density_gained.pdf",   p_gained,   width = 8, height = 4)




ggsave("methylation_A7_C9_separate.pdf", width = 10, height = 5)






