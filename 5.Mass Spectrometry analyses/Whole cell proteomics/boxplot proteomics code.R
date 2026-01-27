#input: normalized proteomics data 
#output: boxplot with LFQ intensities for all three conditions and chosen proteins

library(tidyverse)
library(dplyr)

boxplots_proteomics <- read.csv(".\\proteomics_Jurkat_whole_cell_filtered_logFC_LFQ.csv") %>%
  dplyr::select(Gene_ID:C9_6) %>%
  dplyr::filter(Gene_ID %in% c("CD1A", "CD1B", "CD1C", "CD1D", "NCAM1", "MEF2C")) %>%   #change this for genes of interest
  pivot_longer(-Gene_ID, names_to = "Sample", values_to = "log2LFQ") %>%
  mutate(Condition = substr(Sample, 1, 2)) 

boxplots_proteomics$Condition <- factor(boxplots_proteomics$Condition, levels = c("WT", "A7", "C9"))


boxplots_proteomics$Gene_ID <- factor(boxplots_proteomics$Gene_ID, levels = c("CD1A", "CD1B", "CD1C", "CD1D", "NCAM1", "MEF2C"))  #change this based on genes of interest
p1 <- boxplots_proteomics %>%
  ggplot(aes(x = Condition, y = log2LFQ, fill = Condition)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("black", "#8B008B", "#FFA500")) +
  geom_jitter(width = 0.25, size = 2.8, fill = "black") +
  facet_wrap(~Gene_ID, scales = "free") +
  theme_bw(base_size = 24) +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  labs(y = expression("log"[2] *"(tpm+1)"), legend = NA, title = "Proteomics")
p1
