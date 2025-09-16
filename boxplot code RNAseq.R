#input: normalized tpm from RNAseq 
#note: one sample is mislabelled so it was in added to the beginning, that is why colnames is needed
#output: boxplots with nomrlized tpm across conditions

library(dplyr)
library(readxl)
library(tidyverse)

log2_tpm_jurkat <- read.csv("#path to normalized tpm file") %>%
  dplyr::filter(geneID %in% c(#list of genes))

colnames(log2_tpm_jurkat) <- c("geneID", "A7_1", "WT_1", "WT_2", "WT_3", 
                                          "A7_2", "A7_3", "C9_1", "C9_2", "C9_3")


boxplots <- log2_tpm_jurkat %>%
  pivot_longer(-geneID, names_to = "Sample", values_to = "mRNA") %>%
  mutate(Condition = substr(Sample, 1, 2))
boxplots$Condition <- factor(boxplots$Condition, levels = c("WT", "A7", "C9"))
boxplots$geneID <- factor(boxplots$geneID, levels = c(#list of genes))
p1 <- boxplots %>%
  ggplot(aes(x = Condition, y = mRNA, fill = Condition)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("black", "#8B008B", "#FFA500")) +
  geom_jitter(width = 0.25, size = 2.8, fill = "black") +
  facet_wrap(~geneID, scales = "free") +
  theme_bw(base_size = 24) +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  labs(y = expression("log"[2] *"(tpm+1)"), legend = NA, title = "RNA-Seq")
p1


ggsave("#path to where to save the pdf",
       plot = p1, width = 16, height = 12,
       dpi = 1000, units = "cm", device = cairo_pdf)


