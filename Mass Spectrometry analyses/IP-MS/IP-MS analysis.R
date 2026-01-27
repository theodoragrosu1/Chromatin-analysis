library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase) 
library(Biobase) 
#library(GSVA) 
library(gprofiler2) 
library(clusterProfiler) 
library(msigdbr) 
library(gplots)
library(enrichplot)
library(ggplot2)
library(ggrepel)


difference <- read_xlsx(".\\comparing_stats_analyses_for_FDR_calculation.xlsx", sheet = "volcano-plot") %>%
  select("Gene", "Difference") %>%
  unique() %>%
  as.data.frame()

fdr <- read_xlsx(".\\comparing_stats_analyses_for_FDR_calculation.xlsx", sheet = "permutation-FDR") 
fdr_filtered <- fdr %>%
  dplyr::filter(q.value < 0.05) 
fdr_filtered <- fdr_filtered %>%
  select("Gene", "log.p.value", "FDR") %>%
  unique() %>%
  as.data.frame()


suz12_ip <- left_join(difference, fdr_filtered, by = "Gene") %>%
  drop_na()

suz12_ip$FDR <-  p.adjust(suz12_ip$log.p.value, method = 'hochberg', n = nrow(suz12_ip))



vplot <- ggplot(suz12_ip) +
  aes(y=log.p.value, x=Difference, text = paste("Symbol:", Gene)) +
  geom_point(size=2, alpha = 0.5) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 0.5, xmax = 5.1, ymin = -log10(0.1), ymax = 17, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -0.5, xmax = -5.1, ymin = -log10(0.1), ymax = 17, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -1.1, y = 15.4, label = paste0("Proteins pulled\n in A7: 158"), size = 4, colour = "#2C467A") + #labels are manually added after running summary(results.t.all downstream)
  annotate(geom = "text", x = 1.1, y = 15.4, label = paste0("Proteins pulled\n in WT: 105"), size = 4, colour = "#BE684D") +
  labs(title="Volcano plot",
       subtitle = "SUZ12 IP WT - A7",
       caption="FDR <= 0.1", 
       y = "-log10(FDR)", 
       x = "Difference(MaxLFQ)\n") +
  theme_bw(base_size = 20) + 
  geom_label_repel(data=filter(suz12_ip, abs(Difference)>2 & abs(FDR) > 0.1), aes(label=Gene), max.overlaps = 15)+
  geom_label_repel(data=filter(suz12_ip, Gene %in% c("SMARCC2", "EZH1", "USP7")), aes(label=Gene), max.overlaps = 15)


vplot


###redone analysis of DEGs - try to look for complex members separately 
suz12_all <- read_xlsx(".\\DEG IP-MS.xlsx", sheet = 1) %>%
  select("Gene", "pvalue", "Difference") %>%
  unique() %>%
  as.data.frame()

suz12_all <- suz12_all[-501,]

suz12_up_wt <- suz12_all %>%
  dplyr::filter(Difference > 1) %>%
  dplyr::filter(pvalue <= 0.1)

suz12_up_a7 <- suz12_all %>%
  dplyr::filter(Difference < -1) %>%
  dplyr::filter(pvalue <= 0.1)

vplot <- ggplot(suz12_all) +
  aes(y=-log10(pvalue), x=Difference, text = paste("Symbol:", Gene)) +
  geom_point(size=2, alpha = 0.5) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 0.5, xmax = 5.1, ymin = -log10(0.1), ymax = 12, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -0.5, xmax = -5.1, ymin = -log10(0.1), ymax = 12, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -3.8, y = 11.2, label = paste0("39 proteins\n in EZH2 KO"), size = 6, colour = "#2C467A") + #labels are manually added after running summary(results.t.all downstream)
  annotate(geom = "text", x = 3.8, y = 11.2, label = paste0("102 proteins\n in EZH2 WT"), size = 6, colour = "#BE684D") +
  labs(title="SUZ12 IP WT - A7",
       subtitle = "Volcano Plot",
       caption="FDR <= 0.1\n`Fold change > 1", 
       y = "-log10(p-value)", 
       x = "Fold change\n") +
  theme_bw(base_size = 20) +
  geom_label_repel(data=filter(suz12_all,  Gene %in% c("EZH1", "EZH2", "EED", "AEBP2", "JARID2", "SUZ12", "PALI1", "PALI2", "EPOP", "PHF1", "MTF2", "PHF19", "RBBP4", "RBBP7")), 
                   aes(label=Gene), max.overlaps = 25)


vplot


#EZH2 IP volcano plot
ezh2_all <- read_xlsx(".\\DEG IP-MS.xlsx", sheet = 2) %>%
  select("Gene", "pvalue", "Difference") %>%
  unique() %>%
  as.data.frame()

ezh2_all <- ezh2_all[-200,]

ezh2_up_wt <- ezh2_all %>%
  dplyr::filter(Difference > 1) %>%
  dplyr::filter(pvalue <= 0.1)

ezh2_up_a7 <- ezh2_all %>%
  dplyr::filter(Difference < -1) %>%
  dplyr::filter(pvalue <= 0.1)




vplot <- ggplot(ezh2_all) +
  aes(y=-log10(pvalue), x=Difference, text = paste("Symbol:", Gene)) +
  geom_point(size=2, alpha = 0.5) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 0.5, xmax = 5.1, ymin = -log10(0.1), ymax = 12, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -0.5, xmax = -5.1, ymin = -log10(0.1), ymax = 12, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -4, y = 11, label = paste0("13 proteins\n in EZH2 KO"), size = 5, colour = "#2C467A") + #labels are manually added after running summary(results.t.all downstream)
  annotate(geom = "text", x = 4, y = 11, label = paste0("91 proteins\n in EZH2 WT"), size = 5, colour = "#BE684D") +
  labs(title="EZH2 IP WT - A7",
       subtitle = "Volcano Plot",
       caption="FDR <= 0.1\n`Fold change > 1", 
       y = "-log10(p-value)", 
       x = "Fold change\n") +
  theme_bw(base_size = 20) +
  geom_label_repel(data=filter(ezh2_all, Gene %in% c("EZH1", "EZH2", "EED", "AEBP2", "JARID2", "SUZ12", "PALI1", "PALI2", "EPOP", "PHF1", "MTF2", "PHF19", "RBBP4", "RBBP7")), 
                   aes(label=Gene), max.overlaps = 25)


vplot

