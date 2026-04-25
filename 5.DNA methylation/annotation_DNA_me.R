library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# 1. Load manifest — skip first 7 rows (Illumina header junk)
manifest <- fread("path\\MethylationEPIC v2.0 Files\\EPIC-8v2-0_A2.csv", skip = 7)

# Keep probe ID and hg38 coordinates only
probe_coords <- manifest[, .(IlmnID, CHR, MAPINFO, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island)]
probe_coords <- probe_coords[!is.na(MAPINFO) & CHR != ""]

mean_beta_values <- read.csv("path\\Mean_Beta_values_Sesame_QCDPB_Jurkat_cell_line_930659CpGs(Mean_Beta_values_Sesame_QCDPB_J).csv")  # columns: cg_id, beta
head(mean_beta_values)
dim(mean_beta_values)
colnames(mean_beta_values)

mean_beta_values <- mean_beta_values %>%
  mutate(
    delta_A7 = Jurkat.A7 - Jurkat.WT,
    delta_C9 = Jurkat.C9 - Jurkat.WT
  ) %>%
  na.omit()
###hypermethylated regions
hypermethylated_A7_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_A7 > 0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.A7", "delta_A7"))
row.names(hypermethylated_A7_vs_WT) <- hypermethylated_A7_vs_WT$Name

hypermethylated_C9_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_C9 > 0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.C9", "delta_C9"))
row.names(hypermethylated_C9_vs_WT) <- hypermethylated_C9_vs_WT$Name

hypermethylated_clones_common <- left_join(hypermethylated_C9_vs_WT, hypermethylated_A7_vs_WT, by = "Name") %>%
  na.omit()

head(hypermethylated_clones_common)

##hypomethylated regions
hypomethylated_A7_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_A7 < -0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.A7", "delta_A7"))
row.names(hypomethylated_A7_vs_WT) <- hypomethylated_A7_vs_WT$Name

hypomethylated_C9_vs_WT <- mean_beta_values %>%
  dplyr::filter(delta_C9 < -0.4) %>%
  dplyr::select(c("Name", "Jurkat.WT", "Jurkat.C9", "delta_C9"))
row.names(hypomethylated_C9_vs_WT) <- hypomethylated_C9_vs_WT$Name

hypomethylated_clones_common <- left_join(hypomethylated_C9_vs_WT, hypomethylated_A7_vs_WT, by = "Name") %>%
  na.omit()
head(hypomethylated_clones_common)


###annotate regions based on proximity to genes
probe_coords_clean <- probe_coords %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, UCSC_RefGene_Group, .keep_all = TRUE) 

annotation_genes_clean <- probe_coords_clean %>%
  mutate(
    region_priority = case_when(
      grepl("TSS200", UCSC_RefGene_Group) ~ "TSS200",
      grepl("TSS1500", UCSC_RefGene_Group) ~ "TSS1500",
      grepl("5UTR", UCSC_RefGene_Group) ~ "5UTR",
      grepl("3UTR", UCSC_RefGene_Group) ~ "3UTR",
      grepl("exon", UCSC_RefGene_Group) ~ "Exon",
      TRUE ~ "intron/intergenic"
    )
  ) %>%
  select(Name, region_priority)


mean_beta_values_annotated <- mean_beta_values %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

###Analysis for separate
hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##Anaylsis for common sites
hypermethylated_clones_common_annotated <- hypermethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_clones_common_annotated <- hypomethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##################Visualization for WT vs A7
df1 <- mean_beta_values_annotated$region_priority %>%
  as.data.frame()

df1 <- df1 %>% mutate(Sample = "Total_sites")
colnames(df1) <- c("annotation", "Sample")
df1_summary <- as.data.frame(table(df1$annotation))
colnames(df1_summary) <- c("annotation", "Total_sites")


df2 <- hypermethylated_A7_vs_WT_annotated[,5] %>%
  as.data.frame()
df2 <- df2 %>% mutate(Sample = "Hypermethylated_in_A7")
colnames(df2) <- c("annotation", "Sample")
df2_summary <- as.data.frame(table(df2$annotation))
colnames(df2_summary) <- c("annotation", "Hypermethylated_in_A7")

df3 <- hypomethylated_A7_vs_WT_annotated[,5] %>%
  as.data.frame()
df3 <- df3 %>% mutate(Sample = "Hypomethylated_in_A7")
colnames(df3) <- c("annotation", "Sample")
df3_summary <- as.data.frame(table(df3$annotation))
colnames(df3_summary) <- c("annotation", "Hypomethylated_in_A7")

library(tidyverse)
###For the number plot, it maybe adds nothing to keep the total sites because what does it say? 
df_all_numbers <- left_join(df2_summary, df3_summary, by ="annotation") %>%
  pivot_longer(cols = -annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(annotation != "intron/intergenic")

df_all_numbers$annotation <- factor(df_all_numbers$annotation, 
                                            levels = c("3UTR", "Exon", 
                                                       "TSS200", "TSS1500",  
                                                       "5UTR", "intron/intergenic"))

df_all_numbers$Condition <- factor(df_all_numbers$Condition, 
                                           levels = c("Hypermethylated_in_A7", "Hypomethylated_in_A7"))
#For % plot
df_all <- bind_rows(df1, df2, df3)
table(df_all$Sample)
df_summary <- df_all %>%
  group_by(Sample, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  na.omit()

library(ggplot2)
library(colorspace)

annotated_regions_plot_perc <- df_summary %>%
  ggplot(aes(x = Sample, y = percent, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "TSS200" = darken("#009E73", amount = 0.5),
    "TSS1500" = lighten("#009E73", amount = 0.5),
    "5UTR" = "#0072B2",
    "Exon" = "#D55E00",
    "3UTR" = "#FEB95F",
    "intron/intergenic" = "#CC79A7"
  ),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_classic(base_size = 26) +
  coord_flip() +
  labs(x = "Condition") +
  ylab("Percentage (%)") +
  theme(panel.grid = element_blank())


annotated_regions_peaks_plot <- df_all_numbers %>%
  ggplot(aes(x = Condition, y = value, fill = annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(
    "TSS200" = darken("#009E73", amount = 0.5),
    "TSS1500" = lighten("#009E73", amount = 0.5),
    "5UTR" = "#0072B2",
    "Exon" = "#D55E00",
    "3UTR" = "#FEB95F",
    "intron/intergenic" = "#CC79A7"
  ), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 10000)) +
  labs(x = "Condition") +
  ylab("No. DNAme probes")



#################Visualization for WT vs C9
df1 <- mean_beta_values_annotated$region_priority %>%
  as.data.frame()

df1 <- df1 %>% mutate(Sample = "Total_sites")
colnames(df1) <- c("annotation", "Sample")
df1_summary <- as.data.frame(table(df1$annotation))
colnames(df1_summary) <- c("annotation", "Total_sites")


df2 <- hypermethylated_C9_vs_WT_annotated[,5] %>%
  as.data.frame()
df2 <- df2 %>% mutate(Sample = "Hypermethylated_in_C9")
colnames(df2) <- c("annotation", "Sample")
df2_summary <- as.data.frame(table(df2$annotation))
colnames(df2_summary) <- c("annotation", "Hypermethylated_in_C9")

df3 <- hypomethylated_C9_vs_WT_annotated[,5] %>%
  as.data.frame()
df3 <- df3 %>% mutate(Sample = "Hypomethylated_in_C9")
colnames(df3) <- c("annotation", "Sample")
df3_summary <- as.data.frame(table(df3$annotation))
colnames(df3_summary) <- c("annotation", "Hypomethylated_in_C9")

###For the number plot, it maybe adds nothing to keep the total sites because what does it say? 
df_all_numbers <- left_join(df2_summary, df3_summary, by ="annotation") %>%
  pivot_longer(cols = -annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(annotation != "intron/intergenic")

df_all_numbers$annotation <- factor(df_all_numbers$annotation, 
                                    levels = c("3UTR", "Exon", 
                                               "TSS200", "TSS1500",  
                                               "5UTR", "intron/intergenic"))

df_all_numbers$Condition <- factor(df_all_numbers$Condition, 
                                   levels = c("Hypermethylated_in_C9", "Hypomethylated_in_C9"))
#For % plot
df_all <- bind_rows(df1, df2, df3)
table(df_all$Sample)
df_summary <- df_all %>%
  group_by(Sample, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  na.omit()

library(ggplot2)
library(colorspace)

annotated_regions_plot_perc <- df_summary %>%
  ggplot(aes(x = Sample, y = percent, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "TSS200" = darken("#009E73", amount = 0.5),
    "TSS1500" = lighten("#009E73", amount = 0.5),
    "5UTR" = "#0072B2",
    "Exon" = "#D55E00",
    "3UTR" = "#FEB95F",
    "intron/intergenic" = "#CC79C9"
  ),
  guide = guide_legend(reverse = TRUE)
  ) +
  theme_classic(base_size = 26) +
  coord_flip() +
  labs(x = "Condition") +
  ylab("Percentage (%)") +
  theme(panel.grid = element_blank())


annotated_regions_peaks_plot <- df_all_numbers %>%
  ggplot(aes(x = Condition, y = value, fill = annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(
    "TSS200" = darken("#009E73", amount = 0.5),
    "TSS1500" = lighten("#009E73", amount = 0.5),
    "5UTR" = "#0072B2",
    "Exon" = "#D55E00",
    "3UTR" = "#FEB95F",
    "intron/intergenic" = "#CC79C9"
  ), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 5000)) +
  labs(x = "Condition") +
  ylab("No. DNAme probes")



###annotate regions based on CGIs
probe_coords_clean <- probe_coords %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, Relation_to_UCSC_CpG_Island, .keep_all = TRUE) 

annotation_genes_clean <- probe_coords_clean %>%
  mutate(
    region_priority = case_when(
      grepl("Island", Relation_to_UCSC_CpG_Island) ~ "Island",
      grepl("S_Shelf", Relation_to_UCSC_CpG_Island) ~ "S_Shelf",
      grepl("S_Shore", Relation_to_UCSC_CpG_Island) ~ "S_Shore",
      grepl("N_Shelf", Relation_to_UCSC_CpG_Island) ~ "N_Shelf",
      grepl("N_Shore", Relation_to_UCSC_CpG_Island) ~ "N_Shore",
      TRUE ~ "Open sea"
    )
  ) %>%
  dplyr::select(Name, region_priority)


mean_beta_values_annotated <- mean_beta_values %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

###Analysis for separate
hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 


##Anaylsis for common sites
hypermethylated_clones_common_annotated <- hypermethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

hypomethylated_clones_common_annotated <- hypomethylated_clones_common %>%
  left_join(
    annotation_genes_clean,
    by = "Name"
  ) 

##################Visualization for WT vs A7
df1 <- mean_beta_values_annotated$region_priority %>%
  as.data.frame()

df1 <- df1 %>% mutate(Sample = "Total_sites")
colnames(df1) <- c("annotation", "Sample")
df1_summary <- as.data.frame(table(df1$annotation))
colnames(df1_summary) <- c("annotation", "Total_sites")


df2 <- hypermethylated_A7_vs_WT_annotated[,5] %>%
  as.data.frame()
df2 <- df2 %>% mutate(Sample = "Hypermethylated_in_A7")
colnames(df2) <- c("annotation", "Sample")
df2_summary <- as.data.frame(table(df2$annotation))
colnames(df2_summary) <- c("annotation", "Hypermethylated_in_A7")

df3 <- hypomethylated_A7_vs_WT_annotated[,5] %>%
  as.data.frame()
df3 <- df3 %>% mutate(Sample = "Hypomethylated_in_A7")
colnames(df3) <- c("annotation", "Sample")
df3_summary <- as.data.frame(table(df3$annotation))
colnames(df3_summary) <- c("annotation", "Hypomethylated_in_A7")

###For the number plot, it maybe adds nothing to keep the total sites because what does it say? 
df_all_numbers <- left_join(df2_summary, df3_summary, by ="annotation") %>%
  pivot_longer(cols = -annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(annotation != "Open sea")

df_all_numbers$annotation <- factor(df_all_numbers$annotation, 
                                    levels = c("S_Shelf", "S_Shore", 
                                               "Island", "N_Shelf",  
                                               "N_Shore"))

df_all_numbers$Condition <- factor(df_all_numbers$Condition, 
                                   levels = c("Hypermethylated_in_A7", "Hypomethylated_in_A7"))
#For % plot
df_all <- bind_rows(df1, df2, df3)
table(df_all$Sample)
df_summary <- df_all %>%
  group_by(Sample, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  na.omit()

library(ggplot2)
library(colorspace)

annotated_regions_plot_perc <- df_summary %>%
  ggplot(aes(x = Sample, y = percent, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "N_Shelf" = darken("#009E73", amount = 0.5),
    "N_Shore" = lighten("#009E73", amount = 0.5),
    "S_Shelf" = darken("#CC79A7", amount = 0.5),
    "S_Shore" = lighten("#CC79A7", amount = 0.5),
    "Island" = "#FEB95F",
    "Open sea" = "#0072B2" 
  ),
  guide = guide_legend(reverse = TRUE)
  ) +
  theme_classic(base_size = 26) +
  coord_flip() +
  labs(x = "Condition") +
  ylab("Percentage (%)") +
  theme(panel.grid = element_blank())


annotated_regions_peaks_plot <- df_all_numbers %>%
  ggplot(aes(x = Condition, y = value, fill = annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(
    "N_Shelf" = darken("#009E73", amount = 0.5),
    "N_Shore" = lighten("#009E73", amount = 0.5),
    "S_Shelf" = darken("#CC79A7", amount = 0.5),
    "S_Shore" = lighten("#CC79A7", amount = 0.5),
    "Island" = "#FEB95F" 
  ), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 10000)) +
  labs(x = "Condition") +
  ylab("No. DNAme probes")



#################Visualization for WT vs C9
df1 <- mean_beta_values_annotated$region_priority %>%
  as.data.frame()

df1 <- df1 %>% mutate(Sample = "Total_sites")
colnames(df1) <- c("annotation", "Sample")
df1_summary <- as.data.frame(table(df1$annotation))
colnames(df1_summary) <- c("annotation", "Total_sites")


df2 <- hypermethylated_C9_vs_WT_annotated[,5] %>%
  as.data.frame()
df2 <- df2 %>% mutate(Sample = "Hypermethylated_in_C9")
colnames(df2) <- c("annotation", "Sample")
df2_summary <- as.data.frame(table(df2$annotation))
colnames(df2_summary) <- c("annotation", "Hypermethylated_in_C9")

df3 <- hypomethylated_C9_vs_WT_annotated[,5] %>%
  as.data.frame()
df3 <- df3 %>% mutate(Sample = "Hypomethylated_in_C9")
colnames(df3) <- c("annotation", "Sample")
df3_summary <- as.data.frame(table(df3$annotation))
colnames(df3_summary) <- c("annotation", "Hypomethylated_in_C9")

###For the number plot, it maybe adds nothing to keep the total sites because what does it say? 
df_all_numbers <- left_join(df2_summary, df3_summary, by ="annotation") %>%
  pivot_longer(cols = -annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(annotation != "Open sea")

df_all_numbers$annotation <- factor(df_all_numbers$annotation, 
                                    levels = c("S_Shelf", "S_Shore", 
                                               "Island", "N_Shelf",  
                                               "N_Shore", "Open sea"))

df_all_numbers$Condition <- factor(df_all_numbers$Condition, 
                                   levels = c("Hypermethylated_in_C9", "Hypomethylated_in_C9"))
#For % plot
df_all <- bind_rows(df1, df2, df3)
table(df_all$Sample)
df_summary <- df_all %>%
  group_by(Sample, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  na.omit()

library(ggplot2)
library(colorspace)

annotated_regions_plot_perc <- df_summary %>%
  ggplot(aes(x = Sample, y = percent, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "N_Shelf" = darken("#009E73", amount = 0.5),
    "N_Shore" = lighten("#009E73", amount = 0.5),
    "S_Shelf" = darken("#CC79A7", amount = 0.5),
    "S_Shore" = lighten("#CC79A7", amount = 0.5),
    "Island" = "#FEB95F",
    "Open sea" = "#0072B2" 
  ),
  guide = guide_legend(reverse = TRUE)
  ) +
  theme_classic(base_size = 26) +
  coord_flip() +
  labs(x = "Condition") +
  ylab("Percentage (%)") +
  theme(panel.grid = element_blank())


annotated_regions_peaks_plot <- df_all_numbers %>%
  ggplot(aes(x = Condition, y = value, fill = annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(
    "N_Shelf" = darken("#009E73", amount = 0.5),
    "N_Shore" = lighten("#009E73", amount = 0.5),
    "S_Shelf" = darken("#CC79A7", amount = 0.5),
    "S_Shore" = lighten("#CC79A7", amount = 0.5),
    "Island" = "#FEB95F"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 1000)) +
  labs(x = "Condition") +
  ylab("No. DNAme probes")

#################################annotate the common ones to genes based on the manifest
genes <- manifest[, .(IlmnID, UCSC_RefGene_Name)]
genes_clean <- genes %>%
  mutate(Name = sub("_.*", "", IlmnID)) %>%
  distinct(Name, UCSC_RefGene_Name, .keep_all = TRUE) 

hypomethylated_clones_common_annotated <- hypomethylated_clones_common %>%
  left_join(genes_clean, by = "Name")
hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(genes_clean, by = "Name")
hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(genes_clean, by = "Name")

hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join (genes_clean, by = "Name")
hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join (genes_clean, by = "Name")
hypermethylated_clones_common_annotated <- hypermethylated_clones_common %>%
  left_join(genes_clean, by = "Name")

all_beta_values_annotated <- mean_beta_values %>%
  left_join(genes_clean, by = "Name")

library(WebGestaltR)
hypermethylated_clones_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                 "geneontology_Molecular_Function_noRedundant"),
                              interestGeneType = "genesymbol",
                              referenceGeneType = "genesymbol",
                              interestGene = hypermethylated_clones_common_annotated$UCSC_RefGene_Name,
                              referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                              sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = hypermethylated_clones_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

ggplot(GSEA.chosen, aes(x = `FDR`, 
                       y = fct_reorder(`description`, `FDR`))) +
  geom_point(aes(size = `size`, color = FDR),
             alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(4, 12)) +
  labs(
    title = "Common hypermethylated\n regions in EZH2 KO",
    x = "FDR",
    y = "Term",
    color = "FDR"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15)
  )

hypomethylated_clones_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                            "geneontology_Molecular_Function_noRedundant"),
                                         interestGeneType = "genesymbol",
                                         referenceGeneType = "genesymbol",
                                         interestGene = hypomethylated_clones_common_annotated$UCSC_RefGene_Name,
                                         referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                                         sigMethod = "fdr", fdrThr = 0.1)
#no regions were detected

############################################do separate analysis
hypomethylated_A7_vs_WT_annotated <- hypomethylated_A7_vs_WT %>%
  left_join(genes_clean, by = "Name")

hypermethylated_A7_vs_WT_annotated <- hypermethylated_A7_vs_WT %>%
  left_join(genes_clean, by = "Name")

all_beta_values_annotated <- mean_beta_values %>%
  left_join(genes_clean, by = "Name")

library(WebGestaltR)
hypermethylated_a7_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                            "geneontology_Molecular_Function_noRedundant"),
                                         interestGeneType = "genesymbol",
                                         referenceGeneType = "genesymbol",
                                         interestGene = hypermethylated_A7_vs_WT_annotated$UCSC_RefGene_Name,
                                         referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                                         sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = hypermethylated_a7_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

ggplot(GSEA.chosen[1:20,], aes(x = `FDR`, 
                        y = fct_reorder(`description`, `FDR`))) +
  geom_point(aes(size = `size`, color = FDR),
             alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(4, 12)) +
  labs(
    title = "Hypermethylated regions in A7 vs WT",
    x = "FDR",
    y = "Term",
    color = "FDR"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15)
  )

hypomethylated_a7_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                           "geneontology_Molecular_Function_noRedundant"),
                                        interestGeneType = "genesymbol",
                                        referenceGeneType = "genesymbol",
                                        interestGene = hypomethylated_A7_vs_WT_annotated$UCSC_RefGene_Name,
                                        referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                                        sigMethod = "fdr", fdrThr = 0.1)
#####no significant term enriched

hypomethylated_C9_vs_WT_annotated <- hypomethylated_C9_vs_WT %>%
  left_join(genes_clean, by = "Name")

hypermethylated_C9_vs_WT_annotated <- hypermethylated_C9_vs_WT %>%
  left_join(genes_clean, by = "Name")

all_beta_values_annotated <- mean_beta_values %>%
  left_join(genes_clean, by = "Name")

library(WebGestaltR)
hypermethylated_C9_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                        "geneontology_Molecular_Function_noRedundant"),
                                     interestGeneType = "genesymbol",
                                     referenceGeneType = "genesymbol",
                                     interestGene = hypermethylated_C9_vs_WT_annotated$UCSC_RefGene_Name,
                                     referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                                     sigMethod = "fdr", fdrThr = 0.1)
GSEA.chosen = hypermethylated_C9_go  #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets

ggplot(GSEA.chosen[], aes(x = `FDR`, 
                               y = fct_reorder(`description`, `FDR`))) +
  geom_point(aes(size = `size`, color = FDR),
             alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(4, 12)) +
  labs(
    title = "Hypermethylated regions in C9 vs WT",
    x = "FDR",
    y = "Term",
    color = "FDR"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15)
  )

hypomethylated_c9_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                       "geneontology_Molecular_Function_noRedundant"),
                                    interestGeneType = "genesymbol",
                                    referenceGeneType = "genesymbol",
                                    interestGene = hypomethylated_C9_vs_WT_annotated$UCSC_RefGene_Name,
                                    referenceGene = all_beta_values_annotated$UCSC_RefGene_Name,
                                    sigMethod = "fdr", fdrThr = 0.1)
#one term found: organic acid binding?


######visualize overlaps
library(VennDiagram)
venn.diagram(list("Hypermethylated in A7" = hypermethylated_A7_vs_WT_annotated$Name,
                  "Hypermethylated in C9" = hypermethylated_C9_vs_WT_annotated$Name),
             lwd = 0, cex = 1.6, cat.cex = 1, cat.pos = c(-30,1), cat.dist = c(0.05, 0.147), print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#B66681", "#FF7F92"),
             "path\\DNAmethylation\\overlap_hypermethylated_regions.png", disable.logging = TRUE)


venn.diagram(list("Hypomethylated in A7" = hypomethylated_A7_vs_WT_annotated$Name,
                  "Hypomethylated in C9" = hypomethylated_C9_vs_WT_annotated$Name),
             lwd = 0, cex = 1.6, cat.cex = 1, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), fill = c("#e7c6ff", "#B9C0FF"),
             "path\\DNAmethylation\\overlap_hypomethylated_regions.png", disable.logging = TRUE)


