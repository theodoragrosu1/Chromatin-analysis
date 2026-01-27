#input file: raw counts RNAseq (Polonen et al, 2024)
#input file: clinical data from sequenced patients (Polonen et al, 2024)
#output file: annotated patients from PRC2 alterations
#output file: boxplots for gene expression based on disease status
#output: Wilcoxon test for all comparisons of gene expression between patient groups

library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase)
library(Biobase)
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(gplots)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(tximport)
library(ggrepel)


cnv_samples <- read_xlsx("#path to clinical data", sheet = 20) %>%
  dplyr::filter(chrom %in% c("chr7", "chr11", "chr17"))


IDs <- read_xlsx("#path to clinical data", sheet =2) %>%
  dplyr::select("USI", "ETP.STATUS") %>%
  dplyr::rename("sample"="USI") %>%
  dplyr::rename("ETP_status" = "ETP.STATUS") %>%
  drop_na()

#grch37
#ezh2 chr7:148504464-148581441
#eed chr11:85955806-85989785
#suz12 chr17:30264044-30328057

prc2_cna <- cnv_samples %>%
  dplyr::mutate(prc2_status = case_when(chrom == "chr7" & loc.start<148504464 & loc.end>148581441 & LogRatio < 0.7 ~ "ezh2_loss",
                                        chrom == "chr11" & loc.start<85955806 & loc.end>85989785 & LogRatio < 0.7 ~ "eed_loss",
                                        chrom == "chr17" & loc.start<30264044 & loc.end>30328057 & LogRatio < 0.7 ~ "suz12_loss"))


prc2_lost <- prc2_cna %>%
  dplyr::select(c(sample, prc2_status)) %>%
  dplyr::filter(prc2_status != "NA") %>%
  dplyr::left_join(IDs, by = "sample")


#mutations sheet14
prc2_mutations <- read_xlsx("#path to clinical data", sheet =15) %>%
  dplyr::filter(gene == "EZH2" | gene == "SUZ12" | gene == "EED") %>%
  dplyr::select(sample, gene, mutation_class, ExonicFunc) 


#read patient samples RNAseq raw counts
rna_seq <- read.csv("#path to raw counts", sep = '\t') %>%
  dplyr::rename("geneID" = "X")

rna_seq$geneID <- gsub("\\..*","",rna_seq$geneID)


annots <- select(org.Hs.eg.db, keys=rna_seq$geneID, 
                 columns="SYMBOL", keytype="ENSEMBL")

annots %>%
  drop_na()

annots <- annots %>%
  dplyr::rename("geneID" = "ENSEMBL")

rna_seq_counts <- left_join(rna_seq, annots, by = "geneID")

rna_seq_counts<- rna_seq_counts[,-1]


#R crashes when grouping and summarizing across 6000 genes, so simplify by only summarizing genes that appear twice
summarized_data <- rna_seq_counts %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

summarized_data %>%
  drop_na()
single_row_groups <- rna_seq_counts %>%
  group_by(SYMBOL) %>%
  dplyr::filter(n() == 1) %>%
  ungroup()

rna_seq_counts_filtered <- bind_rows(summarized_data, single_row_groups) %>%
  dplyr::filter(duplicated(SYMBOL))
rownames(rna_seq_counts_filtered) <- NULL
rownames(rna_seq_counts_filtered) <- rna_seq_counts_filtered$'SYMBOL'
rna_seq_counts_filtered <- rna_seq_counts_filtered[,-1] %>%
  as.data.frame()



t_all_rna_seq_samples <- rna_seq_counts_filtered %>%
  colnames()

rna_seq_samples <- IDs %>%
  dplyr::filter(sample %in% t_all_rna_seq_samples) %>%
  dplyr::pull("sample") %>%
  unique()

samples_intersection <- intersect(intersect(cnv_samples, prc2_mutations), rna_seq_samples)


####make study design file
study_design <- IDs %>%
  dplyr::select(sample, ETP_status) %>%
  dplyr:: mutate(prc2_status = case_when(sample %in% unique(c(prc2_mutations$sample, prc2_lost$sample)) ~ "prc2_mut", 
                                         TRUE ~ "prc2_wt")) %>%
  dplyr::mutate(prc2_loss_type = case_when(sample %in% setdiff(prc2_mutations$sample, prc2_lost$sample) ~ "prc2_sequence_mut", 
                                           sample %in% setdiff(prc2_lost$sample, prc2_mutations$sample) ~ "prc2_cn_alteration", 
                                           sample %in% intersect(prc2_lost$sample, prc2_mutations$sample) ~ "prc2_sequence_mut_&_cn_alteration",
                                           TRUE ~ "prc2_wt"))




  CD1_t_all_rna_seq <- rna_seq_counts %>%
    dplyr::filter(SYMBOL %in% c("CD1A", "CD1B", "CD1C", "CD1D", "CD1E")) %>%
    #dplyr::select(-SYMBOL) %>%
    dplyr::select(any_of(c("SYMBOL", study_design$sample))) %>%
    pivot_longer(cols = -SYMBOL, values_to = "RNA_expression", names_to = "patient_id") %>%
    mutate(PRC2_status = case_when(patient_id %in% pull(study_design[study_design[, "prc2_status"] 
                                                                     == "prc2_mut", "sample"], "sample")
                                   
                                   ~ "prc2_mut", 
                                   patient_id %in% pull(study_design[study_design[, "prc2_status"] 
                                                                     == "prc2_wt", "sample"], "sample")
                                   ~ "prc2_wt")) %>%
    mutate(ETP_status = case_when(patient_id %in% pull(study_design[study_design[, "ETP_status"] == "ETP", "sample"], "sample") ~ "ETP",
                                  patient_id %in% pull(study_design[study_design[,"ETP_status"] == "Non-ETP", "sample"], "sample") ~ "Non-ETP",
                                  patient_id %in% pull(study_design[study_design[, "ETP_status"] == "Near-ETP", "sample"], "sample") ~ "Near-ETP"))
  
  
  CD1_t_all_rna_seq %>%
    dplyr::filter(ETP_status != "NA") %>%
    ggplot(aes(x = ETP_status, y = log2(RNA_expression+1), fill = interaction(PRC2_status, ETP_status))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 16, position = position_jitterdodge(), 
               size = 1, alpha = 0.4) +
    scale_fill_manual(values = c("#80A1D4", "#DD2D4A", 
                                 "#80A1D4", "#DD2D4A", 
                                 "#80A1D4", "#DD2D4A")) +
    theme_bw(base_size = 20) +
    facet_wrap(~SYMBOL) +
    labs(title = "CD1 expression in paediatric T-ALL patients") +
    theme(legend.position = "none") +
    labs(x = "ETP status", y = "log2(tpm+1)") 


for (cd in c("CD1A", "CD1B", "CD1C", "CD1D", "CD1E")) {
  CD1_t_all_rna_seq_not_etp <- CD1_t_all_rna_seq %>%
    dplyr::filter(ETP_status %in% c( "ETP", "Non-ETP")) %>%
    dplyr::filter(SYMBOL == cd)
  print(cd)
  test_res_tmp <- wilcox.test(CD1_t_all_rna_seq_not_etp[CD1_t_all_rna_seq_not_etp$PRC2_status == "prc2_wt",]$RNA_expression, 
                              CD1_t_all_rna_seq_not_etp[CD1_t_all_rna_seq_not_etp$PRC2_status == "prc2_mut",]$RNA_expression, 
                              paired=FALSE)
  print(paste("Wilcoxon test p val = ", test_res_tmp[["p.value"]]))
}

  # Extract expression values
  ezh2 <- as.numeric(rna_seq_counts_filtered["EZH2", ])
  cd1a <- as.numeric(rna_seq_counts_filtered["CD1A", ])
  
  # Make scatterplot
  plot(
    x = ezh2,
    y = cd1a,
    xlab = "EZH2 expression",
    ylab = "CD1A expression",
    main = "Scatterplot of EZH2 vs CD1A",
    pch = 19, col = "blue"
  )
  
  # Add regression line (optional)
  abline(lm(cd1a ~ ezh2), col = "red")
  
