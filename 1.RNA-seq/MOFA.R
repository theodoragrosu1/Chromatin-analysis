library(MOFA2)
library(tidyverse)
library(dplyr)

rna_normalised <- read.csv("path\\normalised_expression_jurkat_wt_vs_clones (1).csv", row.names = 1)
colnames(rna_normalised) <- c("A7_1", "WT_1", "WT_2", "WT_3", 
                               "A7_2", "A7_3", "C9_1", "C9_2", "C9_3")
groups <- sub("_.*", "", colnames(rna_normalised))
groups
rna_normalised_avg <- sapply(unique(groups), function(g) {
  rowMeans(rna_normalised[, groups == g])
})

rna_normalised_avg %>%
  as.matrix()
rna_normalised_avg <- rna_normalised_avg[, c("WT", "A7", "C9")]

methylation_beta_values <- read.csv("path\\Mean_Beta_values_Sesame_QCDPB_Jurkat_cell_line_930659CpGs.csv", row.names=1) %>%
  as.matrix()
colnames(methylation_beta_values) <- c("WT", "A7", "C9")

#ATAC peaks need to be first CPM-normalised
library(edgeR)
counts_raw <- read.table("path\\ATAC-peaks\\homer_input\\peak_counts.txt", header=TRUE, row.names=1, skip=1)
counts <- counts_raw[,6:ncol(counts_raw)]  # keep only numeric counts
colnames(counts) <- c("WT", "A7", "C9")

# Create DGEList
dge <- DGEList(counts=counts)

# CPM normalization (log2(CPM+1))
cpm_counts <- cpm(dge, log=TRUE)

# Optional: filter low-count peaks
keep <- rowSums(cpm_counts > 1) >= 2
cpm_filtered <- cpm_counts[keep,]
#Save peak info for biological interpetation
peak_info <- counts_raw[keep, 1:6]

# Save for MOFA
write.table(cpm_filtered, "path\\ATAC-peaks\\homer_input\\peak_counts_cpm_log2_filtered.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(peak_info, "path\\ATAC-peaks\\homer_input\\peak_coordinates_filtered.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Filter the top 2000 genes because MOFA could become unstable otherwise; we have over 13000 in the bgeinning
gene_var <- apply(rna_normalised_avg,1,var)
rna_top <- rna_normalised_avg[order(gene_var, decreasing=TRUE)[1:2000], ]

#Do the same filtering with the DNA methylation probes
cpg_var <- apply(methylation_beta_values,1,var)
meth_top <- methylation_beta_values[order(cpg_var, decreasing=TRUE)[1:2000], ]

#Do the same filtering with the ATAC peaks
atac <- read.table("path\\ATAC-peaks\\homer_input\\peak_counts_cpm_log2_filtered.txt", row.names=1)
atac <- as.matrix(atac)

var_peaks <- apply(atac,1,var)
atac_top <- atac[order(var_peaks,decreasing=TRUE)[1:2000],]

#Build input object for MOFA
data_list <- list(
  RNAseq = rna_top,
  ATAC = atac_top,
  Methylation = meth_top
)

MOFAobject <- create_mofa(data_list)

sample_metadata <- data.frame(
  sample = colnames(rna_normalised_avg),
  genotype = c("WT","A7","C9")
)

samples_metadata(MOFAobject) <- sample_metadata


data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 1   # VERY important: 1 factor for 3 samples

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"   # slow is safer with tiny datasets
train_opts$stochastic <- FALSE
train_opts$maxiter <- 10000

MOFAobject <- prepare_mofa(
  MOFAobject,
  model_options = model_opts,
  training_options = train_opts,
  data_options = data_opts
)

model <- run_mofa(MOFAobject, use_basilisk = TRUE)
plot_variance_explained(model)
plot_weights(model, view="RNAseq", factor=1)
plot_weights(model, view="Methylation", factor=1)
plot_weights(model, view="ATAC", factor=1)

# Get weights for factor 1
weights_list <- get_weights(model)
weights_rna <- weights_list$RNAseq[, 1]
weights_atac <- weights_list$ATAC[, 1]
weights_meth <- weights_list$Methylation[, 1]

top_genes <- names(sort(abs(weights_rna), decreasing = TRUE)[1:50])
top_peaks <- names(sort(abs(weights_atac), decreasing = TRUE)[1:50])
top_cpgs  <- names(sort(abs(weights_meth), decreasing = TRUE)[1:50])
