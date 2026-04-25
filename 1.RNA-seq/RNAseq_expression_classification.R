#input: all omics datasets to cleanup code and extract gene lists for DNA methylation
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase)
library(gprofiler2) 
library(gplots)
library(enrichplot)
library(ggrepel)
library(dplyr)
library(biomaRt)
library(rtracklayer)

#input 1: start with RNAseq: divide genes in differentially expressed, just expressed, and not expressed at all
#for (differentially) expressed genes, these are found in "normalised_expression_wt_vs_a7/c9.csv"
#rule: to define a differentially expressed gene in either WT or clone, logFC > 0.5, the rest are expressed but not differentially

###########differentially expressed genes; divide the conditions from the start in upregulated and downregulated based on the clone you are looking at 
###########otherwise it is going to skew the ATACseq interpretation; peak in promoter should always mean gene expression, no matter if it is differential or not
upregulated_a7 <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_a7_all_genes.csv") %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::select(geneID)

# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Example gene names
gene_names_a7_up <- upregulated_a7

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_a7_up,
  mart = ensembl
)
write.table(results, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_a7.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



upregulated_c9 <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_c9_all_genes.csv") %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::select(geneID)

gene_names_c9_up <- upregulated_c9

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_c9_up,
  mart = ensembl
)


write.table(results, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_c9.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_a7 <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_a7_all_genes.csv") %>%
  dplyr::filter(logFC < -0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::select(geneID)
downregulated_c9 <- read.csv("C:\\Users\\tgrosu\\Downloads\\Jurkat_wt_v_c9_all_genes.csv") %>%
  dplyr::filter(logFC < -0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::select(geneID)

upregulated_wt <- left_join(downregulated_a7, downregulated_c9, by = "geneID")

gene_names_upregulated_wt <- upregulated_wt

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_upregulated_wt,
  mart = ensembl
)
write.table(results, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_wt.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


gtf_file <- "C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\genome_files\\hg38.knownGene.gtf"
gtf_data <- readGFF(gtf_file)

transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_a7.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id

# Step 2: Load full GTF annotation
gtf <- import("C:\\Users\\tgrosu\\OneDrive\\Desktop\\public_datasets\\genome_files\\hg38.knownGene.gtf")  # Replace with your GTF file

# Step 3: Strip version from transcript_id in GTF (if needed)
mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)

# Step 4: Filter GTF by transcript IDs
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]

# Step 5: Export filtered GTF
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_a7.gtf", format = "gtf")

transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_c9.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id
mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_c9.gtf", format = "gtf")

#we can get the upregulated in wt vector, by joining downregulated in a7 and downregulated in c9 joined
transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_wt.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id
mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\upregulated_in_wt.gtf", format = "gtf")


############genes that do not meet the LogFC threshold, are expressed in all conditions but not differentially, so we can keep them 
############select the genes that are not in the upregulated or the downregulated vectors
expr_a7 <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_a7.csv") 
expr_a7 <- anti_join(expr_a7, upregulated_a7, by = "geneID")
expr_a7 <- anti_join(expr_a7, upregulated_wt, by = "geneID")

gene_names_expr_a7 <- expr_a7

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_expr_a7,
  mart = ensembl
)
write.table(results, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_a7.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_a7.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id

mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_a7.gtf", format = "gtf")


expr_c9 <- read.csv("C:\\Users\\tgrosu\\Downloads\\normalised_expression_jurkat_wt_vs_c9.csv") 
expr_c9 <- anti_join(expr_c9, upregulated_c9, by = "geneID")
expr_c9 <- anti_join(expr_c9, upregulated_wt, by = "geneID")
gene_names_expr_c9 <- expr_c9

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_expr_c9,
  mart = ensembl
)

write.table(results, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_c9.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_c9.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id
mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\expr_c9.gtf", format = "gtf")

########the not at all expressed genes according to RNAseq

### read data----
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
design <- read_tsv("C:\\Users\\tgrosu\\OneDrive\\Desktop\\abundance\\study_design_jurkat.txt") 
path <- file.path("C:\\Users\\tgrosu\\OneDrive\\Desktop\\abundance\\", design$sample, "abundance.tsv") # set file paths to mapped data

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) %>% #gene symbols 
  as_tibble() %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select("target_id", "gene_name")

Tx_gene <- tximport(path, #imports the data for all samples
                    type = "kallisto",
                    tx2gene = Tx,
                    txOut = FALSE, #data represented at gene level rather than transcript
                    countsFromAbundance = "lengthScaledTPM", #transcripts per million
                    ignoreTxVersion = TRUE) 
### preprocessing----
sample_labels <- design$sample
DGEList <- DGEList(Tx_gene$counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.df) <- c("geneID", sample_labels)

###tidy data
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                  names_to = "samples", # name of new column
                                  values_to = "expression") # name of new column storing all the data
###filter data
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=3 #user defined - depends on studydesign, eliminates lowly expressed genes
DGEList.filtered <- DGEList[keepers,]
###normalize data
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>% 
  cpm(log=TRUE) %>% 
  as_tibble(rownames = "geneID")

colnames(log2.cpm.filtered.norm) <- c("geneID", sample_labels)

cpm.df <- data.frame(geneID = rownames(cpm), cpm, check.names = FALSE, stringAsFactors = FALSE)
not_expressed <- anti_join(cpm.df, log2.cpm.filtered.norm, by = "geneID")

gene_names_not_expressed <- not_expressed

# Query transcript IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gene_names_not_expressed,
  mart = ensembl
)

write.table(results, file = "C:\\Users\\tgrosu\\OneDrive\\Desktop\\not_expressed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
transcripts <- read.table("C:\\Users\\tgrosu\\OneDrive\\Desktop\\not_expressed.tsv", header = TRUE, stringsAsFactors = FALSE) 
transcript_ids <- transcripts$ensembl_transcript_id
mcols(gtf)$transcript_id_clean <- sub("\\..*$", "", mcols(gtf)$transcript_id)
filtered_gtf <- gtf[mcols(gtf)$transcript_id_clean %in% transcript_ids]
export(filtered_gtf, "C:\\Users\\tgrosu\\OneDrive\\Desktop\\not_expressed.gtf", format = "gtf")

