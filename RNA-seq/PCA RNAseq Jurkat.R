#This script only produces the PCA plot and performs checks on data normalisation (TMM and log2 TPM)
#Input: kallisto abundance for each sample - available on GEO
#Input: study design file available with the code
#Output: This script only produces the PCA plot and performs checks on data normalisation (TMM and log2 TPM)

library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(colorspace)

#function needed later
violin_plot_transcriptomics <- function(df, label) {
  plot <- ggplot(log2.cpm.df.pivot) +
    aes(x = samples, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE, fill = "lightblue") +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 colour = "black", 
                 show.legend = FALSE) +
    labs(y = "log2 expression", x = "sample",
         title = "Log2 Counts per Million (CPM)",
         subtitle = label) +
    theme_bw()
  return(plot)
}


### read data----
targets <- read_tsv("study_design_jurkat.txt") 
path <- file.path("#path to raw data folder", targets$sample, "abundance.tsv") #set file paths to mapped data

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) %>% #gene symbols
  as_tibble() %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select("target_id", "gene_name")

Txi_gene <- tximport(path, #imports the data for all samples
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE, #data represented at gene level rather than transcript
                     countsFromAbundance = "lengthScaledTPM", #transcripts per million
                     ignoreTxVersion = TRUE) 

### preprocessing----

sampleLabels <- targets$sample
DGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

###tidy data
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                  names_to = "samples", # name of new column
                                  values_to = "expression") # name of new column storing all the data


###plot log2 expression (unnormalized and unfiltered)
p1 <- violin_plot_transcriptomics(log2.cpm.df.pivot, "unfiltered, non-normalized")


###filter data
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=3 #user defined - depends on studydesign 
DGEList.filtered <- DGEList[keepers,] #eliminates lowly expressed genes

log2.cpm.filtered <- cpm(DGEList.filtered, log=TRUE) 
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                           names_to = "samples", # name of new column
                                           values_to = "expression") # name of new column storing all the data

###plot log2 expression (filtered and unnormalized)
p2 <- violin_plot_transcriptomics(log2.cpm.filtered.df.pivot, "filtered, non-normalized")

###normalize data 
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(DGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                                names_to = "samples", # name of new column
                                                values_to = "expression") # name of new column storing all the data

###plot normalized and filtered data
p3 <- violin_plot_transcriptomics(log2.cpm.filtered.norm.df.pivot, "filtered, TMM normalized")

###merge the three plots above, to compare
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

#save for PROGENy input 
write.csv(log2.cpm.filtered.norm.df, file = "#path and name to save normalised expression csv", row.names = FALSE)


### multivariate analysis ----
#gets logFC
data.df <- mutate(log2.cpm.filtered.norm.df,
                  jurkat.wt.AVG = (A13 + A11 + A12)/3,
                  jurkat.a7.AVG = (A10 + A14 + A15)/3,
                  jurkat.c9.AVG = (A16 + A17 + A18)/3,
                  #now make columns comparing each of the averages above
                  LogFC = (jurkat.c9.AVG - jurkat.wt.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

group <- targets$group
group <- factor(group)

### pca jurkat ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

#get colours for pca 
dark_okabe <- darken(c("#E69F00", "#56B4E9", "#009E73"), amount = 0.2) 

#plot pca
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, fill = group, colour = group, scale = TRUE) +
  geom_point(size=7, stroke = 1, shape = 21) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  coord_fixed() +
  theme_bw(base_size = 18) +
  scale_fill_manual(labels = c("A7", "C9", "WT"), values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(labels = c("A7", "C9", "WT"), values = dark_okabe) +
  labs(fill = "Group", colour = "Group")

pca.plot