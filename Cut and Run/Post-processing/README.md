This is a folder containing all the scripts I have used for the analysis of the processed fastq files.

* unique and overlapping peaks have been defined using Upset plots
* annotation of peaks was done using ChIPSeeker
* normalization (RPKM) was done using subreadr package in Rstudio and the merged bam files; depending on the hypothesis, we used either the unique peaks (for the WT) or the unique and overlapped peaks between the two EZH2 KO clones 
* scatterplots for gene expression and PTMs were calculated using Pearson correlation regardless of where the peak is 
