#this has to be run on the server if ran at the same time 

#ATAC peaks at H3K27ac signal
computeMatrix reference-point -S ./WT_ac.bigwig ./A7_ac.bigwig ./C9_ac.bigwig \
--referencePoint center -R ./ATAC-peaks/homer_input/WT_unique.bed ./ATAC-peaks/homer_input/A7_unique.bed ./ATAC-peaks/homer_input/C9_unique.bed \
-p max/2 \
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
-o ./heatmap_ac_signal_at_atac_unique_peaks.mat.gz

plotHeatmap -m ./heatmap_ac_signal_at_atac_unique_peaks.mat.gz \
-o ./heatmap_ac_signal_at_atac_peaks.pdf \     
--refPointLabel=TSS \    
--samplesLabel "WT H3K27ac" "A7 H3K27ac"  "C9 H3K27ac" \  
--regionsLabel "WT ATAC peaks" "A7 ATAC peaks" "C9 ATAC peaks" \ 
--missingDataColor=1 \
--xAxisLabel "Distance (kb)" \    
--colorMap Blues Blues Blues \
--zMax 0.3 0.3 0.3  \ 
--dpi=800  \   
--heatmapWidth=8 \
--heatmapHeight=12

#ATAC open peaks at Jurkat WT H3K27me3 signal
computeMatrix reference-point -S ./WT_me3.bigwig \ 
--referencePoint center -R ./ATAC-peaks/homer_input/WT_unique.bed ./ATAC-peaks/homer_input/A7_unique.bed ./ATAC-peaks/homer_input/C9_unique.bed \ 
-p max/2 \ 
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \ 
-bs=100 \ 
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_me3_signal_at_atac_unique_peaks.mat.gz

plotHeatmap -m ./heatmap_me3_signal_at_atac_unique_peaks.mat.gz \ 
-o ./heatmap_me3_signal_at_atac_peaks.pdf  \   
--refPointLabel=TSS  \   
--samplesLabel "WT H3K27me3"  \  
--regionsLabel "WT ATAC peaks" "A7 ATAC peaks" "C9 ATAC peaks" \ 
--missingDataColor=1   \  
--xAxisLabel "Distance (kb)"  \   
--colorMap Greens Greens Greens \ 
--zMax 0.12  \   
--dpi=800   \  
--heatmapWidth=8 \ 
--heatmapHeight=12

#ATAC open peaks at Jurkat H2AK119ub signal
computeMatrix reference-point -S ./WT_ub.bigwig ./A7_ub.bigwig ./C9_ub.bigwig \ 
--referencePoint center -R ./ATAC-peaks/homer_input/WT_unique.bed ./ATAC-peaks/homer_input/A7_unique.bed ./ATAC-peaks/homer_input/C9_unique.bed \ 
-p max/2 \ 
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \ 
-bs=100 \ 
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_ub_signal_at_atac_unique_peaks.mat.gz

plotHeatmap -m ./heatmap_ub_signal_at_atac_unique_peaks.mat.gz \ 
-o ./heatmap_ub_signal_at_atac_peaks.pdf  \   
--refPointLabel=TSS  \   
--samplesLabel "WT H2AK119ub" "A7 H2AK119ub" "C9 H2AK119ub"  \  
--regionsLabel "WT ATAC peaks" "A7 ATAC peaks" "C9 ATAC peaks" \ 
--missingDataColor=1   \  
--xAxisLabel "Distance (kb)"  \   
--colorMap Purples Purples Purples \ 
--zMax 0.3 0.3 0.3  \   
--dpi=800   \  
--heatmapWidth=8 \ 
--heatmapHeight=12

#H3K27ac signal at known genes
computeMatrix reference-point -S ./WT_ac.bigwig ./A7_ac.bigwig ./C9_ac.bigwig \
--referencePoint TSS \
-R ./public_datasets/genome_files/hg38-refGene.gtf \
-p max/2 \
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_ac_signal_at_known_genes.mat.gz

plotHeatmap -m ./heatmap_ac_signal_at_known_genes.mat.gz \    
-out ./heatmap_tss_ac_cut_run.pdf  \   
--refPointLabel=TSS  \   
--samplesLabel "WT H3K27ac" "A7 H3K27ac" "C9 H3K27ac"   \   
--regionsLabel "Genes"   \  
--missingDataColor=1 \    
--legendLocation "best"  \   
--xAxisLabel "Distance (kb)" \    
--colorMap Blues Blues Blues \    
--zMax 0.4 0.4 0.4   \  
--dpi=800   \  
--heatmapWidth=4  \   
--heatmapHeight=8

#H3K27me3 signal at known genes
computeMatrix reference-point -S ./WT_me3.bigwig ./A7_me3.bigwig ./C9_me3.bigwig \
--referencePoint TSS \
-R ./public_datasets/genome_files/hg38-refGene.gtf \
-p max/2 \
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_me3_signal_at_known_genes.mat.gz

plotHeatmap -m ./heatmap_me3_signal_at_known_genes.mat.gz \    
-out ./heatmap_tss_me3_cut_run.pdf  \   
--refPointLabel=TSS  \   
--samplesLabel "WT H3K27me3" "A7 H3K27me3" "C9 H3K27me3"   \   
--regionsLabel "Genes"   \  
--missingDataColor=1 \    
--legendLocation "best"  \   
--xAxisLabel "Distance (kb)" \    
--colorMap Greens Greens Greens \    
--zMax 0.15 0.15 0.15   \  
--dpi=800   \  
--heatmapWidth=4  \   
--heatmapHeight=8

#H2AK119ub signal at known genes
computeMatrix reference-point -S ./WT_ub.bigwig ./A7_ub.bigwig ./C9_ub.bigwig \
--referencePoint TSS \
-R ./public_datasets/genome_files/hg38-refGene.gtf \
-p max/2 \
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_ub_signal_at_known_genes.mat.gz

plotHeatmap -m ./heatmap_ub_signal_at_known_genes.mat.gz \    
-out ./heatmap_tss_ub_cut_run.pdf  \   
--refPointLabel=TSS  \   
--samplesLabel "WT H2AK119ub" "A7 H2AK119ub" "C9 H2AK119ub"   \   
--regionsLabel "Genes"   \  
--missingDataColor=1 \    
--legendLocation "best"  \   
--xAxisLabel "Distance (kb)" \    
--colorMap Purples Purples Purples \    
--zMax 0.4 0.4 0.4   \  
--dpi=800   \  
--heatmapWidth=4  \   
--heatmapHeight=8

#H3K27ac peaks at Myb signal 
computeMatrix reference-point -S ./GSM4420242_TL-1-PJ-Myb_S6.bcov.hg38.sort.split.bw ./GSM4420243_TL-1-EZH2-Myb_S2.bcov.hg38.sort.split.bw \
--referencePoint center -R ./WT_ac_peaks.bed ./A7_ac_peaks.bed ./C9_ac_peaks.bed \
-p max/2 \
-bl ./genome_files/hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
-o ./heatmap_myb_signal_at_ac_peaks.mat.gz

plotHeatmap -m ./heatmap_myb_signal_at_ac_peaks.mat.gz   \ 
-out ./heatmap_myb_signal_ac_peaks.pdf  \   
--refPointLabel=center   \ 
--samplesLabel "WT Myb" "A7 Myb" \ 
--regionsLabel "WT H3K27ac" "A7 H3K27ac" "C9 H3K27ac" \ 
--missingDataColor=1  \   
--legendLocation "best"   \  
--xAxisLabel "Distance (kb)"  \   
--colorMap Reds Reds  \  
--zMax 0.5 0.5  \ 
--dpi=800   \  
--heatmapWidth=4  \   
--heatmapHeight=12

#H3K27ac peaks at public Jurkat H3K27ac signal (to validate our peaks)
computeMatrix reference-point -S ./GSM4851798_Jurkat-H3K27ac_deeptool_normalized.hg38.sort.split.bw \
--referencePoint center -R ./WT_ac_peaks.bed ./A7_ac_peaks.bed ./C9_ac_peaks.bed \
-p max/2 \ 
-bl ./hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_public_ac_signal_at_ac_peaks.mat.gz

plotHeatmap -m ./heatmap_public_ac_signal_at_ac_peaks.mat.gz  \   
-out ./heatmap_public_ac_signal_ac_cut_run.pdf \    
--refPointLabel=Center  \   
--samplesLabel "GSM4851798 Jurkat H3K27ac" \ 
--regionsLabel "WT H3K27ac" "A7 H3K27ac" "C9 H3K27ac"  \   
--missingDataColor=1  \   
--legendLocation "best"  \   
--xAxisLabel "Distance (kb)"  \   
--colorMap Blues Blues Blues  \  
--zMax 3.5 3.5 3.5 \ 
--dpi=800   \
--heatmapWidth=4 \    
--heatmapHeight=8

#H3K27me3 peaks at public Jurkat H3K27me3 signal (to validate our peaks)
computeMatrix reference-point -S ./GSM4851799_Jurkat-H3K27me3_deeptool_normalized.hg38.sort.split.bw \
--referencePoint center -R ./WT_me3_peaks.bed ./A7_me3_peaks.bed ./C9_me3_peaks.bed \
-p max/2 \ 
-bl ./hg38-blacklist.v2.bed \
-bs=100 \
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_public_me3_signal_at_me3_peaks.mat.gz

plotHeatmap -m ./heatmap_public_me3_signal_at_me3_peaks.mat.gz  \   
-out ./heatmap_public_me3_signal_me3_cut_run.pdf \    
--refPointLabel=Center  \   
--samplesLabel "GSM4851799 Jurkat H3K27me3" \ 
--regionsLabel "WT H3K27me3" "A7 H3K27me3" "C9 H3K27me3"  \   
--missingDataColor=1  \   
--legendLocation "best"  \   
--xAxisLabel "Distance (kb)"  \   
--colorMap Greens Greens Greens  \  
--zMax 0.15 0.15 0.15 \ 
--dpi=800   \
--heatmapWidth=4 \    
--heatmapHeight=8

#Polycomb signal (H3K27me3 peaks intersected with H2AK119ub peaks) at H3K27me3 and H2AK119ub signal
computeMatrix reference-point -S ./WT_me3.bigwig ./A7_me3.bigwig ./C9_me3.bigwig \ 
--referencePoint center -R ./WT_PRC1_PRC2_loci.bed \ 
-p max/2 \ 
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \ 
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_PcG_at_me3_signal.mat.gz

computeMatrix reference-point -S ./WT_ub.bigwig ./A7_ub.bigwig ./C9_ub.bigwig \ 
--referencePoint center -R ./WT_PRC1_PRC2_loci.bed \ 
-p max/2 \ 
-bl ./public_datasets/genome_files/hg38-blacklist.v2.bed \
-bs=100 \ 
--beforeRegionStartLength 3000 \ 
--afterRegionStartLength 3000 \ 
-o ./heatmap_PcG_at_ub_signal.mat.gz

plotHeatmap -m ./heatmap_PcG_at_me3_signal.mat.gz  \  
-out ./heatmap_me3_signal_PcG_peaks.pdf  \   
--refPointLabel=center  \  
--samplesLabel "WT H3K27me3" "A7 H3K27me3" "C9 H3K27me3"  \
--regionsLabel "WT PcG loci" \ 
--missingDataColor=1  \   
--legendLocation "best"   \  
--xAxisLabel "Distance (kb)"  \   
--colorMap Greens Greens Greens  \  
--zMax 0.5 0.5 0.5  \  
--dpi=800  \   
--heatmapWidth=8  \   
--heatmapHeight=12

plotHeatmap -m ./heatmap_PcG_at_ub_signal.mat.gz  \  
-out ./heatmap_ub_signal_PcG_peaks.pdf  \   
--refPointLabel=center  \  
--samplesLabel "WT H2AK119ub" "A7 H2AK119ub" "C9 H2AK119ub"  \
--regionsLabel "WT PcG loci" \ 
--missingDataColor=1  \   
--legendLocation "best"   \  
--xAxisLabel "Distance (kb)"  \   
--colorMap Purples Purples Purples  \  
--zMax 0.5 0.5 0.5  \  
--dpi=800  \   
--heatmapWidth=8  \   
--heatmapHeight=12



