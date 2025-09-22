#create bigwigs
bamCoverage -b WT.mRp.clN.sorted.bam -o WT.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b A7.mRp.clN.sorted.bam -o A7.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.mRp.clN.sorted.bam -o C9.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

bamCoverage -b WT.mononucleosomal.sorted.bam -o WT.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b A7.mononucleosomal.sorted.bam -o A7.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.mononucleosomal.sorted.bam -o C9.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

bamCoverage -b WT.open.fragments.sorted.bam -o WT.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b A7.open.fragments.sorted.bam -o A7.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.open.fragments.sorted.bam -o C9.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

#map the transcript IDs of differentially expressed genes to GTF
grep -f ./transcript_ids_upreg_genes_a7.txt ./genome/hg38.knownGene.gtf > ./genome/a7_upreg_genes.knownGene.gtf
grep -f ./transcript_ids_upreg_genes_c9.txt ./genome/hg38.knownGene.gtf > ./genome/c9_upreg_genes.knownGene.gtf
grep -f ./transcript_ids_downreg_genes_a7.txt ./genome/hg38.knownGene.gtf > ./genome/a7_downreg_genes.knownGene.gtf
grep -f ./transcript_ids_downreg_genes_c9.txt ./genome/hg38.knownGene.gtf > ./genome/c9_downreg_genes.knownGene.gtf


#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw A7.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./a7_downreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_a7_downreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_a7_downreg_genes_atac.mat.gz \
    -out ./heatmap_open_fragments_a7_downreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT A7 \
    --regionsLabel "Genes downregulated in A7" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Blues \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10

##############################################


#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw A7.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./a7_upreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_a7_upreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c5_upreg_genes_atac.mat.gz \
    -out ./heatmap_open_fragments_c5_upreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT A7 \
    --regionsLabel "Genes upregulated in A7" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Reds \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10

#############################################



#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C9.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./c9_downreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c9_downreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c9_downreg_genes_atac.mat.gz \
    -out ./heatmap_open_fragments_c9_downreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C9 \
    --regionsLabel "Genes downregulated in C9" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Blues \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10


#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C9.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./c9_upreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c9_upreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c9_upreg_genes_atac.mat.gz \
    -out .//heatmap_open_fragments_c9_upreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C9 \
    --regionsLabel "Genes upregulated in C9" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Reds \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10


#A7 and C9 can also be combined in the same matrix, in that case the reference point for the matrix construction was the common IDs (upregulated and downregulated transcripts) between the two cell lines from the RNAseq, whcih have to be mapped to .gtf first and then matrix can be constructed.