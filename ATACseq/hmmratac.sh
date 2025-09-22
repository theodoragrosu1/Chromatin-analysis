#!/bin/bash

#input files: .bam files after initial QC, alignment to hg38, and duplicate removal
#output files: peak files for each condition without blacklisted regions

cd ./ATAC_Seq/bwa/merged_replicate

samtools view -H WT.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > WT_genome.info 
samtools view -H A7.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > A7_genome.info 
samtools view -H C9.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > C9_genome.info 

echo "Finished genome info files"

mkdir hmmratac_peaks

java -jar HMMRATAC_V1.2.10_exe.jar -b WT.mRp.clN.sorted.bam -i WT.mRp.clN.sorted.bam.bai -g WT_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/WT_peaks
echo "Finished WT peak calling"

java -jar HMMRATAC_V1.2.10_exe.jar -b C5.mRp.clN.sorted.bam -i A7.mRp.clN.sorted.bam.bai -g A7_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/A7_peaks
echo "Finished A7 peak calling"

java -jar HMMRATAC_V1.2.10_exe.jar -b C9.mRp.clN.sorted.bam -i C9.mRp.clN.sorted.bam.bai -g C9_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/C9_peaks
echo "Finished C9 peak calling"