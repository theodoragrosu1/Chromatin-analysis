#!/bin/bash


for file in ./ATACseq/*.mRp.clN.sorted.bam; do
    file_name=$(basename "$file")
    stripped_name="${file_name%.mRp.clN.sorted.bam}"
    TOBIAS ATACorrect --bam "$file" \
    --genome ./genome/genome.fa \
    --peaks ./hmmratac_peaks/"$stripped_name"_ATAC_peaks.bed \
    --blacklist ./hg38-blacklist.v2.bed \
    --outdir ataccorrect_output \
    --cores 22 #this can be changed based on the computer specs
done

echo "ATACorrect step finished"
