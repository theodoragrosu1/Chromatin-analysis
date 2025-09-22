#!/bin/bash

#output: indexed and sorted .bigwigs

cd ./ATAC_Seq/bwa/merged_replicates/

mkdir open_fragments
mkdir mononucleosomal
    
samples=(
    "WT" 
    "A7" 
    "C9"
)

for sample in "${samples[@]}"; do
	echo "Sieving sample $counter: $sample"
    
    #open fragments
    alignmentSieve -b "${sample}.mRp.clN.sorted.bam" \
        --ATACshift \
        -p max/2 \
        --maxFragmentLength 100 \
        -o "./open_fragments/${sample}.open.fragments.bam"
    
    #mononucleosomal fragments
    alignmentSieve -b "${sample}.mRp.clN.sorted.bam" \
        --ATACshift \
        -p max/2 \
        --minFragmentLength 180 \
        --maxFragmentLength 250 \
        -o "./mononucleosomal/${sample}.mononucleosomal.bam"

    ((counter++))
done

echo "Finished sieving"


#deal with open fragments files
cd ./open_fragments

#sorting

for sample in "${samples[@]}"; do
    samtools sort "${sample}.open.fragments.bam" -o "${sample}.open.fragments.sorted.bam"
done

#indexing
for sample in "${samples[@]}"; do
	samtools index "${sample}.open.fragments.sorted.bam"
done
#deal with mononucleosomal files

cd ../mononucleosomal

#sorting
for sample in "${samples[@]}"; do
    samtools sort "${sample}.mononucleosomal.bam" -o "${sample}.mononucleosomal.sorted.bam"
done

cd ./mononucleosomal

#indexing
for sample in "${samples[@]}"; do
	samtools index "${sample}.mononucleosomal.sorted.bam"
done

echo "Finished sorting and indexing"

cd ..

for sample in "${samples[@]}"; do
	echo "Creating bigwig for sample $counter: $sample"
	bamCoverage  -b"./open_fragments/${sample}.open.fragments.sorted.bam" -o "./open_fragments/${sample}.open.fragments.bigwig" \
		--binSize 1 \
		--normalizeUsing RPGC \
		--effectiveGenomeSize 2913022398

	bamCoverage -b "./mononucleosomal/${sample}.mononucleosomal.sorted.bam" -o "./mononucleosomal/${sample}.mononucleosomal.bigwig" \
		--binSize 1 \
		--normalizeUsing RPGC \
		--effectiveGenomeSize 2913022398

    ((counter++))
done

echo "Finished making bigiwgs"