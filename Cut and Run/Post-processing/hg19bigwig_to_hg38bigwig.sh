#!/bin/bash
#this code is used to transform publicly available ChIPseq/Cut and Run datasets from hg19 to hg38 wihtout using UCSC

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_bw_file> <refs_directory>"
    exit 1
fi

# Assign input arguments to variables
input_bw_file=$1
refs_dir=$2

# Extract the base filename without the extension
base_name=$(basename "$input_bw_file" .bw)

# Convert BigWig to BedGraph
bigWigToBedGraph "$input_bw_file" "${base_name}.bedGraph"

# Perform liftover from hg19 to hg38
liftOver "${base_name}.bedGraph" "${refs_dir}/hg19ToHg38.over.chain" "${base_name}.hg38.bedGraph" unMapped

# Sort the resulting BedGraph file
sort -k1,1 -k2,2n "${base_name}.hg38.bedGraph" > "${base_name}.hg38.sort.bedGraph"

# Remove temporary files
rm unMapped
rm -f "${base_name}.hg38.bedGraph"

# Process BedGraph for partitioning and averaging
awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' "${base_name}.hg38.sort.bedGraph" > tmp.bed
bedops --partition tmp.bed | bedmap --echo --mean --delim '\t' - tmp.bed > "${base_name}.hg38.sort.split.bedGraph"

# Clean up temporary files
rm tmp.bed
rm -f "${base_name}.hg38.sort.bedGraph"

# Convert BedGraph back to BigWig
bedGraphToBigWig "${base_name}.hg38.sort.split.bedGraph" "${refs_dir}/chrom.sizes" "${base_name}.hg38.sort.split.bw"

# Final cleanup
rm -f "${base_name}.hg38.sort.split.bedGraph"

echo "Conversion completed: ${base_name}.hg38.sort.split.bw"