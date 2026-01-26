BAM_DIR="./merged_bams"
PEAK_DIR="./macs2_peaks"

echo -e "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP" > ./frip_results.tsv

for bam in "${BAM_DIR}"/*_deduped.bam; do

    # extract sample ID (e.g. A7_me3_1)
    sample=$(basename "$bam" _deduped.bam)

    # construct full path to peak file
    peaks="${PEAK_DIR}/${sample}_IgG_broadpeaks_macs2.bed"

    # sanity check
    if [[ ! -f "$peaks" ]]; then
        echo "WARNING: Peak file not found for $sample"
        continue
    fi

    total=$(samtools view -c -F 4 "$bam")
    in_peaks=$(bedtools intersect -abam "$bam" -b "$peaks" -u | samtools view -c)

    frip=$(echo "scale=4; $in_peaks / $total" | bc)

    echo -e "$sample\t$total\t$in_peaks\t$frip" >> ./frip_results.tsv
done
