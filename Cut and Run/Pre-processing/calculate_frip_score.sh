#!/bin/bash


# --- Paths
BAM_DIR="./deduped"
PEAKS_DIR="./macs2_peaks"
OUTPUT_FILE="./FRiP_scores.csv"

# --- Create output header
echo "Sample,Total_Reads,Reads_in_Peaks,FRiP" > "$OUTPUT_FILE"

# --- Loop through all BAM files (recursively)
find "$BAM_DIR" -type f -name "*.bam" | while read BAM; do

    SAMPLE=$(basename "$BAM" ".bam")

    # Find corresponding peak file
    PEAK_FILE=$(find "$PEAKS_DIR" -type f -name "${SAMPLE}*.bed" | head -n 1)

    if [[ -z "$PEAK_FILE" ]]; then
        echo "⚠️ No peak file found for $SAMPLE — skipping"
        continue
    fi

    echo "=== Calculating FRiP for $SAMPLE ==="

    # Count total mapped reads
    TOTAL_READS=$(samtools view -c -F 4 "$BAM")

    # Count reads overlapping peaks
    READS_IN_PEAKS=$(
        bedtools intersect -u -a "$BAM" -b "$PEAK_FILE" \
        | samtools view -c -
    )

    # Compute FRiP score
    if [[ "$TOTAL_READS" -gt 0 ]]; then
        FRIP_SCORE=$(awk -v a="$READS_IN_PEAKS" -v b="$TOTAL_READS" \
                      'BEGIN{printf "%.6f", a/b}')
    else
        FRIP_SCORE="NA"
    fi

    # Write to CSV
    echo "${SAMPLE},${TOTAL_READS},${READS_IN_PEAKS},${FRIP_SCORE}" >> "$OUTPUT_FILE"

    echo "✅ $SAMPLE → FRiP = $FRIP_SCORE"

done

echo "=== All FRiP scores saved to $OUTPUT_FILE ==="
