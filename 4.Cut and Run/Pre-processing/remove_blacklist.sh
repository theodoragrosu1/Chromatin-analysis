#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Remove hg38 blacklist regions from BAM files
################################################################################

DEDUP_DIR="./deduped"
BLACKLIST="./genome_files/hg38-blacklist.v2.bed"

OUT_DIR="./deduped/blacklist_filtered"

shopt -s nullglob

for BAM in "$DEDUP_DIR"/*_deduped.bam; do
    SAMPLE=$(basename "$BAM" _deduped.bam)
    OUT_BAM="${OUT_DIR}/${SAMPLE}_deduped_blacklist_filtered.bam"

    echo "→ Processing $SAMPLE"
    bedtools intersect -v -abam "$BAM" -b "$BLACKLIST" > "$OUT_BAM"
done

echo "=== Blacklist filtering complete ==="
