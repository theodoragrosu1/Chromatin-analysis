#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Spike-in alignment with Bowtie2 (yeast sacCer); this needs to be indexed first!!!! Use bowtie2-build function on .fa file
###############################################################################

THREADS=8
YEAST_IDX="/mnt/Data/Bond/tgrosu/cutrun/SacCer_index/SacCer_index"

FASTQ_DIR="fastq_trimmed"
SAM_DIR="sams_spikein"
STATS_DIR="stats_spikein"


echo "========================================"
echo "Spike-in alignment (yeast)"
echo "========================================"

for R1 in ${FASTQ_DIR}/*_R1_paired.fq.gz; do

    SAMPLE=$(basename "$R1" _R1_paired.fq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_R2_paired.fq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Missing R2 for $SAMPLE — skipping"
        continue
    fi

    SAM_OUT="${SAM_DIR}/spikein_${SAMPLE}.sam"
    STATS_OUT="${STATS_DIR}/stats_spikein_${SAMPLE}_bowtie2.txt"

    echo ""
    echo "Aligning spike-in for sample: $SAMPLE"
    echo "  R1: $R1"
    echo "  R2: $R2"
    echo "  Output SAM: $SAM_OUT"
    echo "  Stats: $STATS_OUT"

    bowtie2 \
        --end-to-end \
        --dovetail \
        -I 10 \
        -X 700 \
        --no-mixed \
        --no-discordant \
        --threads "$THREADS" \
        -x "$YEAST_IDX" \
        -1 "$R1" \
        -2 "$R2" \
        -S "$SAM_OUT" \
        &> "$STATS_OUT"

done

echo ""
echo "========================================"
echo "Spike-in alignment complete"
echo "========================================"
