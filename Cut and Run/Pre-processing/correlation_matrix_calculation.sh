#!/usr/bin/env bash
#this code can be changed to include just specific PTMs by redefining the initial folder

# Folder containing all BigWig files
BIGWIG_DIR="./merged_bigwigs"
OUTDIR="./correlation"


echo "=== Collecting all BigWigs ==="
mapfile -t BW_FILES < <(find "$BIGWIG_DIR" -maxdepth 1 -type f -name "*.bigwig")

if [ ${#BW_FILES[@]} -eq 0 ]; then
    echo "No BigWig files found in $BIGWIG_DIR"
    exit 1
fi

echo "Found ${#BW_FILES[@]} BigWigs"

# Output matrix file
MATRIX_OUT="${OUTDIR}/all_samples_multibigwig_matrix.npz"
HEATMAP_OUT="${OUTDIR}/all_samples_correlation_heatmap.png"

echo "=== Running multiBigwigSummary ==="
multiBigwigSummary bins \
    -b "${BW_FILES[@]}" \
    -out "$MATRIX_OUT" \
    --binSize 1000

echo "=== Running plotCorrelation ==="
plotCorrelation \
    -in "$MATRIX_OUT" \
    -c pearson \
    -p heatmap \
    -o "$HEATMAP_OUT" \
    --plotTitle "Correlation of All Samples" \
    --whatToPlot heatmap \
    --removeOutliers

echo "=== DONE! ==="
echo "Matrix: $MATRIX_OUT"
echo "Heatmap: $HEATMAP_OUT"
