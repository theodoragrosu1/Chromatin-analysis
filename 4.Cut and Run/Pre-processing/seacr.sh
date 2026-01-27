#!/bin/bash
# ==========================================================
# Run SEACR for replicate 1 and replicate 2 bedGraph pairs
# with bedGraphs stored inside sample-specific subfolders
# ==========================================================

# --- PATHS (edit if needed)
BASE_DIR="/mnt/Data/Bond/tgrosu/cutrun/deduped/merged_bams/merged_bedgraphs"
OUTPUT_DIR="${BASE_DIR}/seacr_stringent_merged"


# --- Define conditions
CELL_LINES=(WT A7 C9)
MARKS=(me3 ub ac)

# --- Function to locate a file anywhere within the directory tree
find_bg () {
    local pattern="$1"
    find "$BASE_DIR" -type f -name "$pattern" | head -n 1
}

# --- Loop through all combinations
for CELL in "${CELL_LINES[@]}"; do
    for MARK in "${MARKS[@]}"; do

            SIGNAL_PATTERN="${CELL}_${MARK}_deduped.sorted.bedgraph"
            CONTROL_PATTERN="${CELL}_IgG_deduped.sorted.bedgraph"

            SIGNAL_BG=$(find_bg "$SIGNAL_PATTERN")
            CONTROL_BG=$(find_bg "$CONTROL_PATTERN")

            if [[ -z "$SIGNAL_BG" || -z "$CONTROL_BG" ]]; then
                echo "⚠️ Skipping ${CELL}_${MARK}: missing file(s)"
                continue
            fi

            OUT_PREFIX="${OUTPUT_DIR}/${CELL}_${MARK}"
            echo "=== Running SEACR for ${CELL}_${MARK} ==="
            echo "   signal: $SIGNAL_BG"
            echo "   control: $CONTROL_BG"

            /mnt/Data/Bond/tgrosu/cutrun/bedgraphs/SEACR_1.3.sh "$SIGNAL_BG" "$CONTROL_BG" norm stringent "$OUT_PREFIX"

            echo "✅ Done: ${OUT_PREFIX}"
        done
    done
done

echo "=== All SEACR replicate runs finished successfully ==="
