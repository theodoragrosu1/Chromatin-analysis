#!/usr/bin/env bash


# Define cell lines and PTMs
CELL_LINES=("WT" "A7" "C9")
PTMS=("me3" "ac" "ub")

GSIZE="2700000000"
OUTDIR="macs2_peaks"

mkdir -p "${OUTDIR}"

for CELL in "${CELL_LINES[@]}"; do
  for PTM in "${PTMS[@]}"; do

    TREATMENT="./${CELL}_${PTM}_deduped.bam"
    CONTROL="./${CELL}_IgG_deduped.bam"
    NAME="${CELL}_${PTM}"

    if [[ ! -f "${TREATMENT}" ]]; then
      echo "Skipping ${NAME}: treatment BAM not found"
      continue
    fi

    if [[ ! -f "${CONTROL}" ]]; then
      echo "Skipping ${NAME}: IgG control BAM not found"
      continue
    fi

    echo "Running MACS2 for ${NAME}"

    macs2 callpeak \
      -t "${TREATMENT}" \
      -c "${CONTROL}" \
      --format BAMPE \
      --gsize "${GSIZE}" \
      --nolambda \
      --to-large \
      --broad --broad-cutoff 0.1 \
      --keep-dup auto \
      --d-min 20 \
      --buffer-size 100000 \
      --qvalue 0.1 \
      --mfold 5 50 \
      --bw 300 \
      --name "${NAME}" \
      --outdir "${OUTDIR}"

  done
done