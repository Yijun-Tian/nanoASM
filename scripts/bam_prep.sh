#!/bin/bash

set -euo pipefail

# ============================
# BAM Preparation Script
# ============================

# Usage:
#   bash bam_prep.sh <bam_list> <output_prefix>

BAM_LIST=${1:-""}
PREFIX=${2:-"output"}
THREADS=${3:-12}

if [[ -z "$BAM_LIST" ]]; then
  echo "Usage: bash bam_prep.sh <bam_list> <output_prefix> [threads]"
  exit 1
fi

OUT_BAM="${PREFIX}.called.bam"
OUT_DIR="${PREFIX}Slice"

echo ">>> Merging BAM files into ${OUT_BAM}"

samtools merge -@ $THREADS -r -c -p --no-PG \
  -o $OUT_BAM \
  $(cat $BAM_LIST)

echo ">>> Indexing merged BAM"
samtools index -@ $THREADS $OUT_BAM

echo ">>> Splitting BAM by chromosome"
mkdir -p $OUT_DIR

for chr in $(samtools idxstats $OUT_BAM | cut -f1 | grep -v "^*"); do
    echo "Processing $chr"
    samtools view -b $OUT_BAM $chr > ${OUT_DIR}/${chr}.bam &
done

wait

echo ">>> Indexing chromosome BAMs"
for bam in ${OUT_DIR}/*.bam; do
    samtools index -@ 2 $bam &
done

wait

echo ">>> Done."
