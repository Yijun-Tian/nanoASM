#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=24G
#SBATCH --output=dorado_%A_%a.out
#SBATCH --error=dorado_%A_%a.err
#SBATCH --gpus-per-task=a30:1
#SBATCH --time=12:00:00

set -euo pipefail

# Usage:
#   sbatch --array=0-N scripts/dorado_basecall_cg.sh <pod5_list.txt> <output_dir>

POD5_LIST=${1:-""}
OUTPUT_DIR=${2:-""}

if [[ -z "$POD5_LIST" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: sbatch --array=0-N dorado_basecall_cg.sh <pod5_list.txt> <output_dir>"
    exit 1
fi

REFERENCE="HG38ALIGNmt.mmi"
BASECALL_MODEL="/home/4470655/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
MODBASE_MODEL="/home/4470655/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3"

mapfile -t INPUT < "$POD5_LIST"
POD5_FILE="${INPUT[$SLURM_ARRAY_TASK_ID]}"

BASENAME=$(basename -s .pod5 "$POD5_FILE")
OUT_BAM="${OUTPUT_DIR}/${BASENAME}.bam"

mkdir -p "$OUTPUT_DIR"

if [[ -s "$OUT_BAM" ]]; then
    echo "Basecalling already finished: $OUT_BAM"
    exit 0
fi

echo "Processing $POD5_FILE"
srun ~/dorado-0.9.1-linux-x64/bin/dorado basecaller \
    --device cuda:0 \
    --reference "$REFERENCE" \
    "$BASECALL_MODEL" \
    "$POD5_FILE" \
    --kit-name SQK-NBD114-24 \
    --batchsize 5120 \
    --modified-bases-threshold 0 \
    --modified-bases-models "$MODBASE_MODEL" \
    > "$OUT_BAM"

echo "Done: $OUT_BAM"
