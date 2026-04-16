#!/bin/bash
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=8G
#SBATCH --output=pp_%A_%a.out
#SBATCH --error=pp_%A_%a.err
#SBATCH --array=0-29
#SBATCH --time=2-00:00:00

set -euo pipefail

module load BCFtools/1.17-GCC-12.2.0
module load SAMtools/1.17-GCC-12.2.0

INPUT_LIST=${1:-""}
INPUT_DIR=${2:-""}
REF=${3:-""}
DBSNP=${4:-""}
THREADS=${5:-32}

if [[ -z "$INPUT_LIST" || -z "$INPUT_DIR" || -z "$REF" || -z "$DBSNP" ]]; then
    echo "Usage: sbatch variant_call_phase.sh <bam_list.txt> <input_dir> <reference.fa> <dbsnp_tab.gz> [threads]"
    exit 1
fi

mapfile -t INPUT < "$INPUT_LIST"
BAM_PATH="${INPUT[$SLURM_ARRAY_TASK_ID]}"

sample=$(basename -s .sorted.bam "$BAM_PATH")
BAM="${INPUT_DIR}/${sample}.sorted.bam"
OUTPUT_DIR="${INPUT_DIR}/${sample}.Phase"
OUTPUT_PREFIX="${sample}"

mkdir -p "$OUTPUT_DIR"

echo ">>> Processing sample: $sample"

echo ">>> Variant calling with PEPPER-Margin-DeepVariant"
singularity exec --bind /usr/lib/locale/ pepper_deepvariant_r0.8.sif \
    run_pepper_margin_deepvariant call_variant \
    -b "$BAM" \
    -f "$REF" \
    -o "$OUTPUT_DIR" \
    -p "$OUTPUT_PREFIX" \
    -t "$THREADS" \
    --pepper_report_snp_above_freq 0.1 \
    --pepper_report_indel_above_freq 0.1 \
    --ont_r10_q20

echo ">>> Annotating variants with dbSNP"
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    "${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz" | \
bcftools annotate \
    -a "$DBSNP" \
    -c CHROM,POS,ID,REF,ALT \
    - \
    -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}.dbSnp153.vcf.gz"

echo ">>> Filtering SNPs"
bcftools view \
    -i 'TYPE="snp" && QUAL>=10' \
    "${OUTPUT_DIR}/${OUTPUT_PREFIX}.dbSnp153.vcf.gz" \
    -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}.refineSNP.vcf.gz"

echo ">>> Phasing with LongPhase"
~/longphase phase \
    --ont \
    -s "${OUTPUT_DIR}/${OUTPUT_PREFIX}.refineSNP.vcf.gz" \
    -b "$BAM" \
    -q 10 \
    -r "$REF" \
    -t 31 \
    -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Hap"

echo ">>> Haplotagging reads"
~/longphase haplotag \
    -s "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Hap.vcf" \
    -b "$BAM" \
    -r "$REF" \
    -t 32 \
    -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Haptag" \
    2> "${OUTPUT_DIR}/${OUTPUT_PREFIX}.proc.txt"

echo ">>> Indexing haplotagged BAM"
samtools index -@ 31 "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Haptag.bam"

echo ">>> Extracting phased variants"
bcftools view \
    -i 'GT ~ "|"' \
    "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Hap.vcf" \
    > "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Hap.phaseV.vcf"

bcftools view -H \
    "${OUTPUT_DIR}/${OUTPUT_PREFIX}.Hap.phaseV.vcf" \
    > "${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.tab"

echo ">>> Done: $sample"
