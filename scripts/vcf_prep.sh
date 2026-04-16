#!/bin/bash

set -euo pipefail

# ============================
# VCF Preparation Script
# ============================

# Usage:
#   bash vcf_prep.sh <vcf_list> <dbsnp_file> <prefix>

VCF_LIST=${1:-""}
DBSNP=${2:-""}
PREFIX=${3:-"cohort"}
THREADS=${4:-8}

if [[ -z "$VCF_LIST" || -z "$DBSNP" ]]; then
  echo "Usage: bash vcf_prep.sh <vcf_list> <dbsnp_file> <prefix> [threads]"
  exit 1
fi

echo ">>> Step 1: Process individual VCFs"

while read vcf; do
    sample=$(basename -s .vcf.gz $vcf)
    echo "Processing $sample"

    # annotate ID + dbSNP
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' $vcf | \
    bcftools annotate -a $DBSNP -c CHROM,POS,ID,REF,ALT - \
    -o ${sample}.dbSnp.vcf.gz

    # filter SNPs
    bcftools filter --threads $THREADS -i 'TYPE="snp"' ${sample}.dbSnp.vcf.gz | \
    bcftools view --threads $THREADS -i 'QUAL>=10 && VAF>=0.1 && ID ~ "^rs"' \
    -o ${sample}.refineSNP.vcf.gz

    bcftools index ${sample}.refineSNP.vcf.gz

done < $VCF_LIST

# ============================
# Step 2: Merge cohort VCF
# ============================

echo ">>> Step 2: Merge VCFs"

ls *.refineSNP.vcf.gz > ${PREFIX}.list.txt

bcftools merge --threads 15 \
    --file-list ${PREFIX}.list.txt \
    --force-samples \
    --merge snps -0 | \
bcftools +fill-tags -- -t AC,AN | \
bcftools view --trim-alt-alleles \
    -o ${PREFIX}.merge.refineSNP.vcf.gz

# ============================
# Step 3: Filter variants
# ============================

echo ">>> Step 3: Filtering variants"

bcftools view --min-ac 3 -m2 -M2 ${PREFIX}.merge.refineSNP.vcf.gz | \
bcftools view -i 'COUNT(FMT/DP >= 10 & FMT/AD[:0] >= 5 & GT="het") >= 3 && QUAL >= 20' \
> ${PREFIX}.filtered.vcf

echo ">>> Done."
