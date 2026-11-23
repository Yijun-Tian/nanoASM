#!/bin/bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# ASM core pipeline
# ============================================================
# Usage:
#   asm_core.sh <input.vcf> <input.bamdir> <output.dmr.txt> \
#               <matrices.dir/> <segmentationfault.log> [na_thr]
#
# Arguments:
#   input.vcf              Chunked VCF summary file
#   input.bamdir           Directory containing chromosome BAMs
#   output.dmr.txt         Output DMR result file
#   matrices.dir/          Directory for matrix output
#   segmentationfault.log  Log for metilene segfaults
#   na_thr                 Optional NA fraction threshold (default: 0.90)
# ============================================================

if [[ "$#" -lt 5 || "$#" -gt 6 ]]; then
  echo "Usage: $0 <input.vcf> <input.bamdir> <output.dmr.txt> <matrices.dir/> <segmentationfault.log> [na_thr]"
  echo "  na_thr: drop samples with estimated NA fraction >= na_thr (default 0.90)"
  exit 1
fi

INPUT_VCF="$1"
INPUT_BAMDIR="$2"
OUTPUT="$3"
MATRICES_DIR="$4"
LOG_FILE="$5"
NA_THR="${6:-0.90}"

# ------------------------------------------------------------
# Validate NA threshold
# ------------------------------------------------------------
awk -v v="$NA_THR" 'BEGIN{
  if (v !~ /^([0-9]*\.?[0-9]+)$/ || v < 0 || v > 1) {
    print "Error: na_thr must be a number between 0 and 1. Got:", v > "/dev/stderr"
    exit 1
  }
}'

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
# Resolve script directory so relative paths work regardless of where the script is invoked from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Configurable paths (can be overridden by environment variables or CLI flags)
GENOME_REF="${GENOME_REF:-${SCRIPT_DIR}/demo/hg38.alignmt.fa}"
METILENE_DIR="${METILENE_DIR:-$(dirname "$(command -v metilene)")}"
LABEL_SCRIPT="${LABEL_SCRIPT:-${SCRIPT_DIR}/label.py}"

MODKIT="${MODKIT:-$(command -v modkit || true)}"
if [[ -z "${MODKIT}" || ! -x "${MODKIT}" ]]; then
  echo "Error: modkit not found on PATH and MODKIT is not set." >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "Error: samtools not found on PATH." >&2
  exit 1
fi

MIN_KEEP=3

mkdir -p "$MATRICES_DIR"
: > "$OUTPUT"

START_TIME=$(date +%s)

while IFS= read -r line; do
  [[ "$line" =~ ^# ]] && continue
  [[ -z "$line" ]] && continue

  IFS=$'\t' read -r CHR POS ID REF ALT DP QUAL GT00 GT01 GT11 <<< "$line" || continue

  ADJUSTED_POS=$((POS - 1))
  IFS=',' read -r -a ALT_ARRAY <<< "$ALT"
  ALT0="${ALT_ARRAY[0]}"

  VARIANT_STRING="$ID:$REF:${ALT_ARRAY[*]}:$CHR:$ADJUSTED_POS"

  BAM_CHR="${INPUT_BAMDIR}/${CHR}.bam"
  OUT_BAM="${ID}.bam"
  OUT_BG_DIR="${ID}.bg"

  python "$LABEL_SCRIPT" -i "$BAM_CHR" -o "$OUT_BAM" -v "$VARIANT_STRING"
  samtools index "$OUT_BAM"

  if [[ "$(samtools view -c "$OUT_BAM")" -eq 0 ]]; then
    rm -f "$OUT_BAM" "$OUT_BAM.bai"
    continue
  fi

  REGION=$(samtools view "$OUT_BAM" | awk '
    {
      c=$6; s=$4; r=0;
      while (match(c, /([0-9]+)([MIDNSHP=X])/)) {
        len=substr(c, RSTART, RLENGTH-1);
        type=substr(c, RSTART+RLENGTH-1, 1);
        if (type ~ /[MDN=X]/) r+=len;
        c=substr(c, RSTART+RLENGTH);
      }
      e=s+r-1;
      if (NR==1) {left=s; right=e; chr=$3;}
      if (s<left) left=s;
      if (e>right) right=e;
    }
    END {print chr ":" left "-" right}
  ')

  "${MODKIT}" pileup \
    -t 4 \
    --region "$REGION" \
    --filter-threshold 0.99609375 \
    --partition-tag VA \
    --partition-tag RG \
    --ref "$GENOME_REF" \
    --cpg "$OUT_BAM" \
    --combine-strands "$OUT_BG_DIR/"

  pushd "$OUT_BG_DIR" >/dev/null

  # ------------------------------------------------------------
  # Step 1: build .bg files and count unique sites per sample
  # ------------------------------------------------------------
  SITE_COUNTS_TSV="${ID}.bg.site_counts.tsv"
  : > "$SITE_COUNTS_TSV"

  for j in *.called.bed; do
    [[ -e "$j" ]] || continue
    PREFIX="${j%.called.bed}"

    awk 'BEGIN{FS=OFS="\t"} $4~/m/ && $5>=5 {print $1,$2,$3,$11}' "$j" \
      | bedtools sort -i - \
      | awk 'BEGIN{FS=OFS="\t"} !seen[$1 FS $2 FS $3]++' \
      > "${PREFIX}.bg"

    if [[ ! -s "${PREFIX}.bg" ]]; then
      rm -f "${PREFIX}.bg"
      echo -e "${PREFIX}\t0" >> "$SITE_COUNTS_TSV"
      continue
    fi

    N_SITES=$(wc -l < "${PREFIX}.bg")
    echo -e "${PREFIX}\t${N_SITES}" >> "$SITE_COUNTS_TSV"
  done

  BG_FILES=( *.bg )
  if [[ "${#BG_FILES[@]}" -eq 0 ]]; then
    popd >/dev/null
    rm -rf "$OUT_BG_DIR/"
    rm -f "$OUT_BAM" "$OUT_BAM.bai"
    continue
  fi

  # ------------------------------------------------------------
  # Step 2: build union CpG set across all .bg files
  # ------------------------------------------------------------
  UNION_SITES_BED="${ID}.bg.union_sites.bed"

  cat ./*.bg \
    | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' \
    | bedtools sort -i - \
    | bedtools merge -i - \
    > "$UNION_SITES_BED"

  TOTAL_SITES=$(wc -l < "$UNION_SITES_BED")

  if [[ "$TOTAL_SITES" -eq 0 ]]; then
    popd >/dev/null
    rm -rf "$OUT_BG_DIR/"
    rm -f "$OUT_BAM" "$OUT_BAM.bai"
    continue
  fi

  # ------------------------------------------------------------
  # Step 3: estimate sparsity and filter samples before merging
  # ------------------------------------------------------------
  KEEP_LIST="${ID}.bg.keep_samples.txt"
  DROP_LIST="${ID}.bg.drop_samples.txt"
  SPARSITY_TSV="${ID}.bg.sparsity.tsv"

  : > "$KEEP_LIST"
  : > "$DROP_LIST"

  awk -v tot="$TOTAL_SITES" -v thr="$NA_THR" 'BEGIN{FS=OFS="\t"}
  {
    sample=$1
    n=$2+0
    present_frac = (tot>0 ? n/tot : 0)
    est_na = 1 - present_frac

    print sample, n, tot, present_frac, est_na

    if (est_na < thr) print sample > "'"$KEEP_LIST"'"
    else              print sample > "'"$DROP_LIST"'"
  }' "$SITE_COUNTS_TSV" > "$SPARSITY_TSV"

  # ------------------------------------------------------------
  # Step 4: update ref/alt groups after filtering
  # ------------------------------------------------------------
  GRP1_ALL=( "${ID}@${REF}_"*.bg )
  GRP2_ALL=( "${ID}@${ALT0}_"*.bg )

  GRP1_FILES=()
  GRP2_FILES=()

  for f in "${GRP1_ALL[@]}"; do
    [[ -e "$f" ]] || continue
    SAMPLE="${f%.bg}"
    if grep -Fxq "$SAMPLE" "$KEEP_LIST"; then
      GRP1_FILES+=( "$f" )
    fi
  done

  for f in "${GRP2_ALL[@]}"; do
    [[ -e "$f" ]] || continue
    SAMPLE="${f%.bg}"
    if grep -Fxq "$SAMPLE" "$KEEP_LIST"; then
      GRP2_FILES+=( "$f" )
    fi
  done

  N1=${#GRP1_FILES[@]}
  N2=${#GRP2_FILES[@]}

  if [[ "$N1" -lt "$MIN_KEEP" || "$N2" -lt "$MIN_KEEP" ]]; then
    echo "[$ID] After sparsity filter (thr=$NA_THR): n1=$N1 n2=$N2 total_sites=$TOTAL_SITES. Skipping."
    popd >/dev/null
    rm -rf "$OUT_BG_DIR/"
    rm -f "$OUT_BAM" "$OUT_BAM.bai"
    continue
  fi

  GRP1=$(IFS=,; echo "${GRP1_FILES[*]}")
  GRP2=$(IFS=,; echo "${GRP2_FILES[*]}")

  # ------------------------------------------------------------
  # Step 5: merge filtered samples only
  # ------------------------------------------------------------
  perl "${METILENE_DIR}/metilene_input.pl" \
    --in1 "$GRP1" \
    --in2 "$GRP2" \
    --out "${ID}.ref.alt.5mCG.input" \
    --h1 ref \
    --h2 alt

  set +e
  "${METILENE_DIR}/metilene_linux64" \
    -t 2 -a ref -b alt -d 10 -M 500 -m 10 \
    "${ID}.ref.alt.5mCG.input" \
    > "${ID}.ref.alt.5mCG.output"
  RC=$?
  set -e

  if [[ "$RC" -eq 139 ]]; then
    {
      echo "===== $ID  $(date) ====="
      echo "Segmentation fault"
      echo
    } >> "$LOG_FILE"

    popd >/dev/null
    rm -rf "$OUT_BG_DIR/"
    rm -f "$OUT_BAM" "$OUT_BAM.bai"
    continue
  fi

  TMP="$(awk -v va="$VARIANT_STRING" -v N1="$N1" -v N2="$N2" '
    BEGIN{FS=OFS="\t"}
    ($7<=0.05 || $8<=0.05){print $0,va,N1,N2}
  ' "${ID}.ref.alt.5mCG.output")"

  if [[ -n "$TMP" ]]; then
    echo "$TMP" >> "$OUTPUT"

    SIG="$(awk -v va="$VARIANT_STRING" -v N1="$N1" -v N2="$N2" '
      BEGIN{FS=OFS="\t"}
      ($4<=1e-2 && ($5*$5)>=400){print $0,va,N1,N2}
    ' "${ID}.ref.alt.5mCG.output")"

    if [[ -n "$SIG" ]]; then
      sed '1d' "${ID}.ref.alt.5mCG.input" \
        | cat <(
            echo "chrom,pos,$GRP1,$GRP2" \
            | sed -E 's/,/\t/g;s/.bg//g;s/rs[0-9]+@//g;s/Canc.//g;s/Norm.//g;s/.CG//g'
          ) - \
        > "${MATRICES_DIR}/${ID}.matrix.txt"
    fi
  fi

  popd >/dev/null
  rm -rf "$OUT_BG_DIR/"
  rm -f "$OUT_BAM" "$OUT_BAM.bai"

done < "$INPUT_VCF"

echo "It took $(( $(date +%s) - START_TIME )) seconds to finish this chunk"
echo "Finished processing VCF"
