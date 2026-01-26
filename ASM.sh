#!/bin/bash
set -euo pipefail
shopt -s nullglob

# before run
# ml load BEDTools/2.30.0-GCC-11.2.0
# ml load Python/3.8.6-GCCcore-10.2.0
# ml load SAMtools/1.17-GCC-12.2.0

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <input.vcf> <input.bamdir> <output.dmr.txt> <matrices.dir/> <segmentationfault.log>"
  exit 1
fi

input_vcf="$1"
input_bamdir="$2"
output="$3"
matrices="$4"
log="$5"

GenomeRef=/share/lab_wangl/Yijun_Tian/WGS_Prostate/hg38.alignmt.fa
metilene_dir=/share/lab_wangl/Yijun_Tian/WGS_Prostate/metilene_v0.2-8
label_script=/share/lab_wangl/Yijun_Tian/WGS_Prostate/labelAlign/label.py

mkdir -p "$matrices"
: > "$output"   # overwrite; change to touch if you want append-only

start_time=$(date +%s)

while IFS= read -r line; do
  [[ "$line" =~ ^# ]] && continue

  IFS=$'\t' read -r chr pos id ref alt dp qual gt00 gt01 gt11 <<< "$line"

  adjusted_pos=$((pos - 1))

  IFS=',' read -r -a alt_array <<< "$alt"
  alt0="${alt_array[0]}"

  variant_string="$id:$ref:${alt_array[*]}:$chr:$adjusted_pos"

  bam_chr="${input_bamdir}/${chr}.bam"
  out_bam="${id}.bam"
  out_bg_dir="${id}.bg"

  python "$label_script" -i "$bam_chr" -o "$out_bam" -v "$variant_string"
  samtools index "$out_bam"

  # Skip if no reads
  if [ "$(samtools view -c "$out_bam")" -eq 0 ]; then
    rm -f "$out_bam" "$out_bam.bai"
    continue
  fi

  region=$(samtools view "$out_bam" | awk '
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

  ~/dist/modkit pileup -t 4 --region "$region" \
    --filter-threshold 0.99609375 \
    --partition-tag VA --partition-tag RG \
    --ref "$GenomeRef" --cpg "$out_bam" --combine-strands "$out_bg_dir/"

  pushd "$out_bg_dir" >/dev/null

  # convert called.bed -> *.bg
  for j in *.called.bed; do
    [ -e "$j" ] || continue
    prefix="${j%.called.bed}"

    awk 'BEGIN{FS=OFS="\t"} $4~/m/ && $5>=5 {print $1,$2,$3,$11}' "$j" \
      | bedtools sort -i - > "${prefix}.bg"

    if [ ! -s "${prefix}.bg" ]; then
      rm -f "${prefix}.bg"
    fi
  done

  grp1_files=( "${id}@${ref}_"*.bg )
  grp2_files=( "${id}@${alt0}_"*.bg )
  n1=${#grp1_files[@]}
  n2=${#grp2_files[@]}

  if [ "$n1" -ge 3 ] && [ "$n2" -ge 3 ]; then
    grp1=$(IFS=,; echo "${grp1_files[*]}")
    grp2=$(IFS=,; echo "${grp2_files[*]}")

    perl "${metilene_dir}/metilene_input.pl" --in1 "$grp1" --in2 "$grp2" \
      --out "${id}.ref.alt.5mCG.input" --h1 ref --h2 alt
    
    # Have to pay attention to the metilene DMR calc process since it estimate missing methylation profiling based on existing values. The only way to adjust this was to increase the -X and -Y option value. Belows are a few exercise to eliminate the segment fault error but have gone in vein.
    # awk 'NR>1{count =0; cont = $0; while (match($0, /NA/)) {count++; $0=substr($0, RSTART+RLENGTH)}} {print count;next;}' <(cut -f 3-$((n1+2)) $id.ref.alt.5mCG.raw.input) | awk '!/^$/' |awk -v N1=$n1 'BEGIN{print "m1"} {print N1-$1}' > m1.txt
    # awk 'NR>1{count =0; cont = $0; while (match($0, /NA/)) {count++; $0=substr($0, RSTART+RLENGTH)}} {print count;next;}' <(cut -f $((n1+3))- $id.ref.alt.5mCG.raw.input) | awk '!/^$/' |awk -v N2=$n2 'BEGIN{print "m2"} {print N2-$1}' > m2.txt
    # paste m1.txt m2.txt $id.ref.alt.5mCG.raw.input | wk '{FS=OFS="\t"}NR==1{print $0;}NR>1{sum = 0; for (i = 5; i <= NF; i++) sum += $i; if (sum > 0) print $0;}' | awk '$1 >= 3 && $2 >=3 {print}' | cut -f 3-  > $id.ref.alt.5mCG.filter.input
    # add -X3 -Y3 to eliminate segmentation fault error (2026-01)
    # add -X3 -Y3 is a bad idea, since for groups with many allelic readout (say 30), if the program sees 3 value, it will try to guess those missing, therefore contirbuting to many false positive hits.
	
	# -M metilene works in two steps, first it pre-segments the whole data into windows so that no large gaps without data are possible. The -M parameter sets this length in nts. The default value of 300 means, that the whole genome is cut whenever a stretch of 300nts or more without data (CpGs) is found. E.g., if you the user does not want to find DMRs with stretches without CpGs longer than 200nt, the option "-M 200" should be used.

    # -m The length parameter -m sets the minimum value of CpGs/data points a DMR need to contain to be reported. As we use a top-down approach, starting with long windows and segmenting them to short significant DMRs, this is also a stop-criteria. Windows that contain a smaller number of CpGs are not considered and skipped.
	
	# -d The option -d sets the minimum mean methylation difference between both groups for a window to be reported as a DMR. This prevents to call regions with very small but significant significant methylation differences. We think that most users do not want to call smaller mean differences than 0.1, as the difference would be too small to term those regions as differentially methylated.

   set +e
   "${metilene_dir}/metilene_linux64" -t 2 -a ref -b alt -d 10 -M 300 -m 10 "${id}.ref.alt.5mCG.input" > "${id}.ref.alt.5mCG.output"
   rc=$?
   set -e

   if [ $rc -eq 139 ]; then
       {
         echo "===== $id  $(date) ====="
         echo "Segmentation fault"
         echo
        } >> "$log"

        # IMPORTANT: you're inside $id.bg
        popd >/dev/null

        # cleanup
        rm -rf "$out_bg_dir"
        rm -f "$out_bam" "$out_bam.bai"

        continue
   fi

    tmp="$(awk -v va="$variant_string" -v N1="$n1" -v N2="$n2" '
      BEGIN{FS=OFS="\t"}
      ($7<=0.05 || $8<=0.05){print $0,va,N1,N2}
    ' "${id}.ref.alt.5mCG.output")"

    if [ -n "$tmp" ]; then
      echo "$tmp" >> "$output"

      sig="$(awk -v va="$variant_string" -v N1="$n1" -v N2="$n2" '
        BEGIN{FS=OFS="\t"}
        ($4<=1e-3 && ($5*$5)>=400){print $0,va,N1,N2}
      ' "${id}.ref.alt.5mCG.output")"

      popd >/dev/null

      if [ -n "$sig" ]; then
        sed '1d' "${out_bg_dir}/${id}.ref.alt.5mCG.input" | \
          cat <(echo "chrom,pos,$grp1,$grp2" | sed -E 's/,/\t/g;s/.bg//g;s/rs[0-9]+@//g;s/Canc.//g;s/Norm.//g;s/.CG//g') - \
          > "${matrices}/${id}.matrix.txt"
          rm -rf "$out_bg_dir/"
          rm -f "$out_bam" "$out_bam.bai"
      else
        rm -rf "$out_bg_dir/"
        rm -f "$out_bam" "$out_bam.bai"
      fi

    else
      # no significant hits -> cleanup to avoid thousands of bg dirs
      popd >/dev/null
      rm -rf "$out_bg_dir/"
      rm -f "$out_bam" "$out_bam.bai"
    fi
  else
    popd >/dev/null
    rm -rf "$out_bg_dir/"
    rm -f "$out_bam" "$out_bam.bai"
    echo "this is a mission impossible"
  fi

done < "$input_vcf"

echo "It took $(( $(date +%s) - start_time )) seconds to finish this chunk"
echo "Finished processing VCF"
