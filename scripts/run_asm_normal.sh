#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --output=ASM_%A_%a.out
#SBATCH --error=ASM_%A_%a.err
#SBATCH --chdir=/share/lab_wangl/Yijun_Tian/RL_ASM/NORMALTMP/
#SBATCH --array=0-9999
#SBATCH --time=12:00:00

set -euo pipefail

SCRIPT="/share/lab_wangl/Yijun_Tian/RL_ASM/ASM.sh"
VCF_CHUNK_LIST="/share/lab_wangl/Yijun_Tian/RL_ASM/Normal.VCFchunklist.txt"
BAMDIR="/share/lab_wangl/Yijun_Tian/RL_ASM/NormalSlice"
LOG="/share/lab_wangl/Yijun_Tian/RL_ASM/ASM.SegmentationFault.Normal.log"

MATOUTPUT="/share/lab_wangl/Yijun_Tian/RL_ASM/NormalMAT"
DMROUTPUT="/share/lab_wangl/Yijun_Tian/RL_ASM/NormalDMR"

mkdir -p "$MATOUTPUT" "$DMROUTPUT"
touch "$LOG"

ml load Python/3.8.6-GCCcore-10.2.0
ml load BEDTools/2.30.0-GCC-11.2.0
ml load SAMtools/1.17-GCC-12.2.0
ml load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0

mapfile -t INLIST < "$VCF_CHUNK_LIST"

func() {
  local chunk_file="$1"
  local id
  local file

  id=$(basename -s .txt "$chunk_file" | cut -d"_" -f 2)
  file="${DMROUTPUT}/${id}.output.txt"

  echo "$file"

  if [[ -s "$file" ]]; then
    echo "this chunk is finished"
  else
    srun "$SCRIPT" "$chunk_file" "$BAMDIR" "$file" "$MATOUTPUT" "$LOG" 0.9
  fi
}

func "${INLIST[$SLURM_ARRAY_TASK_ID]}"
