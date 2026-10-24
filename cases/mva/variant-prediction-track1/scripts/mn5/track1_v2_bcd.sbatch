#!/bin/bash
#SBATCH --job-name=track1_v2_bcd
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:50:00
set -euo pipefail
umask 077
module load intel/2023.2.0
module load samtools/1.23.1
module load singularity/4.1.5
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
B=$MVA_BASE
IN=$B/jobs/track1_v2_bcd_20261016
OUT=$B/out/track1_v2_bcd_20261016_${SLURM_JOB_ID}
mkdir -m 700 "$OUT"
python3 "$IN/track1_v2_bcd.py" "$IN" "$OUT"
printf '%s\n' "$OUT" > "$IN/result_location_${SLURM_JOB_ID}.txt"
printf 'JOB_COMPLETE\n' > "$OUT/completion.status"
