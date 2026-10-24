#!/bin/bash
#SBATCH --job-name=W3_spliceai
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=00:30:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/spliceai_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
test -d "$JOBS/wheels"
python3 "$JOBS/spliceai_keys.py" "$JOBS/spliceai_keys.tsv" "$REF" "$OUT"
test -s "$OUT/spliceai_keys.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
