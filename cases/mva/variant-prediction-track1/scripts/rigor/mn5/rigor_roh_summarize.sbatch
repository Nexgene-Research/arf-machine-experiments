#!/bin/bash
#SBATCH --job-name=RIGOR_rohS
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=acc
#SBATCH --time=00:20:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
SRC=$OUTROOT/roh_47437286/roh.rg.tsv
test -s "$SRC"
OUT=$OUTROOT/roh_sum_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
python3 "$JOBS/rigor_roh.py" "$SRC" "$REF.fai" "$OUT/roh_summary.tsv"
test -s "$OUT/roh_summary.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
