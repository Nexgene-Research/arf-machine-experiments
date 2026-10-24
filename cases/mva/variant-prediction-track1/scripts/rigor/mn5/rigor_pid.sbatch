#!/bin/bash
#SBATCH --job-name=RIGOR_pid
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=08:00:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
OUT=$OUTROOT/pid_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
SCAN=$(cat "$OUTROOT/scan.success")
test -s "$SCAN/keys_impact.tsv"
test -d "$GATK_ROOT"
python3 "$JOBS/rigor_pid.py" "$SCAN/keys_impact.tsv" "$GATK_ROOT" "$OUT"
test -s "$OUT/pid_pair_counts.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
