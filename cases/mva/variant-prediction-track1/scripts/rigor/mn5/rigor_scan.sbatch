#!/bin/bash
#SBATCH --job-name=RIGOR_scan
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=08:00:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
OUT=$OUTROOT/scan_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
echo "${VCF_SHA}  ${VCF}" | sha256sum -c -
python3 --version > "$OUT/python.version"
python3 "$JOBS/rigor_scan.py" "$VCF" "$JOBS/ACMG_SF_v3.2.genes.txt" "$OUT"
test -s "$OUT/vaf_counts.tsv"
test -s "$OUT/keys_impact.tsv"
printf '%s\n' "$OUT" > "$OUTROOT/scan.success"
printf '%s\n' "$OUT" > "$OUT/success"
