#!/bin/bash
#SBATCH --job-name=W3_pileup
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=04:00:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/pileup_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
python3 --version > "$OUT/python.version"
echo "${VCF_SHA}  ${VCF}" | sha256sum -c - > "$OUT/inputs.sha256" || true
python3 "$JOBS/pileup_published.py" "$JOBS/published_sites.vcf" "$BAM" "$REF" "$BCF" "$OUT"
test -s "$OUT/published_site_gt_new.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
