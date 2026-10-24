#!/bin/bash
#SBATCH --job-name=RIGOR_cln
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=04:00:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
OUT=$OUTROOT/clinvar_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
echo "${VCF_SHA}  ${VCF}" | sha256sum -c -
test -s "$CLIN"
python3 "$JOBS/rigor_clinvar.py" "$VCF" "$CLIN" "$JOBS/ACMG_SF_v3.2.genes.txt" "$OUT"
test -s "$OUT/clinvar_plp_summary.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
