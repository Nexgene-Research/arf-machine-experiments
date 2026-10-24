#!/bin/bash
#SBATCH --job-name=RIGOR_roh
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=12:00:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
OUT=$OUTROOT/roh_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
echo "${VCF_SHA}  ${VCF}" | sha256sum -c -
if ! test -s "$VCF.csi" && ! test -s "$VCF.tbi"; then
  bcf index -c "$VCF"
fi
bcf --version > "$OUT/bcftools.version"
# RG lines only. Per-site ST output is not kept.
bcf roh --AF-tag gnomad_AF --AF-dflt 0.4 -G 30 "$VCF" | awk 'BEGIN{OFS="\t"} $1=="RG" {print}' > "$OUT/roh.rg.tsv"
python3 "$JOBS/rigor_roh.py" "$OUT/roh.rg.tsv" "$REF.fai" "$OUT/roh_summary.tsv"
test -s "$OUT/roh_summary.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
