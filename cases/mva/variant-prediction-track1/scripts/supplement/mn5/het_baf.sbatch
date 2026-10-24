#!/bin/bash
#SBATCH --job-name=SUP_baf
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=04:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/het_baf_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
VCF=
for cand in \
  "$BASE/out/critical_annotation_47397018/with_clinvar.vcf.gz" \
  "$BASE/out/critical_annotation_47392756/with_clinvar.vcf.gz" \
  "$BASE/out/fastpath_20261016/full_47392756/with_gnomad.vcf.gz"
do
  if test -s "$cand"; then
    VCF=$cand
    break
  fi
done
test -n "$VCF" && test -s "$VCF"
printf 'vcf=%s\n' "$VCF" > "$OUT/inputs.txt"
python3 "$JOBS/het_baf.py" "$VCF" "$OUT/het_baf.tsv"
printf '%s\n' "$OUT/het_baf.tsv" > "$OUT/success"
