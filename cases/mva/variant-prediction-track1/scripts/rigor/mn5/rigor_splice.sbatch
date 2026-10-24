#!/bin/bash
#SBATCH --job-name=RIGOR_spl
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=acc
#SBATCH --time=04:00:00
source "${MVA_BASE:?}/jobs/track1_rigor/common.sh"
OUT=$OUTROOT/splice_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
echo "${VCF_SHA}  ${VCF}" | sha256sum -c -
test -s "$JOBS/tier1_sixgene_coding_splice20.bed"
test -s "$JOBS/alphagenome_splicing_sixgene_slice.tsv"
test -s "$JOBS/alphagenome_splicing_sixgene_background.tsv"
python3 "$JOBS/rigor_splice.py" \
  "$VCF" \
  "$JOBS/tier1_sixgene_coding_splice20.bed" \
  "$JOBS/alphagenome_splicing_sixgene_slice.tsv" \
  "$JOBS/alphagenome_splicing_sixgene_background.tsv" \
  "$OUT"
test -s "$OUT/splice_sixgene_counts.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
