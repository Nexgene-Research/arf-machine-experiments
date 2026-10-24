#!/bin/bash
#SBATCH --job-name=SUP_hg002
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=acc
#SBATCH --time=02:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/hg002_null_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
VCF=$BASE/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
if ! test -s "$VCF"; then
  VCF=$BASE/annot/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
fi
if ! test -s "$VCF"; then
  VCF=${MVA_BASE:?}/stage/T5/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
fi
test -s "$VCF"
python3 "$JOBS/hg002_null.py" "$VCF" "$JOBS/tier1_sixgene_exon_utr.bed" "$OUT/null_hg002.tsv"
printf '%s\n' "$OUT/null_hg002.tsv" > "$OUT/success"
