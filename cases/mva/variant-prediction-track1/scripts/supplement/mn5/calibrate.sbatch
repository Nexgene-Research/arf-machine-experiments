#!/bin/bash
#SBATCH --job-name=SUP_cal
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=08:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/calibration_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
CLIN=
AM=
REV=
for cand in "$BASE/annot/clinvar.vcf.gz" "$BASE/annot/clinvar_GRCh38.vcf.gz"; do
  test -s "$cand" && CLIN=$cand && break
done
for cand in "$BASE/annot/AlphaMissense_hg38.tsv.gz"; do
  test -s "$cand" && AM=$cand && break
done
for cand in "$BASE/annot/revel-v1.3_all_chromosomes.zip"; do
  test -s "$cand" && REV=$cand && break
done
test -n "$CLIN" && test -n "$AM" && test -n "$REV"
printf 'clinvar=%s\nam=%s\nrevel=%s\n' "$CLIN" "$AM" "$REV" > "$OUT/inputs.txt"
python3 "$JOBS/calibrate_predictors.py" "$CLIN" "$AM" "$REV" "$OUT/predictor_lr.tsv"
printf '%s\n' "$OUT/predictor_lr.tsv" > "$OUT/success"
