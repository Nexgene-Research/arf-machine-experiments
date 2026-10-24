#!/bin/bash
#SBATCH --job-name=W3_mosaic
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=12:00:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/mosaic_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
test -s "$GATK_SIF"
SITES=$JOBS/sites.hg38.vcf.gz
test -s "$SITES" && test -s "$SITES.tbi"
gatk CollectAllelicCounts \
  -L "$SITES" \
  -imr OVERLAPPING_ONLY \
  -I "$BAM" \
  -R "$REF" \
  -O "$OUT/allelicCounts.tsv"
gatk ModelSegments \
  --allelic-counts "$OUT/allelicCounts.tsv" \
  --output "$OUT" \
  --output-prefix mosaic
python3 "$JOBS/mosaic_table.py" "$OUT/mosaic" "$OUT/mosaic_segments.tsv"
test -s "$OUT/mosaic_segments.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
