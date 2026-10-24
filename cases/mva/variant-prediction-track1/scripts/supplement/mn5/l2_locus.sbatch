#!/bin/bash
#SBATCH --job-name=SUP_l2loc
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l2_locus_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
printf 'gene\tchrom\tstart\tend\tstatus\tmean_depth\tnote\n' > "$OUT/l2_locus.tsv"
GIR=$BASE/out/l2_align_full_47356275
ready=1
for n in 1 2 3 4; do
  if ! test -s "$GIR/lane_${n}/success"; then
    ready=0
  fi
done
if test "$ready" -eq 0; then
  printf 'BUB1B\tchr15\t40153277\t40186682\tgiraffe_incomplete\t\t47356275 lanes not all successful; no second copy submitted\n' >> "$OUT/l2_locus.tsv"
  printf 'BUB1\tchr2\t110637528\t110678263\tgiraffe_incomplete\t\t47356275 lanes not all successful; no second copy submitted\n' >> "$OUT/l2_locus.tsv"
  printf '%s\n' "$OUT/l2_locus.tsv" > "$OUT/success"
  exit 0
fi
# Prefer a merged BAM if the compare step already wrote one; else per-lane mean.
MERGED=$GIR/GRCh38.sorted.bam
if test -s "$MERGED"; then
  for rec in "BUB1B chr15:40153277-40186682" "BUB1 chr2:110637528-110678263"; do
    set -- $rec
    gene=$1
    region=$2
    chrom=${region%%:*}
    rest=${region#*:}
    start=${rest%-*}
    end=${rest#*-}
    mean=$(sam coverage -r "$region" "$MERGED" | awk 'NR==2 {print $7}')
    printf '%s\t%s\t%s\t%s\tdone\t%s\tgiraffe merged BAM coverage; not a genotype\n' \
      "$gene" "$chrom" "$start" "$end" "${mean:-}" >> "$OUT/l2_locus.tsv"
  done
else
  printf 'BUB1B\tchr15\t40153277\t40186682\tlanes_unmerged\t\t47356275 succeeded per lane; no merged BAM yet\n' >> "$OUT/l2_locus.tsv"
  printf 'BUB1\tchr2\t110637528\t110678263\tlanes_unmerged\t\t47356275 succeeded per lane; no merged BAM yet\n' >> "$OUT/l2_locus.tsv"
fi
printf '%s\n' "$OUT/l2_locus.tsv" > "$OUT/success"
