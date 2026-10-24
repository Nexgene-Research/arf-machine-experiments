#!/bin/bash
#SBATCH --job-name=SUP_exondp
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=08:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l3_depth_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
IMG=$BASE/containers/mosdepth_0.3.14--h87be163_2.sif
BED=$JOBS/tier1_sixgene_exon_utr.bed
test -s "$IMG" && test -s "$BED"
singularity exec --cleanenv "$IMG" mosdepth -t 8 --mapq 20 --no-per-base --by "$BED" --thresholds 1,10,20,30 "$OUT/exon" "$BAM"
python3 "$JOBS/exon_depth_table.py" "$OUT/exon.regions.bed.gz" "$OUT/exon_depth.tsv"
printf '%s\n' "$OUT/exon_depth.tsv" > "$OUT/success"
