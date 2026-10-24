#!/bin/bash
#SBATCH --job-name=W3_cnv
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=06:00:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/cnv_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
BED=$JOBS/gene_body_2kb.bed
test -s "$MOS" && test -s "$BED"
singularity exec --cleanenv "$MOS" mosdepth -t 8 --mapq 20 --no-per-base --by "$BED" "$OUT/genebody" "$BAM"
python3 "$JOBS/sixgene_cnv.py" "$OUT/genebody.regions.bed.gz" "$OUT/sixgene_cnv.tsv"
test -s "$OUT/sixgene_cnv.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
