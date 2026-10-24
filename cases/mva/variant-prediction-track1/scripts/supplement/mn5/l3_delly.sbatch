#!/bin/bash
#SBATCH --job-name=SUP_delly
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=acc
#SBATCH --time=12:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l3_delly_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
IMG=$BASE/containers/delly_2.7.0--h3752d28_0.sif
BED=$JOBS/tier1_sixgene_exon_utr.bed
test -s "$IMG" && test -s "$BED"
printf 'bam=%s\nref=%s\n' "$BAM" "$REF" > "$OUT/inputs.txt"
singularity exec --cleanenv "$IMG" delly sr -g "$REF" -h 16 -o "$OUT/delly.bcf" "$BAM"
test -s "$OUT/delly.bcf"
bcf view -Oz -o "$OUT/delly.vcf.gz" "$OUT/delly.bcf"
python3 "$JOBS/filter_sv_windows.py" "$OUT/delly.vcf.gz" "$BED" "$OUT/delly_tier1_sv.tsv" delly
printf '%s\n' "$OUT/delly_tier1_sv.tsv" > "$OUT/success"
