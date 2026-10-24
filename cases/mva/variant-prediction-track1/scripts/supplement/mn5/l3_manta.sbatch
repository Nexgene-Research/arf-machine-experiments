#!/bin/bash
#SBATCH --job-name=SUP_manta
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=acc
#SBATCH --time=12:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l3_manta_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
IMG=$BASE/containers/manta_1.6.0--h9ee0642_2.sif
BED=$JOBS/tier1_sixgene_exon_utr.bed
test -s "$IMG" && test -s "$BED"
printf 'bam=%s\nref=%s\n' "$BAM" "$REF" > "$OUT/inputs.txt"
singularity exec --cleanenv "$IMG" configManta.py --bam "$BAM" --referenceFasta "$REF" --runDir "$OUT/workflow"
singularity exec --cleanenv "$IMG" "$OUT/workflow/runWorkflow.py" -m local -j 16
VCF=$OUT/workflow/results/variants/diploidSV.vcf.gz
test -s "$VCF"
python3 "$JOBS/filter_sv_windows.py" "$VCF" "$BED" "$OUT/manta_tier1_sv.tsv" manta
printf '%s\n' "$OUT/manta_tier1_sv.tsv" > "$OUT/success"
