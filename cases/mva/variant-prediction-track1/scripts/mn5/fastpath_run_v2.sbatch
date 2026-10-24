#!/bin/bash
#SBATCH --job-name=FASTPATH_v2
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log
set -euo pipefail
umask 077
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
WORK=$BASE/out/fastpath_20261016
SRC=${MVA_SOURCE_VCF:?}
REF=$BASE/out/l1_20261015/reference_47349911/GRCh38_no_alt.fa
BCF=$BASE/containers/bcftools_1.24--h487d631_1.sif
EV=$BASE/jobs/fastpath_20261016/echtvar
DB=$BASE/echtvar/gnomad_4.1_whole_genome.zip
CODE=$BASE/jobs/fastpath_20261016/fastpath_stream_v2.py
MODE=${1:?pilot or full}
[[ $MODE == pilot || $MODE == full ]] || exit 64
OUT=$WORK/${MODE}_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec 2>"$OUT/stderr.log"
module load singularity/4.1.5
test -s "$SRC" && test -s "$REF" && test -s "$REF.fai" && test -s "$DB" && test -x "$EV"
sha256sum "$SRC" "$REF" "$DB" "$CODE" "$EV" > "$OUT/input_sha256.txt"
singularity exec --cleanenv "$BCF" bcftools --version > "$OUT/bcftools.version"
"$EV" --version > "$OUT/echtvar.version"
LIMIT=0
if [[ $MODE == pilot ]]; then LIMIT=3000; fi
python3 "$CODE" prepare --input "$SRC" --limit "$LIMIT" --reference-fai "$REF.fai" --unsupported "$OUT/unsupported_contigs.json" |
    singularity exec --cleanenv "$BCF" bcftools norm -c e -f "$REF" -m -any -Oz -o "$OUT/normalized.vcf.gz" -
singularity exec --cleanenv "$BCF" tabix -p vcf "$OUT/normalized.vcf.gz"
"$EV" anno -e "$DB" "$OUT/normalized.vcf.gz" "$OUT/with_gnomad.vcf.gz"
singularity exec --cleanenv "$BCF" bcftools view -h "$OUT/with_gnomad.vcf.gz" | grep -E '##INFO=<ID=(FP_|gnomad_)' > "$OUT/annotation_headers.txt"
python3 "$CODE" summarize --input "$OUT/with_gnomad.vcf.gz" --output "$OUT/blind"
sha256sum "$OUT/blind"/* "$OUT/unsupported_contigs.json" > "$OUT/output_sha256.txt"
printf '%s\n' "$OUT/blind/manifest.json" > "$OUT/success"
