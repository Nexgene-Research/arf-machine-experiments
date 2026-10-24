#!/bin/bash
#SBATCH --job-name=SUP_L1cons
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --partition=acc
#SBATCH --array=1-22%4
set -euo pipefail
umask 027
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
N=${SLURM_ARRAY_TASK_ID:?}
CHR=chr${N}
TAG=sup20261017
RUN=$L1
OUT=$OUTROOT/l1_consensus_${TAG}/$CHR
mkdir -p "$OUTROOT/l1_consensus_${TAG}"
if test -s "$OUT/success"; then
  echo "already complete $CHR"
  exit 0
fi
rm -rf "$OUT"
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
DV=$RUN/call_deepvariant_47354612/$CHR
GATK=$RUN/call_gatk_47354613/$CHR
STRELKA=$RUN/call_strelka_47354614/$CHR
for dir in "$DV" "$GATK" "$STRELKA"; do
  test -s "$dir/normalized.vcf.gz"
done
MERGE=$JOBS/l1_consensus_merge.py
test -s "$MERGE"
python3 "$MERGE" \
  --chrom "$CHR" --sample WGS_EX2312012 \
  --deepvariant "$DV/normalized.vcf.gz" \
  --gatk "$GATK/normalized.vcf.gz" \
  --strelka "$STRELKA/normalized.vcf.gz" \
  --consensus "$OUT/consensus.vcf" \
  --discordance "$OUT/caller_discordance.tsv" \
  --summary "$OUT/summary.json"
printf '%s\n' "$OUT/caller_discordance.tsv" > "$OUT/success"
