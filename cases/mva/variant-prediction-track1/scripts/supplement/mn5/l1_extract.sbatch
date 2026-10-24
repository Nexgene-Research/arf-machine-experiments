#!/bin/bash
#SBATCH --job-name=SUP_L1win
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=acc
#SBATCH --time=01:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l1_windows_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
ROOT=$OUTROOT/l1_consensus_sup20261017
python3 "$JOBS/l1_extract_windows.py" "$ROOT" "$JOBS/tier1_sixgene_exon_utr.bed" "$OUT/caller_discordance.tsv"
test -s "$OUT/caller_discordance.tsv"
printf '%s\n' "$OUT/caller_discordance.tsv" > "$OUT/success"
