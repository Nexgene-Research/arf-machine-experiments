#!/bin/bash
#SBATCH --job-name=track1_scores
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:40:00
set -euo pipefail
umask 077
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
B=$MVA_BASE
S=$B/out/critical_annotation_47397018
O=$B/out/track1_47397018
mkdir -p "$O"
exec > >(tee "$O/scores.log") 2>&1
python3 "$B/jobs/track1/track1_scores.py" "$S" "$O" \
  "$B/annot/AlphaMissense_hg38.tsv.gz" \
  "$B/annot/revel-v1.3_all_chromosomes.zip"
python3 "$B/jobs/track1/track1_pairs.py" "$O" "$O"
echo SCORES_PAIRS_COMPLETE > "$O/scores.status"
