#!/bin/bash
#SBATCH --job-name=track1_reconstruct
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
set -euo pipefail
umask 077
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
B=$MVA_BASE
S=$B/out/track1_47397018
O=$B/out/track1_reconstruction_20261016
mkdir -p "$O"
python3 "$B/jobs/track1/37_select_readbacked_reconstruction.py" \
  --ranked-full "$S/ranked_pairs_real.tsv" \
  --read-stage "$S/read_phase_top60.tsv" \
  --hypotheses "$S/hypotheses.tsv" \
  --scored "$S/candidates_scored.tsv" \
  --top60 "$S/top60.tsv" \
  --out "$O/top10_reconstructed.tsv" \
  --counts "$O/primary_counts.json"
sha256sum "$O/top10_reconstructed.tsv" "$O/primary_counts.json"
