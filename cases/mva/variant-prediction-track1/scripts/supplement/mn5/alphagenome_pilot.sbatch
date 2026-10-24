#!/bin/bash
#SBATCH --job-name=SUP_ag
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/alphagenome_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
QUEUE=$JOBS/alphagenome_queue.tsv
python3 "$JOBS/alphagenome_infer.py" \
  "$OUT/indel_inference.tsv" \
  "$OUT/indel_inference.STATUS.json" \
  "$QUEUE"
printf '%s\n' "$OUT/indel_inference.tsv" > "$OUT/success"
