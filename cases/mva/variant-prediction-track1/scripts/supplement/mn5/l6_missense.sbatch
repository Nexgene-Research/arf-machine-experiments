#!/bin/bash
#SBATCH --job-name=SUP_l6
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/l6_missense_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
# Coverage first. ThermoMPNN/ESM run only if a covering chain exists and the
# environment is present; otherwise the TSV records the gap.
PDB=$BASE/stage/T6/pdb
if ! test -d "$PDB"; then
  PDB=$BASE/models/pdb
fi
python3 "$JOBS/l6_missense.py" "$OUT/asn1002lys.tsv" "$PDB"
# Weights are staged; inference code is not. Do not invent ddG or ESM logits.
if grep -q $'\tyes\t' "$OUT/asn1002lys.tsv"; then
  echo 'covering residue present; ThermoMPNN/ESM code not staged; scores left blank'
fi
printf '%s\n' "$OUT/asn1002lys.tsv" > "$OUT/success"
