#!/bin/bash
#SBATCH --job-name=W3_vbid
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=06:00:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/vbid_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
test -s "$VBID"
PREFIX=$JOBS/1000g.phase3.10k.b38.vcf.gz.dat
test -s "${PREFIX}.UD" && test -s "${PREFIX}.bed" && test -s "${PREFIX}.mu"
sha256sum "${PREFIX}.UD" "${PREFIX}.bed" "${PREFIX}.mu" > "$OUT/inputs.sha256"
singularity exec --cleanenv "$VBID" verifybamid2 \
  --BamFile "$BAM" \
  --Reference "$REF" \
  --SVDPrefix "$PREFIX" \
  --Output "$OUT/verify" \
  --NumThread 8
test -s "$OUT/verify.selfSM"
printf '%s\n' "$OUT" > "$OUT/success"
