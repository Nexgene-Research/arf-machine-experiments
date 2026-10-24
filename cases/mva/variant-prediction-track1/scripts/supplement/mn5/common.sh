#!/bin/bash
# Shared MN5 paths for the Track 1 evidence supplement. Source from sbatch.
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
L1=$BASE/out/l1_20261015
JOBS=$BASE/jobs/track1_supplement
OUTROOT=$BASE/out/track1_supplement
mkdir -p "$OUTROOT"
if test -s "$L1/bam.success"; then
  BAM=$(cat "$L1/bam.success")
elif test -s "$L1/merge_47353088/marked.bam"; then
  BAM=$L1/merge_47353088/marked.bam
else
  echo 'BLOCKED: no L1 marked BAM' >&2
  exit 2
fi
if test -s "$L1/reference.success"; then
  REF=$(cat "$L1/reference.success")
else
  REF=$L1/reference_47349911/GRCh38_no_alt.fa
fi
test -s "$BAM" && test -s "$REF" && test -s "$REF.fai"
if test -s "$BAM.bai"; then
  :
elif test -s "${BAM%.bam}.bai"; then
  :
else
  echo 'BLOCKED: BAM index missing' >&2
  exit 2
fi
module load singularity/4.1.5
module load intel/2023.2.0
module load samtools/1.23.1
SAM=$BASE/containers/samtools_1.24--h9dcdb79_1.sif
BCF=$BASE/containers/bcftools_1.24--h487d631_1.sif
test -s "$SAM" && test -s "$BCF"
sam() { singularity exec --cleanenv "$SAM" samtools "$@"; }
bcf() { singularity exec --cleanenv "$BCF" bcftools "$@"; }
