#!/bin/bash
# Shared MN5 paths for Track 1 evidence wave 3. Source from sbatch.
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
L1=$BASE/out/l1_20261015
JOBS=$BASE/jobs/track1_wave3
OUTROOT=$BASE/out/track1_wave3
mkdir -p "$OUTROOT"
if test -s "$L1/merge_47353088/marked.bam"; then
  BAM=$L1/merge_47353088/marked.bam
else
  echo 'BLOCKED: no L1 marked BAM' >&2
  exit 2
fi
REF=$L1/reference_47349911/GRCh38_no_alt.fa
test -s "$BAM" && test -s "$REF" && test -s "$REF.fai"
VCF=
for cand in \
  "$BASE/out/critical_annotation_47397018/with_clinvar.vcf.gz" \
  "$BASE/out/critical_annotation_47392756/with_clinvar.vcf.gz"
do
  if test -s "$cand"; then
    VCF=$cand
    break
  fi
done
test -n "$VCF" && test -s "$VCF"
VCF_SHA=6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d
module load singularity/4.1.5
module load intel/2023.2.0
module load samtools/1.23.1
SAM=$BASE/containers/samtools_1.24--h9dcdb79_1.sif
BCF=$BASE/containers/bcftools_1.24--h487d631_1.sif
GATK_SIF=$BASE/containers/gatk4_4.6.2.0--py310hdfd78af_0.sif
MOS=$BASE/containers/mosdepth_0.3.14--h87be163_2.sif
VBID=$BASE/containers/verifybamid2_2.0.3--hc004090_0.sif
test -s "$SAM" && test -s "$BCF"
sam() { singularity exec --cleanenv "$SAM" samtools "$@"; }
bcf() { singularity exec --cleanenv "$BCF" bcftools "$@"; }
gatk() { singularity exec --cleanenv "$GATK_SIF" gatk "$@"; }
