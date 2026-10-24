#!/bin/bash
# Shared MN5 paths for the Track 1 rigor-pass jobs. Source from sbatch.
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
JOBS=$BASE/jobs/track1_rigor
OUTROOT=$BASE/out/track1_rigor
mkdir -p "$OUTROOT"
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
REF=$BASE/out/l1_20261015/reference_47349911/GRCh38_no_alt.fa
test -s "$REF" && test -s "$REF.fai"
GATK_ROOT=$BASE/out/l1_20261015/call_gatk_47354613
CLIN=$BASE/annot/clinvar.vcf.gz
test -s "$CLIN"
module load singularity/4.1.5
module load intel/2023.2.0
BCF_SIF=$BASE/containers/bcftools_1.24--h487d631_1.sif
test -s "$BCF_SIF"
bcf() { singularity exec --cleanenv "$BCF_SIF" bcftools "$@"; }
# Frozen annotated VCF from lineage.md (job 47397018).
VCF_SHA=6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d
