#!/bin/bash
#SBATCH --job-name=track1_reads
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:50:00
set -euo pipefail
umask 077
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
B=$MVA_BASE
O=$B/out/track1_47397018
BAM=$B/out/l1_20261015/merge_47353088/marked.bam
exec > >(tee -a "$O/reads.log") 2>&1
module load intel
module load samtools/1.23.1
module load singularity/4.1.5
test -s "$BAM"
if ! test -s "$BAM.bai"; then
  if test -s "${BAM%.bam}.bai"; then ln -sf "${BAM%.bam}.bai" "$BAM.bai"
  else samtools index -@ 4 "$BAM"
  fi
fi
python3 "$B/jobs/track1/track1_reads.py" "$O" "$O"
python3 "$B/jobs/track1/track1_epcr.py" "$O" \
  "$B/out/critical_annotation_47397018/with_clinvar.vcf.gz"
echo READS_EPCR_COMPLETE > "$O/reads.status"
