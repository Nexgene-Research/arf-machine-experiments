#!/bin/bash
#SBATCH --job-name=L1_B_merge_markdup
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
RUN=$BASE/out/l1_20261015
ALIGN_ARRAY=${1:?Pass completed alignment array ID}
OUT=$RUN/merge_${SLURM_JOB_ID}
mkdir "$OUT"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; printf "exit_code=%s\n" "$rc" > "$OUT/exit_code.txt"; exit "$rc"' EXIT
INPUTS=()
for lane in 1 2 3 4; do
    d=$RUN/align_${ALIGN_ARRAY}/lane_${lane}
    test -s "$d/success" && test -s "$d/validated.json"
    INPUTS+=("$(cat "$d/success")")
done
module load singularity/4.1.5
module load intel/2023.2.0
module load samtools/1.23.1
SAM=$BASE/containers/samtools_1.24--h9dcdb79_1.sif
GATK=$BASE/containers/gatk4_4.6.2.0--py310hdfd78af_0.sif
sam() { singularity exec --cleanenv "$SAM" samtools "$@"; }
sam quickcheck -v "${INPUTS[@]}"
sam merge -@ 15 -o "$OUT/merged.bam" "${INPUTS[@]}"
sam quickcheck -v "$OUT/merged.bam"
mkdir "$OUT/picard_tmp"
singularity exec --cleanenv "$GATK" gatk --java-options '-Xmx96g -XX:ParallelGCThreads=8' MarkDuplicates -I "$OUT/merged.bam" -O "$OUT/marked.bam" -M "$OUT/duplicate_metrics.txt" --TMP_DIR "$OUT/picard_tmp" --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 --CREATE_INDEX true --VALIDATION_STRINGENCY STRICT
sam quickcheck -v "$OUT/marked.bam"
test -s "$OUT/marked.bai"
sam flagstat -@ 15 -O json "$OUT/marked.bam" > "$OUT/flagstat.json"
sam idxstats "$OUT/marked.bam" > "$OUT/idxstats.tsv"
sam view -H "$OUT/marked.bam" > "$OUT/header.sam"
sam stats -@ 15 "$OUT/marked.bam" > "$OUT/samtools_stats.txt"
python3 - "$OUT" "$RUN/align_${ALIGN_ARRAY}" <<'PY'
import json, pathlib, sys
out, align = map(pathlib.Path, sys.argv[1:])
lanes = [json.loads((align/('lane_%d'%i)/'validated.json').read_text()) for i in range(1,5)]
stat = json.loads((out/'flagstat.json').read_text())
assert stat['QC-passed reads']['primary'] == sum(x['primary_reads'] for x in lanes)
rg = [dict(x.split(':',1) for x in l.split('\t')[1:]) for l in (out/'header.sam').read_text().splitlines() if l.startswith('@RG\t')]
assert len(rg) == 4 and {x['ID'] for x in rg} == {x['read_group']['ID'] for x in lanes}
assert len({x['SM'] for x in rg}) == 1
with (out/'validated.json').open('x') as f:
    json.dump(dict(primary_reads=stat['QC-passed reads']['primary'], read_groups=rg, flagstat=stat), f, indent=2)
print('MERGE_MARKDUP_VALIDATED', stat['QC-passed reads']['primary'])
PY
set -C
printf '%s\n' "$OUT/marked.bam" > "$RUN/bam.success"
