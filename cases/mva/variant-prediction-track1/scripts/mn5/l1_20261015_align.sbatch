#!/bin/bash
#SBATCH --job-name=L1_B_lane
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
BASE=$MVA_BASE
RUN=$BASE/out/l1_20261015
TASK=${SLURM_ARRAY_TASK_ID:?}
test "$TASK" -ge 1 && test "$TASK" -le 4
OUT=$RUN/align_${SLURM_ARRAY_JOB_ID}/lane_${TASK}
mkdir -p "${OUT%/*}"
mkdir "$OUT"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; printf "exit_code=%s\n" "$rc" > "$OUT/exit_code.txt"; exit "$rc"' EXIT
REF=$(cat "$RUN/reference.success")
META=$(cat "$RUN/metadata.success")
for suffix in 0123 amb ann bwt.2bit.64 pac fai; do test -s "$REF.$suffix"; done
IFS=$'\t' read -r task lane r1 r2 rgid sm lb pu pl < <(sed -n "$((TASK+1))p" "$META/lanes.tsv")
test "$task" = "$TASK"
test -s "$r1" && test -s "$r2"
module load singularity/4.1.5
module load intel/2023.2.0
module load samtools/1.23.1
SAM=$BASE/containers/samtools_1.24--h9dcdb79_1.sif
BWA=$BASE/containers/bwa-mem2_2.3--he70b90d_0.sif
sam() { singularity exec --cleanenv "$SAM" samtools "$@"; }
RG="@RG\\tID:${rgid}\\tSM:${sm}\\tLB:${lb}\\tPU:${pu}\\tPL:${pl}"
printf 'task=%s lane=%s RG=%s\n' "$TASK" "$lane" "$RG"
# 14 alignment workers + 4 sort threads + I/O/main threads <= 20 CPU budget.
/usr/bin/time -v singularity exec --cleanenv "$BWA" bwa-mem2 mem -t 14 -K 100000000 -R "$RG" "$REF" "$r1" "$r2" 2> "$OUT/bwa.log" |
    sam sort -@ 3 -m 2G -T "$OUT/sorttmp" -o "$OUT/lane.sorted.bam" -
sam quickcheck -v "$OUT/lane.sorted.bam"
sam index -@ 7 "$OUT/lane.sorted.bam"
sam flagstat -@ 7 -O json "$OUT/lane.sorted.bam" > "$OUT/flagstat.json"
sam idxstats "$OUT/lane.sorted.bam" > "$OUT/idxstats.tsv"
sam view -H "$OUT/lane.sorted.bam" > "$OUT/header.sam"
python3 - "$OUT" "$REF.fai" "$rgid" "$sm" "$pu" <<'PY'
import json, pathlib, sys
out, fai, rgid, sm, pu = sys.argv[1:]
p = pathlib.Path(out)
lines = (p/'header.sam').read_text().splitlines()
sq = [dict(x.split(':', 1) for x in l.split('\t')[1:]) for l in lines if l.startswith('@SQ\t')]
expected = [(x.split('\t')[0], int(x.split('\t')[1])) for x in open(fai)]
assert [(s['SN'], int(s['LN'])) for s in sq] == expected
rg = [dict(x.split(':', 1) for x in l.split('\t')[1:]) for l in lines if l.startswith('@RG\t')]
assert len(rg) == 1 and rg[0]['ID'] == rgid and rg[0]['SM'] == sm and rg[0]['PU'] == pu
stat = json.loads((p/'flagstat.json').read_text())
v = stat['QC-passed reads']; failed = stat['QC-failed reads']
assert v['primary'] > 0 and v['primary'] % 2 == 0 and v['read1'] == v['read2']
assert failed['total'] == 0
assert v['mapped'] > 0
with (p/'validated.json').open('x') as f:
    json.dump(dict(reference_contigs=len(sq), read_group=rg[0], primary_reads=v['primary'], flagstat=stat), f, indent=2)
print('LANE_VALIDATED', rgid, v['primary'])
PY
set -C
printf '%s\n' "$OUT/lane.sorted.bam" > "$OUT/success"
