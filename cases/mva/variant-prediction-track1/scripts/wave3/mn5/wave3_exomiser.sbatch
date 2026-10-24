#!/bin/bash
#SBATCH --job-name=W3_exomiser
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=acc
#SBATCH --time=16:00:00
source "${MVA_BASE:?}/jobs/track1_wave3/common.sh"
OUT=$OUTROOT/exomiser_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
module load java-openjdk/17.0.11+9 || module load java-openjdk/21.0.3+9
DATA=$BASE/annot/exomiser
mkdir -p "$DATA"
fetch() {
  dest=$1
  url=$2
  if test -s "$dest"; then
    return 0
  fi
  curl -L --fail --retry 3 --max-time 600 -o "$dest" "$url" \
    || wget -O "$dest" "$url" \
    || return 1
}
fetch "$DATA/exomiser-cli-14.1.0-distribution.zip" \
  https://github.com/exomiser/exomiser/releases/download/14.1.0/exomiser-cli-14.1.0-distribution.zip || true
fetch "$DATA/2406_hg38.zip" \
  https://data.monarchinitiative.org/exomiser/latest/2406_hg38.zip || true
fetch "$DATA/2406_phenotype.zip" \
  https://data.monarchinitiative.org/exomiser/latest/2406_phenotype.zip || true
if test -s "$DATA/2406_hg38.zip" && ! test -d "$DATA/2406_hg38"; then
  unzip -n "$DATA/2406_hg38.zip" -d "$DATA"
fi
if test -s "$DATA/2406_phenotype.zip" && ! test -d "$DATA/2406_phenotype"; then
  unzip -n "$DATA/2406_phenotype.zip" -d "$DATA"
fi
if ! test -s "$DATA/exomiser-cli-14.1.0-distribution.zip" || ! test -d "$DATA/2406_hg38" || ! test -d "$DATA/2406_phenotype"; then
  echo 'BLOCKED: Exomiser CLI or 2406 data pack missing after download attempt' > "$OUT/caveat.txt"
  python3 "$JOBS/exomiser_top30.py" "$OUT/missing.genes.tsv" "$OUT/exomiser_top30.tsv"
  exit 2
fi
unzip -n "$DATA/exomiser-cli-14.1.0-distribution.zip" -d "$OUT"
CLI_DIR=$(find "$OUT" -maxdepth 1 -type d -name 'exomiser-cli-*' | head -1)
test -n "$CLI_DIR"
python3 "$JOBS/extract_rare_vcf.py" "$VCF" "$JOBS/spliceai_keys.tsv" "$OUT/rare.vcf.gz"
bcf index -t "$OUT/rare.vcf.gz"
python3 - "$OUT" "$JOBS/hpo.txt" <<'PY'
from pathlib import Path
import sys
out = Path(sys.argv[1])
hpo = [ln.strip() for ln in Path(sys.argv[2]).read_text().splitlines() if ln.strip()]
hpo_s = ", ".join(hpo)
(out / "analysis.yml").write_text(
    f"""analysis:
  genomeAssembly: hg38
  vcf: {out}/rare.vcf.gz
  proband: WGS_EX2312012
  hpoIds: [{hpo_s}]
  analysisMode: PASS_ONLY
  inheritanceModes:
    AUTOSOMAL_DOMINANT: 0.1
    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1
    AUTOSOMAL_RECESSIVE_COMP_HET: 0.2
    X_RECESSIVE: 0.1
    X_DOMINANT: 0.1
  frequencySources: [THOUSAND_GENOMES, TOPMED, UK10K, GNOMAD_E, GNOMAD_G]
  pathogenicitySources: [REVEL, MVP, ALPHA_MISSENSE, SPLICE_AI]
  hiPhivePrioritiser: [human]
  frequencyFilter: 0.01
outputOptions:
  outputDirectory: {out}
  outputFileName: exomiser
  outputFormats: [TSV_GENE, TSV_VARIANT, HTML]
"""
)
PY
cat > "$OUT/application.properties" <<EOF
exomiser.phenotype.data-version=2406
exomiser.hg38.data-version=2406
exomiser.data-directory=$DATA
exomiser.hg38.data-directory=$DATA/2406_hg38
exomiser.phenotype.data-directory=$DATA/2406_phenotype
EOF
JAR=$(ls "$CLI_DIR"/exomiser-cli-*.jar | head -1)
java -Xmx12g -jar "$JAR" \
  --analysis "$OUT/analysis.yml" \
  --spring.config.location="$OUT/application.properties" \
  || echo 'Exomiser CLI failed' >> "$OUT/caveat.txt"
GENES=$(find "$OUT" -name '*.genes.tsv' | head -1)
python3 "$JOBS/exomiser_top30.py" "${GENES:-$OUT/missing.genes.tsv}" "$OUT/exomiser_top30.tsv"
test -s "$OUT/exomiser_top30.tsv"
printf '%s\n' "$OUT" > "$OUT/success"
