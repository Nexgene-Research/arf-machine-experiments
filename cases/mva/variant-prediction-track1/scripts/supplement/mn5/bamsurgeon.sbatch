#!/bin/bash
#SBATCH --job-name=SUP_spike
#SBATCH --account=ehpc926
#SBATCH --qos=acc_ehpc
#SBATCH --partition=acc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
source "${MVA_BASE:?}/jobs/track1_supplement/common.sh"
OUT=$OUTROOT/bamsurgeon_${SLURM_JOB_ID:?}
mkdir "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$OUT/exit_code.txt"; exit $rc' EXIT
IMG=$BASE/containers/bamsurgeon_1.4.1--pyhdfd78af_0.sif
test -s "$IMG"
# GRCh38-mapped published alleles only. Others stay unresolved.
# L1 REF is bwa-mem2 indexed; addsnv's default bwa aln needs its own index,
# so each chromosome FASTA is indexed inside the container.
cat > "$OUT/spikes.tsv" <<'EOF'
chr15	40117088	G	A	Ochiai2014_BUB1B
chr2	110677994	A	C	Carvalhal2022_BUB1
chr5	914504	C	T	Yost2017_TRIP13_R354X
chr5	907987	G	C	Yost2017_TRIP13_splice
EOF
printf 'label\tchrom\tpos\tref\talt\tstatus\trank_note\n' > "$OUT/positive_controls.tsv"
while read -r chrom pos ref alt label; do
  region=${chrom}:$((pos-10000))-$((pos+10000))
  slice=$OUT/${label}.bam
  sam view -b "$BAM" "$region" > "$slice"
  sam index "$slice"
  chromfa=$OUT/${chrom}.fa
  if test -s "$OUTROOT/bamsurgeon_47433723/${chrom}.fa.bwt"; then
    chromfa=$OUTROOT/bamsurgeon_47433723/${chrom}.fa
  elif ! test -s "$chromfa.bwt"; then
    sam faidx "$REF" "$chrom" > "$chromfa"
    singularity exec --cleanenv "$IMG" bwa index "$chromfa"
  fi
  varfile=$OUT/${label}.var
  printf '%s\t%s\t%s\t1.0\n' "$chrom" "$pos" "$pos" > "$varfile"
  spiked=$OUT/${label}.spiked.bam
  if singularity exec --cleanenv --pwd "$OUT" "$IMG" addsnv.py \
      -v "$varfile" -f "$slice" -r "$chromfa" -o "$spiked" \
      -p 8 --mindepth 10 -d 0.5; then
    sam index "$spiked"
    bcf mpileup -f "$REF" -r "$region" -Ou "$spiked" | bcf call -m -v -Ov -o "$OUT/${label}.vcf"
    if grep -P "^${chrom}\t${pos}\t" "$OUT/${label}.vcf" | grep -q "$alt"; then
      printf '%s\t%s\t%s\t%s\t%s\trecovered\tregional_bcftools_call\n' \
        "$label" "$chrom" "$pos" "$ref" "$alt" >> "$OUT/positive_controls.tsv"
    else
      printf '%s\t%s\t%s\t%s\t%s\tnot_recovered\tregional_bcftools_call\n' \
        "$label" "$chrom" "$pos" "$ref" "$alt" >> "$OUT/positive_controls.tsv"
    fi
  else
    printf '%s\t%s\t%s\t%s\t%s\tspike_failed\taddsnv_coverdiff_or_remap\n' \
      "$label" "$chrom" "$pos" "$ref" "$alt" >> "$OUT/positive_controls.tsv"
  fi
done < "$OUT/spikes.tsv"
printf 'CEP57\t.\t.\t.\t.\tunresolved\tsource_cDNA_not_mapped_to_GRCh38\n' >> "$OUT/positive_controls.tsv"
printf 'CENATAC\t.\t.\t.\t.\tunresolved\tsource_cDNA_not_mapped_to_GRCh38\n' >> "$OUT/positive_controls.tsv"
printf 'MAD1L1\t.\t.\t.\t.\tunresolved\tsource_cDNA_not_mapped_to_GRCh38\n' >> "$OUT/positive_controls.tsv"
printf '%s\n' "$OUT/positive_controls.tsv" > "$OUT/success"
