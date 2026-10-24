#!/bin/bash
#SBATCH --job-name=CRITICAL_annotate
#SBATCH --account=ehpc926
#SBATCH --qos=acc_debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=01:50:00
set -euo pipefail
umask 077
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${MVA_BASE:?Set MVA_BASE to the MN5 work root}"
B=$MVA_BASE
I=$B/out/fastpath_20261016/full_47392756/with_gnomad.vcf.gz
R=$B/out/l1_20261015/reference_47349911/GRCh38_no_alt.fa
O=$B/out/critical_annotation_${SLURM_JOB_ID:?}
mkdir "$O"
exec > >(tee "$O/run.log") 2>&1
trap 'rc=$?; echo "$rc" > "$O/exit_code.txt"; exit "$rc"' EXIT
module load singularity/4.1.5
BCF=$B/containers/bcftools_1.24--h487d631_1.sif
SNP=$B/containers/snpeff_5.2--hdfd78af_1.sif
DB=$B/annot/snpEff_v5_2_GRCh38.mane.1.0.ensembl.zip
CLIN=$B/annot/clinvar.vcf.gz
test -s "$I" && test -s "$R.fai" && test -s "$BCF" && test -s "$SNP" && test -s "$DB" && test -s "$CLIN.tbi"
echo '4aa21a3477584295a08fdfb22a577b60612d1d7b6cfc88a367d737a89b4eee5c  '"$I" | sha256sum -c -
sha256sum "$I" "$R" "$BCF" "$SNP" "$DB" "$CLIN" "$B/annot/AlphaMissense_hg38.tsv.gz" "$B/annot/revel-v1.3_all_chromosomes.zip" > "$O/inputs.sha256"
bcf() { singularity exec --cleanenv "$BCF" bcftools "$@"; }
bcf --version > "$O/bcftools.version"
if ! test -s "$I.csi" && ! test -s "$I.tbi"; then bcf index -c "$I"; fi
bcf index -n "$I" > "$O/records_in.txt"
test "$(cat "$O/records_in.txt")" = 5035176
mkdir "$O/snpeff"
unzip -q "$DB" -d "$O/snpeff"
printf 'data.dir = %s/snpeff/data\nGRCh38.mane.1.0.ensembl.genome : GRCh38.mane.1.0.ensembl\n' "$O" > "$O/snpeff.config"
singularity exec --cleanenv "$SNP" snpEff -version > "$O/snpeff.version" 2>&1
# Annotate every normalized allele, including non-PASS and common alleles, BEFORE selection.
singularity exec --cleanenv "$SNP" snpEff -c "$O/snpeff.config" -dataDir "$O/snpeff/data" -noStats GRCh38.mane.1.0.ensembl "$I" 2> "$O/snpeff.stderr" |
    bcf view -Oz -o "$O/with_consequence.vcf.gz" -
bcf index -c "$O/with_consequence.vcf.gz"
test "$(bcf index -n "$O/with_consequence.vcf.gz")" = 5035176
bcf annotate -a "$CLIN" -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNHGVS -Oz -o "$O/with_clinvar.vcf.gz" "$O/with_consequence.vcf.gz"
bcf index -c "$O/with_clinvar.vcf.gz"
test "$(bcf index -n "$O/with_clinvar.vcf.gz")" = 5035176
# Full annotated VCF and index are retained on MN5, never shipped to the laptop.
# Selection runs after all 5,035,176 alleles have consequence and ClinVar annotations.
python3 - "$O" "$R.fai" <<'PY'
import csv,gzip,json,sys,collections
from pathlib import Path
o=Path(sys.argv[1]); allowed={x.split('\t')[0] for x in Path(sys.argv[2]).read_text().splitlines()}
stats=collections.Counter(); consequences=collections.Counter()
tiers={'BUB1B','CEP57','TRIP13','CENATAC','MAD1L1','BUB1'}
head=['chrom','pos','ref','alt','raw_record','filter','qual','gt','gq','dp','ad','vaf','gene','transcript','consequence','impact','hgvs','gnomad_af','gnomad_popmax','gnomad_homalt','af_state','clinvar','alphamissense','revel','phase']
def info(s):return dict(z.split('=',1) if '=' in z else (z,'') for z in s.split(';'))
def num(x):
 try:return float(x)
 except (TypeError,ValueError):return None
with gzip.open(o/'with_clinvar.vcf.gz','rt') as f, (o/'candidates_full.tsv').open('w',newline='') as out, (o/'compound_het_candidates.tsv').open('w',newline='') as pairs, (o/'homozygous_candidates.tsv').open('w',newline='') as hom, (o/'singleton_candidates.tsv').open('w',newline='') as single, (o/'candidate_universe.vcf').open('w') as vcf:
 writers=[csv.writer(h,delimiter='\t') for h in (out,pairs,hom,single)]
 for w in writers:w.writerow(head)
 by_gene=collections.defaultdict(list)
 for line in f:
  if line.startswith('#'):
   vcf.write(line);continue
  stats['annotated']+=1
  a=line.rstrip('\n').split('\t'); chrom,pos,_,ref,alt,qual,flt=a[:7]
  if chrom not in allowed:stats['dropped_nonreference_contig']+=1;continue
  if not (chrom in {'chrX','chrY','chrM'} or chrom[3:].isdigit() and 1<=int(chrom[3:])<=22):stats['dropped_nonprimary']+=1;continue
  if alt=='*' or set(ref+alt)-set('ACGT'):stats['symbolic_or_ambiguous']+=1;continue
  inf=info(a[7]); ann=inf.get('ANN','').split(','); best=ann[0].split('|') if ann else []
  gene=best[3] if len(best)>3 else '.'; effect=best[1] if len(best)>1 else '.'
  impact=best[2] if len(best)>2 else '.'; tx=best[6] if len(best)>6 else '.'
  hgvs=best[9] if len(best)>9 else '.'; consequences[effect]+=1
  af=num(inf.get('gnomad_AF')); mx=num(inf.get('gnomad_faf95max'))
  maxaf=max([v for v in (af,mx) if v is not None],default=None)
  state='missing' if maxaf is None else 'observed'
  tier=gene in tiers
  clin=inf.get('CLNSIG','.')
  clinical=any(x in clin.lower() for x in ('pathogenic','uncertain_significance','vus')) and 'benign' not in clin.lower()
  rare=maxaf is not None and maxaf<=.01
  if rare:stats['af_le_0.01']+=1
  if state=='missing':stats['af_missing']+=1
  if tier:stats['tier1_rescue']+=1
  if clinical:stats['clinvar_p_lp_vus']+=1
  if not (rare or state=='missing' or tier or clinical):stats['prefilter_excluded']+=1;continue
  stats['prefilter_universe']+=1
  consequence=impact in ('HIGH','MODERATE') or any(k in effect.lower() for k in ('splice','utr'))
  if not (consequence or clinical or tier):stats['consequence_excluded']+=1;continue
  stats['biologically_filtered']+=1
  fmt=dict(zip(a[8].split(':'),a[9].split(':'))) if len(a)>9 else {}
  ad=fmt.get('AD','.'); adval=[num(x) for x in ad.split(',')]; vaf=adval[1]/sum(adval[:2]) if len(adval)>=2 and None not in adval[:2] and sum(adval[:2])>0 else '.'
  row=[chrom,pos,ref,alt,inf.get('FP_ORIG','.'),flt,qual,fmt.get('GT','.'),fmt.get('GQ','.'),fmt.get('DP','.'),ad,vaf,gene,tx,effect,impact,hgvs,inf.get('gnomad_AF','.'),inf.get('gnomad_faf95max','.'),inf.get('gnomad_HomAlt','.'),state,clin,'.','.', 'unknown']
  writers[0].writerow(row);vcf.write(line)
  gt=fmt.get('GT','.').replace('|','/')
  if gt=='1/1':writers[2].writerow(row)
  elif gt in ('0/1','1/0'):by_gene[gene].append(row)
  elif gt in ('1','.','./.'):writers[3].writerow(row)
 for gene,rows in by_gene.items():
  for row in rows:(writers[1] if len(rows)>=2 else writers[3]).writerow(row)
stats['same_gene_pair_candidate_alleles']=sum(len(r) for r in by_gene.values() if len(r)>=2)
(o/'consequence_counts.tsv').write_text('class\tcount\n'+''.join(f'{k}\t{v}\n' for k,v in sorted(consequences.items())))
(o/'candidate_universe_manifest.json').write_text(json.dumps({'source':str(o/'with_clinvar.vcf.gz'),'no_preconsequence_cap':True,'counts':stats,'predictor_scores':'not_assigned_in_this_stage'},indent=2,sort_keys=True))
PY
bcf view -Oz -o "$O/candidate_universe.vcf.gz" "$O/candidate_universe.vcf"
bcf index -c "$O/candidate_universe.vcf.gz"
sha256sum "$O"/candidates_full.tsv "$O"/compound_het_candidates.tsv "$O"/homozygous_candidates.tsv "$O"/singleton_candidates.tsv "$O"/candidate_universe_manifest.json "$O"/consequence_counts.tsv "$O"/candidate_universe.vcf.gz "$O"/candidate_universe.vcf.gz.csi > "$O/outputs.sha256"
echo ANNOTATION_COMPLETE > "$O/success"
