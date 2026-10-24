# Compound heterozygous *BUB1B* variants p.Leu737\* and p.Asn1002Lys identified by whole-genome sequencing in a child with mosaic variegated aneuploidy

**İlsu Duman, Özge Beğli, Hatice Mutlu, Kutsal Baran Ozkurt\***\
Nexgene AI ([nexgene.ai](https://nexgene.ai)). Correspondence: [kutsal.ozkurt@nexgene.ai](mailto:kutsal.ozkurt@nexgene.ai)

*MVA Hackathon 2026, Track 1 (variant prediction). Team: Nexgene AI. Reference genome: GRCh38 (no-alt analysis set). Report date: 24 October 2026.*

> **Answer.** The child carries two heterozygous variants in *BUB1B*, the gene for autosomal recessive mosaic variegated aneuploidy syndrome 1 (MVA1; OMIM 257300):
>
> - **`chr15:40209701 T>G`**, NM_001211.6:c.2210T>G, p.(Leu737\*). This is a stop-gain in exon 17 of 23, predicted to trigger nonsense-mediated decay. It is reported Pathogenic/Likely pathogenic in ClinVar (VCV000533901, two stars).
> - **`chr15:40220612 T>G`**, NM_001211.6:c.3006T>G, p.(Asn1002Lys). This is a novel missense in the BUBR1 kinase domain. Its AlphaMissense score of 0.923 meets the calibrated PP3_Moderate threshold.
>
> Together they reproduce the genotype class that defines biallelic MVA1: one truncating allele combined with one kinase-domain substitution that destabilises BUBR1 [1, 2]. This pair is the single row of our submission. The automated Track 1 scorer awarded it full rank points, so it agrees with the clinically confirmed result.

![Graphical abstract. Seven HPO terms and blood whole-genome sequencing were reduced from 5,035,176 variant records to one *BUB1B* genotype. The two alleles sit upstream of and within the BUBR1 kinase domain. Six orthogonal evidence lines support the call; physical phase is the remaining confirmatory step.](figures/graphical_abstract.png)

## Abstract

**Background.** Mosaic variegated aneuploidy (MVA) is an ultra-rare recessive disorder. It is characterised by constitutional mosaic gains and losses of whole chromosomes, growth restriction and a high risk of childhood cancer. Most families carry biallelic variants in *BUB1B*, which encodes the spindle assembly checkpoint kinase BUBR1. Rarer forms involve *CEP57*, *TRIP13*, *CENATAC* and *MAD1L1*.

**Methods.** We analysed paired-end blood whole-genome sequencing (median autosomal depth 44.4×) and seven phenotype terms from one affected child. A pre-registered, four-angle strategy ran blinded to candidate outputs: genotype, phenotype, mechanism and published cases. It was combined with an independent re-alignment on the MareNostrum 5 supercomputer, three-caller recall on that alignment, and targeted read-level validation. Read-backed and caller-based phasing, copy-number, structural, splicing and genome-wide aneuploidy screens, and calibrated variant-effect predictors followed. Variants were classified with the ACMG/AMP framework and ClinGen SVI refinements.

**Results.** From 5,035,176 variant records and 246,062 genotype hypotheses, two rare heterozygous *BUB1B* variants were prioritised: p.Leu737\* and p.Asn1002Lys.

- *Authenticity.* Both are supported by high-quality bidirectional reads. They are recalled by an independent alignment and by two additional germline callers (DeepVariant and Strelka2) on that alignment, and have balanced allele fractions (0.55 and 0.45).
- *Rarity and gene context.* They are rare (gnomAD AF 7.4×10⁻⁵ and 5×10⁻⁷). *BUB1B* is the only one of six MVA genes carrying two rare coding alleles.
- *Effect.* The stop-gain meets PVS1, with the premature stop 747 nt upstream of the last exon junction. The missense meets PP3_Moderate (AlphaMissense 0.923) and lies ten residues from the published MVA kinase-domain substitution L1012P.
- *Phenotype.* Phenotype-only ranking placed *BUB1B* first among eight candidate genes, and it stayed first in 83 of 100 weight permutations.
- *Competing explanations.* We found no copy-number change, structural variant or splice-altering second hit in any MVA gene. None of 13 published pathogenic sites was present. Contamination was negligible (FREEMIX 2×10⁻⁹), and there was no long run of homozygosity.
- *Phase.* The alleles lie 10,911 bp apart, beyond short-read phasing range.

**Conclusions.** The compound heterozygous *BUB1B* genotype p.[Leu737\*];[Asn1002Lys] explains the child's MVA. It is concordant with the clinically confirmed answer. Parental genotyping or long-read phasing, together with a BUBR1 abundance assay, would raise p.Asn1002Lys from uncertain significance to likely pathogenic under current ACMG/AMP rules.

**Keywords:** mosaic variegated aneuploidy; *BUB1B*; BUBR1; spindle assembly checkpoint; compound heterozygosity; whole-genome sequencing; rare disease; variant interpretation

## 1. Introduction

### 1.1 Mosaic variegated aneuploidy and the spindle assembly checkpoint

Accurate chromosome segregation depends on the spindle assembly checkpoint (SAC). This surveillance pathway delays anaphase until every kinetochore is attached to spindle microtubules [3]. BUBR1, encoded by *BUB1B*, is a core SAC component. Its N-terminal TPR region and KEN boxes bind CDC20 as part of the mitotic checkpoint complex. Its central region binds BUB3 and kinetochore scaffolds (KNL1). Its C-terminal kinase-like domain (residues 766–1050 in UniProt O60566) is required for protein stability rather than for catalytic activity [2]. Reduced BUBR1 abundance weakens the checkpoint and destabilises kinetochore–microtubule attachments. The result is random gains and losses of whole chromosomes in different cells: "variegated" aneuploidy.

Biallelic *BUB1B* variants were first identified as the cause of MVA by Hanks et al. (2004) [1], establishing MVA1 (OMIM 257300). Affected children typically show:

- intrauterine and postnatal growth restriction;
- microcephaly;
- mosaic aneuploidy of many different chromosomes in lymphocytes and fibroblasts, often with premature chromatid separation;
- a strongly increased risk of embryonal tumours, notably rhabdomyosarcoma, Wilms tumour and leukaemia.

Matsuura et al. (2006) described monoallelic truncating *BUB1B* variants combined with a low-expressing haplotype in Japanese families with premature chromatid separation (PCS) syndrome [4]. Ochiai et al. (2014) later mapped that hypomorphic haplotype to an intergenic variant upstream of *BUB1B* [5]. A homozygous *BUB1B* variant has also been associated with gastrointestinal neoplasia [6].

Other genes cause clinically related phenotypes:

- *CEP57* (MVA2) [7];
- *TRIP13* (MVA3; Wilms tumour with chromosome missegregation) [8];
- the minor-spliceosome component *CENATAC* [9];
- *MAD1L1*, with very high tumour susceptibility [10];
- *BUB1*, a paralogue of *BUB1B*; biallelic *BUB1* variants are associated with microcephaly and developmental delay rather than MVA [11].

Allelic series in mice show that the type and combination of *Bub1b* alleles determine phenotype severity [12].

### 1.2 The molecular rule that guides interpretation

The decisive mechanistic insight for this case comes from Suijkerbuijk et al. (2010) [2]. They studied patient cell lines carrying biallelic *BUB1B* mutations and found that almost every compound heterozygote combined **one truncating allele with one amino-acid substitution located in or near the kinase domain**. Truncating alleles produced no stable transcript. Kinase-domain substitutions (for example R814H, L844F, I909T, Q921H and L1012P) produced protein with a 5–10-fold lower steady-state level through increased turnover. Overall BUBR1 abundance therefore fell below the threshold required for checkpoint function. Re-expressing wild-type BUBR1 rescued the defect.

This "null plus destabilised kinase domain" model predicts exactly the genotype we report here.

### 1.3 The case and the task

The hackathon released blood whole-genome sequencing (FASTQ and a supplied variant call set) from a child with MVA, together with a clinical description. The description was encoded in seven Human Phenotype Ontology (HPO) terms [13]:

- rhabdomyosarcoma (HP:0002859);
- nephrocalcinosis (HP:0000121);
- short stature (HP:0004322);
- failure to thrive (HP:0001508);
- skeletal muscle atrophy (HP:0003202);
- premature birth (HP:0001622);
- small for gestational age (HP:0001518).

A history of recurrent spontaneous abortion (HP:0200067) belongs to the parents and was not scored as a feature of the child.

Track 1 asks for a ranked list of the variant(s) driving the condition. Our aim was to produce a single, defensible genotype rather than a long list, and to document every piece of evidence for and against it.

## 2. Data and methods

### 2.1 Data and computing environment

The input data were:

- paired-end blood WGS FASTQ files;
- the organisers' variant calls (Sentieon pipeline; submitted as `PROBAND01`);
- the clinical description.

The controlled data were processed only on the MareNostrum 5 (MN5) supercomputer at the Barcelona Supercomputing Center and on an access-controlled workstation. Tool versions, container digests, job identifiers and file hashes are recorded in Appendix A.

### 2.2 Pre-registration

Before any candidate record was inspected we fixed the analysis sequence, the Tier-1 gene set (*BUB1B*, *CEP57*, *TRIP13*, *CENATAC*, *MAD1L1*, *BUB1*), phenotype weights, discriminating predictions for seven alternative mechanisms, and the probability policy (the pre-registered analysis plan; Appendix A).

Four evidence angles were kept blinded until convergence:

- (A) genotype, rarity and predicted effect across coding, splice, regulatory, CNV/SV and mosaic classes;
- (B) phenotype only;
- (C) mechanism only;
- (D) published cases only.

Equal scores were resolved by a dated addendum to the analysis plan (Appendix A).

### 2.3 Alignment, variant calling and annotation

We realigned the FASTQs independently on MN5 against GRCh38 no-alt with bwa-mem2 2.3 [14], then marked duplicates. We used bcftools 1.24 [15] for targeted recall and GATK [16] for physical-phasing tags.

On the same bwa-mem2 BAM we ran three germline callers on the autosomes: DeepVariant [42], GATK HaplotypeCaller [16] and Strelka2 [43]. An allele counted as a vote when the caller reported a PASS filter and a diploid ALT genotype (0/1 or 1/1). Raw GATK HaplotypeCaller output has an empty FILTER field because variant-quality score recalibration was not applied, so GATK genotypes were recorded but did not count as a vote. Both submitted alleles reached 2/3 support (DeepVariant and Strelka2 PASS 0/1; GATK 0/1, unfiltered).

The supplied call set (5,035,176 records on primary contigs) was annotated with:

- snpEff 5.2 [17] against MANE Select v1.5 transcripts [18];
- gnomAD allele frequencies [19];
- ClinVar [20];
- AlphaMissense [21] and REVEL [22];
- a precomputed AlphaGenome variant-impact table [23].

All alleles were REF-checked against the project FASTA, left-aligned and decomposed.

A contig-prefix mismatch in the first ClinVar join returned no annotations; the analysis uses the corrected, chromosome-prefixed join.

### 2.4 Prioritisation

Biologically filtered candidates were defined as rare (gnomAD AF ≤ 0.01) protein-altering or splice-region alleles, with rescue of all Tier-1 alleles. They were combined into single-gene genotype hypotheses: homozygous, compound heterozygous, or a single Tier-1 allele.

Each hypothesis received a case-level ranking heuristic, ePCR. Under the pre-registered weakest-link rule, ePCR is bounded by its least-supported link and is capped at 0.1 while a critical link is unverified, such as phase. It is a ranking device, not a calibrated clinical probability. The 0.1 value in the CSV reflects this cap.

The top read-supported hypotheses were adjudicated on five criteria: read-backed authenticity, caller concordance, established gene–disease inheritance, the reviewed molecular effect of the weaker allele, and phenotype fit.

### 2.5 Validation layers

Each layer, its method and its role (file mapping in Appendix A):

| Layer | Method | Role |
|---|---|---|
| Read authentication | MAPQ ≥ 60, base quality ≥ 20 pileups; strand, soft-clip, SA/XA checks | Confirms both alleles in high-quality reads |
| Caller concordance | supplied Sentieon VCF versus independent bwa-mem2 + bcftools recall | Confirms genotype recall on a second alignment |
| Three-caller consensus | DeepVariant, GATK HaplotypeCaller and Strelka2 on the same BAM; 2/3 rule (PASS + diploid ALT) | Independent genotype recall; GATK FILTER empty, not a vote |
| Phasing | WhatsHap 2.8 regional read-backed phasing [24]; GATK PID/PGT tags | Tests whether short reads can place the pair in trans |
| Copy number | mosdepth [25] exon and gene-body depth (MAPQ ≥ 20) | Screens six MVA genes for deletion or duplication |
| Structural variants | Manta [26] and Delly [27] over the six genes | Screens rearrangements as a hidden second hit |
| Splicing | SpliceAI [28] on candidates and 28 rare noncoding keys; AlphaGenome splice table | Screens cryptic splice second hits |
| Missense effect | AlphaMissense and REVEL with ClinGen-calibrated intervals [29, 30]; tail likelihood ratios on ClinVar two-star missense with *BUB1B* held out | Calibrates PP3 for p.Asn1002Lys |
| Aneuploidy | 100 kb arm depth; GATK ModelSegments allele fractions | Bulk dosage screen of all autosomes |
| Sample QC | VerifyBamID2 [31], somalier [32], bcftools/RoH [33] | Contamination, relatedness sketch and runs of homozygosity |
| Published alleles | read-level genotyping of 13 GRCh38-mapped pathogenic sites in MVA genes | Tests known causal alleles |
| Controls | 100 fixed-seed HPO-weight permutations; HG002 truth-set comparison; chr20 HG002/HG003 benchmark [34] | Phenotype-rank stability and pipeline accuracy |
| Knowledge graph | leave-one-out rediscovery of known MVA gene–disease edges in PrimeKG [35] | Independent gene–disease context |

### 2.6 Variant classification

Variants were classified with the ACMG/AMP 2015 standards [36] and the following refinements:

- the ClinGen SVI PVS1 decision tree [37];
- calibrated computational evidence for PP3/BP4 [29, 30];
- splicing guidance [38];
- the Bayesian point scale [39, 40]: Pathogenic ≥ 10 points, Likely pathogenic 6–9, VUS 0–5.

Following the pre-registration, only the strongest calibrated predictor in a correlated family was counted. A predictor without calibration contributed no evidence.

## 3. Results

### 3.1 Sample quality

The genome is technically sound for both germline and dosage analysis:

- Median autosomal depth was 44.4× (100 kb bins, MAPQ ≥ 20).
- The non-pseudoautosomal X ratio was 0.494, consistent with a single X chromosome.
- VerifyBamID2 estimated a contamination fraction (FREEMIX) of 2.1×10⁻⁹ across 10,000 panel SNPs at a mean depth of 42.5×, so there is no detectable contamination or chimerism. Somalier's uncalibrated CHARR screen was 0.008.
- Of 3,079,630 heterozygous calls, 90.2% had allele fractions between 0.25 and 0.75.
- Runs of homozygosity were short. No window was ≥ 1.5 Mb, and F_ROH restricted to long windows was 0. This argues against parental consanguinity or uniparental isodisomy as the route to homozygosity.

### 3.2 From five million variants to one genotype

![**Figure 1. Analysis workflow and candidate funnel.** (a) Stages run on MareNostrum 5 and locally. (b) Unique remaining counts at each layer on a log scale (n = 5,035,176 annotated records to 1 adjudicated genotype). The stages count different units (records, alleles, genotype hypotheses, rows), so the bars are not a single subtractive flow; for example, the 244,650 *MAD1L1* pair combinations are combinatorial hypotheses, not independent diagnoses.](figures/fig1_workflow_funnel.png)

The funnel stages (Figure 1):

1. The 5,035,176 annotated records reduced to a rare prefilter universe of 185,662.
2. Filtering left 2,531 biologically filtered candidate alleles.
3. These formed 246,062 genotype hypotheses: 245,781 compound heterozygous pairs, 280 homozygous genotypes and one Tier-1 singleton.
4. Read checks of the top 60 hypotheses retained 28 SNV-pair hypotheses with clean read support.
5. Ten were carried forward as the initial ranked list of ten genotypes.

Adjudication of that list removed nine rows:

- **BUB1 intronic pairs.** Two intronic *BUB1* alleles have no splice signal (SpliceAI Δ = 0.00) and no established recessive MVA mechanism.
- **BUB1B stop + intron.** The stop plus intron `chr15:40216470 A>G` pair has a second allele with SpliceAI Δ = 0.00 and an AlphaGenome percentile of 47.6.
- **Other pairs.** These had alleles that failed read authentication or exceeded the recessive frequency ceiling.

One genotype remained: *BUB1B* p.Leu737\* with p.Asn1002Lys. This is the submitted CSV row.

### 3.3 The *BUB1B* locus and the predicted consequence of each allele

![**Figure 2. The *BUB1B* locus and the BUBR1 protein.** (a) MANE Select NM_001211.6 gene model, 23 exons on the plus strand of chr15, with both variants. The two alleles are 10,911 bp apart, roughly 20 times the length of a short-read fragment. (b) BUBR1 domains from UniProt O60566 (BUB1 N-terminal/TPR region 62–226; KNL1-binding segment 152–185; D-box 224–232; kinase domain 766–1050). The proband's variants are shown above; MVA-associated substitutions characterised by Suijkerbuijk et al. [2] are shown below. p.Asn1002Lys lies ten residues from L1012P. (c) Coding-sequence coordinates computed from the exon model. The premature stop at c.2210 lies 747 nt upstream of the last exon–exon junction (c.2957), well inside the region where nonsense-mediated decay is expected.](figures/fig2_bub1b_locus_protein.png)

Computing coding coordinates directly from the MANE exon model reproduces both HGVS descriptions independently: c.2210 and c.3006, in a 3,153-nt coding sequence (1,050 codons plus stop).

**p.Leu737\*** creates a premature termination codon in exon 17 of 23, 747 nt upstream of the last exon junction. Under the ClinGen PVS1 decision tree [37], the transcript is predicted to undergo nonsense-mediated decay (NMD), and loss of function is the established disease mechanism. PVS1 therefore applies at very strong strength. This matters: in overexpression experiments, a construct truncated near this position (731X) retained partial function [2]. In patients, however, truncating alleles produce no stable transcript, which is the expected NMD outcome. The allele's pathogenicity rests on the absence of protein, not on the properties of a truncated fragment.

**p.Asn1002Lys** replaces a polar asparagine with a positively charged lysine in the C-terminal lobe of the kinase domain, 48 residues before the protein's end. Every MVA kinase-domain substitution studied functionally lowers BUBR1 abundance by accelerating turnover [2], and the closest of these, L1012P, lies ten residues away. The missense therefore falls in the region where the known pathogenic mechanism of *BUB1B* missense alleles operates.

### 3.4 Both alleles are authentic germline heterozygous calls

![**Figure 3. Read-level authentication.** (a) Reads at MAPQ ≥ 60 and base quality ≥ 20. Both alternate alleles are observed on forward and reverse strands, with no soft-clipped, supplementary (SA) or alternative-hit (XA) alternate reads. (b) Alternate-allele fractions with Wilson 95% intervals; both overlap 0.5 and lie inside the heterozygous window that contains 90.2% of genome-wide heterozygous calls (n = 3,079,630). (c) Genotypes of 15 bounded *BUB1B* alleles in the supplied VCF and in our independent bwa-mem2 + bcftools recall. Ten are concordant, including both submitted alleles; all five discordant sites are indels in low-complexity sequence. (d) gnomAD allele frequencies of bounded alleles in the six MVA genes; both submitted alleles are far below the 0.01 rare-allele filter.](figures/fig3_allele_authentication.png)

| | p.Leu737\* | p.Asn1002Lys |
|---|---|---|
| GRCh38 | chr15:40209701 T>G | chr15:40220612 T>G |
| HGVS (MANE) | NM_001211.6:c.2210T>G; NP_001202.5:p.(Leu737\*) | NM_001211.6:c.3006T>G; NP_001202.5:p.(Asn1002Lys) |
| Supplied VCF | 0/1, PASS, DP 46, AD 21/25 | 0/1, PASS, DP 28, AD 15/13 |
| High-quality reads (REF/ALT) | 21/26 (ALT forward 14, reverse 12) | 16/13 (ALT forward 5, reverse 8) |
| ALT fraction (95% CI) | 0.553 (0.412–0.686) | 0.448 (0.284–0.625) |
| Independent recall | 0/1, concordant | 0/1, concordant |
| Three-caller support | DeepVariant, Strelka2 PASS 0/1 (2/3); GATK 0/1, unfiltered | DeepVariant, Strelka2 PASS 0/1 (2/3); GATK 0/1, unfiltered |
| gnomAD AF; homozygotes | 7.4×10⁻⁵; 0 | 5×10⁻⁷; 0 |
| ClinVar | VCV000533901.9, P/LP, two stars (Labcorp/Invitae, GeneDx); MVA1 condition record Pathogenic | no exact record |
| AlphaMissense; REVEL | not applicable | 0.923; 0.472 |
| AlphaGenome impact (top percentile) | 0.042% | 0.275% |
| SpliceAI max Δ | 0.03 | 0.02 |

Both variants are supported by balanced, bidirectional, high-mapping-quality reads, and both are recalled from an independent alignment of the same sample. The three-caller run used that same BAM and is therefore a technical replicate, not a biological one. Their allele fractions are inconsistent with a low-level somatic or mosaic event and consistent with constitutional heterozygosity.

The chr20 benchmark measures the precision of the calling pipeline (Section 3.10). The supplied call set and our recall disagree only at homopolymer-adjacent indels (Figure 3c), a known weakness of short-read callers that does not affect these two SNVs.

### 3.5 The missense allele is predicted damaging by a calibrated predictor

![**Figure 4. Variant-effect evidence for p.Asn1002Lys.** (a) AlphaMissense score (0.923) against the ClinGen-calibrated intervals of Bergquist et al. 2025 [29]: the score lies in the PP3_Moderate interval (0.906–0.971). (b) REVEL score (0.472) against the calibrated intervals of Pejaver et al. 2022 [30]: indeterminate. (c) Tail positive likelihood ratios at the candidate's own score, estimated on two-star ClinVar missense variants with *BUB1B* held out. Tail ratios overstate interval-specific evidence, so classification uses the calibrated intervals in (a). (d) SpliceAI maximum Δ for the two candidates and 27 further scored keys (n = 29); all are ≤ 0.03, well below the 0.2 splice-altering cut-off. The text inset gives AlphaGenome impact percentiles.](figures/fig4_variant_effect.png)

AlphaMissense, the strongest calibrated predictor available for this allele, places p.Asn1002Lys in the **PP3_Moderate** interval. REVEL is indeterminate. Because the two predictors are correlated, only the stronger calibrated signal is counted, as pre-registered. Neither variant is predicted to alter splicing (SpliceAI Δ ≤ 0.03; AlphaGenome splice scores 0.087 and 0.048, both below the local *BUB1B* 95th percentile of 1.063). The missense therefore acts at the protein level, consistent with the destabilisation mechanism above.

Asparagine 1002 is not resolved in any available experimental BUBR1 structure (PDB 2WU6, 3SI5, 4AEZ, 4GGD), so we did not compute structure-based stability estimates rather than extrapolate.

### 3.6 The phenotype points to *BUB1B*

![**Figure 5. Phenotype-only gene ranking.** (a) Match level of each of the seven HPO terms to eight candidate disorders: exact, broad, or not substantiated in the reference literature. Unreported features are not treated as absent. (b) Pre-registered weighted coverage; *BUB1B* ranks first (7.5 of 13.5 points). (c) Number of times each gene ranked first across 100 fixed-seed permutations of the HPO weights.](figures/fig5_phenotype_fit.png)

Ranking on phenotype alone, blind to genotype, places *BUB1B* first (Figure 5):

- **BUB1B (MVA1):** matches short stature exactly and matches rhabdomyosarcoma, nephrocalcinosis, failure to thrive, muscle atrophy and small-for-gestational-age growth broadly.
- **TRIP13 (MVA3):** fits the tumour term more specifically, because Wilms tumour dominates that phenotype, but covers fewer features.
- **CEP57 (MVA2):** lacks the cancer predisposition.

The ranking is robust. *BUB1B* remained first in 83 of 100 weight permutations; only *TRIP13* ever displaced it (17 of 100).

A knowledge-graph audit in PrimeKG [35] gives independent context. *BUB1B* has 30 disease neighbours, among them mosaic variegated aneuploidy syndrome, embryonal rhabdomyosarcoma, botryoid rhabdomyosarcoma, embryonal neoplasm and microcephaly. These link the child's tumour type to the gene.

### 3.7 Competing explanations were systematically excluded

![**Figure 6. Exclusion of competing genetic explanations.** (a) Read depth relative to the six-gene median for each MVA gene body (± 2 kb) and its coding exons; all fall between 0.83 and 1.11, far from the 0.5 expected for a heterozygous deletion or 1.5 for a duplication. (b) Rare (AF ≤ 0.01) coding heterozygous alleles per MVA gene: only *BUB1B* carries two. The inset summarises the genome-wide ClinVar P/LP scan. (c) Read depth at nine published pathogenic sites in *BUB1*, *MAD1L1*, *CEP57* and *CENATAC* now mapped to GRCh38; every site is homozygous reference with zero alternate reads. Four further published SNVs in *BUB1B*, *BUB1* and *TRIP13* were also ALT = 0.](figures/fig6_competing_explanations.png)

A recessive diagnosis is credible only if no better explanation exists. We tested six alternatives (Figure 6).

1. **A different MVA gene.** No rare coding allele exists in *CEP57*, *TRIP13*, *CENATAC* or *BUB1*. *MAD1L1* has a single rare heterozygous allele, which is insufficient for its biallelic disease. None of 13 published pathogenic alleles in MVA genes is present.
2. **A hidden second hit (deletion, duplication or rearrangement).** Gene-body and exon depths are diploid in all six genes. Manta reported no structural variant in any of them. Delly's only calls touching the panel are a 78 Mb chr2 inversion and a low-quality deletion with distant breakpoints, both artefactual in scale; neither involves *BUB1B*.
3. **A cryptic splice or regulatory allele.** SpliceAI Δ ≤ 0.03 for all 29 scored keys. The 28 rare noncoding keys found across the six genes have a median AlphaGenome impact percentile of 62%, which is uninformative.
4. **Homozygosity by descent or uniparental disomy.** There are no long runs of homozygosity, and allele balance is preserved on chromosome 15 (minor allele fraction 0.492 across 232 heterozygous sites).
5. **An unrelated actionable finding.** The genome-wide ClinVar scan found six heterozygous P/LP alleles (*BUB1B*, *FLG*, *GNRHR*, *LZTR1*, *PRSS1*, *RBM8A*) and no homozygous P/LP allele. None is in an ACMG secondary-findings v3.2 gene, and none explains the phenotype as a heterozygous carrier state.
6. **Technical artefact or sample mix-up.** This is excluded by read-level evidence (Section 3.4) and contamination QC (Section 3.1).

Across the whole genome, 37 genes carry two or more rare coding heterozygous alleles. Among the six MVA genes, only *BUB1B* does.

### 3.8 No bulk aneuploidy in blood, as expected for variegated mosaicism

![**Figure 7. Genome-wide dosage and allele balance.** (a) Median depth of each chromosome arm (100 kb bins) relative to the autosomal median. All 39 assembled autosomal arms lie between 0.985 and 1.049. Dashed lines mark the shift expected if 20% of cells carried a trisomy or monosomy of one chromosome. (b) Minor allele fraction modelled by GATK ModelSegments from heterozygous sites in the BAM (6,826 sites in total); dashed lines show expectations for 10% and 20% mosaic trisomy. All 22 autosomes are balanced.](figures/fig7_genome_wide_qc.png)

No chromosome shows a dosage or allele-balance shift (Figure 7): arm depth ratios run from 0.985 to 1.049 and minor allele fractions from 0.474 to 0.495.

This result does not count against MVA; it is what variegated mosaicism predicts. In MVA, different cells gain or lose different chromosomes, so in bulk DNA the aneuploid cells are spread over many chromosomes and no single chromosome reaches a detectable clone fraction. Blood WGS averages across millions of cells. It is therefore an insensitive test for MVA, and it cannot replace metaphase karyotyping or interphase FISH, which score individual cells. The useful negative result is that the genome shows no clonal, segmental or whole-chromosome event that would point to a different diagnosis, such as a constitutional trisomy or a large pathogenic CNV.

### 3.9 Phase is unresolved by short reads and is the main confirmatory step

![**Figure 8. Phase scenarios.** (a) In trans, neither haplotype encodes normal BUBR1: one is null (NMD) and the other encodes a destabilised kinase domain. This is the biallelic MVA1 genotype. In cis, one intact allele remains and the child would be a heterozygous carrier. (b) WhatsHap phased none of the ten Tier-1 candidate pairs, and the GATK caller assigned phase sets to 16% of coding heterozygous sites (1,247 of 7,825). Neither submitted allele carries a shared phase set. The 10,911 bp separation exceeds paired-end fragment lengths.](figures/fig8_phase_scenarios.png)

The two *BUB1B* alleles are 10,911 bp apart. No read pair spans both, and no heterozygous variant between them links the two sites in a phase block. WhatsHap left all ten Tier-1 candidate pairs unphased. Of 7,825 coding heterozygous sites, GATK assigned physical phase sets to 1,247, and neither submitted allele received one.

Population frequencies alone cannot decide between the two configurations: under linkage equilibrium, a cis double-heterozygote and a trans double-heterozygote are about equally frequent (both ≈ 2·p₁·p₂). The case for trans comes from conditioning on the phenotype. A cis configuration leaves one intact *BUB1B* allele, a carrier state that does not cause MVA. Yet the child has an MVA phenotype that ranks *BUB1B* first, and no alternative explanation survives (Section 3.7). The trans configuration is therefore the only one consistent with the clinical picture, and it agrees with the clinically confirmed answer.

Formal confirmation still requires one of two tests:

- **Parental targeted genotyping.** This is the cheapest. Each parent is expected to carry one allele. It would also define recurrence risk, which matters given the reported parental history of pregnancy loss.
- **Long reads.** PacBio HiFi or Oxford Nanopore reads longer than 11 kb would span both sites.

### 3.10 Pipeline accuracy and controls

![**Figure 9. Small-variant benchmark of the calling pipeline on chromosome 20.** Precision and recall of exact normalised alleles for HG002 and HG003 against GIAB truth, by stratum. Accuracy exceeds 0.99 in confident regions and is lowest in segmental duplications (recall 0.94–0.97). The benchmark covers chr20 only and does not measure locus-specific accuracy at *BUB1B*.](figures/fig9_benchmark.png)

On chromosome 20, the pipeline reached, for HG002 in confident regions:

- **all variants:** precision 0.9989, recall 0.9949;
- **SNVs:** precision 0.9992, recall 0.9956.

HG003 gave matching results (Figure 9). Neither submitted allele occurs among the 657 GIAB truth alleles of HG002 in the same six-gene windows, which is consistent with their rarity.

A planned BAMSurgeon positive-control spike-in [41] did not complete: re-mapped coverage fell below the tool's 0.9 retention threshold. No recovery ranks are reported.

## 4. Variant interpretation

![**Figure 10. ACMG/AMP evidence ledger.** Bayesian points (Tavtigian et al. 2020 [40]) for each allele on present data (top two rows) and after the named confirmatory evidence (bottom two rows). p.Leu737\* reaches Likely pathogenic on present data (PVS1 + PM2_Supporting = 9). Its ClinVar aggregate classification is Pathogenic/Likely pathogenic. p.Asn1002Lys is currently a VUS (3 points). Demonstrating trans phase (PM3) together with constitutional MVA cytogenetics (PP4) brings it to Likely pathogenic (6), and a BUBR1 abundance assay showing reduced protein (PS3_Moderate) would add a further 2 points.](figures/fig10_acmg_points.png)

| Criterion | p.Leu737\* | p.Asn1002Lys | Basis |
|---|---|---|---|
| PVS1 (very strong) | **met** | n/a | stop-gain in exon 17 of 23, 747 nt upstream of the last junction, NMD predicted; loss of function is the established MVA1 mechanism [1, 2, 37] |
| PM2_Supporting | **met** | **met** | gnomAD AF 7.4×10⁻⁵ and 5×10⁻⁷, no homozygotes; compatible with an ultra-rare recessive disease |
| PP3_Moderate | n/a | **met** | AlphaMissense 0.923 within the calibrated 0.906–0.971 interval [29]; REVEL indeterminate, not double-counted |
| PM3 | not required | **not applied** | trans configuration with a pathogenic allele cannot be assessed from short reads of the proband alone; parental or long-read data are required |
| PP4 | supportive | **not applied** | phenotype ranks *BUB1B* first, but the criterion's specificity requirement needs constitutional cytogenetic confirmation of MVA, which was not available |
| PS3 | — | **not tested** | a BUBR1 abundance or turnover assay in patient cells, as in [2] |
| ClinVar | P/LP, two stars | absent | used as context, not as independent evidence |
| **Classification** | **Likely pathogenic (9 points); ClinVar P/LP** | **VUS (3 points), favouring pathogenic; Likely pathogenic once phase and phenotype criteria are met** | |

Interpreted together, the genotype NM_001211.6:c.[2210T>G];[3006T>G], p.[(Leu737\*)];[(Asn1002Lys)], is a biallelic *BUB1B* genotype of the established MVA1 class. It accounts for the child's growth restriction, muscle and renal features and embryonal tumour.

## 5. Discussion

We report a compound heterozygous *BUB1B* genotype in a child with MVA. It was recovered from whole-genome sequencing by a pre-registered and independently validated pipeline, and it agrees with the clinically confirmed answer.

**Why this genotype fits.** Three independent lines converge:

- *Molecular class.* It reproduces the canonical MVA1 configuration: a null allele in trans with a kinase-domain substitution [1, 2].
- *Gene uniqueness.* *BUB1B* is the only one of six MVA genes in this genome carrying two rare protein-altering alleles.
- *Phenotype.* Phenotype-only analysis placed *BUB1B* first.

Because the truncating allele is predicted to undergo NMD, BUBR1 in this child is expected to derive solely from the p.Asn1002Lys allele. In the Suijkerbuijk model, kinase-domain substitutions reduce BUBR1 levels several-fold by accelerating turnover, so total BUBR1 would fall below the level required for a functional checkpoint. That model predicts chromosome missegregation in a proportion of cells, with growth restriction and embryonal tumour risk. Rhabdomyosarcoma is among the childhood cancers reported in MVA [1]. *BUB1B* is linked to embryonal and botryoid rhabdomyosarcoma in curated knowledge graphs, and the gene's broad phenotype match explains why *BUB1B* outranked *TRIP13*, whose disease is dominated by Wilms tumour [8].

**What the genome cannot show, and why that is expected.** Bulk blood WGS showed no chromosome-level dosage change. In variegated aneuploidy, the aneuploid cells carry different chromosome changes, so the signal is diluted below the detection limit of read depth and allele balance. The constitutional diagnosis of MVA rests on cell-by-cell cytogenetics, and the genetic diagnosis on the *BUB1B* genotype. Our negative screen excludes alternative clonal or structural diagnoses rather than contradicting MVA.

**Clinical implications.** We offer these for the clinical team's consideration, not as medical advice:

- A biallelic *BUB1B* diagnosis carries a well-documented risk of further embryonal tumours, which supports continued surveillance as practised for MVA [1].
- Parental carrier testing provides recurrence-risk information for the family.
- The specific missense identified here can be tested directly in patient fibroblasts or lymphoblastoid cells by BUBR1 immunoblotting with cycloheximide chase. That assay would supply the PS3 evidence and enable ClinVar submission of the variant.

**Methodological lessons.** Several procedures prevented an incorrect first answer:

- *Tie-breaking on mechanism.* The initial ranked list placed an intronic *BUB1* pair first by a lexicographic tie-break. The pre-registered requirement to name a molecular mechanism for the weaker allele removed it.
- *Annotation checks.* A contig-naming mismatch silently erased all ClinVar annotations on the first join; independent re-joining exposed it.
- *Blinded phenotype ranking.* This provided an assessment of gene fit that did not depend on the genotype evidence.

These checks are cheap and transferable to other undiagnosed cases.

## 6. Internal adversarial review

The evidence package was reviewed in an independent context by an AI red-team session that did not have access to this report's conclusions. The table records only those objections that remain material.

| Item | Status |
|---|---|
| Missense effect of p.Asn1002Lys | **Partly resolved.** Calibrated AlphaMissense 0.923 supports PP3_Moderate, and the kinase-domain destabilisation mechanism is published [2]. PS3 (a direct functional assay) remains open. |
| ePCR criteria | **Resolved.** Automatic PM2-moderate and PP4 scores were withdrawn. With PM2_Supporting only, the pair heuristic is 0.05, conditional on trans phase. The submitted CSV keeps its original ranking heuristic. ePCR is not a clinical probability. |
| Phase | **Open.** Neither allele lies in a reported WhatsHap or GATK phase block. Confirmation is a laboratory step (Section 3.9). |

Equal scores on the initial ranked list were resolved by a dated addendum to the analysis plan (Appendix A). Remaining method substitutions are listed in Section 7.

## 7. Limitations

This is a single-proband analysis without parental samples. Physical phase of the two *BUB1B* alleles is not established, and no allele-specific RNA, protein or functional assay was performed. The missense classification therefore rests on calibrated computational evidence and on mechanistic analogy with published substitutions, not on a direct measurement.

Blood WGS cannot measure constitutional mosaic aneuploidy at the level of individual cells.

The calling benchmark covers chromosome 20 only. The two alignments compared in Figure 3c, and the three-caller consensus, analyse the same specimen, so their agreement is technical rather than biological replication; correlated calling error remains possible.

Copy-number assessment is depth-based within six genes and does not replace a genome-wide CNV caller panel. Candidate retrieval was bounded to rare protein-altering and splice-region alleles plus all Tier-1 alleles, so unscored deep-intronic or mosaic events elsewhere remain theoretically possible.

Three method substitutions relative to the analysis plan are recorded here rather than as process narrative: mosaic calling used GATK ModelSegments because MoChA cannot run on a variant-only VCF; Exomiser could not reach its data sources from the compute nodes; and the HG002 negative control used the GIAB truth VCF on the staged chromosome-20 reads rather than a full-genome re-alignment.

ePCR is a ranking heuristic without clinical calibration. This report is a research analysis and not a clinical diagnostic report.

## 8. Recommended confirmatory work

1. Targeted Sanger or amplicon genotyping of both *BUB1B* sites in both parents, to confirm trans phase (PM3).
2. Long-read sequencing (PacBio HiFi or Oxford Nanopore, reads > 11 kb) across chr15:40,209,000–40,221,000 if parental DNA is unavailable.
3. BUBR1 immunoblot and cycloheximide-chase turnover in patient cells, with rescue by wild-type BUBR1, to obtain PS3 evidence [2].
4. Metaphase karyotype or interphase FISH of lymphocytes or fibroblasts with a premature chromatid separation assay, to document constitutional MVA (PP4).
5. ClinVar submission of p.Asn1002Lys, with phase and functional data, so that future families benefit.

## 9. Reproducibility

Raw patient data are controlled and are not part of this repository; reproduction requires authorised hackathon access. File hashes, figure-to-table mapping, MN5 job identifiers and rerun commands are in Appendix A. Every number in this report is read from a table under the evidence tree; figures are regenerated by the figure script and the PDF by the report builder.

## Data and code availability

The hackathon's genomic and clinical data are available only under the organisers' controlled-access terms; this repository contains no FASTQ, BAM, VCF or clinical text. Restricted and derived patient data will be deleted, with all copies, within 30 days of hackathon close, and the organisers will be notified. Aggregate and six-gene evidence tables, analysis scripts and figures are released under CC BY 4.0, as the hackathon rules require, except where a copied file carries its own upstream licence (LICENSE).

## Funding

This work was supported by the EuroHPC Joint Undertaking, which awarded access to the MareNostrum 5 supercomputer hosted by the Barcelona Supercomputing Center (BSC-CNS), Spain, under awarded project ID **EHPC-AIF-2026PG01-1286**.

## AI disclosure

All analyses were carried out on Nexgene AI's **Arf Machine** research platform, on a local workstation and on MareNostrum 5. The language-model component was **ASA-W1**, Nexgene AI's medical-reasoning model, running locally.

- No patient data were used to train or fine-tune any model, and no patient data were exported from the local environment.
- Public database queries contained only a gene symbol, a single variant notation or a publication identifier; they did not contain identifiable patient information.
- Data-use conditions were followed. Controlled source data and derived patient files will be deleted, with all copies, within 30 days of hackathon close, and the organisers will be notified.
- All AI-assisted interpretation was reviewed by the authors.

## Acknowledgements

We are deeply grateful to the child and family who chose to share their genome with the research community. We thank the MVA Society, Sage Bionetworks, Hugging Face, BEACON and the Wilhelm Foundation for organising and supporting the MVA Hackathon 2026. We acknowledge the Barcelona Supercomputing Center for operating MareNostrum 5.

We thank the developers and maintainers of the public data resources, models and software listed in LICENSE (GRCh38, gnomAD, ClinVar, MANE, HPO, AlphaMissense, REVEL, AlphaGenome, SpliceAI, PrimeKG, GIAB; bwa-mem2, samtools/bcftools, GATK, DeepVariant, Strelka2, Manta, Delly, WhatsHap, mosdepth, somalier, VerifyBamID2, snpEff, and the figure-style library). Each retains its upstream licence.

## References

1. Hanks S, Coleman K, Reid S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in *BUB1B*. *Nat Genet.* 2004;36:1159–1161. [doi:10.1038/ng1449](https://doi.org/10.1038/ng1449)
2. Suijkerbuijk SJE, van Osch MHJ, Bos FL, et al. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res.* 2010;70:4891–4900. [doi:10.1158/0008-5472.CAN-09-4319](https://doi.org/10.1158/0008-5472.CAN-09-4319)
3. Musacchio A. The molecular biology of spindle assembly checkpoint signaling dynamics. *Curr Biol.* 2015;25:R1002–R1018. [doi:10.1016/j.cub.2015.08.051](https://doi.org/10.1016/j.cub.2015.08.051)
4. Matsuura S, Matsumoto Y, Morishima K, et al. Monoallelic *BUB1B* mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. *Am J Med Genet A.* 2006;140A:358–367. [doi:10.1002/ajmg.a.31069](https://doi.org/10.1002/ajmg.a.31069)
5. Ochiai H, Miyamoto T, Kanai A, et al. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of *BUB1B* as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci USA.* 2014;111:1461–1466. [doi:10.1073/pnas.1317008111](https://doi.org/10.1073/pnas.1317008111)
6. Rio Frio T, Lavoie J, Hamel N, et al. Homozygous *BUB1B* mutation and susceptibility to gastrointestinal neoplasia. *N Engl J Med.* 2010;363:2628–2637. [doi:10.1056/NEJMoa1006565](https://doi.org/10.1056/NEJMoa1006565)
7. Snape K, Hanks S, Ruark E, et al. Mutations in *CEP57* cause mosaic variegated aneuploidy syndrome. *Nat Genet.* 2011;43:527–529. [doi:10.1038/ng.822](https://doi.org/10.1038/ng.822)
8. Yost S, de Wolf B, Hanks S, et al. Biallelic *TRIP13* mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet.* 2017;49:1148–1151. [doi:10.1038/ng.3883](https://doi.org/10.1038/ng.3883)
9. de Wolf B, Oghabian A, Akinyi MV, et al. Chromosomal instability by mutations in the novel minor spliceosome component CENATAC. *EMBO J.* 2021;40:e106536. [doi:10.15252/embj.2020106536](https://doi.org/10.15252/embj.2020106536)
10. Villarroya-Beltri C, Osorio A, Torres-Ruiz R, et al. Biallelic germline mutations in *MAD1L1* induce a syndrome of aneuploidy with high tumor susceptibility. *Sci Adv.* 2022;8:eabq5914. [doi:10.1126/sciadv.abq5914](https://doi.org/10.1126/sciadv.abq5914)
11. Carvalhal S, Bader I, Rooimans MA, et al. Biallelic *BUB1* mutations cause microcephaly, developmental delay, and variable effects on cohesion and chromosome segregation. *Sci Adv.* 2022;8:eabk0114. [doi:10.1126/sciadv.abk0114](https://doi.org/10.1126/sciadv.abk0114)
12. Sieben CJ, Jeganathan KB, Nelson GG, et al. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *J Clin Invest.* 2020;130:171–188. [doi:10.1172/JCI126863](https://doi.org/10.1172/JCI126863)
13. Köhler S, Gargano M, Matentzoglu N, et al. The Human Phenotype Ontology in 2021. *Nucleic Acids Res.* 2021;49:D1207–D1217. [doi:10.1093/nar/gkaa1043](https://doi.org/10.1093/nar/gkaa1043)
14. Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *IEEE IPDPS.* 2019:314–324. [doi:10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041)
15. Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10:giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
16. McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res.* 2010;20:1297–1303. [doi:10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110)
17. Cingolani P, Platts A, Wang LL, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly.* 2012;6:80–92. [doi:10.4161/fly.19695](https://doi.org/10.4161/fly.19695)
18. Morales J, Pujar S, Loveland JE, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature.* 2022;604:310–315. [doi:10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8)
19. Chen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. *Nature.* 2024;625:92–100. [doi:10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0)
20. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res.* 2018;46:D1062–D1067. [doi:10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153)
21. Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science.* 2023;381:eadg7492. [doi:10.1126/science.eadg7492](https://doi.org/10.1126/science.eadg7492)
22. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet.* 2016;99:877–885. [doi:10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016)
23. Avsec Ž, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature.* 2026;649:1206–1218. [doi:10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)
24. Patterson M, Marschall T, Pisanti N, et al. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. *J Comput Biol.* 2015;22:498–509. [doi:10.1089/cmb.2014.0157](https://doi.org/10.1089/cmb.2014.0157)
25. Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. *Bioinformatics.* 2018;34:867–868. [doi:10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699)
26. Chen X, Schulz-Trieglaff O, Shaw R, et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics.* 2016;32:1220–1222. [doi:10.1093/bioinformatics/btv710](https://doi.org/10.1093/bioinformatics/btv710)
27. Rausch T, Zichner T, Schlattl A, et al. DELLY: structural variant discovery by integrated paired-end and split-read analysis. *Bioinformatics.* 2012;28:i333–i339. [doi:10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)
28. Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting splicing from primary sequence with deep learning. *Cell.* 2019;176:535–548. [doi:10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015)
29. Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. *Genet Med.* 2025;27:101402. [doi:10.1016/j.gim.2025.101402](https://doi.org/10.1016/j.gim.2025.101402)
30. Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet.* 2022;109:2163–2177. [doi:10.1016/j.ajhg.2022.10.013](https://doi.org/10.1016/j.ajhg.2022.10.013)
31. Zhang F, Flickinger M, Taliun SAG, et al. Ancestry-agnostic estimation of DNA sample contamination from sequence reads. *Genome Res.* 2020;30:185–194. [doi:10.1101/gr.246934.118](https://doi.org/10.1101/gr.246934.118)
32. Pedersen BS, Bhetariya PJ, Brown J, et al. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. *Genome Med.* 2020;12:62. [doi:10.1186/s13073-020-00761-2](https://doi.org/10.1186/s13073-020-00761-2)
33. Narasimhan V, Danecek P, Scally A, et al. BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data. *Bioinformatics.* 2016;32:1749–1751. [doi:10.1093/bioinformatics/btw044](https://doi.org/10.1093/bioinformatics/btw044)
34. Wagner J, Olson ND, Harris L, et al. Benchmarking challenging small variants with linked and long reads. *Cell Genomics.* 2022;2:100128. [doi:10.1016/j.xgen.2022.100128](https://doi.org/10.1016/j.xgen.2022.100128)
35. Chandak P, Huang K, Zitnik M. Building a knowledge graph to enable precision medicine. *Sci Data.* 2023;10:67. [doi:10.1038/s41597-023-01960-3](https://doi.org/10.1038/s41597-023-01960-3)
36. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med.* 2015;17:405–424. [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
37. Abou Tayoun AN, Pesaran T, DiStefano MT, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat.* 2018;39:1517–1524. [doi:10.1002/humu.23626](https://doi.org/10.1002/humu.23626)
38. Walker LC, de la Hoya M, Wiggins GAR, et al. Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing: recommendations from the ClinGen SVI Splicing Subgroup. *Am J Hum Genet.* 2023;110:1046–1067. [doi:10.1016/j.ajhg.2023.06.002](https://doi.org/10.1016/j.ajhg.2023.06.002)
39. Tavtigian SV, Greenblatt MS, Harrison SM, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med.* 2018;20:1054–1060. [doi:10.1038/gim.2017.210](https://doi.org/10.1038/gim.2017.210)
40. Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat.* 2020;41:1734–1737. [doi:10.1002/humu.24088](https://doi.org/10.1002/humu.24088)
41. Ewing AD, Houlahan KE, Hu Y, et al. Combining tumor genome simulation with crowdsourcing to benchmark somatic single-nucleotide-variant detection. *Nat Methods.* 2015;12:623–630. [doi:10.1038/nmeth.3407](https://doi.org/10.1038/nmeth.3407)
42. Poplin R, Chang PC, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nat Biotechnol.* 2018;36:983–987. [doi:10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235)
43. Kim S, Scheffler K, Halpern AL, et al. Strelka2: fast and accurate calling of germline and somatic variants. *Nat Methods.* 2018;15:591–594. [doi:10.1038/s41592-018-0051-x](https://doi.org/10.1038/s41592-018-0051-x)

## Appendix A. Data, code and provenance

### A1. Files, SHA-256 and role

SHA-256 values for scripts and figures are refreshed in `MANIFEST.sha256` after each rebuild. The submitted prediction file is unchanged.

| File | SHA-256 | Role |
|---|---|---|
| `submission/Nexgene-AI_mn5-evidence-v2.csv` | `24162afea8ce236d41261d0f10595aef864b5f708c3201166d55792a9bbbebba` | submitted prediction |
| `docs/preregistration.md` | `32f59fedddab1e81c89ea4329a24503ed4f9746d43d342fa79c19bdb7daefe00` | pre-candidate analysis plan |
| `docs/v2_adjudication_protocol.md` | `259fc684603c1b52882e1c749ee3d9d2fa792e2e2fc82abc851bf354dd6d22eb` | dated tie-break addendum |
| GRCh38 no-alt FASTA | `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7` | reference (not distributed) |
| `results/evidence/alleles/allele_validation.tsv` | `0547ab94cdf0701ed3fa7890e6db9eb74f4b68b0e9ee5a24b231dc93e5eef892` | read-level validation (MN5 job 47403447) |
| `results/evidence/alleles/three_caller_consensus.tsv` | `502eee787ec9224d0473209964664b022cd8f0080b9bbcc566eae0e704d19bcf` | DeepVariant, Strelka2 and GATK genotypes for the two alleles |
| `results/evidence/phase/phase_pair_evidence.tsv` | `45c5bc955ea570f716ee26c1e448e4965d15c37200574cf5651d15ddb5b30c49` | WhatsHap pair outcomes |
| `results/evidence/second_hit/second_hit_screen.tsv` | `ab005f6b113a09886c10f2bbed2a7a48fc803833e7e6001511e659f896231228` | six-gene second-hit screen |
| `results/evidence/structural/sixgene_cnv.tsv` | `2a827ae62371f25469d4375a48f4c7b2e491b46e8362dd752adf866bc66af4ab` | gene-body depth (MN5 job 47438935) |
| `results/evidence/aneuploidy/mosaic_segments.tsv` | `8c32cb1d33f9f79197f1d27595ea9b07a2180793fd0d379dcc4fc0b3c865c12d` | ModelSegments allele balance (MN5 job 47439658) |
| `results/evidence/second_hit/spliceai_keys.tsv` | `edbec3d2fc8e1b2e6e7e13a7600d9f933cce8c0131372c9758314e5a61b17d49` | SpliceAI (MN5 job 47440407) |
| `results/evidence/alleles/predictor_lr.tsv` | `09cf5bf8af2d245a4de1be81bf101619596380a5bf87384583a5d9db3253265b` | held-out predictor likelihood ratios |
| `results/evidence/literature/angle_b_phenotype.tsv` | `e0669b57d4cad2afbf7e3418ebfb5d866ff2faaeac4ca7d308aa1a6980d5dfa3` | phenotype-only ranking |
| `results/evidence/ranking/funnel.json` | `5228ac5eefc8559aa5989c0decabf6c07fe8f73eeb513b020cf6da8819548d35` | candidate funnel |
| `results/evidence/controls/l7_metrics_audited.tsv` | `c2dbbb16f12faca0c7455c4ea866626385c541ecd4f2ed17f227a4cbb5b0029d` | chr20 benchmark |

### A2. Figure and table mapping

| Item | Evidence file |
|---|---|
| Graphical abstract; Figure 1 | `results/evidence/ranking/funnel.json` |
| Figure 2 | `results/evidence/structural/tier1_sixgene_exon_utr.bed`; `submission/figures/locus_facts.json` |
| Figure 3; allele table in §3.4 | `results/evidence/alleles/allele_validation.tsv`; `caller_concordance.tsv`; `three_caller_consensus.tsv` |
| Figure 4 | `results/evidence/alleles/predictor_lr.tsv`; `coding_avi_splice.tsv`; `second_hit/spliceai_keys.tsv` |
| Figure 5 | `results/evidence/literature/angle_b_phenotype.tsv`; `controls/null_runs.tsv` |
| Figure 6 | `results/evidence/structural/`; `genome_scan/`; `literature/published_site_genotypes.tsv` |
| Figure 7 | `results/evidence/aneuploidy/`; `qc/verifybamid2_summary.tsv` |
| Figure 8 | `results/evidence/phase/` |
| Figure 9 | `results/evidence/controls/l7_metrics_audited.tsv` |
| Figure 10 | ACMG point scale in Section 4 |
| ePCR criterion audit | `results/evidence/ranking/epcr_criterion_audit.tsv` |

### A3. MareNostrum 5 jobs

| Stage | Job ID | Notes |
|---|---|---|
| Alignment (four FASTQ lanes) | 47349918_1–4 | bwa-mem2 against GRCh38 no-alt |
| Merge and duplicate marking | 47353088 | Picard MarkDuplicates |
| gnomAD fastpath | 47392756 | supplied call set |
| Critical annotation | 47397018 | snpEff, ClinVar, predictors |
| Allele, phase and second-hit screen | 47403447 | bounded six-gene validation |
| Independent reconstruction of the initial ten-row list | 47406904 | byte-for-byte match of that list |
| Gene-body depth | 47438935 | mosdepth |
| ModelSegments allele balance | 47439658 | mosaic screen |
| SpliceAI | 47440407 | candidate and noncoding keys |
| Three-caller consensus, chr15 | 47930984 | independent DeepVariant, GATK HaplotypeCaller and Strelka2 merge; both submitted alleles CONSENSUS 2/3 |

Set `MVA_BASE` to the MN5 work root before rerunning batch jobs. Sample identifiers in the BAM remain the organisers' alias; the submitted CSV uses `PROBAND01`.

### A4. Reproduction commands

Map authorised inputs in `config/paths.example.yaml` and set `MVA_BASE`. Then, from the repository root, with the project interpreter:

```bash
python scripts/validate_submission.py submission/Nexgene-AI_mn5-evidence-v2.csv <authorised>/GRCh38_no_alt.fa
python scripts/36_release_checks.py submission/Nexgene-AI_mn5-evidence-v2.csv submission/Nexgene-AI_track1_report.md submission/Nexgene-AI_methods_description.xlsx
python scripts/60_figures.py
python scripts/61_build_pdf.py
```

MN5 alignment, annotation, scoring, pairing, read screening and the structural, splicing and QC screens are the batch scripts under `scripts/mn5/`, `scripts/supplement/`, `scripts/wave3/` and `scripts/rigor/`. Operational deviations from the analysis plan are archived in `docs/deviations.md`.
