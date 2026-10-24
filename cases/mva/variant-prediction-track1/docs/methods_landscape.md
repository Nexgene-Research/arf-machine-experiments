# Bioinformatics tools, databases, and AI models for variant detection and prioritization in rare and ultra-rare disease

**Literature synthesis (English)**  
Search window closed **17 October 2026**.  
This file is an evidence synthesis. It does not claim that every record matching the queries was read. A “who used what” row appears only when a paper or an official pipeline description names both the tool and the cohort. Preprints are labelled. Absence of a citation means the fact was not found in the sources below, not that a laboratory does not use the method.

---

## 1. Question and definitions

Rare-disease genomics fails when three jobs are collapsed into one score.

| Layer | Question | Typical output |
| ----- | -------- | -------------- |
| **Detection** | Is this allele present in the reads? | VCF/gVCF of SNVs, indels, SVs, CNVs, STRs, mosaic sites, mtDNA |
| **Prioritization** | Given this phenotype and pedigree, which of the millions of alleles is the cause? | Ranked gene/variant list (Exomiser, AI-MARRVEL, LIRICAL, …) |
| **Classification** | Does the evidence meet a clinical standard for this variant–disease pair? | ACMG/AMP five-class call (Richards et al., 2015) |

A caller with the best GIAB F1 can still miss a private clinical allele. A ranker trained on solved genes can bury an n-of-1 gene. A pathogenic ClinVar star does not prove that the allele explains *this* patient (Mao et al., 2024).

**Rare disease.** EU / Orphanet convention: a condition affecting no more than 1 in 2,000 people (EU Orphan Medicinal Products Regulation; Nguengang Wakap et al., 2020). Some reviews use 1 in 2,500 (Ferreira, 2019, as cited by Del Gobbo and Boycott, 2025).

**Ultra-rare.** NICE and several European usages treat diseases affecting fewer than about 1 in 50,000 as ultra-rare. Orphanet epidemiology shows that 84.5% of annotated rare diseases have point prevalence below 1 in 1,000,000 (Nguengang Wakap et al., 2020). In this report, ultra-rare also covers n-of-1 presentations and “unsolvable” syndromes with no known locus (Steyaert et al., 2025).

**Scope.** Germline (and constitutional mosaic) methods used on rare or ultra-rare Mendelian disease: short-variant calling, SV/CNV calling, repeat expansions, mosaic and mitochondrial calling, RNA-outlier detection, annotation databases, variant-effect predictors, and case-level AI rankers. Somatic oncology callers appear only where a rare-disease methods paper reused them for mosaicism (Barbitoff et al., 2024, citing Mutect2 in mosaic best practices).

**Excluded.** Single-patient case reports that do not describe a method; wet-lab protocol recipes; pathogen or toxin work; HGMD record contents (licensed; named only as a resource).

---

## 2. How the corpus was assembled

### 2.1 Sources and date

Queries were run on **17 October 2026** against:

- Europe PMC REST (`/search`, JSON `hitCount`)
- NCBI PubMed E-utilities (`esearch`, `rettype=count`)

Deep reading covered systematic reviews, GIAB-style benchmarks, consortium pipeline papers, and primary tool papers. PubMed counts are closer to a biomedical title/abstract filter. Europe PMC counts are larger and are inflated when a query term is a common English word or a short acronym (DROP, FRASER, DDD, n-of-1). Tightened follow-up queries (Q2b, Q7b, Q10b) are the ones used for interpretation.

### 2.2 Fixed queries and hit counts

| ID | Query | Europe PMC | PubMed | Note |
| -- | ----- | ---------- | ------ | ---- |
| Q1 | `("rare disease" OR "rare diseases" OR Mendelian OR "undiagnosed disease") AND ("variant calling" OR "variant detection" OR "variant discovery")` | 9,136 | 255 | Core detection literature |
| Q2 | `("ultra-rare" OR ultrarare OR "n-of-1") AND (variant OR mutation) AND (genome OR exome OR sequencing)` | 48,608 | 962 | Inflated by “n-of-1” |
| Q2b | `("ultra-rare disease" OR "ultra rare disease" OR ultrarare OR "ultra-rare disorders") AND (variant OR mutation) AND (exome OR genome OR sequencing)` | 1,487 | 889 | Preferred ultra-rare count |
| Q3 | `(Exomiser OR Genomiser OR "AI-MARRVEL" OR DeepVariant OR "HaplotypeCaller") AND ("rare disease" OR Mendelian)` | 2,222 | 42 | Named callers/rankers |
| Q4 | `(gnomAD OR ClinVar OR PanelApp OR DECIPHER) AND ("rare disease") AND (interpretation OR prioritization OR annotation)` | 3,884 | 80 | Database use |
| Q5 | `("long-read" OR HiFi OR nanopore) AND ("rare disease") AND ("structural variant" OR SV OR "copy number")` | 799 | 36 | LRS SV |
| Q6 | `(SpliceAI OR AlphaMissense OR REVEL OR CADD) AND ("rare disease" OR Mendelian) AND (prioritization OR pathogenicity)` | 4,094 | 63 | VEPs |
| Q7 | `(DROP OR OUTRIDER OR FRASER) AND ("rare disease" OR Mendelian OR "RNA-seq")` | 40,774 | 535 | Inflated by DROP/FRASER |
| Q7b | `(OUTRIDER OR "FRASER2" OR "aberrant splicing" OR "RNA outlier") AND ("rare disease" OR Mendelian OR undiagnosed)` | 2,083 | 126 | Preferred RNA-outlier count |
| Q7c | `("Detection of RNA Outliers" OR "DROP pipeline") AND (RNA)` | 53 | 267 | DROP by name |
| Q8 | `(ExpansionHunter OR TRGT OR "repeat expansion") AND ("rare disease" OR Mendelian)` | 2,217 | 123 | Repeats |
| Q9 | `("100000 Genomes" OR "100,000 Genomes" OR "Solve-RD" OR "Undiagnosed Diseases Network" OR "Deciphering Developmental Disorders") AND (pipeline OR "variant calling" OR Exomiser)` | 1,412 | 78 | Consortia |
| Q10 | `(mosaic OR mosaicism OR mtDNA OR mitochondrial) AND ("rare disease" OR Mendelian) AND ("variant calling" OR detection)` | 21,235 | 308 | Inflated |
| Q10b | `(mosaic OR mosaicism) AND ("rare disease" OR Mendelian) AND ("variant calling" OR "mosaic variant")` | 1,164 | 17 | Preferred mosaic count |
| Q11 | `("Genome Aggregation Database" OR gnomAD) AND ("v4" OR "version 4")` | 1,715 | 87 | gnomAD v4 |
| Q12 | `("100,000 Genomes Project" OR "100000 Genomes Project") AND (rare) AND (diagnostic OR yield OR pipeline)` | 963 | 10 | 100kGP |
| Q13 | `("Deciphering Developmental Disorders" OR DDD study) AND (exome OR diagnostic)` | 2,787,318 | 119 | Europe PMC “DDD” explosion; use PubMed |
| Q14 | `("Undiagnosed Diseases Network") AND (exome OR genome) AND (diagnostic)` | 447 | 224 | UDN |

**What the counts mean.** Q1 + Q2b + Q5 + Q7b + Q8 already exceed several thousand records. The field is larger than any one review. This synthesis therefore reports *methods that papers actually name*, not an exhaustive paper-by-paper census.

### 2.3 Inclusion for deep reading

A paper entered the tool/cohort tables if it (i) described a method, benchmark, or production pipeline, and (ii) was used for rare/ultra-rare or undiagnosed Mendelian disease, or is the primary method paper of a tool that a rare-disease pipeline later names. Case reports that only mention “exome sequencing” without a caller were not used for adoption claims.

### 2.4 Anchor reviews and pipeline papers

| Paper | Why it is an anchor |
| ----- | ------------------- |
| Barbitoff et al., *Brief. Bioinform.* 2024 ([doi:10.1093/bib/bbad508](https://doi.org/10.1093/bib/bbad508), PMC10810331) | Short-read germline discovery for rare-disease diagnostics: aligners, callers, SV/CNV, mosaic, annotation, remaining challenges |
| Barbitoff et al., *BMC Genomics* 2022 ([doi:10.1186/s12864-022-08365-3](https://doi.org/10.1186/s12864-022-08365-3)) | 14 GIAB WES/WGS datasets × 4 aligners × 9 caller/filter combinations |
| Cooperstein et al., *Genome Med.* 2025 ([doi:10.1186/s13073-025-01546-1](https://doi.org/10.1186/s13073-025-01546-1), PMC12539062) | Exomiser/Genomiser parameter optimization on 386 diagnosed UDN probands |
| Mao et al., *NEJM AI* 2024 ([doi:10.1056/aioa2300009](https://doi.org/10.1056/aioa2300009), PMC11221788) | AI-MARRVEL: random forest trained on diagnosed cases; benchmarks vs Exomiser, LIRICAL, PhenIX, Xrare |
| Steyaert et al., *Genome Res.* 2025 ([doi:10.1101/gr.279414.124](https://doi.org/10.1101/gr.279414.124), PMID 40138663) | Solve-RD 10× PacBio HiFi in 114 undiagnosed families; names the full caller stack |
| Del Gobbo and Boycott, *Genome Res.* 2025 (PMC12047273) | Added diagnostic yield of long-read sequencing after negative short-read testing |
| Smedley et al., *N. Engl. J. Med.* 2021 ([doi:10.1056/NEJMoa2035790](https://doi.org/10.1056/NEJMoa2035790), PMID 34758253) | 100,000 Genomes rare-disease pilot: panels, Exomiser, yield |
| Dutta et al., *Eur. J. Hum. Genet.* / medRxiv 2025 ([doi:10.1101/2025.08.19.25333674](https://doi.org/10.1101/2025.08.19.25333674); journal [doi:10.1038/s41431-026-02210-x](https://doi.org/10.1038/s41431-026-02210-x)) | Scottish Genomes Partnership ONT SV calling after 100kGP short-read analysis |
| Yépez et al., *Genome Med.* 2022 ([doi:10.1186/s13073-022-01019-9](https://doi.org/10.1186/s13073-022-01019-9)) | DROP in 303 fibroblast RNA-seq cases after WES |
| Cohen et al., *Genet. Med.* 2022 ([doi:10.1016/j.gim.2022.02.007](https://doi.org/10.1016/j.gim.2022.02.007)) | Genomic Answers for Kids: ES, srGS, HiFi-GS callers |
| Deciphering Developmental Disorders Study, *Nature* 2017 ([doi:10.1038/nature21062](https://doi.org/10.1038/nature21062), PMC6016744) | DDD de novo architecture: BWA, GATK, DeNovoGear |

---

## 3. Pipeline map

```text
FASTQ
  -> QC (FastQC / MultiQC) and optional adapter trim
  -> align (BWA-MEM, Isaac/iSAAC, pbmm2, minimap2, …)
  -> mark duplicates / (optional) BQSR
  -> call
       SNVs/indels   |  SVs/CNVs  |  STRs  |  mosaic/mtDNA
  -> annotate (VEP/SnpEff/ANNOVAR + gnomAD, ClinVar, constraint, splice scores)
  -> phenotype-aware rank (Exomiser, Genomiser, LIRICAL, AI-MARRVEL, …)
  -> optional RNA outliers (DROP: OUTRIDER, FRASER, MAE)
  -> ACMG/AMP review, orthogonal validation, matchmaking
```

```mermaid
flowchart LR
  reads[Aligned reads] --> callers[Callers]
  callers --> annot[Databases and annotation]
  annot --> ai[VEPs and case-level models]
  ai --> review[Manual and ACMG review]
```

Barbitoff et al. (2024) separate *identification of all variants present* from *interpretation that selects the causal allele*. Wright, FitzPatrick and Firth (2018) review diagnostic yields of 8–70% depending on disease group and inclusion; Barbitoff et al. (2024) cite a mean trio-WGS diagnosis rate of 42% across disease groups.

---

## 4. Short-variant callers

Two algorithm classes dominate germline SNV/indel discovery (Barbitoff et al., 2024):

1. **Pileup-based.** Allele counts (and, in deep-learning callers, pileup images or tensors) at each site. SAMtools/BCFtools and VarScan are the classical form. DeepVariant and Clair3 are the high-accuracy form.
2. **Haplotype-based.** Local assembly, realignment, and genotyping. GATK HaplotypeCaller, FreeBayes, Strelka2, Octopus.

Joint (cohort) calling reports homozygous-reference genotypes and raises power. GATK natively joint-genotypes large cohorts. GLnexus joint-genotypes GATK and DeepVariant GVCFs (Yun et al., 2021; used by Solve-RD HiFi, Steyaert et al., 2025).

### 4.1 GIAB-style accuracy is not clinical coverage

Barbitoff et al. (2022) ran 45 aligner × caller × filter combinations on 14 GIAB WES/WGS datasets (HG001 and the Ashkenazi and Han trios). Findings that matter for rare disease:

- **The caller dominates the aligner**, except that Bowtie2 should not be used for medical calling.
- **DeepVariant** had the best F1 for SNPs and indels, for every aligner they tested.
- Clair3, Octopus (standard filters), and Strelka2 were the closest runners-up.
- FreeBayes lagged on indels (F1 at least 5.7% lower than other methods).
- GATK HaplotypeCaller with **2D CNN filtering** and Octopus **random-forest filtering** over-filtered WES (high coverage), not WGS.
- Differences were mostly **recall**, not precision. Filtering can remove true diagnostic alleles; Barbitoff et al. (2024) advise keeping filtered calls visible for reanalysis.
- Concordance on six non-GIAB African and Russian samples did not show obvious GIAB overfitting for DeepVariant, Strelka2, or Octopus with standard filters.

The same group (Barbitoff et al., 2024) then quantified the clinical hole in GIAB high-confidence regions:

- Shared GIAB v4.2 high-confidence sequence: **2.37 Gbp (75.2%)** of the genome and **30.4 Mbp (86.6%)** of coding sequence.
- **4.7 Mbp of CDS** remain hard to call (sex chromosomes, MHC, segmental duplications, SVs). Chromosome X accounts for 1.3 of those 4.7 Mbp.
- Unique calls and unique misses are **up to 100-fold more common** outside high-confidence regions.
- **17,607 (12.9%)** ClinVar pathogenic variants fall outside high-confidence regions (more than expected from interval length, driven in part by MHC).
- **10,604** ClinVar pathogenic variants (v. 2023-03-26) sit on X or Y (**7.8%** of pathogenic variants vs 4.0% expected).
- An individual exome has **more than 4,000** variants outside high-confidence regions (HG001 estimate).
- GRCh37/38 still carry **reference minor alleles**, including known disease alleles (example: *SCARF2* rs587606297 / Van den Ende–Gupta). Coding RMA sites: 12,709 (GRCh37), 11,897 (GRCh38), 12,587 (T2T), with only 3,003 shared between GRCh38 and T2T.

PrecisionFDA Truth Challenge pipelines dropped F1 by up to ~10% in difficult-to-map regions (Olson et al., cited by Barbitoff et al., 2024). A caller that is best on SNVs in high-confidence CDS is therefore not automatically best on private, low-frequency, sex-chromosome, or paralog alleles.

### 4.2 Tool catalog (short variants)

| Tool | Class | Primary paper | Rare-disease use named in sources |
| ---- | ----- | ------------- | --------------------------------- |
| **GATK HaplotypeCaller** | Local assembly + pair-HMM genotyping; optional CNN filter | McKenna et al., 2010; DePristo et al., 2011; Poplin et al., 2018 (GATK4 scaling) | “Most widely used germline caller” (Barbitoff et al., 2024). DDD: GATK 3.1.1 BAM improvement + HaplotypeCaller 3.2.2 multisample + VQSR (DDD Study, 2017). SeqOne GermlineVar includes GATK (Ruzicka et al., 2025, preprint). Mosaic discovery option: HC with ploidy 50 (Wang et al., cited by Barbitoff et al., 2024) |
| **DeepVariant** | CNN on pileup images | Poplin et al., 2018, *Nat. Biotechnol.* ([doi:10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235)) | Best GIAB CDS F1 in Barbitoff et al., 2022. Solve-RD HiFi: DeepVariant 1.1.0 (Steyaert et al., 2025). GA4K HiFi: DeepVariant v1.0 PacBio model (Cohen et al., 2022). PEPPER-Margin-DeepVariant for ONT/HiFi (Shafin et al., 2021) |
| **Strelka2** | Haplotype caller with a paired-sample model; germline mode | Kim et al., 2018, *Nat. Methods* 15:591–594 ([doi:10.1038/s41592-018-0051-x](https://doi.org/10.1038/s41592-018-0051-x)) | Strong GIAB runner-up (Barbitoff et al., 2022). Fast relative to DeepVariant/Clair3 in that benchmark |
| **FreeBayes** | Bayesian haplotype caller | Garrison and Marth, 2012 | Weaker indel F1 (Barbitoff et al., 2022). Still used in some clinical stacks: SeqOne GermlineVar runs FreeBayes with AF ≥10% on exomes and ≥5% on panels (Ruzicka et al., 2025, preprint) |
| **Octopus** | Haplotype caller; optional random-forest filter | Cooke et al., 2021 | Competitive with Strelka2/Clair3 on default filters; RF filter not recommended on WES (Barbitoff et al., 2022) |
| **Clair3** | Pileup + full-alignment deep-learning; designed for long reads, also tested on short reads | Zheng et al., 2022 (HKU-BAL) | GIAB short-read runner-up, but adapter contamination and very high WES coverage can hurt pileup models (Barbitoff et al., 2024). Production use in this corpus is long-read |
| **Sentieon** | Accelerated GATK-compatible implementation | Freed et al. / Sentieon technical reports; UDN methods cite Sentieon | UDN GS/ES joint calling of SNVs/indels with Sentieon after CGAP per-sample GVCFs (Cooperstein et al., 2025) |
| **DRAGEN** | FPGA haplotype caller (Illumina) | Vendor evaluations; GIAB/PrecisionFDA comparisons exist | **Not found** as a named production caller in the consortium method papers read for this review (100kGP pilot and SGP SRS used iSAAC; DDD used GATK; Solve-RD HiFi used DeepVariant). Listed here because clinical laboratories use it; parameters were not extracted from a rare-disease methods paper in this search |
| **Isaac / iSAAC** | Illumina aligner (and historically a caller) | Raczy et al., 2013 | 100kGP / SGP SRS: 150 bp HiSeqX reads aligned with iSAAC (Dutta et al., 2025, citing the GEL pipeline) |
| **BWA-MEM** | Burrows–Wheeler aligner | Li, 2013 | Barbitoff “gold standard” short-read aligner. DDD: BWA 0.59 backtrack on hs37d5 (DDD Study, 2017). Yépez et al., 2022: BWA 0.7.5a to hg19 for diagnostic WES |
| **pbmm2** | minimap2 wrapper for PacBio | Li, 2018/2021; PacBio workflow | Solve-RD: pbmm2 1.4.0 to GRCh38 (Steyaert et al., 2025). GA4K: pbmm2 v1.4.0 (Cohen et al., 2022) |
| **minimap2** | Long-read (and short-read) aligner | Li, 2018 | SGP ONT: minimap2 v2.23 `-x map-ont` (Dutta et al., 2025) |
| **PEPPER-Margin-DeepVariant** | Haplotype-aware long-read small-variant calling | Shafin et al., 2021, *Nat. Methods* | ONT SNV calling that can outperform short-read SNVs in segmental duplications; phases 85–92% of genes. Used in Stanford ultrarapid nanopore ICU genomes (Gorzynski et al., 2022, *N. Engl. J. Med.*) |
| **GLnexus** | Joint genotyping of GVCFs | Yun et al., 2021 | Solve-RD trios: GLnexus 1.3.1 on DeepVariant (Steyaert et al., 2025) |

European WGS recommendations cited by Barbitoff et al. (2024) ask laboratories to validate pipelines on monozygotic twins and previously solved cases, not only GIAB (Weiss et al. / European recommendations, 2021, as cited).

**Aligner note.** Barbitoff et al. (2022) found Isaac and Novoalign slightly ahead of BWA-MEM on F1, with Bowtie2 clearly behind. Rare-disease production papers in this corpus still name BWA or iSAAC for short reads and pbmm2/minimap2 for long reads.

---

## 5. Structural variants, copy number, repeats, mosaic, and mtDNA

### 5.1 Why this layer matters for ultra-rare disease

Short-read WGS diagnostic yield in large rare-disease programmes is typically **25–41%** (Dutta et al., 2025, citing 100kGP and related studies). Most of those diagnoses are coding SNVs/indels. SVs (≥50 bp), CNVs, inversions, repeat expansions, and mosaic alleles are exactly the classes that remain after a “negative exome.” Del Gobbo and Boycott (2025) put the added yield of long-read sequencing after negative short-read genome sequencing at **7–17%** in emerging ~100-patient cohorts, driven mainly by SV resolution, STR characterisation, and phasing.

### 5.2 Short-read SV and CNV callers

Barbitoff et al. (2024) name DELLY, Manta, and Lumpy as the most cited short-read SV callers (Google Scholar), and Manta as a component of GATK-SV, which built gnomAD-SV. They recommend callers that combine discordant pairs, split reads, depth, and local assembly; GRIDSS and Manta are strong candidates, with GRIDSS breakend notation harder to interpret. Ensemble merging (Parliament2) is not clearly better than one state-of-the-art caller (citation to Kosugi / Mahmoud reviews). CNV calling is mostly read-depth; CNVkit and GATK gCNV are the commonly named examples. Panel CNV tools (DECoN, panelcn.MOPS) outperform WES/WGS CNV callers in the benchmarks they cite.

| Tool | What it does | Rare-disease use named in sources |
| ---- | ------------ | --------------------------------- |
| **Manta** | Assembly-based SV caller (Illumina) | 100kGP / SGP SRS SV calling (Dutta et al., 2025). SVRare aggregated Manta+Canvas calls from 71,408 100kGP rare-disease participants (Yu et al., 2021, preprint) |
| **Canvas** | Read-depth CNV caller (Illumina) | Same 100kGP/SGP stack |
| **DELLY** | Paired-end + split-read SV caller | Named as widely used (Barbitoff et al., 2024); no consortium version pinned in the papers read here |
| **Lumpy** | Probabilistic SV evidence integration | Same |
| **GRIDSS** | Assembly + breakend SV caller | SeqOne GermlineVar (Ruzicka et al., 2025, preprint) |
| **CNVkit** | Hybrid capture / WGS depth CNV | Named in Barbitoff et al., 2024 |
| **GATK gCNV** | Cohort read-depth CNV | Named in Barbitoff et al., 2024 |
| **ExomeDepth** | WES/panel CNV vs reference set | Cacheiro et al., 2017: 30-gene HSP panel, 83 patients; ExomeDepth 5/5 multi-exon deletions vs panelcn.MOPS 4/5 and CNVPanelizer 3/5 |
| **ClinSV** | Clinical-grade short-read SV | Barbitoff et al., 2024 cite ClinSV as claiming clinical-grade short-read SV calling (Minoche et al., 2021) |
| **GATK-SV** | Multi-algorithm integration | Used to construct gnomAD-SV (Collins et al. / Karczewski resources; Barbitoff et al., 2024) |
| **SVRare** | Cohort aggregation + rarity + phenotype for SVs | Yu et al., 2021 preprint: 554,060,126 Manta/Canvas events in 71,408 100kGP rare-disease participants; pilot of 4,313 families found 36 novel protein-coding disrupting SVs on diagnostic-grade genes; estimated ≥4-fold increase in SV-based yield; recall 6/7 previously made SV diagnoses excluding two miscalls (85.7%). SGP used SVRare retrospectively on SRS (Dutta et al., 2025) |
| **SVAFotate** | Annotate SVs with population AF from CCDG, gnomAD-SV, 1000G, TopMed | SGP LRS AF annotation, v0.1.0 (Dutta et al., 2025), gnomAD-SV v4.2.1 |
| **AnnotSV** | SV annotation | Solve-RD: AnnotSV 3.1.1 (Steyaert et al., 2025) |
| **Jasmine** | Population-scale SV merge | Solve-RD: 251,672 unique SVs in 166 unrelated individuals (Steyaert et al., 2025) |

### 5.3 Long-read SV callers

| Tool | Platform | Rare-disease use named in sources |
| ---- | -------- | --------------------------------- |
| **pbsv** | PacBio HiFi | Solve-RD: pbsv 2.4.0, min 20 bp (Steyaert et al., 2025). GA4K: pbsv 2.4.0 (Cohen et al., 2022) |
| **Sniffles2** | ONT and PacBio | SGP family calling: Sniffles v2.0.6 multisample, `--minsvlen 50`, `--tandem-repeats`, long-event caps 1 Mb (Dutta et al., 2025). Best Mendelian consistency among the family modes they tested |
| **cuteSV** | ONT and PacBio | SGP single-sample: cuteSV v1.0.13 with minimap2; highest F1 on GIAB HG002 ONT (F1 0.93 at 33×, 0.79 at 10×) (Dutta et al., 2025) |
| **SURVIVOR** | SV VCF merge | SGP benchmarking consensus (Dutta et al., 2025) |
| **HiPhase** | PacBio phasing of SNV+SV+TR | Holt et al., 2024; N50 0.48 Mb, 88% of genes fully phased (cited by Del Gobbo and Boycott, 2025) |
| **Paraphase** | HiFi haplotypes of paralogs (*SMN1/2* and 160 long segmental duplications) | Chen et al., 2023, *Am. J. Hum. Genet.*; Genomics England Research Consortium among authors |

Solve-RD HiFi (mean depth **9.8×**, one SMRT cell) still called a mean of **55,658 SVs ≥20 bp** (23,385 ≥50 bp) and **4,700,505 SNVs** per person, with 13,481 SVs and 43,172 SNVs private to a family (Steyaert et al., 2025). Candidate de novo SVs fell from 32.6 to 8.5 after QC and to **0.81 per trio** after IGV, illustrating how noisy modest-coverage SV calling is without visual and orthogonal checks.

SGP ONT (74 people, 24 families) retained **60,022** filtered SVs; **23,024–25,009** per family (median 23,814). After quality, rarity, PanelApp green-gene, and inheritance filters: **392** candidate SVs genome-wide and **8** in diagnostic-grade genes. Pathogenic or likely pathogenic de novo SVs in **3 of 24** families, including an *AUTS2* inversion and events at *DLX5/6* and *FN1* (Dutta et al., 2025). Benchmarking: **minimap2 + cuteSV** for single samples; **minimap2 + Sniffles2** for families.

### 5.4 Repeat expansions

| Tool | Input | Rare-disease use named in sources |
| ---- | ----- | --------------------------------- |
| **ExpansionHunter** | Short-read WGS; sequence graphs | Dolzhenko et al., 2017/2019. First STR caller that consistently genotypes repeats longer than the read. Used for *HTT* interruption analysis (Wright et al., 2019, cited therein). Stranger annotates ExpansionHunter VCFs with pathogenic ranges (Clinical Genomics, Sweden) |
| **GangSTR, TREDPARSE, STRetch** | Short-read STR | Named as subsequent expansion callers (Dolzhenko et al., 2019) |
| **TRGT** (Tandem Repeat Genotyper) | PacBio HiFi | Dolzhenko et al., 2024, *Nat. Biotechnol.* Solve-RD: TRGT 0.3.3 at **56** known disease STR loci (Steyaert et al., 2025). Diagnoses: *NOP56* SCA36, *DAB1* SCA37, *RFC1* CANVAS |
| **Stranger** | Annotates ExpansionHunter/TRGT with pathologic copy-number ranges | Clinical Genomics (Sweden) software; TRGT mode parses `GT.MC` |

Del Gobbo and Boycott (2025) review LRS as a single assay for size, motif interruptions, and methylation of pathogenic STRs (*FMR1*, *HTT*, *RFC1*, *DMPK*, *FGF14* SCA27B, *NOTCH2NLC*, *SAMD12*, and others). Genome-wide novel STR discovery is still limited by the size of LRS population catalogs.

### 5.5 Mosaic variants

Mosaic SNVs and SVs cause a non-trivial fraction of congenital neurological and developmental disease (Barbitoff et al., 2024). Unlike tumour calling, there is usually **no matched normal**. Wang et al. (cited by Barbitoff) recommend:

- candidate discovery with **Mutect2** or **GATK HaplotypeCaller ploidy 50**
- heavy filtering on VAF, population frequency, and accessible regions
- sensitivity still ~**65%** even at ~250× WGS

Cooperstein et al. (2025) excluded mosaic UDN diagnoses from Exomiser benchmarking because those cases were often solved by methods other than standard ES/GS, and because VAF filters of 15–85% would drop mosaic heterozygotes.

No rare-disease consortium methods paper in this corpus published a production mosaic-caller version string beyond those citations.

### 5.6 Mitochondrial genome

Smedley et al. (2021) state that 14% of 100kGP pilot diagnoses required the combination of research and automated analysis, and that this combination was critical for etiologic **noncoding, structural, and mitochondrial** variants and for coding variants poorly covered by exome sequencing. SeqOne GermlineVar includes the **GATK mitochondrial pipeline** (Ruzicka et al., 2025, preprint). Dedicated mtDNA callers (e.g. MToolBox, Mutserve, mitoCaller) were **not named with versions** in the consortium pipeline papers read for this review.

---

## 6. RNA and splicing evidence

RNA-seq does not replace DNA calling. It supplies functional evidence for VUS, finds splicing/expression defects that point back to a DNA allele, and can reveal variants outside the capture kit.

| Study | Tissue / N | Method | Yield |
| ----- | ---------- | ------ | ----- |
| Kremer et al., 2017, *Nat. Commun.* | Fibroblasts, 48 undiagnosed mitochondriopathy | Aberrant expression, splicing, mono-allelic expression | **5/48 (10%)** diagnoses; median 1 expression outlier, 5 splicing events, 6 MAE rare variants |
| Frésard et al., 2019, *Nat. Med.* | Blood + large controls | Outlier expression | Rare-disease gene identification from blood transcriptome |
| Cummings et al., 2017, *Sci. Transl. Med.* | Muscle, Mendelian neuromuscular | Transcriptome vs DNA | Splice-defect diagnoses |
| Yépez et al., 2022, *Genome Med.* | Fibroblasts, 303 (197 WES-inconclusive) | **DROP**: OUTRIDER, FRASER, MAE, plus RNA variant calling | **32/197 = 16%** (95% CI 11–22%) of WES-inconclusive cases; 25/32 (78%) were expression outliers; 28% of 46 causal variants already P/LP; SpliceAI recovered 12/21 (57%) pathogenic splice variants at Δ>0.5 |
| Gonorazky et al., 2019; Lee et al., 2020 | Neuromuscular / mixed Mendelian | Clinical RNA-seq | Cited by Del Gobbo and Boycott as part of the 8–36% RNA-seq added-yield range |

**DROP** (Detection of RNA Outliers Pipeline; Yépez et al., *Nat. Protoc.*) wraps:

- **OUTRIDER** — autoencoder denoising of expression; FDR outliers
- **FRASER / FRASER 2** — aberrant splicing (ψ5/ψ3, θ)
- **MAE** — GATK ASEReadCounter + DESeq2/tMAE
- optional RNA-seq variant calling with gnomAD annotation

nf-core/drop documents that outlier detection needs a **sufficient cohort (authors recommend >30 samples)**. Mixing local BAMs with external counts is supported; splicing junction identity must match exactly.

**SpliceAI** (Jaganathan et al., 2019, *Cell*) is a 32-layer residual CNN that predicts acceptor/donor gain or loss from sequence. The authors estimate that **9–11%** of pathogenic mutations in rare genetic disease are cryptic splice variants; de novo predicted splice-altering variants are enriched in autism and intellectual disability and validated in 21/28 RNA-seq patients. Yépez et al. (2022) found SpliceAI more sensitive than MMSplice on their pathogenic splice set, and still weak on deep-intronic and coding splice variants — the reason RNA-seq exists.

Cooperstein et al. (2025) added SpliceAI to Exomiser’s pathogenicity sources specifically to recover splice diagnoses without using CADD’s compressed 0–1 scale.

Long-read RNA-seq (Iso-Seq) resolves full-length isoforms, intron retention, and phasing of splice effects (examples in *MFN2*, *CLPB*; Del Gobbo and Boycott, 2025). It is complementary evidence, not a DNA caller.

---

## 7. Databases

Each row is the question the resource answers, the version current in the papers/sites used for this synthesis, and the ultra-rare caveat. HGMD is named only as a licensed resource; no HGMD records were retrieved.

| Resource | Question it answers | Version / scale in sources | Ultra-rare caveat |
| -------- | ------------------- | -------------------------- | ----------------- |
| **gnomAD** | Population AF, popmax, gene constraint (LOEUF, pLI, missense Z), coverage | v4.1: **730,947 exomes + 76,215 genomes**, GRCh38 ([gnomAD downloads](https://gnomad.broadinstitute.org/downloads)). Genome constraint map: 76,156 genomes (Chen et al., 2024, *Nature*). v4.1 regenerated constraint; non-UKB exome subset exposed | Absence is **not** pathogenicity. Biobanks still contain rare-disease individuals (ClinGen SVI, June 2025 gnomAD v4 note: use calculated BA1/BS1/PM2 frequencies, not strict absence). Ancestry-specific AF still required (Barbitoff et al., 2024) |
| **gnomAD-SV / gnomAD CNV** | SV/CNV frequency | gnomAD-SV used by GATK-SV; SGP annotated with gnomAD-SV v4.2.1 via SVAFotate | Short-read SV catalogues miss many insertions and repeats that LRS sees |
| **ClinVar** | Submitted clinical classifications | Landrum et al.; Cooperstein et al. disabled ClinVar whitelist for benchmarking because UDN deposits solved variants | One-third of AI-MARRVEL diagnostic variants were **not** ClinVar pathogenic; only 8% of 1,586 ClinVar-P alleles in those cases were the true diagnosis (Mao et al., 2024). Conflicting/single-submitter records need care |
| **OMIM** | Gene–disease validity and phenotype | Continuously updated; >6,400 phenotypes with known basis (OMIM gene map, cited by Del Gobbo and Boycott, 2025) | Lag for genes reported in the last year; Cooperstein et al. recovered some such cases only with multi-species hiPHIVE |
| **Orphanet** | Epidemiology, ORPHAcodes, prevalence bins | Nguengang Wakap et al., 2020 | Prevalence bins are sparse for ultra-rare / n-of-1 |
| **HPO** | Phenotype encoding for rankers | Köhler et al.; UDN median 20 HPO terms (range 1–117) in Cooperstein et al. | Noisy or missing HPO drops Exomiser/AIM; random HPO still inflates known-gene ranks and must not be used on unsolved cases |
| **Gene2Phenotype / DDG2P** | Developmental-disorder gene–disease confidence | DDD; PanelApp R27 paediatric superpanel remapped G2P to GenCC terms | G2P is stronger for paediatric DD than for adult ultra-rare metabolic disease |
| **PanelApp** | Diagnostic-grade (“green”) virtual panels and MOI | Genomics England; 100kGP and NHS GMS tiering (Martin et al.; PanelApp site) | Green genes miss new associations. SGP and 100kGP tier only green genes as Tier 1/2 |
| **DECIPHER** | Patient genotype–phenotype sharing, CNV tracks, protein browser | Firth et al., 2009; DECIPHER 2023 review | Consent-limited; not a population AF source |
| **GenCC** | Harmonised gene–disease validity (ClinGen, PanelApp, G2P, OMIM, Orphanet, …) | DiStefano et al., 2022 | Validity ≠ this patient’s diagnosis |
| **ClinGen** | Expert panels, VCEPs, ACMG specifications, dosage sensitivity | SVI June 2025 gnomAD v4 letter | Gene-specific rules lag new genes |
| **DGV / dbVar** | Historical SV catalogues | MacDonald et al. / NCBI | Mixed clinical-grade; prefer gnomAD-SV / CoLoRSdb for rarity |
| **CoLoRSdb** | PacBio LRS SV/STR frequency | >1,400 PacBio genomes (HPRC, HGSVC, GA4K, …) ([colorsdb.org](https://colorsdb.org/)) | Still small vs gnomAD; ancestry thin |
| **1000 Genomes ONT** | LRS population SVs | 1,019 samples at 16.9× (Schloissnig et al., 2024); 100 samples ≥30× (Gustafson et al., 2024) | Intermediate coverage is for SVs, not small-variant truth |
| **All of Us LRS** | Planned population HiFi | Mahmoud et al., 2024; first 1,000 on Researcher Workbench | Registered access |
| **LOVD / LSDBs** | Locus-specific variants | Fokkema et al. | Patchy gene coverage |
| **Matchmaker Exchange** | n-of-1 gene matching across nodes (GeneMatcher, DECIPHER, MyGene2, RD-Connect GPAP, …) | Philippakis et al., 2015. Solve-RD submitted *PSMA3* via MME and GPAP (Steyaert et al., 2025) | The method of record when no second case exists |
| **RD-Connect GPAP** | Genome-phenome analysis platform | Laurie et al., 2022 | Solve-RD data interpretation |
| **MANE** | Clinical transcript set (RefSeq + Ensembl) | Morales et al., 2022 | Annotation discordance otherwise reaches 40%+ for pLoF (Barbitoff et al., 2024, citing McCarthy/Frankish literature) |
| **VEP / SnpEff / ANNOVAR** | Molecular consequence | McLaren et al., 2016; Cingolani et al.; Wang et al. | Tool × transcript set discordance, especially indels |
| **LOFTEE / pext** | High- vs low-confidence pLoF | Karczewski et al., 2020; Cummings et al. | Still miss population-specific false pLoF (Barbitoff et al., 2024) |
| **dbNSFP** | Precomputed missense scores | Liu et al. | Score vintage and transcript must match the pipeline |
| **HGMD** | Licensed mutation database | Stenson et al. | Not audited here. AIM used HGMD as a feature; AIM-NDG drops it |

ClinGen SVI (June 2025) asked VCEPs to move to **gnomAD v4**, keep calculated frequency thresholds, and **not** require a “non-cancer” gnomAD slice (TCGA was already removed in v4).

---

## 8. AI models, split by job

One model does not do every job. The tables below keep callers, variant-effect predictors (VEPs), splice models, and case-level rankers apart. Phenotype-only face-matching is listed last so it is not described as a variant caller.

### 8.1 Callers that are neural networks

Covered in §4: **DeepVariant**, **Clair3**, **PEPPER-Margin-DeepVariant**, GATK **CNNScoreVariants**. Barbitoff et al. (2022) found DeepVariant robust on GIAB and non-GIAB samples; other ML filters (GATK 2D CNN, Octopus RF) were coverage-sensitive.

### 8.2 Missense and generic VEPs

| Model | Supervision | Rare-disease use named in sources |
| ----- | ----------- | --------------------------------- |
| **REVEL** | Ensemble of missense tools; ClinVar-heavy training (Ioannidis et al., 2016) | Exomiser default pathogenicity source with MVP (Cooperstein et al., 2025). Cooperstein warn that ClinVar overlap can inflate ranks on solved UDN variants |
| **MVP** | Missense pathogenicity (Qi et al.) | Exomiser default with REVEL |
| **AlphaMissense** | AlphaFold-based; trained on population-frequency weak labels, **not** ClinVar (Cheng et al., 2023, *Science*). Classifies 89% of proteome missense as likely benign or likely pathogenic | Added to Exomiser with SpliceAI; modest cohort-wide gain, larger for missense (Cooperstein et al., 2025). NHS Genomic AI log: Guy’s and St Thomas’ / Synnovis |
| **CADD** | Combined annotation, genome-wide (Kircher / Rentzsch) | Useful beside Exomiser; **do not** feed Exomiser-normalised CADD as a max-score source (almost all coding variants compress to ≥0.98) (Cooperstein et al., 2025) |
| **SpliceAI** | See §6 / §8.3 | Exomiser + AIM feature |
| **SIFT / PolyPhen-2 / MutationTaster** | Classical | Overestimate pathogenicity and dominate Exomiser’s max-score if included (Cooperstein et al., 2025) |
| **CAPICE** | Gradient-boosting of CADD-like features for rare variants (Li et al., 2020) | Named in NHS Genomic AI tool log as a rare-exome ensemble; not re-benchmarked in Cooperstein or Mao |
| **BayesDel** | Deleteriousness with/without AF (Feng et al., 2017) | CAGI/clinical VEP comparisons (Livesey et al. / CAGI assessments, 2025) |
| **PrimateAI-3D** | 3D structure + primate constraint, not ClinVar-supervised (Gao et al., 2023) | Independent VEP class; not in the 2025 Exomiser recommended set |
| **EVE** | Deep generative model of MSA, no disease labels (Frazer et al., 2021). 36M variants / 3,219 disease genes | Independent evidence for VUS; gene coverage limited to modelled proteins |
| **ESM-1b** | Protein language model (Brandes et al., 2023; Rives et al., 2021) | Same independent class as EVE/AlphaMissense |
| **MutPred2 / MutPred-LoF** | Supervised structural/sequence | Named by Barbitoff et al., 2024 for pLoF/missense |
| **M-CAP, ClinPred, VEST4, Eigen, FATHMM-XF** | Various supervised ensembles | Common in clinical score panels (CAGI 6 assessment) |

Cooperstein et al. (2025) **recommended pathogenicity sources for Exomiser**: REVEL + MVP + AlphaMissense + SpliceAI. They **excluded** CADD and SIFT/PolyPhen/MutationTaster from the max-score pool.

### 8.3 Splicing models

| Model | Role | Named use |
| ----- | ---- | --------- |
| **SpliceAI** | Sequence-only splice disruption (Jaganathan et al., 2019) | Exomiser, Genomiser, AIM, SGP compound-het SNV filter (VEP + masked SpliceAI), Yépez 2022 RNA-seq comparison, NHS labs |
| **MMSplice** | Modular splicing (Cheng et al., 2019) | Weaker than SpliceAI on Yépez pathogenic set (8/21 vs 12/21) |
| **Pangolin** | Deep splice model | Best with SpliceAI vs massively parallel splicing assays (Wagner et al. / Genome Biol. splice benchmark, 2023) |
| **ReMM** | Regulatory Mendelian mutation score for noncoding | Genomiser default; Cooperstein et al. found ReMM **hurt** their small noncoding UDN set (splice-heavy) and still recommend keeping ReMM + SpliceAI for promoter/enhancer cases |

### 8.4 Case-level rankers (phenotype + genotype)

These tools consume a VCF, HPO terms, and often a PED file. They output a ranked gene/variant list. They are **not** SV-complete unless explicitly stated.

| Tool | Algorithm | Training / knowledge | Named rare-disease use and numbers |
| ---- | --------- | -------------------- | ---------------------------------- |
| **Exomiser** | Frequency × pathogenicity × phenotype (PhenIX / PHIVE / hiPHIVE) | Knowledge-based; regular annotation releases | 100kGP pilot: 77% of diagnoses as top candidate, 86% top 3, 88% top 5 (Smedley et al., 2021, n=490 diagnoses). Later 100kGP analysis: 82.6 / 91.3 / 92.4 / 93.6% at top / 3 / 5 / 10 of **4,877** diagnoses; singletons similar (Smedley et al., 2022, *Hum. Mutat.*). UDN optimization (Cooperstein et al., 2025, v14.0.0, dataset 2406): GS top-10 **49.7% → 85.5%** (raw default → optimized); filtered default already 66.6%; ES **67.3% → 88.2%**. Independent 17 newly diagnosed UDN GS probands: 22/23 variants in top 30. Implemented in UDN Mosaic platform |
| **Genomiser** | Exomiser + noncoding + ReMM | Same group (Smedley et al., 2016) | UDN GS noncoding cohort: top-10 **15.0% → 40.0%** (Cooperstein et al., 2025). Required for coding+noncoding compound hets; not a replacement for Exomiser (~30% fewer coding diagnoses in top 30, ~3× runtime) |
| **LIRICAL** | Likelihood ratio of HPO vs genomic findings (Robinson et al., 2020) | Knowledge-based | 100kGP singletons: 85.2% of diagnoses in top 5 vs 94.3% Exomiser (Smedley et al., 2022). AIM benchmark comparator (Mao et al., 2024). Yuan et al. 2022 DDD/in-house benchmark: LIRICAL and AMELIE complementary |
| **PhenIX** | Human-only phenotype + variant (Robinson / Exomiser) | Human clinical associations only | AIM comparator. Cooperstein: PhenIX 62.8% GS top-10 vs hiPHIVE default 66.6% vs human-only hiPHIVE 82.8% |
| **Xrare** | ML jointly modelling phenotype and ACMG-like genetic features (Li et al., 2019, *Genet. Med.*) | Supervised | AIM comparator. Jacobsen et al. / *Trends Genet.* review: similar top-1 to Exomiser/LIRICAL on one 134-patient set |
| **AMELIE** | Matches patient to primary literature (Birgmeier et al., 2020, *Sci. Transl. Med.*) | Literature NLP | No local programmatic access in AIM benchmark; similar DDD top-1 to AIM, AIM +7.6% at top 5 (Mao et al., 2024, citing Birgmeier) |
| **AI-MARRVEL (AIM)** | Random forest; 56 raw + 47 expert-engineered features; VEP + MARRVEL DBs (DGV, DECIPHER, ClinVar, OMIM); SpliceAI | Trained on **>3.5 million variants** from diagnosed cases. DiagLab: 1,044 train / 58 test; also 75 UDN and 200 DDD as independent tests. 80/20 split of remaining train for RF hyperparameters. Weighted impurity for 7,791:1 negative:diagnostic ratio | Outperformed Exomiser, LIRICAL, PhenIX, Xrare on three cohorts (Mao et al., 2024). Confidence metric: 98% precision, identified 57% of diagnosable cases among 871. Recessive-specific and trio-specific models. AIM-NDG (no OMIM/ClinVar/HGMD) recovered UDN novel genes *MYCBP2* and *TMEM161B*. **Does not call/rank SVs or CNVs.** Trained mostly on exomes/coding; MAF filter <0.01 with ClinVar/HGMD rescue |
| **DiagAI** | Linear model on UP² (XGBoost ACMG-class probabilities, 74 features, ClinVar 04-2024) + PhenoGenius + PanelApp/OMIM inheritance | **Preprint** (Ruzicka et al., 2025, medRxiv [doi:10.1101/2025.02.04.25321641](https://doi.org/10.1101/2025.02.04.25321641)). 196 diagnosed nephrology exomes of 966 sequenced | 94.9% of causal variants shortlisted with HPO (median shortlist 12); 90.8% without (median 9); 74% vs 42% top-ranked. Authors report outperforming Exomiser v13 and AIM on this cohort. Calling stack: FreeBayes, GATK, GRIDSS, AluMEI, GATK-mitochondrial |
| **GEM, MOON, Emedgene** | Commercial phenotype-aware rankers | Named in Mao et al. Table S1 | Not independently re-run in the open benchmarks above |
| **VAAST, Phevor, Phen-Gen, Phenolyzer** | Older gene rankers | Named in Mao et al. | Superseded in the UDN/100kGP papers by Exomiser-class tools |

**Cooperstein et al. (2025) Exomiser/Genomiser recipe (UDN, implemented in Mosaic):**

1. Joint family VCF; filter heterozygous VAF 15–85%, GQ ≥20, drop `*` alleles.
2. Pathogenicity sources: REVEL + MVP + AlphaMissense + SpliceAI (Genomiser: also ReMM).
3. **Human-only hiPHIVE** (not default multi-species+PPI).
4. Pedigree + inheritance filters; **enable ClinVar whitelist** on unsolved cases (disabled during benchmarking).
5. Review top 30 contributing variants; flag genes that hit top 30 at p≤0.3 in ≥5% of the local cohort (86 such genes in 231 GS cases) rather than hard-blacklisting.
6. If empty: check pedigree/penetrance; rerun proband-only; for GS, run Genomiser.
7. Optional p≤0.3 to widen the list for RNA/SV integration (drops 1.7% of GS diagnostic variants).

Human-only hiPHIVE raised GS top-10 from 66.6% to 82.8%. Adding AM+SpliceAI added 2.7% more. The cost: genes with disease associations **younger than the OMIM/Orphanet freeze** can lose phenotype score (seen in 4 of 17 new UDN diagnoses). Multi-species hiPHIVE remains a fallback.

**Mao et al. (2024) AIM facts that matter for ultra-rare disease:**

- Average person carries **15.5** ClinVar-pathogenic variants; ~5.6 of those are recessive-only. 94% of cases were still diagnosed with one or two variants.
- Removing ClinVar features (AIM-withoutVarDB) only slightly dropped accuracy.
- Irrelevant HPO: −11% top-1, but top-5 still 78%.
- Feature ablation: MAF, database status, and phenotype have **similar** effect sizes; none is sufficient alone.
- AIM learned that genes with ≥3 coding variants in one person are less likely rare-disease genes.
- AIM-Recessive recovered 26/41 remaining recessive pairs as top-1.
- AIM-Trio top-1 40% → 45% (16% of diagnostic variants were inherited dominant).
- Limitation stated by the authors: no SV/CNV, exome-trained, coding-heavy.

These two papers do **not** declare a universal winner. Cooperstein selected Exomiser because it is the most widely deployed open tool with GS/noncoding support. Mao reports AIM ahead of Exomiser on their three test sets. Different truth sets, different Exomiser versions, and AIM’s lack of SVs explain the disagreement. Ultra-rare n-of-1 genes punish any ranker that depends on OMIM/ClinVar; AIM-NDG and multi-species hiPHIVE are the documented workarounds, not a solved problem.

### 8.5 Phenotype-only tools (not variant callers)

**GestaltMatcher / Face2Gene** match facial photographs to syndromes (FDNA). The NHS Genomic AI Network tool log lists Face2Gene Gestalt for patient identification. They do not call or rank genomic variants. **ClinPhen**, **PhenoTips**, **PhenoBCBERT / PhenoGPT** extract HPO from notes (AIM uses ClinPhen; UDN uses PhenoTips). Bad HPO is a first-class failure mode for every ranker above.

---

## 9. Who used which stack

Rows require a method paper (or official pipeline write-up) that names the tool and the cohort. Diagnostic yields are those papers’ yields, not a meta-analytic average.

| Programme | Cohort (as published) | Detection stack | Prioritization / interpretation | Yield / added yield |
| --------- | --------------------- | --------------- | ------------------------------- | ------------------- |
| **100,000 Genomes Project rare-disease pilot** (Smedley et al., 2021) | 4,660 participants, 2,183 families, 161 disorder categories | Illumina WGS; applied virtual panels from **PanelApp**; automated filter to rare, segregating, predicted-damaging coding variants; research analysis for noncoding/SV/mtDNA | **PanelApp** green-gene **tiering** (Tier 1 LoF/de novo protein-altering on applied panels; Tier 2 other variants on those genes; Tier 3 the rest). **Exomiser** phenotype search across all genes. HPO collected at recruitment | **25%** of probands; **35%** if the disorder was likely monogenic; 40–55% for ID, hearing, vision. 14% of diagnoses needed research+automated (noncoding, SV, mtDNA, poor exome coverage). Exomiser: 77/86/88% at top/3/5 of 490 diagnoses. 25% of diagnoses changed clinical management. Burden testing in 57,000 genomes: 3 new genes, 19 new associations |
| **100kGP follow-on Exomiser/LIRICAL** (Smedley et al., 2022) | 4,877 diagnoses | WGS VCF GRCh37/38 | Exomiser 82.6–93.6% top through top 10; LIRICAL 85.2% top 5 on singletons. Failures: 1% low-quality VCF flags, 2% high AF, 3% incomplete penetrance. Exomiser ranked 74% of SV diagnoses top, 89% top 5 in that analysis | Same programme, larger diagnosed set |
| **100kGP SVRare** (Yu et al., 2021, preprint) | 71,408 rare-disease participants; 4,313-family pilot | **Manta + Canvas** aggregated by **SVRare** (554,060,126 events) | Rarity + coding disruption + phenotype/panels | 36 novel diagnostic-grade-gene SVs in the pilot; estimated ≥4-fold SV-yield increase |
| **Scottish Genomes Partnership SRS** (cited in Dutta et al., 2025) | 999 genomes, 394 families | HiSeqX 150 bp, **iSAAC** to GRCh37/38, GEL SNV/indel pipeline | PanelApp Tiers 1–2 then Exomiser / Tier 3 | **23%** families |
| **SGP LRS** (Dutta et al., 2025) | 24 unsolved trios/quads, 74 people | ONT R9.4.1 ± R10.4.1; **minimap2** 2.23; **Sniffles2** 2.0.6 family calling; **cuteSV** 1.0.13 single-sample; VEP 110.1; SVAFotate (CCDG, gnomAD-SV 4.2.1, 1000G, TopMed); **SVRare** on historical Manta/Canvas | PanelApp green genes; slivar de novo/AR/XL; compound het = LRS SV + Illumina SNV (gnomAD v4.1 AF<0.01, SpliceAI, CADD) | **3/24** families with P/LP de novo SVs (*AUTS2*, *DLX5/6*, *FN1*) |
| **Solve-RD HiFi** (Steyaert et al., 2025) | 293 people, 114 families (93 unsolved neuro/NMD/epilepsy; 21 “unsolvable” syndromes). ERNs RND, EURO-NMD, ITHACA, EpiCARE | PacBio Sequel IIe, 1 SMRT cell, mean **9.8×** HiFi; **pbmm2** 1.4.0 → hg38; **DeepVariant** 1.1.0; **GLnexus** 1.3.1; **WhatsHap** 1.1.0; **TRGT** 0.3.3 (56 loci); **pbsv** 2.4.0; **AnnotSV** 3.1.1; Jasmine merge | ERN gene lists for inherited SVs <100 kb and inherited SNVs; genome-wide de novo; visual IGV; orthogonal PCR/Sanger/OGM/RNA-seq | **11/93 (11.8%)** unsolved families diagnosed; **5.4%** extra candidate SVs; **1/21** unsolvable (de novo *TUBA1A* misdiagnosed as Aicardi). Variant classes: DMD/NHS inversions, *REEP1* exon deletion, *MYOT* tandem dup, STR expansions (*NOP56*, *DAB1*, *RFC1*), deep-intronic *DMD*, splice *SPAST*, phased *TTN* compound het |
| **Solve-RD reanalysis of short-read resource** (Laurie et al., 2025, *Nat. Med.*) | Pan-European RD resource | Reanalysis (paper named in Steyaert discussion) | New diagnoses from existing data | Cited as the short-read companion; details not re-extracted here beyond the citation |
| **Undiagnosed Diseases Network** (Splinter et al., 2018, *N. Engl. J. Med.*) | 1,519 referred; 601 accepted; 382 complete evaluations | Clinical + research ES/GS (two sequencing cores) | Multidisciplinary; later computational HPO (PhenoTips in Cooperstein et al., 2025) | **132/382 (35%)**. 74% of diagnoses by ES/GS (98); ES through UDN 55/195 (28%); GS 32/165 (19%). 31 new syndromes |
| **UDN Exomiser deployment** (Cooperstein et al., 2025) | 5,412 GS samples / 1,772 families (Nov 2023); 1,454 ES samples / 314 families (Jun 2024); benchmark 386 diagnosed probands | **CGAP** alignment to GRCh38 with decoys+ALTs; per-sample GVCF; **Sentieon** joint SNV/indel (no SVs in this VCF) | Exomiser/Genomiser v14.0.0, dataset 2406, CADD v1.7 and ReMM v0.4 on disk; BCFtools 1.16 filters; Mosaic UI | Optimized GS top-10 85.5%; ES 88.2%; Genomiser noncoding top-10 40%. **SVs and mosaics not in this VCF** |
| **UDN / Baylor AIM** (Mao et al., 2024) | 1,102 DiagLab (1,044+58), 75 UDN, 200 DDD | Clinical ES VCFs | AIM vs Exomiser, LIRICAL, PhenIX, Xrare | AIM higher gene rank on those three sets; AIM-NDG flagged *MYCBP2*, *TMEM161B* |
| **Deciphering Developmental Disorders** (Wright / DDD Study, 2015–2017) | 1,133 then 4,293–4,294 children | **BWA** 0.59, hs37d5; **GATK** 3.1.1 realign+BQSR; **HaplotypeCaller** 3.2.2 multisample; **VQSR**; **DeNovoGear** 0.54; CNVs: CNsolidate, CoNVex, CIFER; UPD: UPDio; mosaicism: triPOD; midsized indels: Indelible | DDG2P; DECIPHER reporting; ExAC/1000G/internal parental AF | 2015: **27%** of 1,133. 2017: 8,361 candidate DNMs (mean 1.95/proband); 95 genome-wide significant genes; 42.5% estimated to carry pathogenic de novo coding SNVs/indels. Reanalysis of 1,133 (Wright et al., 2018, *Genet. Med.*): **40%** (454/1,133), 78% of diagnoses de novo |
| **Genomic Answers for Kids** (Cohen et al., 2022) | 960 families / 1,080 patients | srES on all; srGS (Illumina PCR-free, MGI, 10x Linked-Read — 10x later dropped); PacBio **HiFi** on a subset (n=274/472 in platform table; ~550 HiFi in contemporaneous reporting). HiFi: **pbmm2** 1.4.0, **pbsv** 2.4.0, **DeepVariant** v1.0 PacBio, **hifiasm** 0.9 | Machine-learning variant prioritization; PhenoTips; dbGaP sharing | **11%** if prior negative testing; **34.5%** if testing-naive. SVs from GS added up to **13%** of new diagnoses in previously unsolved. HiFi: **>4-fold** more rare coding SVs than srGS; ~5% more SNVs at 25× HiFi vs 35× Illumina; nearly double SV discovery. VUS/GUS remain 58% of nondiagnostic cases. Rural follow-up (2024): >4,700 probands sequenced, >1,370 diagnosed by end of 2023 |
| **GA4K methylation / 5-base HiFi** (Cheung et al., 2023) | 276 people, 152 undiagnosed paediatric families | PacBio 5-base HiFi | Rare hypermethylation outliers | *DIP2B* repeat-expansion hypermethylation diagnosis |
| **GA4K pangenome SVs** (Groza et al., 2024) | 574 assemblies from paediatric RD + 94 controls | Graph pangenome | Improved rare SV prioritisation | Novel diagnostic SV in *KMT2E* |
| **Hiatt et al. NDD LRS** (2021, 2024) | 6 then 96 undiagnosed NDD/MCA/myopathy | PacBio LRS | De novo and gene-based SV | 2024: 16/96 with disease-relevant variants; **7/96 (7.3%)** exclusively LRS |
| **Sanford Kobayashi et al., 2022** | 30 patients, 26 paediatric RD families | PacBio HiFi | Challenging-region SNVs | One extra diagnosis after SR-GS: hemizygous *IKBKG* stop-loss; 98% of SRS “dark” regions covered |
| **AlAbdi et al., 2023** | 34 AR families, ES-negative | PacBio LRS, ROH analysis | AR missing-variant hunt | 13/34 diagnosed; **8/34 (23.5%)** not detectable by SR-ES |
| **DROP clinical RNA-seq** (Yépez et al., 2022) | 303 fibroblasts (197 WES-inconclusive), mostly mitochondrial | WES: Agilent v5/v6, BWA 0.7.5a, hg19. RNA: DROP (OUTRIDER, FRASER, MAE, RNA calling) | Manual OMIM/phenotype overlay; some WGS follow-up | **16%** of WES-inconclusive diagnosed |
| **SeqOne / nephrology lab** (Ruzicka et al., 2025, **preprint**) | 966 adult nephropathy exomes, 196 diagnosed | Twist exome, NovaSeq, BWA, GermlineVar (FreeBayes, GATK, GRIDSS, AluMEI, GATK-mito) | **DiagAI** + UP² + PhenoGenius + PanelApp 240807 | 24% diagnostic rate in the lab; DiagAI shortlist metrics in §8 |
| **NHS Genomic AI Network tool log** (public PDF, Jan 2025) | NHS deployments, not a cohort paper | N/A | Lists **Exomiser**, **AI-MARRVEL**, **AlphaMissense**, **SpliceAI**, **REVEL**, Face2Gene Gestalt | Deployment log, not a yield study. Cited only as named use |

**NHS Genomic Medicine Service.** PanelApp documentation states that version 1+ green panels are used for GMS rare-disease **tiering** (Tier 1/2). This review did not re-fetch a GMS methods paper that names the current aligner/caller (historically the 100kGP iSAAC/Manta/Canvas path). Absence of that paper here is not evidence that GMS still uses iSAAC.

---

## 10. What is specific to ultra-rare disease

1. **Private alleles.** Solve-RD saw 13,481 family-private SVs and 43,172 family-private SNVs at 10× HiFi. Population AF from gnomAD is necessary but not sufficient; a private SV still needs a local LRS frequency panel (CoLoRSdb, Jasmine callset, SGP AF).

2. **Reference-allele and ancestry bias.** Thousands of coding RMAs sit in GRCh37/38/T2T (Barbitoff et al., 2024). Two SNVs in an HSP panel were missed by every caller because the reference carried the rare allele (Cacheiro et al., 2017). popmax and ancestry-matched AF reduce false “novel” calls in underrepresented groups.

3. **GIAB high-confidence regions miss clinical space.** 12.9% of known pathogenic variants and 4.7 Mbp of CDS sit outside the regions that produce 99% F1 scores. Sex-chromosome and paralog genes (*DMD*, *SMN1*, *PKD1*, *IKBKG*, *NEB*, *TTN*) are enriched here. Paraphase, HiFi, and graph genomes are the documented technical answers.

4. **SV frequency databases are still weak.** SVRare’s value in 100kGP was exactly that Manta/Canvas calls had no usable AF. LRS finds 3–5× more SVs than SRS (Del Gobbo and Boycott, 2025); without CoLoRSdb-scale AF, most of those SVs remain uninterpretable. Solve-RD say this is currently the main limit on LRS yield.

5. **Repeats, inversions, and noncoding splice.** Short-read panels and even SR-GS miss copy-neutral inversions (*DMD*, *AUTS2* in Solve-RD and SGP), motif-interrupted STRs (*NOP56*, *DAB1*, *RFC1*, *FGF14*), and deep-intronic splice alleles (*DMD* c.31+36947G>A). TRGT, ExpansionHunter, SpliceAI, and RNA-seq are complementary, not redundant.

6. **Supervised rankers fail on n-of-1 genes.** Exomiser human-only hiPHIVE and AIM both rely on human gene–disease maps. Cooperstein et al. lost rank on associations younger than one year. AIM-NDG drops OMIM/ClinVar/HGMD and still only matches other tools that *keep* those databases. For a true novel gene, **Matchmaker Exchange** (Philippakis et al., 2015) is the method, as used for Solve-RD *PSMA3*.

7. **Pedigree error is a hard filter.** Cooperstein et al. excluded 24 GS diagnostic variants because of incomplete penetrance or misphenotyped relatives; 21/24 came back on proband-only VCFs. Ultra-rare dominant disease with mildly affected parents will look like a de novo pipeline failure.

8. **“Unsolvable” syndromes are a different sampling frame.** Solve-RD found no shared locus in Aicardi, Hallermann–Streiff, Gomez–Lopez-Hernandez, Pai, or OAVS families, and only one mis-syndrome SNV. Boycott et al. (2018, cited therein) list mosaicism, methylation, oligogenic, and non-genetic causes as first-class hypotheses — not more of the same SNV caller.

9. **Reanalysis is a method.** Tan et al. (cited by Barbitoff et al., 2024): median **15%** new diagnoses on reanalysis. DDD 27% → 40% on the same 1,133 families. AIM’s confidence metric exists because laboratories cannot manually reread every unsolved ES every year.

10. **Phenotype quality is part of the bioinformatics.** UDN median 20 HPO terms; dropping all HPO collapses Exomiser; random HPO falsely helps known genes. Ultra-rare adult presentations are often under-HPO’d compared with paediatric DD.

---

## 11. Limits of this synthesis

- **Search window:** queries run 17 October 2026. Later tool versions are not in the tables.
- **Not a census of thousands of full texts.** Europe PMC Q1 alone is 9,136 records. Deep reading was reviews, benchmarks, consortium methods, and primary tool papers. Adoption is reported only when those papers name the tool.
- **Index bias.** English-language PubMed/Europe PMC. Panel software from a single hospital that never published the stack is invisible here.
- **Acronym inflation.** Q7, Q10, Q13 Europe PMC counts are not biological. Tightened queries are the ones interpreted.
- **Licensed databases.** HGMD contents were not used. Commercial pipelines (DRAGEN, Emedgene, MOON, GEM, Franklin) appear only when a peer-reviewed or preprint evaluation named them.
- **Preprints.** DiagAI (Ruzicka et al., 2025) and SVRare (Yu et al., 2021) are labelled. Numbers may change in journal versions.
- **Conflicting benchmarks.** AIM vs Exomiser vs DiagAI vs LIRICAL results are cohort- and version-specific. This report does not pick a winner.
- **GIAB F1 is not a clinical SLA.** Highest CDS F1 (DeepVariant in Barbitoff et al., 2022) is the right *benchmark* result. Clinical callers still need solved-case validation, sex-chromosome and paralog tests, and a policy for filtered-but-visible alleles (Barbitoff et al., 2024; European WGS recommendations cited therein).
- **SVs in UDN Exomiser paper.** Cooperstein et al. jointly called only SNVs/indels. UDN SV practice is not documented by that paper.
- **No patient-level reanalysis of VCFs** was performed for this note.

---

## 12. Practical map (what to run for which hole)

This is a reading of the papers, not a laboratory SOP.

| If the unsolved case looks like… | Methods the literature actually used |
| -------------------------------- | ----------------------------------- |
| Negative panel/ES, likely coding SNV | WGS or re-call ES with DeepVariant or GATK-HC; Exomiser with human-only hiPHIVE + REVEL/MVP/AM/SpliceAI; or AIM if ES+HPO |
| Negative SR-GS, suspected SV | Short-read: Manta+Canvas+SVRare (100kGP). Long-read: pbsv (HiFi) or Sniffles2/cuteSV (ONT); AnnotSV; LRS AF (CoLoRSdb/gnomAD-SV) |
| Ataxia / anticipation / known STR | ExpansionHunter on SR-GS; TRGT on HiFi; motif and methylation, not just length |
| Missing second AR allele | LRS phasing (WhatsHap/HiPhase/PEPPER); targeted LRS (Miller et al., Werner syndrome 8/9 second hits, cited by Del Gobbo and Boycott) |
| Splice VUS or muscle/mitochondrial phenotype | DROP (OUTRIDER/FRASER/MAE) on a ≥~30-sample batch; SpliceAI is adjunct |
| Paralog / dark gene (*SMN1*, *IKBKG*, *PKD1*) | HiFi + Paraphase or gene-specific LRS |
| Mosaic malformation | Deep WGS; Mutect2 or HC ploidy 50; do not use 15–85% VAF Exomiser filters |
| Novel gene, n-of-1 | AIM-NDG or multi-species hiPHIVE as hypothesis generators; **Matchmaker Exchange** as the test |
| Unsolvable named syndrome | Do not expect another coding SNV; consider methylation (5-base HiFi), mosaicism, non-genetic |

---

## Bibliography

DOIs and stable URLs. Papers are cited from the text that was read, not from titles alone.

### Reviews, definitions, guidelines

1. Barbitoff YA, et al. Bioinformatics of germline variant discovery for rare disease diagnostics: current approaches and remaining challenges. *Brief Bioinform.* 2024;25(2):bbad508. [doi:10.1093/bib/bbad508](https://doi.org/10.1093/bib/bbad508). PMC10810331.
2. Barbitoff YA, et al. Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery. *BMC Genomics.* 2022;23:155. [doi:10.1186/s12864-022-08365-3](https://doi.org/10.1186/s12864-022-08365-3).
3. Del Gobbo GF, Boycott KM. The additional diagnostic yield of long-read sequencing in undiagnosed rare diseases. *Genome Res.* 2025. PMC12047273.
4. Nguengang Wakap S, et al. Estimating cumulative point prevalence of rare diseases: analysis of the Orphanet database. *Eur J Hum Genet.* 2020;28:165–173. [doi:10.1038/s41431-019-0508-0](https://doi.org/10.1038/s41431-019-0508-0).
5. Ferreira CR. The burden of rare diseases. *Am J Med Genet A.* 2019;179:885–892. [doi:10.1002/ajmg.a.61124](https://doi.org/10.1002/ajmg.a.61124).
6. Richards S, et al. Standards and guidelines for the interpretation of sequence variants: ACMG/AMP. *Genet Med.* 2015;17:405–424. [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30).
7. Wright CF, FitzPatrick DR, Firth HV. Paediatric genomics: diagnosing rare disease in children. *Nat Rev Genet.* 2018;19:253–268. [doi:10.1038/nrg.2017.116](https://doi.org/10.1038/nrg.2017.116).
8. Kernohan KD, Boycott KM. The expanding diagnostic toolbox for rare genetic diseases. *Nat Rev Genet.* 2024;25:401–415. [doi:10.1038/s41576-023-00683-w](https://doi.org/10.1038/s41576-023-00683-w).
9. ClinGen SVI VCEP Review Committee. Communication to ClinGen VCEPs on gnomAD v4. June 2025. [PDF](https://clinicalgenome.org/site/assets/files/9445/june_2025_communication_to_clingen_vceps_from_clingen_vcep_review_committee.pdf).
10. NICE. Strategic principles for rare diseases (ultra-rare <1 in 50,000). [NICE PMG46 resource](https://www.nice.org.uk/process/pmg46/resources/nice-strategic-principles-a-complementary-approach-to-public-health-social-care-and-rare-disease-topics-13430355949/chapter/strategic-principles-for-rare-diseases).

### Consortia and cohort pipelines

11. Smedley D, et al. 100,000 Genomes Pilot on Rare-Disease Diagnosis in Health Care — Preliminary Report. *N Engl J Med.* 2021;385:1868–1880. [doi:10.1056/NEJMoa2035790](https://doi.org/10.1056/NEJMoa2035790). PMID 34758253.
12. Smedley D, et al. Phenotype-driven approaches to enhance variant prioritization and diagnosis of rare disease (100kGP Exomiser/LIRICAL). *Hum Mutat.* 2022. PMC9288531.
13. Yu J, et al. SVRare: discovering disease-causing structural variants in the 100K Genomes Project. medRxiv. 2021. [doi:10.1101/2021.10.15.21265069](https://doi.org/10.1101/2021.10.15.21265069). **Preprint.**
14. Dutta P, et al. Detecting pathogenic structural variation in families with undiagnosed rare disease in a national genome project. *Eur J Hum Genet.* 2026. [doi:10.1038/s41431-026-02210-x](https://doi.org/10.1038/s41431-026-02210-x); medRxiv [doi:10.1101/2025.08.19.25333674](https://doi.org/10.1101/2025.08.19.25333674).
15. Steyaert W, et al. Unraveling undiagnosed rare disease cases by HiFi long-read genome sequencing. *Genome Res.* 2025;35:755–768. [doi:10.1101/gr.279414.124](https://doi.org/10.1101/gr.279414.124). PMID 40138663.
16. Laurie S, et al. Genomic reanalysis of a pan-European rare-disease resource yields new diagnoses. *Nat Med.* 2025;31:478–489. [doi:10.1038/s41591-024-03420-w](https://doi.org/10.1038/s41591-024-03420-w).
17. Splinter K, et al. Effect of genetic diagnosis on patients with previously undiagnosed disease (UDN). *N Engl J Med.* 2018;379:2131–2139. [doi:10.1056/NEJMoa1714458](https://doi.org/10.1056/NEJMoa1714458). PMID 30304647.
18. Cooperstein IB, et al. An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser. *Genome Med.* 2025. [doi:10.1186/s13073-025-01546-1](https://doi.org/10.1186/s13073-025-01546-1). PMC12539062.
19. Deciphering Developmental Disorders Study. Large-scale discovery of novel genetic causes of developmental disorders. *Nature.* 2015;519:223–228. [doi:10.1038/nature14135](https://doi.org/10.1038/nature14135).
20. Deciphering Developmental Disorders Study. Prevalence and architecture of de novo mutations in developmental disorders. *Nature.* 2017;542:433–438. [doi:10.1038/nature21062](https://doi.org/10.1038/nature21062). PMC6016744.
21. Wright CF, et al. Making new genetic diagnoses with old data: iterative reanalysis in 1,133 DDD families. *Genet Med.* 2018. [doi:10.1038/gim.2017.246](https://doi.org/10.1038/gim.2017.246).
22. Cohen ASA, et al. Genomic answers for children: dynamic analyses of >1000 pediatric rare disease genomes. *Genet Med.* 2022;24:1336–1348. [doi:10.1016/j.gim.2022.02.007](https://doi.org/10.1016/j.gim.2022.02.007).
23. Hiatt SM, et al. Long-read genome sequencing for the molecular diagnosis of neurodevelopmental disorders. *HGG Adv.* 2021;2:100023. [doi:10.1016/j.xhgg.2021.100023](https://doi.org/10.1016/j.xhgg.2021.100023).
24. Hiatt SM, et al. Long-read genome sequencing and variant reanalysis increase diagnostic yield in neurodevelopmental disorders. *Genome Res.* 2024;34:1747–1762. [doi:10.1101/gr.279227.124](https://doi.org/10.1101/gr.279227.124).
25. AlAbdi L, et al. Beyond the exome: utility of long-read whole genome sequencing in exome-negative autosomal recessive diseases. *Genome Med.* 2023;15:114. [doi:10.1186/s13073-023-01270-8](https://doi.org/10.1186/s13073-023-01270-8).
26. Cheung WA, et al. Direct haplotype-resolved 5-base HiFi sequencing for genome-wide profiling of hypermethylation outliers in a rare disease cohort. *Nat Commun.* 2023;14:3090. [doi:10.1038/s41467-023-38782-1](https://doi.org/10.1038/s41467-023-38782-1).
27. Groza C, et al. Pangenome graphs improve the analysis of structural variants in rare genetic diseases. *Nat Commun.* 2024;15:657. [doi:10.1038/s41467-024-44980-2](https://doi.org/10.1038/s41467-024-44980-2).
28. Gorzynski JE, et al. Ultrarapid nanopore genome sequencing in a critical care setting. *N Engl J Med.* 2022;386:700–702. [doi:10.1056/NEJMc2112090](https://doi.org/10.1056/NEJMc2112090).
29. Laurie S, et al. The RD-Connect genome-phenome analysis platform. *Hum Mutat.* 2022;43:717–733. [doi:10.1002/humu.24353](https://doi.org/10.1002/humu.24353).
30. NHS Genomic AI Network. Genomic AI Tool Log V3. January 2025. [PDF](https://genomicainetwork.nhs.uk/wp-content/uploads/2025/01/AI-Tools-NHS-Deployments-Log-V3.pdf).

### Callers, aligners, SV/STR tools

31. Poplin R, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nat Biotechnol.* 2018;36:983–987. [doi:10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235). PMID 30247488.
32. McKenna A, et al. The Genome Analysis Toolkit. *Genome Res.* 2010;20:1297–1303. [doi:10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110).
33. Yun T, et al. Accurate, scalable cohort variant calls using DeepVariant and GLnexus. *Bioinformatics.* 2021;36:5582–5589. [doi:10.1093/bioinformatics/btaa1081](https://doi.org/10.1093/bioinformatics/btaa1081).
34. Shafin K, et al. Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads. *Nat Methods.* 2021. [doi:10.1038/s41592-021-01299-w](https://doi.org/10.1038/s41592-021-01299-w).
35. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics.* 2018;34:3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191).
36. Chen X, et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics.* 2016;32:1220–1222. [doi:10.1093/bioinformatics/btv710](https://doi.org/10.1093/bioinformatics/btv710).
37. Roller E, et al. Canvas: versatile and scalable detection of copy number variants. *Bioinformatics.* 2016;32:2375–2377. [doi:10.1093/bioinformatics/btw163](https://doi.org/10.1093/bioinformatics/btw163).
38. Dolzhenko E, et al. ExpansionHunter: a sequence-graph-based tool to analyze variation in short tandem repeat regions. *Bioinformatics.* 2019;35:4754–4756. PMC6853681.
39. Dolzhenko E, et al. Characterization and visualization of tandem repeats at genome scale (TRGT). *Nat Biotechnol.* 2024;42:1606–1614. [doi:10.1038/s41587-023-02057-3](https://doi.org/10.1038/s41587-023-02057-3).
40. Chen X, et al. Comprehensive SMN1 and SMN2 profiling for spinal muscular atrophy analysis using long-read PacBio HiFi sequencing (Paraphase). *Am J Hum Genet.* 2023;110:240–250. [doi:10.1016/j.ajhg.2023.01.001](https://doi.org/10.1016/j.ajhg.2023.01.001).
41. Holt JM, et al. HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing. *Bioinformatics.* 2024;40:btae042. [doi:10.1093/bioinformatics/btae042](https://doi.org/10.1093/bioinformatics/btae042).
42. Kim S, et al. Strelka2: fast and accurate calling of germline and somatic variants. *Nat Methods.* 2018;15:591–594. [doi:10.1038/s41592-018-0051-x](https://doi.org/10.1038/s41592-018-0051-x). PMID 30013048.
43. Kosugi S, et al. Comprehensive evaluation of structural variation detection algorithms for whole genome sequencing. *Genome Biol.* 2019;20:117. [doi:10.1186/s13059-019-1720-5](https://doi.org/10.1186/s13059-019-1720-5).

### RNA-seq

44. Kremer LS, et al. Genetic diagnosis of Mendelian disorders via RNA sequencing. *Nat Commun.* 2017;8:15824. [doi:10.1038/ncomms15824](https://doi.org/10.1038/ncomms15824).
45. Yépez VA, et al. Clinical implementation of RNA sequencing for Mendelian disease diagnostics. *Genome Med.* 2022;14:38. [doi:10.1186/s13073-022-01019-9](https://doi.org/10.1186/s13073-022-01019-9).
46. Frésard L, et al. Identification of rare-disease genes using blood transcriptome sequencing and large control cohorts. *Nat Med.* 2019;25:911–919. [doi:10.1038/s41591-019-0457-8](https://doi.org/10.1038/s41591-019-0457-8).
47. Cummings BB, et al. Improving genetic diagnosis in Mendelian disease with transcriptome sequencing. *Sci Transl Med.* 2017;9:eaal5209. [doi:10.1126/scitranslmed.aal5209](https://doi.org/10.1126/scitranslmed.aal5209).
48. Jaganathan K, et al. Predicting splicing from primary sequence with deep learning. *Cell.* 2019;176:535–548.e24. [doi:10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015).

### Rankers and VEPs

49. Mao D, et al. AI-MARRVEL — a knowledge-driven AI system for diagnosing Mendelian disorders. *NEJM AI.* 2024;1(5). [doi:10.1056/aioa2300009](https://doi.org/10.1056/aioa2300009). PMC11221788.
50. Ruzicka J, et al. Clinical evaluation of an AI system for streamlined variant interpretation in genetic testing (DiagAI). medRxiv. 2025. [doi:10.1101/2025.02.04.25321641](https://doi.org/10.1101/2025.02.04.25321641). **Preprint.**
51. Smedley D, et al. Next-generation diagnostics and disease-gene discovery with the Exomiser. *Nat Protoc.* 2015;10:2004–2015.
52. Smedley D, et al. A whole-genome analysis framework for effective identification of pathogenic regulatory variants: Genomiser. *Am J Hum Genet.* 2016;99:595–606.
53. Robinson PN, et al. Interpretable clinical genomics with a likelihood ratio paradigm (LIRICAL). *Am J Hum Genet.* 2020;107:403–417. [doi:10.1016/j.ajhg.2020.08.017](https://doi.org/10.1016/j.ajhg.2020.08.017).
54. Birgmeier J, et al. AMELIE speeds Mendelian diagnosis by matching patient phenotype and genotype to primary literature. *Sci Transl Med.* 2020;12:eaau9113. [doi:10.1126/scitranslmed.aau9113](https://doi.org/10.1126/scitranslmed.aau9113).
55. Li Q, et al. Xrare: a machine learning method jointly modeling phenotypes and genetic evidence for rare disease diagnosis. *Genet Med.* 2019;21:2126–2134. [doi:10.1038/s41436-019-0439-8](https://doi.org/10.1038/s41436-019-0439-8).
56. Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science.* 2023;381:eadg7492. [doi:10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492). PMID 37733863.
57. Ioannidis NM, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet.* 2016;99:877–885. [doi:10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016).
58. Frazer J, et al. Disease variant prediction with deep generative models of evolutionary data (EVE). *Nature.* 2021. [doi:10.1038/s41586-021-04043-8](https://doi.org/10.1038/s41586-021-04043-8).
59. Li S, et al. CAPICE: a computational method for Consequence-Agnostic Pathogenicity Interpretation of Clinical Exome variations. *Genome Med.* 2020.
60. Yuan X, et al. Evaluation of phenotype-driven gene prioritization methods for Mendelian diseases (LIRICAL/AMELIE complementarity on DDD). *Brief Bioinform.* 2022. [doi:10.1093/bib/bbac019](https://doi.org/10.1093/bib/bbac019).

### Population and knowledge resources

61. Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans (gnomAD). *Nature.* 2020;581:434–443. [doi:10.1038/s41586-020-2308-7](https://doi.org/10.1038/s41586-020-2308-7).
62. Chen S, et al. A genomic mutational constraint map using variation in 76,156 human genomes. *Nature.* 2024;625:92–100. [doi:10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0). PMID 38057664.
63. gnomAD v4.1 data. [https://gnomad.broadinstitute.org/downloads](https://gnomad.broadinstitute.org/downloads) (730,947 exomes; 76,215 genomes).
64. Landrum MJ, et al. ClinVar. *Nucleic Acids Res.* (ongoing releases). [https://www.ncbi.nlm.nih.gov/clinvar/](https://www.ncbi.nlm.nih.gov/clinvar/).
65. Firth HV, et al. DECIPHER. *Am J Hum Genet.* 2009;84:524–533.
66. Philippakis AA, et al. The Matchmaker Exchange: a platform for rare disease gene discovery. *Hum Mutat.* 2015;36:915–921. [doi:10.1002/humu.22858](https://doi.org/10.1002/humu.22858).
67. Morales J, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature.* 2022.
68. McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biol.* 2016;17:122. [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4).
69. Liao W-W, et al. A draft human pangenome reference. *Nature.* 2023;617:312–324. [doi:10.1038/s41586-023-05896-x](https://doi.org/10.1038/s41586-023-05896-x).
70. Nurk S, et al. The complete sequence of a human genome (T2T-CHM13). *Science.* 2022;376:44–53. [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987).
71. Mahmoud M, et al. Utility of long-read sequencing for All of Us. *Nat Commun.* 2024;15:837. [doi:10.1038/s41467-024-44804-3](https://doi.org/10.1038/s41467-024-44804-3).
72. Gustafson JA, et al. High-coverage nanopore sequencing of samples from the 1000 Genomes Project. *Genome Res.* 2024;34:2061–2073. [doi:10.1101/gr.279273.124](https://doi.org/10.1101/gr.279273.124).
73. CoLoRSdb. [https://colorsdb.org/](https://colorsdb.org/).
74. Genomics England PanelApp. [https://panelapp.genomicsengland.co.uk/](https://panelapp.genomicsengland.co.uk/).

### Additional method papers cited via the anchors

75. Poplin R, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv. 2018 (GATK4 HaplotypeCaller).
76. DePristo MA, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nat Genet.* 2011;43:491–498.
77. Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics.* 2009/2013 (BWA-MEM).
78. Rausch T, et al. DELLY: structural variant discovery by integrated paired-end and split-read analysis. *Bioinformatics.* 2012;28:i333–i339.
79. Layer RM, et al. LUMPY: a probabilistic framework for structural variant discovery. *Genome Biol.* 2014;15:R84.
80. Cameron DL, et al. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. *Genome Res.* 2017;27:2050–2060.
81. Talevich E, et al. CNVkit: genome-wide copy number from targeted DNA sequencing. *PLoS Comput Biol.* 2016;12:e1004873.
82. Minoche AE, et al. ClinSV: clinical grade structural and copy number variant detection from whole genome sequencing data. *Genome Med.* 2021.
83. Smolka M, et al. Detection of mosaic and population-level structural variants with Sniffles2. *Nat Biotechnol.* 2024.
84. Jiang T, et al. Long-read-based human genomic structural variation detection with cuteSV. *Genome Biol.* 2020;21:189.
85. Yépez VA, et al. Detection of aberrant gene expression events in RNA sequencing data (DROP protocol). *Nat Protoc.* 2021.
86. Bonder MJ / Mertes C, et al. OUTRIDER and FRASER method papers as required by DROP citation rules (Gagneur lab).
87. Cibulskis K, et al. Sensitive detection of somatic point mutations (MuTect). *Nat Biotechnol.* 2013;31:213–219.
88. Boycott KM, Dyment DA, Innes AM. Unsolved recognizable patterns of human malformation. *Am J Med Genet C.* 2018;178:382–386.
89. Wright CF, et al. (DDD 2015 diagnostic workflow commentary). *Lancet.* 2015 (27% yield on 1,133).
90. Sanford Kobayashi E, et al. PacBio HiFi in undiagnosed paediatric rare disease (IKBKG). 2022 (cited by Del Gobbo and Boycott, 2025).
91. Miller DE, et al. Targeted ONT for missing second variants (Werner syndrome 8/9). 2022 (cited by Del Gobbo and Boycott, 2025).
92. Wenger AM, et al. Accurate circular consensus long-read sequencing (PacBio HiFi). *Nat Biotechnol.* 2019.
93. Kirsche M, et al. Jasmine and Iris: population-scale structural variant comparison. *Nat Methods.* 2023;20:408–417.
94. Geoffroy V, et al. AnnotSV. *Bioinformatics.* 2018;34:3572–3574.
95. Martin M, Ebert P, Marschall T. WhatsHap. *Methods Mol Biol.* 2023;2590:127–138.
96. Köhler S, et al. The Human Phenotype Ontology. *Nucleic Acids Res.* (ongoing).
97. DiStefano MT, et al. The Gene Curation Coalition (GenCC). *Genet Med.* 2022.
98. Liu X, et al. dbNSFP. *Hum Mutat.* (ongoing).
99. Rentzsch P, et al. CADD: predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Res.* 2019.
100. Gao H, et al. The landscape of tolerated genetic variation in humans and primates (PrimateAI-3D). *Science.* 2023.
101. Brandes N, et al. Genome-wide prediction of disease variants with a deep protein language model (ESM-1b). *Nat Genet.* 2023.
102. Feng BJ. PERCH: a unified framework for disease gene prioritization (BayesDel). *Hum Mutat.* 2017.
103. Wagner N, et al. Benchmarking splice variant prediction algorithms using massively parallel splicing assays. *Genome Biol.* 2023. [doi:10.1186/s13059-023-03144-z](https://doi.org/10.1186/s13059-023-03144-z).
104. Ebbert MTW, et al. Systematic analysis of dark and camouflaged genes. *Genome Biol.* 2019;20:97.
105. Chaisson MJP, et al. Multi-platform discovery of haplotype-resolved structural variation in human genomes. *Nat Commun.* 2019;10:1784.
106. Logsdon GA, Vollger MR, Eichler EE. Long-read human genome sequencing and its applications. *Nat Rev Genet.* 2020;21:597–614.
107. Cacheiro P, et al. Evaluating the calling performance of a rare disease NGS panel for single nucleotide and copy number variants. *Mol Diagn Ther.* 2017;21:303–313. [doi:10.1007/s40291-017-0268-x](https://doi.org/10.1007/s40291-017-0268-x). PMID 28290094.

---

## Search reproducibility

Queries, Europe PMC hit counts, and PubMed hit counts in §2 were obtained on 17 October 2026 via:

- `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=…&format=json&pageSize=1&resultType=idlist`
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=…&retmode=json&rettype=count`

No VCFs, no patient data, and no licensed mutation-database dumps were downloaded for this note.
