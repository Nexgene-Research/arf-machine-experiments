# Variant-Aware Drug Repurposing in BUB1B-Associated Mosaic Variegated Aneuploidy: A Mechanism-Grounded Hypothesis Portfolio

**İlsu Duman, Özge Beğli, Hatice Mutlu and Kutsal Baran Ozkurt\***

Nexgene AI, https://nexgene.ai

\*Correspondence: kutsal.ozkurt@nexgene.ai

**MVA Hackathon 2026, Track 2 (Drug Repurposing).** Assessment date: 22 October 2026.

> **Disclaimer.** The candidates presented are repurposing *hypotheses for follow-up*, not clinical recommendations or evidence of efficacy, and this work does not constitute medical advice. All clinical features and protein labels are user-supplied observations, not independently verified medical findings. No patient genomic files, cells, or clinical records were analysed or uploaded. No medicine below should be given to any patient on the basis of this document.

---

## Abstract

**Background.** Mosaic variegated aneuploidy (MVA, MIM 257300) is an ultra-rare autosomal-recessive disorder in which biallelic loss of the spindle-assembly-checkpoint (SAC) protein BUBR1 (gene *BUB1B*) produces constitutional chromosome mis-segregation, growth failure and cancer predisposition. We were asked to characterise the disrupted biology of one child and to nominate approved medicines plausibly acting on that biology.

**Case.** Reported biallelic *BUB1B* variants `p.Leu737*` (nonsense) and `p.Asn1002Lys` (missense), with embryonal rhabdomyosarcoma (RMS), low muscle mass and failure to thrive; additional reported features include nephrocalcinosis, short stature, prematurity or small-for-gestational-age birth, and parental recurrent miscarriage.

**Methods.** We performed a bounded, bidirectional evidence synthesis across the literature, variant databases and regulatory labels, graded every claim on six independent axes, and assembled a curated directed causal graph linking variants, protein states, mitotic machinery, cellular phenotypes and candidate drugs. A deterministic, standard-library pipeline exports the machine-readable ledgers and independently validates schema, referential integrity, graph continuity and file hashes.

**Results.** We frame four independent, mechanism-based strategies and a portfolio of approved candidates: **(1) nonsense suppression and BUBR1 restoration** by premature-stop read-through of the `p.Leu737*` allele (azithromycin, amikacin, gentamicin); **(2) proteostasis and metabolic rescue** to raise or stabilise functional BUBR1 (arimoclomol, metformin, trametinib, prescription niacin ER, hydroxychloroquine); **(3) aneuploidy-associated cellular stress modulation through PP2A reactivation** (fingolimod); and **(4) tumour-directed options for RMS** (larotrectinib, vincristine, olaparib). No candidate completes a demonstrated exact-genotype causal chain from drug to faithful cell division; each entry carries its weakest link, decisive test and kill criterion.

**Conclusion.** No approved medicine meets the bar for molecular correction or constitutional disease modification of the reported genotype, and we report that honestly: Tier 1 is empty. The contribution is therefore methodological as much as clinical: a transparently graded, fully falsifiable hypothesis portfolio, an auditable causal graph, and an allele-resolved experimental program whose results can promote, demote or kill any candidate. In doing so we convert an open-ended question of "which drug?" into a finite, prioritised set of decisive experiments.

**Keywords:** mosaic variegated aneuploidy; *BUB1B* / BubR1; spindle-assembly checkpoint; drug repurposing; nonsense read-through; PP2A-B56.

---

## Plain-language summary

A child has two changes in both copies of the *BUB1B* gene. This gene makes BUBR1, a protein that acts as a quality-control officer during cell division: it holds division until every chromosome is correctly gripped, and it recruits a counter-enzyme (PP2A) to steady those grips. With BUBR1 broken, cells distribute their chromosomes unevenly, which explains the growth problems and the high cancer risk seen in this disease, including the muscle tumour (rhabdomyosarcoma) this child has.

We asked a focused question: is there an **already-approved medicine** that could plausibly act on this broken machinery? Our honest answer is that **no approved drug has been shown to fix this exact problem**, so we do not claim a cure. Instead, we lay out four separate, testable repair ideas, each tied to a specific broken step, and for every medicine we state its single most important weakness, the one experiment that would prove or kill it, and the safety problems that would matter for a fragile child. One change in the gene (a premature "stop") could in principle be "read through" by certain antibiotics; the other (a single letter swap) could, if the protein is merely unstable, be supported by drugs that slow its breakdown. A separate track covers the tumour itself. The point of this report is not to say "use these drugs." It is to show the reasoning transparently, show the evidence behind every claim, flag the dangers, and hand over a clear experimental plan that could turn a hypothesis into an answer.

---

## 1. Introduction

BUBR1 is a core component of the spindle-assembly checkpoint: it builds the mitotic checkpoint complex (MCC: BUBR1, BUB3, MAD2 and CDC20) that restrains the anaphase-promoting complex (APC/C^CDC20) until all kinetochores are correctly attached, and it recruits the PP2A-B56 phosphatase to kinetochores to stabilise microtubule attachments. Biallelic *BUB1B* mutations cause MVA1, a constitutional aneuploidy syndrome with growth failure and a markedly elevated risk of childhood cancers including rhabdomyosarcoma and Wilms tumour. ClinGen records the *BUB1B* and MVA relationship as *definitive* and autosomal recessive [E05], and the founding clinical series already linked truncating and missense *BUB1B* variants to embryonal RMS [E06].

The therapeutic question for Track 2 is deliberately narrow: given the disrupted biology, does an *existing, approved* medicine plausibly act on the problem? We treat this as a hypothesis-generation task under strict evidential discipline. Two principles recur throughout. First, **more BUBR1 protein is not the goal; more *functional* BUBR1 at the right place and time is**: raising total abundance could increase useful protein, inert protein, or an actively interfering product [E07, E08]. Second, **a graph connection or an expression change is not a rescue**: measured effects in model systems must be distinguished from inferred and hypothesised links, and from any claim about this child.

Our contribution is to make that discipline operational and to answer the question as fully as the evidence honestly allows. It is threefold: a residue-level account of *how* each allele breaks BUBR1 function; a portfolio of approved medicines mapped, by mechanism rather than by keyword or indication, to the specific broken step each could plausibly address; and a reproducible, auditable apparatus (a typed causal graph, six-axis evidence grades and a deterministic validator) that lets any reader retrace, challenge or extend every inference we draw. The result is not a single recommended drug but something more useful to a clinician-scientist: a ranked, falsifiable map from genotype to candidate intervention, with each edge's strength and each candidate's decisive experiment made explicit.

## 2. Case and variant characterisation

**Reference framework.** MANE Select transcript `NM_001211.6` and `ENST00000287598.11`; protein `NP_001202.5`, 1,050 aa. The DNA descriptions below are *public matching records*, not a reconstruction of the child's genotype; the reference calculation is reproducible with the reference-normalisation script (Annex B) [E03].

| Allele | Public evidence | Mechanistic inference | Critical unknown |
|---|---|---|---|
| Reported `p.Leu737*` (nonsense) | ClinVar VCV000533901 maps a matching allele to `c.2210T>G`; aggregate pathogenic or likely-pathogenic assertions [E01] | Exon-17 termination consistent with nonsense-mediated decay (NMD); if expressed, residues 1 to 736 remain (N-terminal CDC20, KNL1 and BUB3 regions, Lys250, KARD, Lys668) but the C-terminal kinase-like domain (residues 766 to 1050) is lost [E03, E04] | Exact patient nucleotide; RNA decay; product abundance; residual binding or checkpoint activity |
| Reported `p.Asn1002Lys` (missense) | ClinVar VCV004600147 maps *one possible* allele to `c.3006T>A`; single-submitter variant of uncertain significance [E02] | Asn1002 lies in the C-terminal kinase-like domain; possible effect on structure, stability, interactions or splicing [E04] | Exact nucleotide identity, phase, abundance, folding, localisation and useful function |

No formal pathogenicity classification is assigned from protein labels, and a public missense record must not be silently substituted for another nucleotide allele giving the same consequence. Phase (cis versus trans) and any additional allele or copy-number change remain unresolved; these determine whether the alleles are truly compound heterozygous. Population annotations are retained in the supplementary records (Annex B); current gnomAD counts were not independently verified, and missing frequency fields are not evidence of absence.

**Allele to strategy map.** The two alleles motivate distinct, non-exclusive strategies. The nonsense `p.Leu737*` allele is the substrate for read-through (Strategy 1). The missense `p.Asn1002Lys` allele, if unstable but recoverable, motivates proteostasis and metabolic support (Strategy 2) and localisation-specific PP2A reactivation (Strategy 3). The tumour branch (Strategy 4) is separate and conditional on tumour-specific biomarkers, not on constitutional *BUB1B* status.

## 3. Methods

**Evidence curation and grading.** Claims were extracted manually from primary literature, variant databases and regulatory labels, with access recorded as record-read, abstract-only, or full-text-sections (a read-scope label, not a quality score). Each evidence item carries six independent grades (rigour, causal directness, variant and model relevance, exposure feasibility, reproducibility and clinical-context relevance), plus a study-family identifier so that duplicate papers, corrections and registry records are not counted as independent replication. Unknown values carry explicit reasons; no values are fabricated. The full ledger is provided as a machine-readable evidence table (Annex B), and every cited item is reproduced with its grade in the References and evidence ledger below.

**How to read a candidate (evidence legend).** Every candidate below is reported in the same disciplined shape so the reasoning is inspectable rather than asserted:

- **Proposed mechanism.** The specific broken step the drug is meant to act on, not a vague claim to "improve cell health".
- **Evidence.** What was actually shown, and crucially in *what system*. Graph links and this report use four honesty states: *measured* (shown directly in a model), *inferred* (a reasonable deduction from measured facts), *hypothesised* (a proposed but untested link), and *contradicted* (evidence points the other way). A reported number is given only with its model, exposure and species.
- **Key safety note.** The toxicity or feasibility problem that would matter for a fragile child.
- **Weakest link** and **kill criterion** (in the ledgers). The single assumption most likely to be wrong, and the result that would end the hypothesis.

Two inferences we refuse to make are that *abundance equals function*, and that a *short path in the graph equals a rescue*. A conservative complete-chain check therefore counts every non-measured edge as unproven and counts each context transfer (species, tissue, genotype) separately; net biological sign is never computed by multiplying edge signs.

**Causal graph.** We represent the biology as a directed graph of typed nodes (variant, allele, transcript, protein state or function, complex or interaction, process, cell or tissue phenotype, drug, metabolite, target, safety outcome, model, evidence source) and typed edges carrying an evidence state, species, tissue, genotype, exposure and explicit transfer assumptions. The node, edge and path tables and a readable graph rendering are provided as supplementary files (Annex B).

**Reproducibility.** The reference-normalisation script recomputes allele sequence context from the public records; the build script exports the ledgers deterministically; the validation script independently re-reads the final files and checks schemas, identifiers, referential integrity, graph continuity, counts and SHA-256 hashes (validation record, Annex B). The pipeline uses only the Python 3.10+ standard library, with no network, credentials, patient data or stochastic estimators. Validation covers package consistency; it cannot validate patient diagnosis, phase, tissue expression, drug interactions or efficacy.

**Tiering.** Tier 1 requires relevant functional rescue, plausible exposure and no invalidating safety gap. Tier 2 requires bounded, decisive conditions (here, the biomarker-gated tumour branch). Tier 3 is exploratory. Tier X preserves rejected or deprioritised mechanisms. No weighted efficacy probability is computed, and order within a tier is a research sequence, not a probability of benefit.

## 4. Results: four mechanism-based strategies

Each candidate below is an approved medicine mapped, by mechanism rather than by indication or keyword, to a specific disrupted step of BUBR1 biology. For every drug we give its proposed mechanism, the strength *and system* of its evidence, the decisive safety note, and, in the ledgers, the single experiment that would confirm or kill it. The four strategies are independent and non-competing: a clinician could pursue any one without committing to the others, and each closes with its weakest link and a falsifying criterion so that the portfolio is read as a set of testable bets, not a list of endorsements. The candidate table (Annex B) carries every requested column, and consolidated safety appears in Section 7.

### 4.1 Strategy 1: Nonsense suppression and BUBR1 restoration

*Rationale.* If the `p.Leu737*` transcript escapes NMD, aminoglycoside and macrolide antibiotics can promote ribosomal read-through of the premature stop, potentially yielding a longer, partly functional BUBR1. Read-through efficiency is strongly dependent on the exact stop codon and flanking sequence, and a longer protein is not itself proof of useful BUBR1 [E18].

| Drug | Proposed mechanism | Evidence (system and state) | Key safety note |
|---|---|---|---|
| **Gentamicin** | Premature-stop read-through | Dose-dependent read-through in mammalian reporters, strongly stop-context-dependent (mouse NIH3T3; *measured*, not *BUB1B*) [E18] | Nephrotoxicity and irreversible auditory or vestibular toxicity; monitor; not a feasible chronic exposure [E30] |
| **Amikacin** | Premature-stop read-through, possibly lower nephrotoxicity | Aminoglycoside class effect; preclinical only (*inferred* from class, no *BUB1B* data) [E18] | Screen `MT-RNR1` (m.1555A>G, m.1494C>T, m.1095T>C) before use, because carriers are at high risk of aminoglycoside-induced hearing loss (CPIC; *measured* pharmacogenomic risk) [E53] |
| **Azithromycin** | Read-through at roughly 100-fold lower concentration than gentamicin | Macrolide read-through in reporter and *Apc^Min/+* models (*measured* in non-*BUB1B* systems; microbiome caveat) [E46] | Unconfirmed for *BUB1B* `p.Leu737*`; QT prolongation, hepatic and GI effects [E50] |

**Evidence strength:** preclinical reporter and animal data only; *no study used the `p.Leu737*` allele*, and the incorporated amino acid and its function are unmeasured. **Weakest link:** antibacterial exposure does not establish cytosolic ribosomal read-through *in this gene*; the macrolide window (azithromycin's low effective concentration) is the most translationally attractive but must first be confirmed for this allele. **Kill criterion:** no full-length or functional BUBR1 produced from the patient allele, or no checkpoint restoration despite a longer product.

### 4.2 Strategy 2: Proteostasis and metabolic rescue

*Rationale.* Agents reported to raise BUBR1 levels or support folding and half-life of the `p.Asn1002Lys` product. This strategy is justified only if an expressed allele is unstable but *functionally recoverable*; raising abundance of a defective product is not a goal.

| Drug | Proposed mechanism | Evidence (system and state) | Key safety note |
|---|---|---|---|
| **Arimoclomol** | Amplifies the heat-shock response (HSF1 to HSP70 and HSP90); the chaperone rise could refold or stabilise the misfolded `p.Asn1002Lys` BUBR1 | Chaperone induction shown in NPC models (*measured* in NPC, not *BUB1B*); a **newly approved pediatric orphan drug** for Niemann-Pick type C with a low adverse-effect burden [E47, E51]. Caveats: it is a heat-shock-response *co-inducer* rather than an HSP90 inhibitor, its autophagy and HSP70 contribution is uncertain, and a second NPC study found no HSP70 induction (*contradicted*) [E48] | Pediatric approval with a low side-effect profile; creatinine may rise without GFR decline; OCT2 interactions [E51] |
| **Metformin** | Spindle-assembly-checkpoint modulation; the motivating signal is a **dose- and cell-type-dependent rise** in BUBR1 | Direction is not uniform: at lower concentrations BUBR1 rose (markedly in some cell types), whereas higher concentrations and other cell types lowered it, within a glucose-dependent genotoxic-stress context (*measured*, bidirectional) [E23] | Pediatric use acceptable; lactic-acidosis risk with renal impairment or hypoxia (relevant given nephrocalcinosis); B12 depletion [E34] |
| **Trametinib** | MEK inhibition may raise functional BUBR1 (theoretical inversion) | BubR1 *deficiency* enhances MEK and ERK under microtubule stress, motivating the inverse hypothesis (*hypothesised*; not shown for trametinib or these alleles) [E19, E39] | Cardiomyopathy, ocular, haemorrhagic, skin, GI and lung toxicity; developmental signalling effects [E39] |
| **Prescription niacin ER** *(complementary)* | NAD⁺ precursor; NAD-driven, SIRT2-dependent deacetylation (K668) opposes BUBR1 degradation | NMN raised NAD and BUBR1 abundance in mice and MEFs (*measured*); niacin replenished NAD in mitochondrial myopathy (*measured*, adults); niacin is not NMN or NR [E09, E10, E49] | Flushing, hepatic, glucose and urate effects; pediatric efficacy (age 16 and under) not established; not interchangeable with supplements [E52] |
| **Hydroxychloroquine** | Blocks autophagic and lysosomal BUBR1 degradation, raising BUBR1 levels | BUBR1 is an autophagy substrate; HCQ raised BUBR1 **only in senescent cells** by blocking degradation, while the autophagy inducer everolimus lowered it (*measured*, senescence-specific) [E54] | Caution: neuromuscular weakness overlaps with the child's low muscle mass; retinal and cardiac toxicity |

**Evidence strength:** mixed; the SIRT2 and NAD stability axis is mechanistically the best-grounded (*measured* stability data [E09, E10]) but is not demonstrated for either patient allele, and arimoclomol's chaperone effect is directly *contradicted* by a second NPC study [E48]. **Weakest link:** for every agent, whether raising or stabilising BUBR1 yields *functional* checkpoint activity for this genotype rather than more inert or interfering protein; the SIRT2 axis also deacetylates K250 (checkpoint silencing), so direction and timing matter [E11]. **Kill criterion:** abundance rises without functional checkpoint rescue, localisation is wrong, or the stabilised product sequesters partners.

### 4.3 Strategy 3: Aneuploidy-associated cellular stress modulation through PP2A reactivation

*Rationale.* BUBR1 recruits PP2A-B56 to kinetochores to stabilise attachments; in MVA-derived cells, artificially targeting PP2A-B56 to kinetochores rescued chromosome alignment [E12, E14]. A drug that activates PP2A could, *if* its activity were localised via B56 at the kinetochore, improve alignment.

| Drug | Proposed mechanism | Evidence (system and state) | Key safety note |
|---|---|---|---|
| **Fingolimod** | Binds and activates PP2A (via SET or I2PP2A); if localised via B56 at the kinetochore, may improve chromosome alignment | Bulk PP2A activation (*measured*); a kinetochore B56 tether rescued alignment in MVA cells in a *separate* genetic experiment (*measured*); the drug-to-localisation link is *hypothesised* [E12, E36] | Approved from age 10 (MS); first-dose bradyarrhythmia, macular oedema, infection and lymphocyte effects, malignancy concern; additive immunosuppression with antineoplastics [E36] |

**Evidence strength:** the *mechanism model* is the strongest in the portfolio (localised B56 rescue is *measured*), but the *drug link* is the weakest (fingolimod acts in bulk). **Weakest link:** bulk PP2A activation is not kinetochore-localised B56 recruitment. **Kill criterion:** PP2A activation without B56-dependent kinetochore localisation, or without fewer segregation errors in patient-genotype cells.

### 4.4 Strategy 4: Tumour-directed options for RMS

*Rationale.* Standard-of-care plus biomarker-matched candidates for the reported rhabdomyosarcoma. These are tumour-branch hypotheses, conditional on tumour molecular findings, not constitutional MVA rescue.

| Drug | Proposed mechanism | Evidence (system and state) | Key safety note |
|---|---|---|---|
| **Larotrectinib** | Inhibits NTRK fusion kinases; conditional on a somatic NTRK fusion | FDA-approved, tumour-agnostic for confirmed NTRK-fusion solid tumours including pediatric sarcomas (*measured* clinical; conditional on fusion) [E40] | Fusion testing required; CNS, hepatic, skeletal-fracture and haematologic effects; approval restricted to metastatic or no-alternative settings [E40] |
| **Vincristine** | Standard RMS chemotherapy (vinca alkaloid); benefit reported in RMS regimens | Established RMS backbone (clinical standard-of-care; *measured* at population level, not MVA-specific) | Neurotoxicity; dose adjustment important in MVA given chromosomal fragility |
| **Olaparib** | PARP inhibition; hypothesis that BUBR1 loss or chromosomal instability raises reliance on DNA repair | PARP inhibitors are HRR or BRCA-selected; a BUB1B-high tumour context was linked to DNA-repair-associated resistance (*hypothesised* for RMS; no indication) [E24, E55] | Myelosuppression, MDS or AML, pneumonitis, VTE; no RMS or *BUB1B* indication; chromosomal instability is not homologous-recombination deficiency [E55] |

**Evidence strength:** actionable *only* on a matching biomarker. Larotrectinib is well supported *given a confirmed fusion*; vincristine is supportive standard-of-care rather than mechanism-specific; olaparib's synthetic-lethality rationale is *unproven* in RMS. **Weakest link:** none is established by the supplied features: larotrectinib needs a documented fusion and dependency; olaparib needs a demonstrated repair dependency. **Kill criterion:** absent target or fusion, no tumour dependency, or equal toxicity to affected normal tissue.

## 5. Causal knowledge graph

The machine-readable graph (Annex B) makes every intermediate step and transfer inspectable.

- **Best-supported short model path:** kinetochore B56 tether, leading to localised PP2A-B56 and then improved alignment (measured in published MVA models; the extension to *fewer errors per completed division in this genotype* is the decisive missing link) [E12, E14].
- **Conditional longer paths:** read-through to longer BUBR1 to function (Strategy 1); chaperone or metabolic stabilisation to abundance to function (Strategy 2); bulk PP2A activation to localised recruitment (Strategy 3). Every such path exposes hypothesised links and context transfers.
- **Rejected or interrupted paths:** bulk PP2A activation assumed to equal localised B56 recruitment; preserved aneuploid-cell counts assumed to equal restored tissue output; generic mTOR inhibition assumed to improve RMS outcome. A short graph distance does not rescue these claims.

## 6. Proposed experimental program

The program is allele-first and endpoint-disciplined; the primary endpoint is the **fraction of completed divisions producing new segregation errors**, with death, arrest and proliferation reported separately.

1. **Confirm exact alleles and phase**, then map residual function: allele-specific RNA and NMD, N-terminal and C-terminal protein detection, and matched-expression complementation (wild-type versus each allele) with live-cell fidelity.
2. **Discriminate checkpoint failure from attachment failure:** MCC, CDC20 and APC/C measures and PP2A-B56 localisation; genetic B56 tethering as a positive control.
3. **Only then test the compatible medication branch** (Strategies 1, 2 and 3), with exposure-response, target-dependence and division-rate controls, comparing normal and tumour cells.
4. **Tumour branch (Strategy 4)** proceeds on its own biomarkers: confirm an NTRK fusion before larotrectinib; establish a repair dependency before olaparib.

## 7. Safety considerations

Safety is not a footnote in this disease: the patient is a growing child with chromosomal fragility and an active cancer, so an intervention that is routine in adults can be disqualifying here. Five cross-cutting safety matters govern the whole portfolio.

1. **The tumour-versus-normal window is the overriding hazard.** Several interventions that would help a dividing normal cell could equally protect the rhabdomyosarcoma. Restoring mitotic fitness rescued tumour-like cells under MAPK stress [E19], and raising BUBR1 was tumour-promoting in a bladder-cancer context [E24]. **No candidate should advance without a demonstrated normal-versus-tumour differential.**
2. **Pharmacogenomic screening before aminoglycosides (Strategy 1).** `MT-RNR1` variants (m.1555A>G, m.1494C>T, m.1095T>C) confer high risk of permanent aminoglycoside-induced hearing loss; CPIC recommends avoiding aminoglycosides in carriers, and there is *no* evidence that amikacin is intra-class safer [E53]. Genotyping is mandatory, and gentamicin carries irreversible oto- and nephrotoxicity at the exposures discussed [E30].
3. **Pediatric approval status is not pediatric efficacy.** Several labels restrict or omit the child's age group: prescription niacin ER (age 16 and under not established) [E52], and temsirolimus and others (pediatric effectiveness not established) [E32]. Approval for one indication does not transfer an exposure or a safety margin to constitutional MVA.
4. **Organ-specific overlaps with the phenotype.** Hydroxychloroquine's neuromuscular toxicity overlaps the child's low muscle mass; fingolimod causes first-dose bradyarrhythmia and macular oedema and adds immunosuppression on top of chemotherapy [E36]; metformin risks lactic acidosis with any renal compromise (relevant given nephrocalcinosis) [E34]; somatropin is contraindicated with active malignancy [E38].
5. **"Do no harm" dominates a speculative benefit.** Where a mechanism is only *hypothesised* (trametinib's inversion, olaparib's synthetic-lethality rationale) or directly *contradicted* (arimoclomol's chaperone induction [E48]), the known toxicity outweighs the unproven benefit until a decisive test is passed. Per-drug warnings are carried in the candidate table (Annex B) and in the regulatory labels cited in the References and evidence ledger below.

## 8. Significance and alignment with the evaluation criteria

The deliverable was built to be judged on the four Track 2 criteria. The table maps each criterion to concrete, verifiable artefacts; the paragraphs explain the approach behind them.

| Criterion (weight) | What we deliver | Where to verify |
|---|---|---|
| **Scientific rigour (35%)** | MANE-anchored, reproducible variant characterisation; every candidate tied to the mechanism by a typed causal chain, six-axis grades, a stated weakest link and a falsifying kill criterion; measured, inferred, hypothesised and contradicted kept distinct; Tier 1 left empty rather than overclaimed | Sections 2 to 4; evidence ledger, claim audit and validation record (Annex B) |
| **Potential impact (25%)** | An allele-resolved program with one disciplined primary endpoint that can promote or kill each candidate; an immediately actionable diagnostic safety gate (`MT-RNR1`) and one named decisive experiment | Sections 6 and 7; path table (Annex B) |
| **Innovation (25%)** | A typed causal graph with conservative complete-chain accounting; an allele-class to drug-class mapping; a read-through hypothesis coupled to a pharmacogenomic safety gate | Section 5; node and edge tables and the appendices (Annex B) |
| **Scalability (15%)** | A deterministic, standard-library pipeline with an independent validator; a reusable grading schema and allele-class framework | Scripts and validation record (Annex B) |

**Scientific rigour (35%).** The mechanism is characterised, not asserted. The reference framework is MANE-anchored and recomputed from public records (Section 2); the truncation's consequence is given at residue resolution (loss of the kinase-like domain, residues 766 to 1050; retention of the KARD, Lys250 and Lys668); and the missense allele is reported honestly as a single-submitter VUS with unresolved phase. Each repurposing candidate is tied to that mechanism through a typed causal chain and graded on six independent axes, with its single weakest link and a falsifying kill criterion stated in the open. We separate what is *measured* in a model from what is *inferred*, *hypothesised*, or actively *contradicted*, and we correct our own leads where the literature demands it: the bidirectional metformin effect [E23], arimoclomol as a heat-shock co-inducer with a contradicting NPC study [E48], trametinib's reliance on inverting a BubR1-loss finding [E19], and fingolimod's bulk-versus-localised PP2A gap. Leaving Tier 1 empty is itself a rigour statement: we decline to label any approved drug a disease modifier on present evidence.

**Potential impact (25%).** The output is designed to change understanding and care if independently validated. The experimental program (Section 6) is allele-first and fixes one disciplined primary endpoint, the fraction of completed divisions that produce *new* segregation errors, so a single experiment can promote or kill a candidate for this child and for other MVA families. Two contributions are actionable now: an `MT-RNR1` pharmacogenomic screen that prevents permanent aminoglycoside ototoxicity before any read-through attempt [E53], and the identification of the one decisive missing experiment (whether localised PP2A-B56 recruitment reduces errors per division in this genotype [E12, E14]) whose result would move the whole field, not just this case.

**Innovation (25%).** The method is the novelty. We encode the biology as a typed causal graph whose edges carry honesty states and explicit context transfers, and we score complete chains conservatively: every non-measured edge counts as unproven and every species, tissue or genotype transfer is counted separately, never multiplied into a false net sign. On top of this sits an allele-class to drug-class mapping (nonsense to read-through; unstable missense to proteostasis and metabolic support; localisation defect to PP2A reactivation) and a read-through hypothesis coupled to a pharmacogenomic safety gate, so that an efficacy idea and its safety test are proposed as a single unit rather than bolted on afterwards.

**Scalability (15%).** The package is a deterministic, standard-library pipeline: JSON curation inputs build the CSV ledgers and the graph, and an independent validator re-checks schema, referential integrity, graph continuity and SHA-256 hashes with no network or credentials. The six-axis grading schema, the study-family provenance model and the allele-class framework are gene-agnostic and disease-agnostic, so the same machinery extends to other variants, other conditions and larger evidence corpora (for example, DepMap or PRISM sensitivity as a discovery aid), not just this one case.

## 9. Discussion

The strongest surviving hypothesis is that a demonstrated residual-scaffold or attachment deficit may be experimentally rescuable, including by localisation-specific PP2A-B56 bypass, but no approved drug presently completes that chain for the exact alleles. The portfolio's value is that each strategy is independently falsifiable and tied to a specific disrupted step, so a single well-designed experiment can promote or kill a candidate. The read-through strategy (Strategy 1) is attractive because it targets the nonsense allele directly and because azithromycin's low effective concentration offers a plausible therapeutic window; the PP2A strategy (Strategy 3) has the strongest mechanistic model evidence but the weakest drug-to-localisation link; the proteostasis and metabolic strategy (Strategy 2) has the best-grounded stability biology but the least certain direction of effect. The tumour branch (Strategy 4) is the most clinically actionable *if and only if* a matching biomarker is found.

## 10. Limitations

This is a bounded, bidirectional synthesis, not an exhaustive systematic review. Important primary papers were read in accessible sections; some sources were abstract-only and some corrections or supplements remain inaccessible. No raw omics reanalysis, patient testing, exhaustive patent search, or complete approved-drug-library screen was performed. Negative search results are coverage-limited, not evidence of no effect. All clinical observations remain user-reported; public database alleles with matching protein consequences are reconciliation candidates, not proven patient genotypes. The most consequential gaps are the exact-allele functional gap, the absence of a normal-versus-tumour safety window, and tissue and exposure transfers; additional citations to the same first mechanistic step would not close them (claim-scope audit, Annex B).

## 11. Conclusion

No approved medicine is supported for molecular correction or constitutional disease modification of the reported `p.Leu737*` and `p.Asn1002Lys` genotype; **Tier 1 is empty.** We regard this disciplined negative as a result in its own right, and we pair it with the constructive deliverable the question actually needs: a transparently graded portfolio of four mechanism-based repurposing strategies, an auditable genotype-to-candidate causal map, and an allele-resolved experimental program whose next result could promote a candidate out of hypothesis and into evidence. The bar we hold for that promotion is explicit: reproducible exact-genotype rescue of faithful divisions at achievable exposure, with demonstrated target dependence, no selection artefact and a convincing normal-versus-tumour safety window. Meeting it would benefit not only this child but every family living with MVA, and the method that gets there is designed to generalise well beyond this one case.

## Funding

This work was supported by the EuroHPC Joint Undertaking, which awarded access to the MareNostrum 5 supercomputer hosted by the Barcelona Supercomputing Center (BSC-CNS), Spain, under awarded project ID EHPC-AIF-2026PG01-1286.

## AI disclosure

All analyses were carried out on Nexgene AI's Arf Machine research platform, on a local workstation and on MareNostrum 5. The language-model component was ASA-W1, Nexgene AI's medical-reasoning model, running locally.

- No patient data were used to train or fine-tune any model, and no patient data were exported from the local environment.
- Public database queries contained only a gene symbol, a single variant notation or a publication identifier; they did not contain identifiable patient information.
- Data-use conditions were followed. Controlled source data and derived patient files will be deleted, with all copies, within 30 days of hackathon close, and the organisers will be notified.
- All AI-assisted interpretation was reviewed by the authors.

## Acknowledgements

We are deeply grateful to the child and family who chose to share their genome with the research community. We thank the MVA Society, Sage Bionetworks, Hugging Face, BEACON and the Wilhelm Foundation for organising and supporting the MVA Hackathon 2026. We acknowledge the Barcelona Supercomputing Center for operating MareNostrum 5.

## References and evidence ledger

Every claim above is traceable to a graded evidence item. Grades are shown as **rig** (rigour), **dir** (causal directness) and **rel** (variant and model relevance); access is the read-scope (record for a database or label record read, abstract for abstract only, full for full-text sections). The complete machine-readable ledger with all six axes, transfer assumptions and conflicts is provided as the evidence table (Annex B).

| ID | Source and identifier (DOI or accession) | Type and access | Grade (rig/dir/rel) | What it shows or caution |
|---|---|---|---|---|
| E01 | ClinVar VCV000533901 | variant DB, record | n.a. / n.a. / mod | Public `c.2210T>G` (`p.Leu737Ter`) record; does **not** verify the patient allele or phase. |
| E02 | ClinVar VCV004600147 | variant DB, record | n.a. / n.a. / mod | `c.3006T>A` single-submitter VUS; matching protein notation is not a proven nucleotide allele. |
| E03 | MANE NM_001211.6, ENST00000287598.11, NP_001202.5 | reference computation, record | high / n.a. / mod | 1,050-aa reference; places public alleles in exons 17 and 23 (reproducible). |
| E04 | UniProt O60566 | protein annotation, record | n.a. / n.a. / mod | Kinase-like domain, residues 766 to 1050; a 1 to 736 product keeps the N-terminus, Lys250 and Lys668 and loses the C-terminus. |
| E05 | ClinGen BUB1B and MVA1 (CGGV assertion) | gene-disease curation, record | n.a. / n.a. / low | *Definitive*, autosomal-recessive gene-disease validity (not a diagnosis). |
| E06 | Hanks et al., *Nat Genet* 2004 (DOI:10.1038/ng1449) | primary observational, abstract | mod / mod / mod | Founding MVA series; 2 of 5 families had embryonal RMS. |
| E07 | Suijkerbuijk et al., *Cancer Res* 2010 (DOI:10.1158/0008-5472.CAN-09-4319) | primary experiment, full | mod / high / mod | MVA BUBR1 mutants reduce stability and checkpoint; **abundance is not function**. |
| E08 | *J Clin Invest* 2019 (DOI:10.1172/JCI126863; corr. 10.1172/JCI144781) | primary experiment, full | mod / high / mod | BubR1 allele combinations give heterogeneous phenotypes; not scaled by total protein. |
| E09 | North et al., *EMBO J* 2014 (DOI:10.15252/embj.201386907) | primary experiment, full | mod / high / low | SIRT2 deacetylates BUBR1 Lys668, opposing degradation (stability axis). |
| E10 | North et al., *EMBO J* 2014 (same study family) | primary experiment, full | mod / high / low | NMN raised NAD and BUBR1 abundance in mice and MEFs; NMN is not an approved drug. |
| E11 | *Biochem Biophys Res Commun* 2014 (DOI:10.1016/j.bbrc.2014.09.128) | primary experiment, abstract | n.a. / mod / low | SIRT2 also deacetylates Lys250 (checkpoint silencing); direction and timing matter. |
| E12 | Xu et al., *Biol Open* 2013 (DOI:10.1242/bio.20134051) | primary experiment, full | mod / high / mod | Kinetochore-targeted PP2A-B56 rescued alignment in MVA-derived cells (key model path). |
| E14 | Kruse et al., *J Cell Sci* 2013 (DOI:10.1242/jcs.122481) | primary experiment, abstract | n.a. / mod / mod | BUBR1 and PP2A-B56 support stable attachments; artificial targeting bypasses the deficit. |
| E15 | Baker et al., *Nature* 2011 (DOI:10.1038/nature10600) | primary experiment, full | mod / high / mod | INK-ATTAC senescent-cell clearance delayed some tissues; genetic, not dasatinib. |
| E16 | *Cell Rep* 2013 (DOI:10.1016/j.celrep.2013.03.028) | primary and unpublished, full | low / low / mod | Preliminary rapamycin in BubR1 mice gave phenotype rates similar to control (weak negative). |
| E17 | Zhu et al., *Aging Cell* 2015 (DOI:10.1111/acel.12344) | primary experiment, abstract | n.a. / mod / low | Dasatinib and quercetin senolysis is cell-type dependent; no *BUB1B* window shown. |
| E18 | Floquet et al., *PLoS Genet* 2012 (DOI:10.1371/journal.pgen.1002608) | primary experiment, full | mod / high / low | Gentamicin read-through is strongly stop-context dependent; longer is not functional. |
| E19 | Herman et al., *PNAS* 2022 (DOI:10.1073/pnas.2208255119) | primary experiment, full | mod / high / low | RAS and MAPK raise BUBR1 dependence; MEK inhibition rescued alignment after BUBR1 loss (tumour-fitness caution). |
| E21 | ARST0921 (DOI:10.1200/JCO.19.00576) | randomised phase 2, full | high / mod / low | Temsirolimus-containing relapse arm favoured, but against a bevacizumab arm, not chemotherapy alone. |
| E22 | ARST1431 (DOI:10.1016/S1470-2045(24)00255-9) | randomised phase 3, full | high / high / low | Adding temsirolimus did **not** improve EFS in intermediate-risk RMS (66.8% versus 64.8%). |
| E23 | *Food Chem Toxicol* 2022 (DOI:10.1016/j.fct.2022.113129) | primary experiment, abstract | n.a. / mod / low | High-dose, glucose-dependent genotoxic stress; BUBR1 **bidirectional** across cell types. |
| E24 | *Oncogene* 2021 (DOI:10.1038/s41388-021-02021-y) | primary experiment, full | mod / mod / low | High BUB1B drove DNA-repair-associated resistance (tumour-promoting direction); CIN is not HRD. |
| E28 | EMA EPAR Translarna (non-renewal 2025) | regulatory, record | n.a. / n.a. / low | Ataluren EU authorisation not renewed; not currently EU-approved. |
| E30 | Gentamicin label (DailyMed 63413469…) | regulatory, full | n.a. / n.a. / low | Nephrotoxicity, irreversible oto- and vestibular toxicity, neuromuscular blockade. |
| E31 | Sirolimus (Rapamune) label (DailyMed 3275b824…) | regulatory, full | n.a. / n.a. / low | Infection and malignancy boxed warning; no congenital-MVA benefit. |
| E32 | Temsirolimus (Torisel) label (DailyMed 95b7dc92…) | regulatory, full | n.a. / n.a. / low | RCC only; pediatric effectiveness not established; not RMS-approved. |
| E33 | Dasatinib (Sprycel) label (DailyMed 4764f37b…) | regulatory, full | n.a. / n.a. / low | No senolysis indication; marrow suppression, PAH, growth and bone effects. |
| E34 | Metformin IR label (DailyMed 2387a618…) | regulatory, full | n.a. / n.a. / low | Lactic-acidosis risk with renal impairment or hypoxia; B12 depletion. |
| E35 | Potassium citrate ER (Urocit-K) label (DailyMed 72cdea1b…) | regulatory, full | n.a. / mod / low | Pivotal trials non-randomised; hyperkalaemia and GI risk; no congenital-nephrocalcinosis benefit. |
| E36 | Fingolimod (Gilenya) label (DailyMed cc9e1c8c…) | regulatory, full | n.a. / n.a. / low | S1P modulator; bradyarrhythmia, macular oedema, malignancy; no kinetochore PP2A action shown. |
| E37 | Bortezomib (Velcade) label (DailyMed 1521d321…) | regulatory, full | n.a. / n.a. / low | Broad proteasome block is not selective BUBR1 stabilisation; neuropathy, cytopenias. |
| E38 | Somatropin (Norditropin) label (DailyMed 1058e17c…) | regulatory, full | n.a. / n.a. / low | Active-malignancy contraindication; second-neoplasm risk in survivors. |
| E39 | Trametinib (Mekinist) label (DailyMed 0002ad27…) | regulatory, full | n.a. / n.a. / low | BRAF indication only; cardiomyopathy, ocular, haemorrhage; no MVA or RAS-RMS approval. |
| E40 | Larotrectinib FDA label (210861s012, 2025) | regulatory, full | n.a. / mod / low | Tumour-agnostic for confirmed NTRK fusions; approval is eligibility-restricted. |
| E43 | *Nat Commun* 2026 (DOI:10.1038/s41467-026-70521-0) | primary experiment, full | mod / high / low | Rapamycin preserved fly neuroblasts in *bub3* and *rod* loss (not *BUB1B*; no error reduction). |
| E46 | Caspi 2016 (DOI:10.1007/s00109-015-1364-1); Semba 2023 (DOI:10.1016/j.biopha.2023.114968) | primary experiments, abstract | n.a. / mod / low | Macrolide read-through in reporters and *Apc^Min/+*; not *BUB1B*; microbiome caveat. |
| E47 | Kirkegaard et al., *Sci Transl Med* 2016 (DOI:10.1126/scitranslmed.aad9823) | primary experiment, full | mod / mod / low | Arimoclomol altered HSF1 and HSP70 in NPC; clinical mechanism remains unknown. |
| E48 | Pipalia et al., *J Lipid Res* 2021 (DOI:10.1016/j.jlr.2021.100114) | primary experiment, full | mod / high / low | **Contradiction:** no HSP70 induction or storage correction with arimoclomol up to 0.5 mM. |
| E49 | Pirinen et al., *Cell Metab* 2020 (DOI:10.1016/j.cmet.2020.04.008) | clinical intervention, abstract | low / mod / low | Niacin replenished NAD in adult mitochondrial myopathy; not *BUB1B*, not pediatric. |
| E50 | Azithromycin FDA label (050784s040, 2026) | regulatory, full | n.a. / n.a. / low | QT and torsades, hepatic, GI; infant pyloric-stenosis warning; no read-through approval. |
| E51 | Arimoclomol label (DailyMed 5feffc0e…) | regulatory, full | n.a. / n.a. / low | Approved for NPC with miglustat from age 2; creatinine rise without GFR decline; OCT2. |
| E52 | Prescription niacin ER label (DailyMed 52ad4a6a…) | regulatory, full | n.a. / n.a. / low | Lipid indication; age 16 and under not established; not interchangeable with NMN or NR supplements. |
| E53 | CPIC MT-RNR1 and aminoglycosides 2022 (DOI:10.1002/cpt.2309) | guideline, full | high / high / low | Avoid aminoglycosides in MT-RNR1 carriers; **no** intra-class safety ranking (screen first). |
| E54 | Myrianthopoulos et al., *Redox Biol* 2023 (DOI:10.1016/j.redox.2023.102701) | primary experiment, full | mod / high / low | BUBR1 is an autophagy substrate; HCQ raised it only in senescent cells (degradation-blocking). |
| E55 | Olaparib (Lynparza) label (DailyMed 741ff3e3…) | regulatory, full | n.a. / low / low | HRR or BRCA-selected; MDS or AML, pneumonitis, VTE; no RMS or *BUB1B* indication. |

## Annex A. Considered and deprioritised candidates

Six candidates were considered and deprioritised; they are retained in the rejection ledger (Annex B) with explicit reopening rules. Potassium citrate ER is limited to a conditional renal-biochemistry branch [E35]. Sirolimus rests on a fly cell-preservation signal [E43] with no human *BUB1B* tissue-output benefit and a preliminary negative mouse report [E16, E31]. Dasatinib has no demonstrated selective senolytic window in this genotype [E15, E17, E33]. Temsirolimus was negative in a later randomised RMS trial [E21, E22, E32]. Bortezomib provides global proteasome inhibition rather than selective BUBR1 stabilisation [E37]. Somatropin has no qualifying endocrine indication and is contraindicated by active malignancy [E38]. Separate unapproved research probes (NMN, genetic B56 tethering, research MEK and Aurora inhibitors, INK-ATTAC) are listed among the research tools (Annex B); ataluren is outside the current EU-approved shortlist [E28].

## Annex B. Supplementary files and reproducibility

All quantitative claims, the causal graph and the candidate decisions are provided as machine-readable supplementary files, so that any reader can retrace or extend the analysis. No repository has been published; the rare-feature report itself merits privacy review before any public release, and full-text PDFs, copyrighted article text and restricted datasets are deliberately not redistributed.

**Supplementary files.**

| File | Contents |
|---|---|
| `data/evidence.csv` | Claim-level evidence ledger with identifiers, access status, six-axis grades and study-family provenance |
| `data/nodes.csv`, `data/edges.csv`, `data/paths.csv`, `data/graph.md` | The directed causal graph (typed nodes and edges) and its readable rendering |
| `data/candidates.csv`, `data/rejections.csv` | Full candidate columns with strategy grouping; the deprioritised candidates with reopening rules |
| `data/research_tools.csv` | Separate list of unapproved or genetic research probes |
| `docs/claim_audit.md` | Load-bearing claim-scope audit |
| `docs/appendices.md` | Phenotype audit, variant detail, decision axes and the full decision matrix |
| `data/search_log.csv`, `data/source_coverage.csv`, `data/datasets.csv` | Discovery queries, source retrievals and explicitly unsearched resources |
| `data/reference_context.json`, `data/population_annotations.json` | Reproducible allele sequence context and population annotations |
| `reference/` | Public RefSeq, Ensembl, UniProt, ClinVar and HPO snapshots |
| `curation/` | Manual extraction JSON that generates the ledgers |
| `scripts/` | Normalisation, build, independent validation and PDF rendering |
| `data/validation.json`, `data/package_summary.json`, `docs/DELIVERABLES.md`, `environment.txt` | Verification record, package summary, deliverable contract and environment |

**Reproduction.** With Python 3.10+ and its standard library, run `python scripts/normalize_references.py`, then `python scripts/build_package.py`, then `python scripts/validate_package.py`. The CSV ledgers and `data/graph.md` are generated from `curation/`, so corrections belong in the JSON inputs, not in the CSV files. All exports are hash-verified in `data/validation.json`. No network, credentials or patient data are required.
