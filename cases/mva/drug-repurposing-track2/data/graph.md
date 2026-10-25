# Readable directed causal graph

Generated from the exported graph. Bracket labels are edge IDs, evidence state and claim IDs. 
Measured means reported in the stated model, not demonstrated in this child. Inferred, hypothesized and contradicted links remain visible.
Conservative unproven-link count includes every non-measured edge; model-to-case transfers are counted separately.

## P01: PP2A localization bypass rescues alignment in published models

Objective B; short_supported_model; 2 links; 0 unproven links; 2 specified transfer events.

```text
Kinetochore-targeted B56 construct --RECRUITS_TO_KINETOCHORE [G19; measured; E12;E14]--> Kinetochore PP2A-B56
Kinetochore PP2A-B56 --IMPROVES_CHROMOSOME_ALIGNMENT [G20; measured; E12;E14]--> Chromosome alignment
```

**Weakest link:** No approved localization drug; exact-allele fidelity still untested
**Decisive test:** Tether B56 in authorized/isogenic cells, then score completed divisions
**Competing effects:** Recruitment must remain kinetochore- and mitosis-specific; excess stabilization may preserve wrong attachments
**Transfers:** genotype:other_MVA_to_case;tissue:published_lines_to_case

## P02: PP2A recruitment to lower micronuclei through fewer new errors

Objective B; long_conditional; 4 links; 2 unproven links; 3 specified transfer events.

```text
Kinetochore-targeted B56 construct --RECRUITS_TO_KINETOCHORE [G19; measured; E12;E14]--> Kinetochore PP2A-B56
Kinetochore PP2A-B56 --IMPROVES_CHROMOSOME_ALIGNMENT [G20; measured; E12;E14]--> Chromosome alignment
Chromosome alignment --REDUCES_SEGREGATION_ERRORS [G21; hypothesized; E12;E13]--> New chromosome-segregation errors
New chromosome-segregation errors --CAN_GENERATE_MICRONUCLEI [G54; inferred; E19;E20]--> Micronuclei
```

**Weakest link:** G21: alignment to faithful division in this genotype
**Decisive test:** Blinded live-cell lineage tracking plus chromosome measurement
**Competing effects:** Do not use reduced micronuclei from slowed proliferation or cell loss as rescue
**Transfers:** genotype:models_to_case;genotype:other_MVA_to_case;tissue:published_lines_to_case

## P03: SIRT2-stability-functional-rescue chain

Objective A; long_conditional; 5 links; 2 unproven links; 3 specified transfer events.

```text
SIRT2 deacetylase --DEACETYLATES [G08; measured; E09]--> BUBR1 K668 deacetylation
BUBR1 K668 deacetylation --INCREASES_PROTEIN_STABILITY [G09; measured; E09]--> BUBR1 protein stability
BUBR1 protein stability --INCREASES_ABUNDANCE [G10; measured; E09;E10]--> Total BUBR1 abundance
Total BUBR1 abundance --INCREASES_FUNCTIONAL_ACTIVITY [G11; hypothesized; E07;E09;E20]--> Functional BUBR1 checkpoint activity
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
```

**Weakest link:** G11: increased mutant abundance to useful checkpoint function
**Decisive test:** Allele-calibrated expression and SIRT2-dependence matrix
**Competing effects:** K250 checkpoint-silencing effect competes; tumor fitness may increase
**Transfers:** genotype:WT_to_case;species:mouse_to_human;tissue:cell_lines_to_case

## P04: Gentamicin readthrough to functional fidelity

Objective A; long_conditional; 4 links; 3 unproven links; 3 specified transfer events.

```text
Gentamicin injection --PROMOTES_READTHROUGH [G15; measured; E18]--> Premature-stop readthrough
Premature-stop readthrough --PRODUCES_LONGER_BUBR1 [G16; hypothesized; E01;E03;E18]--> Readthrough-derived longer BUBR1
Readthrough-derived longer BUBR1 --RESTORES_FUNCTIONAL_ACTIVITY [G17; hypothesized; E07;E18]--> Functional BUBR1 checkpoint activity
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
```

**Weakest link:** G16/G17 plus exposure and toxicity
**Decisive test:** Native-context transcript assay, proteomics of incorporated residue and function
**Competing effects:** NMD limits substrate; longer protein can remain defective; Tier X clinical translation
**Transfers:** genotype:WT_to_case;genotype:reporters_to_BUB1B;species:mouse_to_human

## P05: Expressed truncated scaffold to attachment rescue

Objective B; short_conditional; 3 links; 2 unproven links; 3 specified transfer events.

```text
BUBR1 residues 1-736 --RETAINS_USEFUL_RECRUITMENT [G22; hypothesized; E04;E12]--> Kinetochore PP2A-B56
Kinetochore PP2A-B56 --IMPROVES_CHROMOSOME_ALIGNMENT [G20; measured; E12;E14]--> Chromosome alignment
Chromosome alignment --REDUCES_SEGREGATION_ERRORS [G21; hypothesized; E12;E13]--> New chromosome-segregation errors
```

**Weakest link:** G22: endogenous p.Leu737Ter product recruits B56
**Decisive test:** N/C terminal immunoblots, co-IP/localization and matched-expression rescue
**Competing effects:** Do not infer exact truncation behaves like overexpressed 1-730 construct
**Transfers:** genotype:engineered_fragment_to_case;genotype:other_MVA_to_case;tissue:published_lines_to_case

## P06: Genetic senescent-cell clearance to tissue improvement

Objective C; short_conditional; 2 links; 0 unproven links; 3 specified transfer events.

```text
INK-ATTAC/AP20187 clearance --REMOVES_MARKED_CELLS [G30; measured; E15]--> p16-high senescent-cell burden
p16-high senescent-cell burden --CONTRIBUTES_TO_DYSFUNCTION [G31; measured; E15;E16]--> Tissue dysfunction
```

**Weakest link:** Mouse-to-child and transgene-to-drug transfer
**Decisive test:** Measure actual senescence burden and tissue function before drug screen
**Competing effects:** Cell removal may damage progenitors; protective senescence and tumor immunity matter
**Transfers:** genotype:hypomorph_to_case;species:mouse_to_human;tissue:aged_mouse_to_child

## P07: Dasatinib senolysis to tissue benefit

Objective C; short_conditional; 2 links; 1 unproven links; 5 specified transfer events.

```text
Dasatinib --SELECTIVELY_REMOVES_CASE_SENESCENT_CELLS [G32; hypothesized; E15;E17]--> p16-high senescent-cell burden
p16-high senescent-cell burden --CONTRIBUTES_TO_DYSFUNCTION [G31; measured; E15;E16]--> Tissue dysfunction
```

**Weakest link:** G32: selective killing at a tolerable exposure
**Decisive test:** Paired senescent/nonsenescent case-cell survival and tissue assays
**Competing effects:** Mixed species path explicit; one mouse genetic result cannot complete a human drug chain
**Transfers:** genotype:hypomorph_to_case;genotype:other_senescence_to_case;species:mouse_to_human;tissue:aged_mouse_to_child;tissue:aging_models_to_child

## P08: Sirolimus mTOR-SASP-tissue chain

Objective C; short_conditional; 3 links; 2 unproven links; 4 specified transfer events.

```text
Sirolimus --INHIBITS_TARGET [G34; measured; E31]--> mTOR signaling
mTOR signaling --DRIVES_CASE_INFLAMMATORY_OUTPUT [G35; hypothesized; E15;E16]--> Senescence-associated inflammatory output
Senescence-associated inflammatory output --CAUSES_CASE_TISSUE_DYSFUNCTION [G36; hypothesized; E15;E16]--> Tissue dysfunction
```

**Weakest link:** G35/G36 absent in case; preliminary negative E16
**Decisive test:** Verify mTOR-dependent SASP then functional organoid rescue
**Competing effects:** Growth inhibition, immune harm and tumor effects may outweigh reduced inflammatory markers
**Transfers:** genotype:general_patients_to_case;genotype:hypomorph_to_case;species:mouse_to_human;tissue:aging_models_to_child

## P09: Biochemically selected renal citrate support

Objective C; short_conditional; 3 links; 2 unproven links; 1 specified transfer events.

```text
Potassium citrate ER --INCREASES_URINARY_CITRATE_AND_PH [G41; measured; E35]--> Urinary citrate/alkali response
Urinary citrate/alkali response --REDUCES_CALCIUM_OXALATE_SUPERSATURATION [G42; inferred; E35]--> Calcium-crystal supersaturation
Calcium-crystal supersaturation --CONTRIBUTES_TO_MINERAL_DEPOSITION [G43; hypothesized; E35]--> Nephrocalcinosis
```

**Weakest link:** Correct nephrocalcinosis cause and mineral chemistry
**Decisive test:** Urine/serum metabolic evaluation; ex vivo crystallization model
**Competing effects:** Alkalinization can favor calcium-phosphate precipitation; potassium retention risk
**Transfers:** genotype:biochemical_subgroup_to_case

## P10: Fusion-selected tumor dependency

Objective E; short_conditional; 2 links; 1 unproven links; 2 specified transfer events.

```text
Larotrectinib --INHIBITS_TARGET [G45; measured; E40]--> Oncogenic NTRK fusion signaling
Oncogenic NTRK fusion signaling --DRIVES_TUMOR_VIABILITY_IF_ONCOGENIC_FUSION [G46; inferred; E40]--> Rhabdomyosarcoma viability/growth
```

**Weakest link:** NTRK fusion/dependency absent from supplied data
**Decisive test:** Tumor RNA fusion confirmation and matched normal/tumor response
**Competing effects:** No BUB1B-specific window; cannot infer from TRK expression or copy gain
**Transfers:** genotype:fusion_positive_to_untyped_RMS;tissue:pan_tumor_to_RMS

## P11: Trametinib-MAPK-Hec1-attachment-fidelity bypass

Objective B; long_conditional; 4 links; 2 unproven links; 5 specified transfer events.

```text
Trametinib --INHIBITS_TARGET [G52; measured; E39]--> MEK-ERK signaling
MEK-ERK signaling --PROMOTES_KINETOCHORE_HYPERPHOSPHORYLATION [G50; measured; E19]--> Kinetochore Hec1 hyperphosphorylation
Kinetochore Hec1 hyperphosphorylation --IMPAIRS_ATTACHMENT_AND_ALIGNMENT [G51; inferred; E19]--> Chromosome alignment
Chromosome alignment --REDUCES_SEGREGATION_ERRORS [G21; hypothesized; E12;E13]--> New chromosome-segregation errors
```

**Weakest link:** Drug substitution, active MAPK in normal cells and functional fidelity
**Decisive test:** Compare trametinib with U0126 and genetic MEK suppression in isogenic cells
**Competing effects:** Feedback/nonlinear attachment balance; reducing mitotic stress can protect tumor
**Transfers:** genotype:BRAF_to_case;genotype:RAS_KD_to_case;genotype:other_MVA_to_case;tissue:ARPE19_to_case;tissue:approved_tumors_to_case

## P12: Fingolimod to localized B56 rescue

Objective B; short_rejected; 3 links; 2 unproven links; 3 specified transfer events.

```text
Fingolimod --RESTORES_LOCALIZED_COMPLEX [G26; hypothesized; E12;E36]--> Kinetochore PP2A-B56
Kinetochore PP2A-B56 --IMPROVES_CHROMOSOME_ALIGNMENT [G20; measured; E12;E14]--> Chromosome alignment
Chromosome alignment --REDUCES_SEGREGATION_ERRORS [G21; hypothesized; E12;E13]--> New chromosome-segregation errors
```

**Weakest link:** G26 searched but no supporting localized rescue found
**Decisive test:** Do not prioritize without localization-specific epistasis evidence
**Competing effects:** S1P immunosuppression is real while the desired bridge is missing
**Transfers:** genotype:other_MVA_to_case;genotype:unrelated_pharmacology_to_case;tissue:published_lines_to_case

## P13: Temsirolimus as generic intermediate-risk RMS enhancer

Objective E; short_rejected; 1 links; 1 unproven links; 1 specified transfer events.

```text
Temsirolimus --IMPROVES_UNSELECTED_INTERMEDIATE_RISK_RMS_EFS [G48; contradicted; E22]--> RMS event-free survival
```

**Weakest link:** Randomized clinical result contradicts generalized benefit claim
**Decisive test:** Only a new biomarker-defined dependency could reopen research
**Competing effects:** Earlier relapse trial and later frontline trial are distinct; no pooled monotherapy claim
**Transfers:** genotype:trial_population_to_case

## P14: Checkpoint-to-micronuclear-inflammation-tissue path

Objective C; long_conditional; 5 links; 5 unproven links; 5 specified transfer events.

```text
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
New chromosome-segregation errors --CAN_GENERATE_MICRONUCLEI [G54; inferred; E19;E20]--> Micronuclei
Micronuclei --ACTIVATES_CASE_CGAS_STING [G55; hypothesized; E15;E19]--> cGAS-STING inflammatory signaling
cGAS-STING inflammatory signaling --DRIVES_CASE_INFLAMMATION [G56; hypothesized; E15;E16]--> Senescence-associated inflammatory output
Senescence-associated inflammatory output --CAUSES_CASE_TISSUE_DYSFUNCTION [G36; hypothesized; E15;E16]--> Tissue dysfunction
```

**Weakest link:** cGAS-STING and tissue-causality links unverified
**Decisive test:** Measure rupture/cGAS/STING and perturb genetically only if activated
**Competing effects:** Antitumor immune suppression and developmental effects; too many gaps for approved-drug tiering
**Transfers:** genotype:WT_to_case;genotype:hypomorph_to_case;genotype:models_to_case;species:mouse_to_human;tissue:aging_models_to_child

## P15: Azithromycin readthrough to functional BUBR1

Objective A; long_conditional; 4 links; 3 unproven links; 5 specified transfer events.

```text
Azithromycin --PROMOTES_READTHROUGH [G62; measured; E46]--> Premature-stop readthrough
Premature-stop readthrough --PRODUCES_LONGER_BUBR1 [G16; hypothesized; E01;E03;E18]--> Readthrough-derived longer BUBR1
Readthrough-derived longer BUBR1 --RESTORES_FUNCTIONAL_ACTIVITY [G17; hypothesized; E07;E18]--> Functional BUBR1 checkpoint activity
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
```

**Weakest link:** Exact allele readthrough and safe intracellular exposure
**Decisive test:** Native-transcript allele-specific rescue assay
**Competing effects:** Antibacterial/microbiome effects cannot count as readthrough; QT and tumor interaction screen
**Transfers:** genotype:APC_to_BUB1B;genotype:WT_to_case;genotype:reporters_to_BUB1B;species:mouse_to_human;tissue:intestine_to_case

## P16: Arimoclomol-chaperone-stability-fidelity chain

Objective A; long_conditional; 5 links; 3 unproven links; 4 specified transfer events.

```text
Arimoclomol --AMPLIFIES_STRESS_RESPONSE [G63; measured; E47;E48;E51]--> HSF1/HSP response
HSF1/HSP response --STABILIZES_USEFUL_MUTANT_BUBR1 [G64; hypothesized; E07;E47;E48]--> BUBR1 protein stability
BUBR1 protein stability --INCREASES_ABUNDANCE [G10; measured; E09;E10]--> Total BUBR1 abundance
Total BUBR1 abundance --INCREASES_FUNCTIONAL_ACTIVITY [G11; hypothesized; E07;E09;E20]--> Functional BUBR1 checkpoint activity
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
```

**Weakest link:** G64 and G11: mutant-specific folding/function
**Decisive test:** Matched-expression allele rescue plus HSF dependence
**Competing effects:** Chaperone induction inconsistent; can preserve malignant cells
**Transfers:** genotype:NPC1_to_BUB1B;genotype:WT_to_case;species:mouse_to_human;tissue:lysosomal_to_mitotic

## P17: Prescription niacin to NAD-SIRT2-functional rescue

Objective A; long_conditional; 7 links; 4 unproven links; 5 specified transfer events.

```text
Prescription niacin ER --INCREASES_METABOLITE [G65; inferred; E49;E52]--> NAD availability
NAD availability --SUPPORTS_PATHWAY_ACTIVITY [G07; inferred; E10]--> SIRT2 deacetylase
SIRT2 deacetylase --DEACETYLATES [G08; measured; E09]--> BUBR1 K668 deacetylation
BUBR1 K668 deacetylation --INCREASES_PROTEIN_STABILITY [G09; measured; E09]--> BUBR1 protein stability
BUBR1 protein stability --INCREASES_ABUNDANCE [G10; measured; E09;E10]--> Total BUBR1 abundance
Total BUBR1 abundance --INCREASES_FUNCTIONAL_ACTIVITY [G11; hypothesized; E07;E09;E20]--> Functional BUBR1 checkpoint activity
Functional BUBR1 checkpoint activity --REDUCES_SEGREGATION_ERRORS [G12; hypothesized; E07;E20]--> New chromosome-segregation errors
```

**Weakest link:** Selective SIRT2 engagement and functional mutant rescue
**Decisive test:** NAD then SIRT2-dependent modification, protein and fidelity
**Competing effects:** Seven links retained because approved precursor is distinct from NMN; no sign multiplication; metabolic and K250 effects compete
**Transfers:** genotype:WT_to_case;genotype:myopathy_to_case;species:mouse_to_human;tissue:adult_muscle_to_child;tissue:cell_lines_to_case

## P18: Sirolimus-related preservation of aneuploid neuroblast counts

Objective C; short_supported_model; 2 links; 1 unproven links; 3 specified transfer events.

```text
Sirolimus --MODULATES_AUTOPHAGY_VIA_TOR [G66; inferred; E43]--> Autophagy/mitophagy balance
Autophagy/mitophagy balance --PRESERVES_NEUROBLAST_COUNT [G67; measured; E43]--> Surviving aneuploid neuroblasts
```

**Weakest link:** Species/gene/exposure transfer; count is not tissue function
**Decisive test:** Human relevant tissue flux and productive progeny assay
**Competing effects:** G68 explicitly blocks assuming preserved cell count restores tissue
**Transfers:** genotype:bub3_rod_to_BUB1B;species:fly_to_human;tissue:brain_to_case

## P19: Genetic mitochondrial support restores fly brain size

Objective C; short_supported_model; 1 links; 0 unproven links; 3 specified transfer events.

```text
Hsp60/Hsp60c or Sod2/GTPx-1 expression --RESTORES_BRAIN_SIZE [G69; measured; E44]--> Brain size/progeny output
```

**Weakest link:** No approved equivalent and no case tissue match
**Decisive test:** Demonstrate mitochondrial causal lesion in justified tissue model
**Competing effects:** Could protect tumor or abnormal cells; no chromosome correction
**Transfers:** genotype:bub3_rod_to_BUB1B;species:fly_to_human;tissue:brain_to_case

## Complete edge adjacency

| Edge | From | Relationship / sign | To | State | Evidence |
|---|---|---|---|---|---|
| G01 | V737 | ANNOTATED_ON / neutral | TX | inferred | E01;E03 |
| G02 | V737 | PREDICTS_NMD / positive | NMD | inferred | E01;E03 |
| G03 | NMD | REDUCES_TRANSCRIPT / negative | RNA737 | hypothesized | E01;E03 |
| G04 | V737 | ENCODES_IF_EXPRESSED / positive | TRUNC | inferred | E03;E04 |
| G05 | V1002 | MAY_ALTER_FUNCTION_OR_STABILITY / unknown | MISS | hypothesized | E02;E07 |
| G06 | NMN | INCREASES_METABOLITE / positive | NAD | measured | E10 |
| G07 | NAD | SUPPORTS_PATHWAY_ACTIVITY / positive | SIRT2 | inferred | E10 |
| G08 | SIRT2 | DEACETYLATES / positive | K668 | measured | E09 |
| G09 | K668 | INCREASES_PROTEIN_STABILITY / positive | STABLE | measured | E09 |
| G10 | STABLE | INCREASES_ABUNDANCE / positive | ABUND | measured | E09;E10 |
| G11 | ABUND | INCREASES_FUNCTIONAL_ACTIVITY / unknown | FUNC | hypothesized | E07;E09;E20 |
| G12 | FUNC | REDUCES_SEGREGATION_ERRORS / negative | ERROR | hypothesized | E07;E20 |
| G13 | SIRT2 | DEACETYLATES / positive | K250 | measured | E11 |
| G14 | K250 | PROMOTES_CHECKPOINT_SILENCING / negative | FUNC | measured | E11 |
| G15 | GENT | PROMOTES_READTHROUGH / positive | READ | measured | E18 |
| G16 | READ | PRODUCES_LONGER_BUBR1 / positive | LONG | hypothesized | E01;E03;E18 |
| G17 | LONG | RESTORES_FUNCTIONAL_ACTIVITY / unknown | FUNC | hypothesized | E07;E18 |
| G18 | GENT | CAUSES_TOXICITY / positive | RENALTOX | measured | E30 |
| G19 | TETHER | RECRUITS_TO_KINETOCHORE / positive | PP2A | measured | E12;E14 |
| G20 | PP2A | IMPROVES_CHROMOSOME_ALIGNMENT / positive | ALIGN | measured | E12;E14 |
| G21 | ALIGN | REDUCES_SEGREGATION_ERRORS / negative | ERROR | hypothesized | E12;E13 |
| G22 | TRUNC | RETAINS_USEFUL_RECRUITMENT / unknown | PP2A | hypothesized | E04;E12 |
| G23 | ZM | INHIBITS_TARGET / negative | AURORA | measured | E13 |
| G24 | AURORA | OPPOSES_ATTACHMENT_STABILITY_IF_EXCESSIVE / negative | ALIGN | inferred | E13 |
| G25 | FING | MODULATES_TARGET / mixed | S1P | measured | E36 |
| G26 | FING | RESTORES_LOCALIZED_COMPLEX / unknown | PP2A | hypothesized | E12;E36 |
| G27 | BORTEZ | INHIBITS_TARGET / negative | PROTEASOME | measured | E37 |
| G28 | PROTEASOME | REDUCES_PROTEIN_STABILITY / negative | STABLE | measured | E09 |
| G29 | BORTEZ | CAUSES_MARROW_TOXICITY / positive | IMMUNE | measured | E37 |
| G30 | CLEAR | REMOVES_MARKED_CELLS / negative | SEN | measured | E15 |
| G31 | SEN | CONTRIBUTES_TO_DYSFUNCTION / positive | TISSUE | measured | E15;E16 |
| G32 | DASA | SELECTIVELY_REMOVES_CASE_SENESCENT_CELLS / unknown | SEN | hypothesized | E15;E17 |
| G33 | DASA | CAUSES_MARROW_TOXICITY / positive | IMMUNE | measured | E33 |
| G34 | SIRO | INHIBITS_TARGET / negative | MTOR | measured | E31 |
| G35 | MTOR | DRIVES_CASE_INFLAMMATORY_OUTPUT / unknown | SASP | hypothesized | E15;E16 |
| G36 | SASP | CAUSES_CASE_TISSUE_DYSFUNCTION / unknown | TISSUE | hypothesized | E15;E16 |
| G37 | SIRO | IMPROVES_BUBR1_HYPOMORPH_PHENOTYPES / negative | TISSUE | contradicted | E16 |
| G38 | SIRO | CAUSES_IMMUNOSUPPRESSION / positive | IMMUNE | measured | E31 |
| G39 | MET | REDUCES_GENOME_INJURY / unknown | STRESS | hypothesized | E23;E34 |
| G40 | MET | INCREASES_GENOTOXIC_STRESS / positive | STRESS | measured | E23 |
| G41 | CITRATE | INCREASES_URINARY_CITRATE_AND_PH / positive | URINE | measured | E35 |
| G42 | URINE | REDUCES_CALCIUM_OXALATE_SUPERSATURATION / negative | CRYSTAL | inferred | E35 |
| G43 | CRYSTAL | CONTRIBUTES_TO_MINERAL_DEPOSITION / positive | NEPHRO | hypothesized | E35 |
| G44 | CITRATE | CAUSES_TOXICITY_IF_POTASSIUM_RETENTION / positive | HYPERK | measured | E35 |
| G45 | LARO | INHIBITS_TARGET / negative | NTRK | measured | E40 |
| G46 | NTRK | DRIVES_TUMOR_VIABILITY_IF_ONCOGENIC_FUSION / positive | TUMOR | inferred | E40 |
| G47 | TEM | INHIBITS_TARGET / negative | MTOR | measured | E32 |
| G48 | TEM | IMPROVES_UNSELECTED_INTERMEDIATE_RISK_RMS_EFS / positive | EFS | contradicted | E22 |
| G49 | U0126 | INHIBITS_TARGET / negative | MAPK | measured | E19 |
| G50 | MAPK | PROMOTES_KINETOCHORE_HYPERPHOSPHORYLATION / positive | HEC1 | measured | E19 |
| G51 | HEC1 | IMPAIRS_ATTACHMENT_AND_ALIGNMENT / negative | ALIGN | inferred | E19 |
| G52 | TRAM | INHIBITS_TARGET / negative | MAPK | measured | E39 |
| G53 | GH | MAY_PROMOTE_ACTIVE_MALIGNANCY / positive | TUMOR | inferred | E38 |
| G54 | ERROR | CAN_GENERATE_MICRONUCLEI / positive | MICRO | inferred | E19;E20 |
| G55 | MICRO | ACTIVATES_CASE_CGAS_STING / unknown | CGAS | hypothesized | E15;E19 |
| G56 | CGAS | DRIVES_CASE_INFLAMMATION / unknown | SASP | hypothesized | E15;E16 |
| G57 | ERROR | MAY_INCREASE_TRANSFORMATION_RISK / positive | CANCERRISK | inferred | E06;E20 |
| G58 | MCC | INHIBITS_APC_CDC20 / negative | APC | measured | E07;E20 |
| G59 | TRUNC | RETAINS_USEFUL_CHECKPOINT_SCAFFOLD / unknown | MCC | hypothesized | E04;E07;E12;E20 |
| G60 | LOC | ENABLES_LOCALIZED_RECRUITMENT / positive | PP2A | measured | E12 |
| G61 | ABUND | MAY_SUPPORT_TUMOR_FITNESS / positive | TUMOR | hypothesized | E24 |
| G62 | AZI | PROMOTES_READTHROUGH / positive | READ | measured | E46 |
| G63 | ARIMO | AMPLIFIES_STRESS_RESPONSE / positive | HSF | measured | E47;E48;E51 |
| G64 | HSF | STABILIZES_USEFUL_MUTANT_BUBR1 / unknown | STABLE | hypothesized | E07;E47;E48 |
| G65 | NIACIN | INCREASES_METABOLITE / positive | NAD | inferred | E49;E52 |
| G66 | SIRO | MODULATES_AUTOPHAGY_VIA_TOR / positive | AUTO | inferred | E43 |
| G67 | AUTO | PRESERVES_NEUROBLAST_COUNT / positive | NB | measured | E43 |
| G68 | NB | RESTORES_TISSUE_OUTPUT / positive | BRAIN | contradicted | E43 |
| G69 | MITOPROBE | RESTORES_BRAIN_SIZE / positive | BRAIN | measured | E44 |
| G70 | MITOPROBE | REDUCES_MITOCHONDRIAL_STRESS_CONSEQUENCES / negative | STRESS | inferred | E44 |
| G71 | AMIK | PROMOTES_READTHROUGH / positive | READ | inferred | E18;E53 |
| G72 | AMIK | CAUSES_TOXICITY / positive | RENALTOX | measured | E53 |
| G73 | HCQ | INHIBITS_AUTOPHAGIC_DEGRADATION / negative | AUTO | measured | E54 |
| G74 | HCQ | INCREASES_ABUNDANCE_IN_SENESCENCE / positive | ABUND | measured | E54 |
| G75 | VINC | REDUCES_RMS_VIABILITY / negative | TUMOR | measured | E22 |
| G76 | OLAP | INHIBITS_TARGET / negative | PARP | measured | E55 |
| G77 | PARP | MAY_REDUCE_TUMOR_SURVIVAL_IF_REPAIR_DEPENDENT / negative | TUMOR | hypothesized | E24;E55 |
