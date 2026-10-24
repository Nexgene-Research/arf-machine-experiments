# Dated Track 1 v2 adjudication amendment

Frozen **2026-10-16 22:10 local clock time**, before submission of any new B+C+D analysis. The v1 ranking and its variants are already known; this is a dated amendment, not a claim of blinded preregistration. Publish this protocol beside v1 and the rank-change table. Source: `docs/preregistration.md` and the verbatim mission in `results/track1/v2_mission.md`.

## Eligibility: exclusion from final CSV if failed

1. GRCh38 normalized REF-matching allele supported by reads on the independent alignment; no obvious mapping or paralog artifact. If authentication cannot be completed, retain the unmodified v1 fallback, label v2 not completed.
2. For a recessive pair, each required allele has gnomAD AF ≤0.01 unless a documented hypomorphic founder exception is supported by a primary source.
3. Gene–disease validity alone does not establish a damaging allele; a pair requires a defensible damaging second-allele mechanism. Cis-proven pairs are excluded as compound heterozygotes. Unknown phase remains unknown, never trans.
4. A single allele may be proposed if an established inheritance model fits. Incidental findings alone may be labelled secondary.

## Score and ordered ties

Retain the exact preregistered `scripts/mn5/track1_epcr.py` and `scripts/32_epcr.py` case-level heuristic: per-allele `P_path = logistic(logit(0.1) + sum of unique justified log LRs)`; unknown-phase pair `q_variant=min(P_a,P_b) × 0.5` (sensitivity 0.25–0.75); case ePCR = weakest-link result times preregistered angle multiplier, capped at **0.1 whenever phase, second-allele mechanism, genotype authenticity or another critical link is unverified**. Computational predictors without empirical calibration (AlphaMissense, REVEL, AlphaGenome, phenotype and mechanism scores) have LR=1. ePCR is an unvalidated heuristic, never a clinical probability. Do not create fractional tie-break scores.

Equal ePCR values remain exactly equal; rank file order lexicographically by these *evidence* tiers, in order:
1. Read-backed **trans** phase (no evidence ≠ trans).
2. Genotype confirmed by separate alignment/caller path (same biological sample, not independent sample replication).
3. Established gene–disease validity and matching inheritance.
4. Reviewed molecular effect: LoF/splice > damaging missense > supported regulatory > unsupported intron/upstream. For recessive pairs assess the weaker allele as the limiting mechanism; do not promote a LoF + unsupported intron pair on the LoF alone.
5. Phenotype fit.
6. Number of independent supportive evidence *dimensions* in Phase L convergence table (one study/case counted once across multiple reports).
7. Population rarity among eligible alleles.
8. Lexicographic normalized GRCh38 allele-pair key.

An existing preregistered ACMG criterion changes `P_path` only when its original requirements apply exactly: e.g. two-star ClinVar P/LP meeting PS1 equivalence or a published same-allele case *confirmed in trans* meeting ClinGen SVI PM3 strength. Log primary source, criterion, affected allele and change; no retroactive score manipulation.

## Adjudication records

For each eligible candidate, produce an evidence card documenting normalized alleles; MAPQ/BQ, strand and caller checks; phase/block/direct versus chained evidence; AF/HomAlt; consequence, transcript, NMD, protein context, gene–disease validity and mode; phenotype; primary literature and independent-source counts; competing mechanism; strongest red-team objection; open critical links. The v1→v2 rank-change table explicitly gives each removal/reordering reason. A screen signal (depth/splits) never becomes a confirmed SV/CNV. If analysis fails, retain v1 intact and report the gap rather than invent evidence.
