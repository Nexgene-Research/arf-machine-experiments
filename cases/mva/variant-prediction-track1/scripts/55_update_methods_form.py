"""Update the filled Track 1 methods workbook in place, keeping the official layout.

Usage: python scripts/55_update_methods_form.py
The blank organiser form is not part of this release, so 43_methods_form.py
cannot regenerate the workbook; this script rewrites only the answers below.
"""

from __future__ import annotations

from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "submission/Nexgene-AI_methods_description.xlsx"

ANSWERS = {
    7: "Nexgene AI",
    8: "1 — mn5-evidence-v2: compound heterozygous BUB1B p.Leu737* / p.Asn1002Lys (single ranked row)",
    9: (
        "GRCh38 blood WGS from one affected child. Independent bwa-mem2 realignment and duplicate marking on the "
        "MareNostrum 5 supercomputer, compared with the supplied Sentieon call set. Three germline callers "
        "(DeepVariant, GATK HaplotypeCaller, Strelka2) were run on that BAM; an allele counted as a vote when PASS "
        "and diploid ALT (raw GATK FILTER empty, not a vote). Annotation: snpEff/MANE v1.5, gnomAD, ClinVar, "
        "AlphaMissense, REVEL and AlphaGenome. A pre-registered four-angle design (genotype, phenotype-only, "
        "mechanism-only, published cases) over six Tier-1 MVA genes and the genome, building homozygous and compound "
        "heterozygous hypotheses (246,062) from 2,531 rare protein-altering/splice alleles. Validation: MAPQ60/BQ20 "
        "bidirectional read support, independent caller concordance, three-caller consensus, WhatsHap and GATK phasing, "
        "six-gene depth CNV, Manta/Delly SV, SpliceAI, genome-wide arm depth and GATK ModelSegments allele balance, "
        "VerifyBamID2, ROH, genotyping of 13 published pathogenic sites, HPO-weight permutations and a chr20 "
        "HG002/HG003 benchmark. Variants classified with ACMG/AMP plus ClinGen SVI (PVS1 decision tree; Bergquist 2025 "
        "/ Pejaver 2022 calibrated PP3). Full methods, figures and references: submission/Nexgene-AI_track1_report.md/.pdf."
    ),
    10: (
        "Nexgene AI Arf Machine research platform, a local workstation and MareNostrum 5. Language model ASA-W1 "
        "(Nexgene AI medical-reasoning model) ran locally. No patient data were exported, used to train or fine-tune "
        "any model. External database queries contained only public gene symbols, single public variant notations or "
        "publication identifiers. Data-use conditions were followed; controlled source and derived patient files will "
        "be deleted, with all copies, within 30 days of hackathon close, and the organisers notified. All AI-assisted "
        "interpretation was reviewed by the authors."
    ),
    13: (
        "Restricted hackathon patient sequence and clinical data (authorised research access), plus public resources only: "
        "GRCh38 reference, gnomAD, ClinVar, MANE, UniProt, HPO, PrimeKG, GIAB truth sets and published literature. "
        "No proprietary third-party datasets were used. Third-party licences are listed in LICENSE."
    ),
    14: (
        "GRCh38 no-alt reference, gnomAD allele frequencies, ClinVar assertions/review status, MANE/snpEff consequences, "
        "HPO/gene-disease sources, independent published case series; query identifiers and dates are in "
        "results/evidence/literature/literature_queries.tsv and the report reference list."
    ),
    15: (
        "Only the gated single-patient FASTQ, supplied VCF and clinical description provided by the organisers. "
        "These files are not included in the release repository. Restricted and derived patient data will be deleted "
        "within 30 days of hackathon close."
    ),
    18: (
        "Alignment, three-caller genotyping, annotation, phasing, CNV/SV, SpliceAI and mosaic screens ran as MareNostrum 5 "
        "Slurm jobs (Appendix A3 of submission/Nexgene-AI_track1_report.md; job identifiers also in docs/lineage.md). "
        "Ranking, ACMG/AMP adjudication and the report were completed on a local workstation. No commercial generative-AI "
        "API was used (ASA-W1 ran locally), so that spend is zero. MN5 allocation was the EuroHPC award "
        "EHPC-AIF-2026PG01-1286; a single wall-clock total was not recorded."
    ),
    19: (
        "We identify the genetic cause of mosaic variegated aneuploidy in one child from blood whole-genome sequencing. "
        "A pre-registered strategy kept genotype, phenotype-only, mechanism-only and published-case evidence blinded until "
        "convergence. From 5,035,176 variant records we built 246,062 single-gene genotype hypotheses, verified the top "
        "hypotheses at read level on an independent MareNostrum 5 alignment, confirmed both alleles by three-caller "
        "consensus (DeepVariant and Strelka2 PASS 0/1; GATK 0/1, unfiltered), and adjudicated them on read authenticity, "
        "caller concordance, gene-disease inheritance, the molecular effect of the weaker allele and phenotype fit. The "
        "single submitted genotype is compound heterozygous BUB1B NM_001211.6:c.2210T>G p.(Leu737*) and c.3006T>G "
        "p.(Asn1002Lys). Both alleles have balanced bidirectional high-quality reads and independent recall; both are rare "
        "(gnomAD AF 7.4e-5, 5e-7). The stop-gain lies 747 nt upstream of the last exon junction (PVS1; ClinVar P/LP). The "
        "missense lies in the BUBR1 kinase domain ten residues from the published MVA substitution L1012P and reaches "
        "PP3_Moderate with calibrated AlphaMissense (0.923); this matches the truncating-plus-destabilised-kinase-domain "
        "genotype class of biallelic MVA1. BUB1B ranked first on phenotype alone (83/100 weight permutations) and is the "
        "only MVA gene with two rare coding alleles; CNV, SV, splice and published-allele screens found no competing "
        "explanation, and contamination and ROH were negligible. Strengths: pre-registration, blinded angles, independent "
        "realignment, three-caller consensus, calibrated predictors, explicit exclusion of alternatives and hash-level "
        "provenance. Limitations: no parental samples or long reads, so the 10,911 bp pair is unphased; no functional "
        "BUBR1 assay; bulk blood WGS cannot measure cell-level mosaic aneuploidy; benchmark accuracy measured on chr20 "
        "only; three-caller agreement is a technical replicate on one BAM. ePCR is a ranking heuristic, not a calibrated "
        "clinical probability."
    ),
}


def main() -> None:
    wb = load_workbook(BOOK)
    ws = wb["Track 1 methods"]
    for row, answer in ANSWERS.items():
        if not ws[f"A{row}"].value:
            raise ValueError(f"question cell A{row} is empty; layout changed")
        ws[f"B{row}"] = answer
    wb.save(BOOK)
    print(f"updated {len(ANSWERS)} Track 1 answers in {BOOK}")


if __name__ == "__main__":
    main()
