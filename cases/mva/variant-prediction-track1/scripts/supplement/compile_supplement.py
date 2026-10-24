"""Copy bounded, already-computed Track 1 evidence into the supplement.

Does not write into release/Nexgene-AI-track1. Refuses BAM/CRAM/VCF/FASTQ
/docx and files over 50 MB. Each copy gets a SOURCE.txt with hash and claim.
"""

from __future__ import annotations

import hashlib
import pathlib
import shutil

DEST = pathlib.Path(__file__).resolve().parents[1]
PROJECT = DEST.parent.parent
FROZEN = PROJECT / "release" / "Nexgene-AI-track1"
LIMIT = 50_000_000
BLOCKED = {".bam", ".cram", ".vcf", ".fastq", ".fq", ".docx"}
PROHIBITED = {
    "scripts/mn5/track1_v2_bcd_collect.py",
    "scripts/mn5/track1_v2_bcd_wait_collect.py",
}

# dest relative to DEST/compiled, source relative to PROJECT, claim, integrates_as
ITEMS = [
    (
        "alleles/allele_validation.tsv",
        "results/track1/allele_validation.tsv",
        "GRCh38 REF-checked alleles with independent bwa-mem2/bcftools recall",
        "results/track1/allele_validation.tsv",
    ),
    (
        "alleles/caller_concordance.tsv",
        "results/track1/caller_concordance.tsv",
        "Supplied VCF versus separate-alignment caller concordance for bounded alleles",
        "results/track1/caller_concordance.tsv",
    ),
    (
        "alleles/variant_crosswalk.tsv",
        "results/track1/variant_crosswalk.tsv",
        "MANE Select HGVS and protein crosswalk for the v2 pair",
        "results/track1/variant_crosswalk.tsv",
    ),
    (
        "alleles/clinvar_exact_audit.tsv",
        "results/track1/clinvar_exact_audit.tsv",
        "Exact-allele ClinVar audit; missing missense is not a negative",
        "results/track1/clinvar_exact_audit.tsv",
    ),
    (
        "alleles/alphagenome_queue.tsv",
        "results/alphagenome_queue.tsv",
        "Five BUB1B indel records queued for local AlphaGenome inference",
        "results/alphagenome_queue.tsv",
    ),
    (
        "phase/phase_pair_evidence.tsv",
        "results/track1/phase_pair_evidence.tsv",
        "WhatsHap pair outcomes; ten of ten unknown",
        "results/track1/phase_pair_evidence.tsv",
    ),
    (
        "phase/phase_blocks_summary.tsv",
        "results/track1/phase_blocks_summary.tsv",
        "Regional WhatsHap blocks over six Tier-1 genes",
        "results/track1/phase_blocks_summary.tsv",
    ),
    (
        "second_hit/second_hit_screen.tsv",
        "results/track1/second_hit_screen.tsv",
        "Six-gene exon/UTR depth and nonspecific split-read screen",
        "results/track1/second_hit_screen.tsv",
    ),
    (
        "second_hit/tier1_evidence_matrix.tsv",
        "results/track1/tier1_evidence_matrix.tsv",
        "Per-gene rare-allele, phase-pair and lowest exon/flank depth ratio",
        "results/track1/tier1_evidence_matrix.tsv",
    ),
    (
        "second_hit/rare_noncoding_keys.tsv",
        "results/track1/rare_noncoding_keys.tsv",
        "Twenty-eight rare noncoding SNV keys from the bounded screen",
        "results/track1/rare_noncoding_keys.tsv",
    ),
    (
        "second_hit/noncoding_avi.tsv",
        "results/track1/noncoding_avi.tsv",
        "Tabix AVI for those twenty-eight keys; uncalibrated LR=1",
        "results/track1/noncoding_avi.tsv",
    ),
    (
        "literature/literature_evidence.tsv",
        "results/track1/literature_evidence.tsv",
        "Primary-source reading extent and direction",
        "results/track1/literature_evidence.tsv",
    ),
    (
        "literature/literature_queries.tsv",
        "results/track1/literature_queries.tsv",
        "Public gene or single-allele queries only",
        "results/track1/literature_queries.tsv",
    ),
    (
        "literature/evidence_convergence.tsv",
        "results/track1/evidence_convergence.tsv",
        "Eight-dimension literature matrix; gene-level is not allele-level",
        "results/track1/evidence_convergence.tsv",
    ),
    (
        "literature/literature_cases.tsv",
        "results/literature_cases.tsv",
        "Partial published MVA/PCS catalogue with GRCh38 where mapped",
        "results/literature_cases.tsv",
    ),
    (
        "literature/genotype_queue.tsv",
        "results/genotype_queue.tsv",
        "Published positions queued for read-level genotyping",
        "results/genotype_queue.tsv",
    ),
    (
        "literature/angle_b_phenotype.tsv",
        "results/angle_b_phenotype.tsv",
        "Blind phenotype-only weighted coverage; not Exomiser",
        "results/angle_b_phenotype.tsv",
    ),
    (
        "literature/angle_c_mechanism.tsv",
        "results/angle_c_mechanism.tsv",
        "Blind mechanism chains from primary papers",
        "results/angle_c_mechanism.tsv",
    ),
    (
        "redteam/internal_review_v2.md",
        "results/redteam/internal_review_v2.md",
        "Single internal review v2 review; blocking eligibility objection",
        "results/redteam/internal_review_v2.md",
    ),
    (
        "redteam/responses_v2.md",
        "results/redteam/responses_v2.md",
        "Lead response sustaining the internal review eligibility gate",
        "results/redteam/responses_v2.md",
    ),
    (
        "redteam/summary.tsv",
        "results/redteam/summary.tsv",
        "v1-row red-team verdicts",
        "results/redteam/summary.tsv",
    ),
    (
        "redteam/responses.md",
        "results/redteam/responses.md",
        "Lead checks against the v1 summary",
        "results/redteam/responses.md",
    ),
    (
        "ranking/funnel.json",
        "results/track1/funnel.json",
        "Saved MN5 funnel counts; units are not a single disjoint flow",
        "results/track1/funnel.json",
    ),
    (
        "ranking/top10_reconstructed.tsv",
        "results/track1/top10_reconstructed.tsv",
        "Independent reconstruction of the ten v1 rows",
        "results/track1/top10_reconstructed.tsv",
    ),
    (
        "ranking/v1_v2_rank_changes.tsv",
        "results/track1/v1_v2_rank_changes.tsv",
        "v1 to provisional v2 rank movement",
        "results/track1/v1_v2_rank_changes.tsv",
    ),
    (
        "ranking/redteam_reconciliation.tsv",
        "results/track1/redteam_reconciliation.tsv",
        "Conditional ePCR 0.05 after stripping automatic PM2/PP4",
        "results/track1/redteam_reconciliation.tsv",
    ),
    (
        "ranking/v2_evidence_cards.md",
        "results/track1/v2_evidence_cards.md",
        "Narrative card for the exploratory BUB1B pair",
        "results/track1/v2_evidence_cards.md",
    ),
    (
        "ranking/candidates_full.tsv",
        "results/candidates_full.tsv",
        "Critical-annotation candidate table returned from MN5",
        "results/candidates_full.tsv",
    ),
    (
        "ranking/Nexgene-Research_wgs-discovery.csv",
        "sub/Nexgene-Research_wgs-discovery.csv",
        "Format-valid v1 ten-row CSV; operational fallback, not a cleared diagnosis",
        "sub/Nexgene-Research_wgs-discovery.csv",
    ),
    (
        "ranking/exclusions.tsv",
        "results/exclusions.tsv",
        "Considered and dropped genes or alleles with reasons",
        "results/exclusions.tsv",
    ),
    (
        "ranking/anomalies.tsv",
        "results/anomalies.tsv",
        "Unresolved observations",
        "results/anomalies.tsv",
    ),
    (
        "qc/qc_full_interim.md",
        "results/qc_full_interim.md",
        "mosdepth and somalier interim QC; not a mosaic karyotype",
        "results/qc_full_interim.md",
    ),
    (
        "qc/l4_full_mosdepth_collection.json",
        "results/l4_full_mosdepth_collection.json",
        "Receipt for 100 kb mosdepth aggregates transferred from MN5",
        "results/l4_full_mosdepth_collection.json",
    ),
    (
        "qc/l4_full_somalier_summary.json",
        "results/l4_full_somalier_summary.json",
        "Somalier fingerprint summary receipt",
        "results/l4_full_somalier_summary.json",
    ),
    (
        "qc/l4_full_qc_reconciliation.json",
        "results/l4_full_qc_reconciliation.json",
        "Checksum reconciliation of L4 QC aggregates",
        "results/l4_full_qc_reconciliation.json",
    ),
    (
        "benchmark/l7_metrics_audited.tsv",
        "results/l7_metrics_audited.tsv",
        "chr20-only HG002/HG003 exact-allele precision and recall",
        "results/l7_metrics_audited.tsv",
    ),
    (
        "provenance/provenance.tsv",
        "results/provenance.tsv",
        "Existing project provenance ledger",
        "results/provenance.tsv",
    ),
]


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_item(dest_rel: str, source_rel: str, claim: str, integrates: str) -> dict:
    if source_rel.replace("\\", "/") in PROHIBITED:
        raise ValueError(f"prohibited: {source_rel}")
    src = PROJECT / source_rel
    if not src.is_file():
        raise FileNotFoundError(src)
    if src.stat().st_size > LIMIT:
        raise ValueError(f"over 50 MB: {src}")
    if src.suffix.lower() in BLOCKED or src.name.lower().endswith(".vcf.gz"):
        raise ValueError(f"blocked suffix: {src}")
    dest = DEST / "compiled" / dest_rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if digest(dest) != digest(src):
            raise FileExistsError(f"existing copy differs: {dest}")
    else:
        shutil.copy2(src, dest)
    sha = digest(dest)
    source_note = dest.with_name(dest.name + ".SOURCE.txt")
    source_note.write_text(
        f"source\t{src.resolve()}\n"
        f"sha256\t{sha}\n"
        f"bytes\t{dest.stat().st_size}\n"
        f"claim\t{claim}\n"
        f"integrates_as\t{integrates}\n",
        encoding="utf-8",
    )
    return {
        "compiled": f"compiled/{dest_rel}",
        "source": source_rel.replace("\\", "/"),
        "sha256": sha,
        "bytes": dest.stat().st_size,
        "claim": claim,
        "integrates_as": integrates,
    }


def write_integration(rows: list[dict]) -> None:
    planned = [
        (
            "computed/phase_span/span.tsv",
            "results/track1/phase_span.tsv",
            "Physical distance between the two BUB1B sites; short-read phase impossible",
        ),
        (
            "computed/l5_alphagenome/coding_snv_scores.tsv",
            "results/track1/coding_avi_splice.tsv",
            "AVI and splicing table hits for the stop and missense SNVs",
        ),
        (
            "computed/graph/evidence_graph.json",
            "results/evidence_graph.json",
            "Auditable nodes and edges from ranked pairs and blinded B/C/D tables",
        ),
        (
            "computed/dosage_baf/arm_depth.tsv",
            "results/track1/arm_depth.tsv",
            "Chromosome-arm depth from 100 kb mosdepth; not a clone fraction",
        ),
        (
            "computed/l1_consensus/caller_discordance.tsv",
            "results/caller_discordance.tsv",
            "Tier-1 locus discordance from finished DeepVariant/GATK/Strelka VCFs",
        ),
        (
            "computed/l3_sv/tier1_sv.tsv",
            "results/track1/tier1_sv.tsv",
            "Manta/Delly records overlapping the six genes",
        ),
        (
            "computed/l3_sv/exon_depth.tsv",
            "results/track1/exon_depth.tsv",
            "mosdepth exon means for the six-gene BED",
        ),
        (
            "computed/dosage_baf/het_baf.tsv",
            "results/track1/het_baf.tsv",
            "Per-chromosome heterozygous allele balance from the annotated VCF",
        ),
        (
            "computed/l6_missense/asn1002lys.tsv",
            "results/track1/asn1002lys_structure.tsv",
            "ThermoMPNN/ESM at p.Asn1002Lys if the structure covers the residue",
        ),
        (
            "computed/l5_alphagenome/indel_inference.tsv",
            "results/alphagenome_scores.tsv",
            "Local AlphaGenome scores for queued BUB1B indels after a parity pilot",
        ),
        (
            "computed/calibration/predictor_lr.tsv",
            "results/track1/predictor_lr.tsv",
            "ClinVar two-star missense LR for AM/REVEL with BUB1B held out",
        ),
        (
            "computed/primekg/kg_rediscovery.tsv",
            "results/kg_rediscovery.tsv",
            "Leave-one-out rediscovery ranks of five known MVA gene-disease links",
        ),
        (
            "computed/null/null_runs.tsv",
            "results/null_runs.tsv",
            "HPO-permutation and HG002-truth false-candidate counts",
        ),
        (
            "computed/positive_controls/positive_controls.tsv",
            "results/positive_controls.tsv",
            "BAMSurgeon recovery ranks of GRCh38-mapped published alleles",
        ),
        (
            "computed/not_available.md",
            "docs/not_available.md",
            "Claims that cannot be computed from present data, without invented numbers",
        ),
    ]
    lines = [
        "# Integration map",
        "",
        "This supplement does not copy into `release/Nexgene-AI-track1`.",
        "A later merge copies a row only when its sha256 still matches this ledger.",
        "",
        "| Supplement path | Frozen-package destination | Role |",
        "|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['compiled']}` | `{row['integrates_as']}` | {row['claim']} |"
        )
    for path, dest, claim in planned:
        lines.append(f"| `{path}` | `{dest}` | {claim} |")
    (DEST / "INTEGRATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not FROZEN.is_dir():
        raise FileNotFoundError(FROZEN)
    for name in ("compiled", "computed", "scripts", "logs"):
        (DEST / name).mkdir(parents=True, exist_ok=True)
    for name in (
        "l1_consensus",
        "l3_sv",
        "dosage_baf",
        "l5_alphagenome",
        "l6_missense",
        "calibration",
        "phase_span",
        "graph",
        "primekg",
        "null",
        "positive_controls",
    ):
        (DEST / "computed" / name).mkdir(parents=True, exist_ok=True)
    rows = [copy_item(*item) for item in ITEMS]
    write_integration(rows)
    manifest = DEST / "MANIFEST.sha256"
    lines = [f"{row['sha256']}  {row['compiled']}" for row in rows]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"copied={len(rows)} dest={DEST}")


if __name__ == "__main__":
    main()
