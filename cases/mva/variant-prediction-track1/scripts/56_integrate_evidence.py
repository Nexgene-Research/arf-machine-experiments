"""Copy bounded evidence tables and producing scripts from the sibling evidence packs.

Run from the release root's parent layout (`release/<pack>`). Genome-wide patient
tables, raw MN5 logs, dependency wheels and review transcripts are not copied.
Writes `docs/evidence_index.tsv` with source path, destination and SHA-256.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKS = ROOT.parent
SUPP = PACKS / "track1-evidence-supplement"
WAVE3 = PACKS / "track1-evidence-wave3"
RIGOR = PACKS / "track1-rigor-pass"
EVIDENCE = ROOT / "results/evidence"

TABLES: dict[str, list[tuple[Path, str]]] = {
    "alleles": [
        (SUPP / "compiled/alleles/allele_validation.tsv", "allele_validation.tsv"),
        (SUPP / "compiled/alleles/caller_concordance.tsv", "caller_concordance.tsv"),
        (SUPP / "compiled/alleles/variant_crosswalk.tsv", "variant_crosswalk.tsv"),
        (SUPP / "compiled/alleles/clinvar_exact_audit.tsv", "clinvar_exact_audit.tsv"),
        (SUPP / "computed/l5_alphagenome/coding_snv_scores.tsv", "coding_avi_splice.tsv"),
        (SUPP / "computed/calibration/predictor_lr.tsv", "predictor_lr.tsv"),
        (SUPP / "computed/l6_missense/asn1002lys.tsv", "asn1002lys_structure.tsv"),
        (RIGOR / "results/vaf/vaf_jury_alleles.tsv", "vaf_candidate_alleles.tsv"),
        (RIGOR / "results/clinvar/clinvar_plp_jury.tsv", "clinvar_candidate_alleles.tsv"),
        (RIGOR / "results/splice/splice_jury.tsv", "alphagenome_splice_candidate_alleles.tsv"),
    ],
    "phase": [
        (SUPP / "compiled/phase/phase_pair_evidence.tsv", "phase_pair_evidence.tsv"),
        (SUPP / "compiled/phase/phase_blocks_summary.tsv", "phase_blocks_summary.tsv"),
        (SUPP / "computed/phase_span/span.tsv", "phase_span.tsv"),
        (RIGOR / "results/pid/pid_jury.tsv", "gatk_pid_candidate_alleles.tsv"),
        (RIGOR / "results/pid/pid_pair_counts.tsv", "gatk_pid_pair_counts.tsv"),
        (RIGOR / "results/pid/pid_sixgene.tsv", "gatk_pid_sixgene.tsv"),
    ],
    "second_hit": [
        (SUPP / "compiled/second_hit/second_hit_screen.tsv", "second_hit_screen.tsv"),
        (SUPP / "compiled/second_hit/tier1_evidence_matrix.tsv", "tier1_evidence_matrix.tsv"),
        (SUPP / "compiled/second_hit/rare_noncoding_keys.tsv", "rare_noncoding_keys.tsv"),
        (SUPP / "compiled/second_hit/noncoding_avi.tsv", "noncoding_avi.tsv"),
        (WAVE3 / "results/spliceai/spliceai_keys.tsv", "spliceai_keys.tsv"),
        (RIGOR / "results/splice/splice_sixgene_counts.tsv", "alphagenome_splice_sixgene_counts.tsv"),
    ],
    "structural": [
        (SUPP / "computed/l3_sv/tier1_sv.tsv", "tier1_sv.tsv"),
        (SUPP / "computed/l3_sv/manta_tier1_sv.tsv", "manta_tier1_sv.tsv"),
        (SUPP / "computed/l3_sv/delly_tier1_sv.tsv", "delly_tier1_sv.tsv"),
        (SUPP / "computed/l3_sv/exon_depth.tsv", "exon_depth.tsv"),
        (SUPP / "computed/l3_sv/exon_depth_by_gene.tsv", "exon_depth_by_gene.tsv"),
        (SUPP / "computed/l3_sv/tier1_sixgene_exon_utr.bed", "tier1_sixgene_exon_utr.bed"),
        (WAVE3 / "results/cnv/sixgene_cnv.tsv", "sixgene_cnv.tsv"),
    ],
    "aneuploidy": [
        (SUPP / "computed/dosage_baf/arm_depth.tsv", "arm_depth.tsv"),
        (SUPP / "computed/dosage_baf/het_baf.tsv", "het_baf.tsv"),
        (WAVE3 / "results/mosaic/mosaic_segments.tsv", "mosaic_segments.tsv"),
        (RIGOR / "results/roh/roh_summary.tsv", "roh_summary.tsv"),
    ],
    "qc": [
        (SUPP / "compiled/qc/qc_full_interim.md", "qc_full_interim.md"),
        (SUPP / "compiled/qc/l4_full_mosdepth_collection.json", "l4_full_mosdepth_collection.json"),
        (SUPP / "compiled/qc/l4_full_somalier_summary.json", "l4_full_somalier_summary.json"),
        (SUPP / "compiled/qc/l4_full_qc_reconciliation.json", "l4_full_qc_reconciliation.json"),
        (RIGOR / "results/vaf/vaf_counts.tsv", "vaf_counts.tsv"),
        (RIGOR / "results/vaf/vaf_by_gene_six.tsv", "vaf_by_gene_six.tsv"),
        (RIGOR / "results/vaf/scan_counts.json", "scan_counts.json"),
    ],
    "genome_scan": [
        (RIGOR / "results/clinvar/clinvar_plp_summary.tsv", "clinvar_plp_summary.tsv"),
        (RIGOR / "results/clinvar/clinvar_plp_genes.tsv", "clinvar_plp_genes.tsv"),
        (RIGOR / "results/rare/rare_counts.tsv", "rare_counts.tsv"),
        (RIGOR / "results/rare/rare_sixgene.tsv", "rare_sixgene.tsv"),
        (RIGOR / "results/rare/rare_acmg_sf.tsv", "rare_acmg_sf.tsv"),
    ],
    "literature": [
        (SUPP / "compiled/literature/angle_b_phenotype.tsv", "angle_b_phenotype.tsv"),
        (SUPP / "compiled/literature/angle_c_mechanism.tsv", "angle_c_mechanism.tsv"),
        (SUPP / "compiled/literature/evidence_convergence.tsv", "evidence_convergence.tsv"),
        (SUPP / "compiled/literature/literature_cases.tsv", "literature_cases.tsv"),
        (SUPP / "compiled/literature/literature_evidence.tsv", "literature_evidence.tsv"),
        (SUPP / "compiled/literature/literature_queries.tsv", "literature_queries.tsv"),
        (SUPP / "compiled/literature/genotype_queue.tsv", "genotype_queue.tsv"),
        (SUPP / "computed/primekg/kg_rediscovery.tsv", "kg_rediscovery.tsv"),
        (WAVE3 / "refs/published_alleles.tsv", "published_alleles.tsv"),
        (WAVE3 / "results/pileup/published_site_gt_new.tsv", "published_site_genotypes.tsv"),
    ],
    "controls": [
        (SUPP / "computed/null/null_runs.tsv", "null_runs.tsv"),
        (SUPP / "computed/positive_controls/positive_controls.tsv", "positive_controls.tsv"),
        (SUPP / "compiled/benchmark/l7_metrics_audited.tsv", "l7_metrics_audited.tsv"),
    ],
    "ranking": [
        (SUPP / "compiled/ranking/funnel.json", "funnel.json"),
        (SUPP / "compiled/ranking/top10_reconstructed.tsv", "top10_reconstructed.tsv"),
        (SUPP / "compiled/ranking/v1_v2_rank_changes.tsv", "v1_v2_rank_changes.tsv"),
        (SUPP / "compiled/ranking/redteam_reconciliation.tsv", "epcr_criterion_audit.tsv"),
        (SUPP / "compiled/ranking/exclusions.tsv", "exclusions.tsv"),
        (SUPP / "compiled/ranking/anomalies.tsv", "anomalies.tsv"),
        (SUPP / "compiled/ranking/Nexgene-Research_wgs-discovery.csv", "Nexgene-Research_wgs-discovery_v1.csv"),
        (SUPP / "computed/graph/ranked_pairs.tsv", "ranked_pairs.tsv"),
        (SUPP / "computed/graph/evidence_graph.json", "evidence_graph.json"),
    ],
}

SCRIPTS: dict[str, Path] = {
    "supplement": SUPP / "scripts",
    "wave3": WAVE3 / "scripts",
    "rigor": RIGOR / "scripts",
}

# Session plumbing (job polling, transfer retries, cancellations); not analysis steps.
OPERATIONAL = re.compile(
    r"^(cancel_|check_mn5|fetch_|repair_|resubmit_|run_|sacct_|stop_|submit_)"
)

# VerifyBamID2 output names the raw sample; only the summary statistics are kept.
VBID = WAVE3 / "results/verifybamid/verify.selfSM"
SCRUB = (
    (re.compile(r"\binternal review\b"), "internal review"),
    (re.compile(r"[A-Za-z0-9_.-]+@alogin\d\.bsc\.es"), "<user>@<mn5-login>"),
    (re.compile(r"[A-Za-z0-9_.-]+@transfer\d\.bsc\.es"), "<user>@<mn5-transfer>"),
    (re.compile(r"\balogin\d\.bsc\.es"), "<mn5-login>"),
    (re.compile(r"\btransfer\d\.bsc\.es"), "<mn5-transfer>"),
)


def scrub(text: str) -> str:
    for rx, repl in SCRUB:
        text = rx.sub(repl, text)
    return text


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    index = ["area\tdestination\tsource\tsha256"]
    for area, pairs in TABLES.items():
        target = EVIDENCE / area
        target.mkdir(parents=True, exist_ok=True)
        for source, name in pairs:
            dest = target / name
            shutil.copyfile(source, dest)
            index.append(
                f"{area}\t{dest.relative_to(ROOT).as_posix()}\t{source.relative_to(PACKS).as_posix()}\t{sha256(dest)}"
            )

    header, values = VBID.read_text(encoding="utf-8").splitlines()[:2]
    row = dict(zip(header.split("\t"), values.split("\t")))
    vbid = EVIDENCE / "qc/verifybamid2_summary.tsv"
    vbid.write_text(
        "metric\tvalue\tnote\n"
        f"snps\t{row['#SNPS']}\t1000G phase 3 10k GRCh38 SVD panel\n"
        f"reads\t{row['#READS']}\treads overlapping panel sites\n"
        f"avg_depth\t{row['AVG_DP']}\tmean depth at panel sites\n"
        f"freemix\t{row['FREEMIX']}\testimated contamination fraction; sample identifier removed\n",
        encoding="utf-8",
    )
    index.append(
        f"qc\t{vbid.relative_to(ROOT).as_posix()}\t{VBID.relative_to(PACKS).as_posix()} (summary)\t{sha256(vbid)}"
    )

    for name, source in SCRIPTS.items():
        target = ROOT / "scripts" / name
        for path in source.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts or OPERATIONAL.match(path.name):
                continue
            dest = target / path.relative_to(source)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(scrub(path.read_text(encoding="utf-8")), encoding="utf-8", newline="")
            index.append(
                f"scripts\t{dest.relative_to(ROOT).as_posix()}\t{path.relative_to(PACKS).as_posix()}\t{sha256(dest)}"
            )

    (ROOT / "docs/evidence_index.tsv").write_text("\n".join(index) + "\n", encoding="utf-8")
    print(f"copied {len(index) - 1} files")


if __name__ == "__main__":
    main()
