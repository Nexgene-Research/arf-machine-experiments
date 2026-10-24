"""Append missing final v2 file hashes without altering historical provenance rows.

Runs once after red-team/report edits, before assembling the release snapshot.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROV = ROOT / "results/provenance.tsv"
ITEMS = [
    "results/track1/v2_mission.md",
    "results/track1/official_requirements_audit.md",
    "results/track1/v2_adjudication_protocol.md",
    "results/track1/bcd_47403447_receipt.json",
    "results/track1/bcd_47403447_audit.md",
    "results/track1/allele_validation.tsv",
    "results/track1/caller_concordance.tsv",
    "results/track1/phase_blocks_summary.tsv",
    "results/track1/phase_pair_evidence.tsv",
    "results/track1/second_hit_screen.tsv",
    "results/track1/tier1_evidence_matrix.tsv",
    "results/track1/rare_noncoding_keys.tsv",
    "results/track1/noncoding_avi.tsv",
    "results/track1/phase_aggregate.png",
    "results/track1/phase_aggregate.pdf",
    "results/track1/literature_queries.tsv",
    "results/track1/literature_evidence.tsv",
    "results/track1/evidence_convergence.tsv",
    "results/track1/variant_crosswalk.tsv",
    "results/track1/clinvar_exact_audit.tsv",
    "results/track1/v2_evidence_cards.md",
    "results/track1/v2_adjudication.tsv",
    "results/track1/v1_v2_rank_changes.tsv",
    "results/track1/v1_v2_rank_changes_redteam.tsv",
    "results/track1/v2_exclusions.tsv",
    "results/track1/redteam_reconciliation.tsv",
    "results/track1/top10_reconstructed.tsv",
    "results/track1/primary_counts_reconstructed.json",
    "results/track1/funnel_primary_reconstructed.json",
    "results/track1/exclusions_primary_reconstructed_v2.tsv",
    "results/track1/legacy_reconstruction_audit.md",
    "results/track1/legacy_reconstruction_verification_v2.json",
    "results/redteam/internal_review_v2.md",
    "results/redteam/responses_v2.md",
    "sub/Nexgene-AI_mn5-evidence-v2.csv",
    "sub/Nexgene-AI_methods_description.xlsx",
    "results/Nexgene-AI_track1_report.md",
    "docs/lineage.md",
]


def main() -> None:
    with PROV.open(newline="", encoding="utf-8") as fh:
        prior = list(csv.DictReader(fh, delimiter="\t"))
    mapped = {row["item"]: row["sha256"] for row in prior}
    output = []
    for rel in ITEMS:
        path = ROOT / rel
        if not path.is_file():
            raise FileNotFoundError(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if rel in mapped:
            if digest != mapped[rel]:
                raise ValueError(
                    f"previous provenance row mismatches final file: {rel}"
                )
            continue
        output.append(
            (
                rel,
                "v2-frozen-2026-10-16",
                "saved script / bounded MN5 receipt / public source log",
                digest,
                "2026-10-16",
                "mission final snapshot; restricted source genomes not included",
            )
        )
    with PROV.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerows(output)
    print(f"Recorded {len(output)} final v2 SHA-256 rows in results/provenance.tsv")


if __name__ == "__main__":
    main()
