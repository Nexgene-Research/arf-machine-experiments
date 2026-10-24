"""Regenerate the independently supported subset of the legacy funnel.

Never edit legacy funnel.json. Fields requiring MN5 full hypothesis/read tables
are deliberately omitted instead of being copied from the legacy artifact.
"""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
TARGET = RESULTS / "track1/funnel_reconstructed.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-counts", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.primary_counts is not None:
        if args.out is None or args.out.exists():
            raise ValueError("new output path required and must not already exist")
        manifest = json.loads(
            (RESULTS / "candidate_universe_manifest.json").read_text()
        )
        primary = json.loads(args.primary_counts.read_text())
        h = primary["hypothesis_class"]
        r = primary["ranked"]
        s = primary["scored"]
        d = primary["read_stage"]
        original = json.loads((RESULTS / "track1/funnel.json").read_text())
        recovered = {
            **manifest["counts"],
            "dropped_gq_or_alt_reads": json.loads(
                (RESULTS / "track1/hypothesis_summary.json").read_text()
            )["dropped_gq_or_alt_reads"],
            "scored_candidates": s["scored_candidates"],
            "am_hit": s["am_hit"],
            "revel_hit": s["revel_hit"],
            "hypotheses": r["hypotheses"],
            "compound_het": h["compound_het"],
            "homozygous": h["homozygous"],
            "tier1_single": h["tier1_single"],
            "het_genes_with_pairs": primary["pair_genes"],
            "mad1l1_compound_het_pairs": primary["tier1_pair_counts"]["MAD1L1"],
            "bub1_compound_het_pairs": primary["tier1_pair_counts"]["BUB1"],
            "bub1b_compound_het_pairs": primary["tier1_pair_counts"]["BUB1B"],
            "cenatac_compound_het_pairs": primary["tier1_pair_counts"]["CENATAC"],
            "cep57_compound_het_pairs": primary["tier1_pair_counts"]["CEP57"],
            "top60_read_checked": primary["top60_read_checked"],
            "top60_both_snv": d["top60_both_snv"],
            "whatshap_spanning_reads_pair0": primary["whatshap_spanning_reads_pair0"],
            "whatshap_spanning_reads_pair1": primary["whatshap_spanning_reads_pair1"],
            "published_sites_pileup_alt": d["published_sites_pileup_alt"],
            "reads_ok_hypotheses": r["reads_ok_hypotheses"],
            "submitted_rows": primary["selected_rows"],
            "submitted_epcr_cap": primary["selected_epcr_cap"],
            "candidates_scored_sha256": primary["source_sha256"][
                "candidates_scored.tsv"
            ],
            "hypotheses_sha256": primary["source_sha256"]["hypotheses.tsv"],
            "ranked_pairs_sha256": primary["source_sha256"]["ranked_pairs_real.tsv"],
            "top10_readbacked_sha256": primary["output_sha256"],
        }
        missing = sorted(set(original) - set(recovered))
        mismatch = {
            key: {"reconstructed": value, "legacy": original.get(key)}
            for key, value in recovered.items()
            if value != original.get(key)
        }
        args.out.write_text(
            json.dumps(
                {
                    "reconstructed": recovered,
                    "legacy_mismatches": mismatch,
                    "unreconstructed_fields": missing,
                    "upstream_sources": [
                        "original MN5 full ranking, hypotheses, scored and read-stage tables",
                        "candidate_universe_manifest.json",
                        "hypothesis_summary.json",
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(
            f"primary_funnel_fields={len(recovered)} mismatches={len(mismatch)} unreconstructed={len(missing)}"
        )
        return
    if TARGET.exists():
        raise FileExistsError(TARGET)
    manifest = json.loads((RESULTS / "candidate_universe_manifest.json").read_text())
    hypotheses = json.loads((RESULTS / "track1/hypothesis_summary.json").read_text())
    epcr = json.loads((RESULTS / "track1/epcr_summary.json").read_text())
    with (RESULTS / "track1_top10_ranked.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {"annotation_job": "47397018", **manifest["counts"]}
    out.update(
        {
            "dropped_gq_or_alt_reads": hypotheses["dropped_gq_or_alt_reads"],
            "hypotheses": hypotheses["hypotheses"],
            "compound_het": hypotheses["compound_het"],
            "homozygous": hypotheses["homozygous"],
            "tier1_single": hypotheses["tier1_single"],
            "het_genes_with_pairs": hypotheses["het_genes_with_pairs"],
            "top60_read_checked": hypotheses["top60"],
            "submitted_rows": len(rows),
            "submitted_epcr_cap": max(float(row["epcr"]) for row in rows),
            "top10_readbacked_sha256": digest(RESULTS / "track1_top10_ranked.tsv"),
        }
    )
    source = RESULTS / "track1/funnel.json"
    original = json.loads(source.read_text())
    comparison = {
        key: [value, original.get(key)]
        for key, value in out.items()
        if value != original.get(key)
    }
    out["legacy_comparison_mismatches"] = comparison
    out["not_reconstructed_from_primary_data"] = [
        "reads_ok_hypotheses",
        "whatshap_spanning_reads_pair0",
        "whatshap_spanning_reads_pair1",
        "published_sites_pileup_alt",
        "top60_both_snv",
        "clinvar_annotation_sha256",
    ]
    TARGET.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        f"Reconstructed {len(out) - 2} fields; {len(comparison)} disagree with legacy; missing raw read counts explicit"
    )


if __name__ == "__main__":
    main()
