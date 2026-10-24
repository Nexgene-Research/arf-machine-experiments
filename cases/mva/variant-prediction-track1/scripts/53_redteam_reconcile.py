"""Recompute red-team disputed phase and capped-heuristic arithmetic.

Retains original provisional CSV unchanged as an audit artifact. No output is
promoted to a valid causal submission when frozen eligibility fails.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/track1/redteam_reconciliation.tsv"
PAIR = ROOT / "results/track1/phase_pair_evidence.tsv"
BLOCKS = ROOT / "results/track1/phase_blocks_summary.tsv"
CALLS = ROOT / "results/track1_top10_ranked.tsv"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    v1 = next(
        r
        for r in rows(CALLS)
        if r["gene"] == "BUB1B"
        and r["pos_1"] == "40209701"
        and r["pos_2"] == "40220612"
    )
    phased = next(
        r
        for r in rows(PAIR)
        if r["allele_1"] == "chr15-40209701-T-G"
        and r["allele_2"] == "chr15-40220612-T-G"
    )
    blocks = [r for r in rows(BLOCKS) if r["gene"] == "BUB1B"]
    if phased["phase"] != "unknown" or phased["same_block"] != "false":
        raise ValueError("source pair state changed")
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "allele",
                "reported_blocks_overlapping_allele",
                "v1_p_path",
                "v1_criteria",
                "lr1_only_p_path",
                "unknown_phase_assumption",
                "lr1_only_pair_heuristic",
                "v1_pair_q",
                "v1_csv_epcr",
                "interpretation",
            )
        )
        for index in (1, 2):
            pos = int(v1[f"pos_{index}"])
            block_count = sum(
                int(r["start_1based"]) <= pos <= int(r["end_1based"]) for r in blocks
            )
            if block_count:
                raise ValueError("reported block overlap changed; audit red-team B3")
            # With no justified independent missense LR, min(P_stop, P_missense)
            # <= prior 0.1 regardless of stop PVS1; 0.1 ceiling never raises it.
            prior, phase_prior = 0.1, float(v1["p_trans"])
            conditional = min(prior, float(v1["p_path_1"])) * phase_prior
            writer.writerow(
                (
                    f"chr15:{pos}",
                    block_count,
                    v1[f"p_path_{index}"],
                    v1[f"criteria_{index}"],
                    f"{prior:.6g}" if index == 2 else "stop consequence review pending",
                    phase_prior,
                    f"{conditional:.6g}",
                    v1["q_variant"],
                    v1["epcr"],
                    "provisional CSV has unsupported PM2-moderate and PP4; eligibility fails regardless of score",
                )
            )
    print(
        f"No reported block overlaps either allele; default LR1-only limiting-allele heuristic={conditional:.6g}; v1={v1['epcr']}; eligibility unresolved"
    )


if __name__ == "__main__":
    main()
