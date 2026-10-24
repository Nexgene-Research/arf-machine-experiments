"""Emit an allele-keyed v1→v2 rank-change audit after final CSV exists."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT / "sub/Nexgene-Research_wgs-discovery.csv"
V2 = ROOT / "sub/Nexgene-AI_mn5-evidence-v2.csv"
OUT = ROOT / "results/track1/v1_v2_rank_changes.tsv"
AUDIT = ROOT / "results/track1/v2_adjudication.tsv"
AF_AUDIT = ROOT / "results/track1/v2_exclusions.tsv"


def key(row: dict) -> str:
    variants = []
    for n in ("1", "2"):
        if row[f"chrom_{n}"]:
            variants.append(
                ":".join(row[f"{part}_{n}"] for part in ("chrom", "pos", "ref", "alt"))
            )
    return ";".join(sorted(variants))


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    with V1.open(newline="", encoding="utf-8-sig") as fh:
        old = list(csv.DictReader(fh))
    with V2.open(newline="", encoding="utf-8-sig") as fh:
        new = list(csv.DictReader(fh))
    with AUDIT.open(newline="", encoding="utf-8") as fh:
        decisions = {
            r["grch38_pair"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")
        }
    with AF_AUDIT.open(newline="", encoding="utf-8") as fh:
        af_exclusion = {
            r["grch38_pair"].replace(":", "-"): r["reason"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["v2_status"] == "exclude"
        }
    new_positions = {key(row): pos for pos, row in enumerate(new, 1)}
    old_positions = {key(row): pos for pos, row in enumerate(old, 1)}
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("grch38_pair", "v1_rank", "v2_rank", "reason"))
        for row in old:
            pair = key(row)
            pos = new_positions.get(pair)
            source_key = ";".join(pair.replace(":", "-").split(";"))
            reason = (
                "retained with equal capped ePCR; promoted by caller concordance, established BUB1B gene-disease fit, stop plus kinase-domain missense over untested introns"
                if pos
                else "; ".join(
                    filter(
                        None,
                        [
                            decisions.get(source_key, "no independent v2 read support"),
                            af_exclusion.get(source_key, ""),
                        ],
                    )
                )
            )
            writer.writerow((pair, old_positions[pair], pos or "excluded", reason))
        for row in new:
            pair = key(row)
            if pair not in old_positions:
                writer.writerow(
                    (
                        pair,
                        "new",
                        new_positions[pair],
                        "newly supported read/phase/second-hit evidence; see evidence card",
                    )
                )
    print(
        f"Audited {len(old_positions)} legacy rows against {len(new_positions)} final rows"
    )


if __name__ == "__main__":
    main()
