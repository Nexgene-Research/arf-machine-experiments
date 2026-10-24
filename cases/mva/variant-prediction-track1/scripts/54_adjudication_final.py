"""Preserve pre-red-team ranks and add final frozen-eligibility disposition."""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "results/track1/v1_v2_rank_changes.tsv"
OUT = ROOT / "results/track1/v1_v2_rank_changes_redteam.tsv"
STOP_MISSENSE = "chr15:40209701:T:G;chr15:40220612:T:G"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    if "**BLOCKING" not in (ROOT / "results/redteam/internal_review_v2.md").read_text(
        encoding="utf-8"
    ):
        raise ValueError("no recorded blocking internal-review objection")
    with PRIOR.open(newline="", encoding="utf-8") as fh:
        before = list(csv.DictReader(fh, delimiter="\t"))
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "grch38_pair",
                "v1_rank",
                "draft_v2_rank",
                "final_v2_eligible_rank",
                "reason",
            )
        )
        for row in before:
            pair = row["grch38_pair"]
            reason = (
                "internal review B1: p.Asn1002Lys domain context alone is not a damaging second-allele mechanism; B2: 0.1 score unreconciled"
                if pair == STOP_MISSENSE
                else row["reason"]
            )
            writer.writerow((pair, row["v1_rank"], row["v2_rank"], "excluded", reason))
    print(f"Final frozen-eligible v2 rows=0; reviewed historical rows={len(before)}")


if __name__ == "__main__":
    main()
