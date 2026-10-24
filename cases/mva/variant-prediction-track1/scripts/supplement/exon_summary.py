"""Per-gene min/median exon depth from the collected mosdepth table."""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
SRC = DEST / "computed" / "l3_sv" / "exon_depth.tsv"
OUT = DEST / "computed" / "l3_sv" / "exon_depth_by_gene.tsv"


def main() -> None:
    by_gene: dict[str, list[float]] = {}
    with SRC.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["interval"].split("|", 1)[0]
            by_gene.setdefault(gene, []).append(float(row["depth"]))
    all_vals = [v for vs in by_gene.values() for v in vs]
    auto = statistics.median(all_vals)
    lines = ["gene\tn_intervals\tmin_depth\tmedian_depth\tratio_to_sixgene_median\tnote\n"]
    for gene, vals in sorted(by_gene.items()):
        lines.append(
            f"{gene}\t{len(vals)}\t{min(vals):.2f}\t{statistics.median(vals):.2f}\t"
            f"{statistics.median(vals)/auto:.3f}\t"
            "mosdepth MAPQ20 exon mean; not a confirmed CNV\n"
        )
    OUT.write_text("".join(lines), encoding="utf-8")
    print(f"genes={len(by_gene)} median={auto:.2f}")


if __name__ == "__main__":
    main()
