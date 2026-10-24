#!/usr/bin/env python3
"""Single-sample gene-body depth ratios. Not Canvas panel gCNV, not a karyotype."""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path


def load_regions(path: Path) -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".bed.gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.split()
            fields = rest[-1] if rest else "."
            # mosdepth regions: chrom start end ... mean
            mean = float(fields) if rest else float("nan")
            if path.name.endswith(".regions.bed.gz") or "regions" in path.name:
                # chrom start end name mean
                name = rest[0] if rest else f"{chrom}:{start}-{end}"
                mean = float(rest[1]) if len(rest) > 1 else float(rest[0])
            else:
                name = rest[0] if rest else f"{chrom}:{start}-{end}"
            rows.append((name, chrom, mean))
    return rows


def main() -> None:
    regions_bedgz, out_tsv = Path(sys.argv[1]), Path(sys.argv[2])
    rows = []
    with gzip.open(regions_bedgz, "rt", encoding="utf-8") as fh:
        for line in fh:
            chrom, start, end, name, mean_s = line.split()[:5]
            rows.append(
                {
                    "gene": name,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "mean_depth": float(mean_s),
                }
            )
    depths = [r["mean_depth"] for r in rows]
    median = sorted(depths)[len(depths) // 2] if depths else float("nan")
    out_rows = []
    for r in rows:
        ratio = r["mean_depth"] / median if median else float("nan")
        if ratio < 0.65:
            call = "loss_candidate"
        elif ratio > 1.35:
            call = "gain_candidate"
        else:
            call = "diploid_compatible"
        out_rows.append(
            {
                **{k: r[k] for k in ("gene", "chrom", "start", "end")},
                "mean_depth": f"{r['mean_depth']:.4f}",
                "ratio_to_sixgene_median": f"{ratio:.4f}",
                "call": call,
                "note": "single-sample mosdepth MAPQ20 gene-body+/-2kb; not Canvas/gCNV panel; not a causal claim",
            }
        )
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "gene",
                "chrom",
                "start",
                "end",
                "mean_depth",
                "ratio_to_sixgene_median",
                "call",
                "note",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(out_rows)


if __name__ == "__main__":
    main()
