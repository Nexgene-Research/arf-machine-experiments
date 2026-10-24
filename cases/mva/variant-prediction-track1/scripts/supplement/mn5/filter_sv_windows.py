#!/usr/bin/env python3
"""Emit SV records whose span overlaps a six-gene exon BED."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path


def windows(bed: Path) -> dict[str, list[tuple[int, int, str]]]:
    out: dict[str, list[tuple[int, int, str]]] = {}
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene = name.split("|", 1)[0]
            out.setdefault(chrom, []).append((int(start), int(end), gene))
    return out


def overlaps(chrom: str, start: int, end: int, index: dict) -> str:
    genes = []
    for lo, hi, gene in index.get(chrom, ()):
        if start < hi and end > lo:
            genes.append(gene)
    return ",".join(sorted(set(genes)))


def main() -> None:
    vcf, bed, dest, caller = sys.argv[1:5]
    index = windows(Path(bed))
    opener = gzip.open if vcf.endswith(".gz") else open
    rows = ["caller\tchrom\tpos\tref\talt\tend\tsvtype\tfilter\tgenes\n"]
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _qual, flt, info = fields[:8]
            end = pos
            svtype = "."
            for item in info.split(";"):
                if item.startswith("END="):
                    end = item.split("=", 1)[1]
                if item.startswith("SVTYPE="):
                    svtype = item.split("=", 1)[1]
            try:
                start_i, end_i = int(pos), int(end)
            except ValueError:
                continue
            if end_i < start_i:
                start_i, end_i = end_i, start_i
            genes = overlaps(chrom, start_i - 1, end_i, index)
            if not genes:
                continue
            rows.append(
                f"{caller}\t{chrom}\t{pos}\t{ref}\t{alt}\t{end}\t{svtype}\t{flt}\t{genes}\n"
            )
    Path(dest).write_text("".join(rows), encoding="utf-8")
    print(f"sv_rows={len(rows)-1}")


if __name__ == "__main__":
    main()
