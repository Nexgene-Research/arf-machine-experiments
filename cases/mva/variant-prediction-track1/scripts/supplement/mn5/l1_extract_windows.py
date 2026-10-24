#!/usr/bin/env python3
"""Keep only six-gene ±1 Mb rows from per-chromosome discordance tables."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

def load_windows(bed: Path) -> dict[str, list[tuple[int, int]]]:
    pad = 1_000_000
    extents: dict[str, list[int]] = {}
    keep = {"BUB1B", "BUB1", "CEP57", "TRIP13", "CENATAC", "MAD1L1"}
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene = name.split("|", 1)[0]
            if gene not in keep:
                continue
            extents.setdefault(chrom + "\t" + gene, []).extend([int(start), int(end)])
    windows: dict[str, list[tuple[int, int]]] = {}
    for key, values in extents.items():
        chrom, _gene = key.split("\t")
        lo, hi = max(0, min(values) - pad), max(values) + pad
        windows.setdefault(chrom, []).append((lo, hi))
    return windows


def keep_row(chrom: str, pos: int, windows: dict[str, list[tuple[int, int]]]) -> bool:
    for lo, hi in windows.get(chrom, ()):
        if lo <= pos <= hi:
            return True
    return False


def main() -> None:
    root = Path(sys.argv[1])
    bed = Path(sys.argv[2])
    dest = Path(sys.argv[3])
    windows = load_windows(bed)
    header = None
    kept = []
    for chrom_dir in sorted(root.glob("chr*")):
        table = chrom_dir / "caller_discordance.tsv"
        if not table.is_file():
            continue
        with table.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if header is None:
                header = reader.fieldnames
            for row in reader:
                pos = int(row["pos"])
                if keep_row(row["chrom"], pos, windows):
                    kept.append(row)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)
    print(f"discordance_rows={len(kept)}")


if __name__ == "__main__":
    main()
