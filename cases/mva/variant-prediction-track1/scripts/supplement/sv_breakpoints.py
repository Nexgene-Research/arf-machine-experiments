"""Keep SVs whose breakpoint falls in a six-gene window. Span-only hits go aside."""

from __future__ import annotations

import csv
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
BED = DEST / "computed" / "l3_sv" / "tier1_sixgene_exon_utr.bed"
PAD = 50_000


def windows() -> dict[str, list[tuple[int, int, str]]]:
    extents: dict[str, list[int]] = {}
    with BED.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene = name.split("|", 1)[0]
            extents.setdefault(chrom + "\t" + gene, []).extend([int(start), int(end)])
    out: dict[str, list[tuple[int, int, str]]] = {}
    for key, vals in extents.items():
        chrom, gene = key.split("\t")
        out.setdefault(chrom, []).append((min(vals) - PAD, max(vals) + PAD, gene))
    return out


def genes_at(index, chrom: str, pos: int) -> set[str]:
    hits = set()
    for lo, hi, gene in index.get(chrom, ()):
        if lo <= pos <= hi:
            hits.add(gene)
    return hits


def main() -> None:
    index = windows()
    rows = []
    for name in ("manta_tier1_sv.tsv", "delly_tier1_sv.tsv"):
        path = DEST / "computed" / "l3_sv" / name
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    kept = []
    for row in rows:
        try:
            pos, end = int(row["pos"]), int(row["end"])
        except ValueError:
            continue
        hits = genes_at(index, row["chrom"], pos) | genes_at(index, row["chrom"], end)
        if not hits:
            continue
        row = dict(row)
        row["breakpoint_genes"] = ",".join(sorted(hits))
        kept.append(row)
    dest = DEST / "computed" / "l3_sv" / "tier1_sv.tsv"
    fields = list(rows[0].keys()) + ["breakpoint_genes"] if rows else []
    # unique fields
    seen = []
    for key in fields:
        if key not in seen:
            seen.append(key)
    with dest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=seen, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)
    print(f"span_overlap={len(rows)} breakpoint_in_window={len(kept)}")


if __name__ == "__main__":
    main()
