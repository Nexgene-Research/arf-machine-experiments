#!/usr/bin/env python3
"""Join six-gene coding/splice-neighbour SNVs with a local AlphaGenome slice.

The 104 GB genome-wide table is not copied. Missing scores are not benign.
This table does not establish causality.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

TIER1 = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1")
JURY = {
    ("chr15", "40209701", "T", "G"): "BUB1B_p.Leu737Ter",
    ("chr15", "40220612", "T", "G"): "BUB1B_p.Asn1002Lys",
}


def prefix(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else "chr" + chrom


def load_windows(bed: Path) -> dict[str, list[tuple[int, int, str]]]:
    out: dict[str, list[tuple[int, int, str]]] = {}
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene = line.rstrip("\n").split("\t")[:4]
            out.setdefault(prefix(chrom), []).append((int(start), int(end), gene))
    return out


def gene_at(windows: dict[str, list[tuple[int, int, str]]], chrom: str, pos: int) -> str | None:
    for lo, hi, gene in windows.get(chrom, ()):
        if lo < pos <= hi:
            return gene
    return None


def load_slice(path: Path) -> dict[tuple[str, str, str, str], tuple[str, str]]:
    out: dict[tuple[str, str, str, str], tuple[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, pos, ref, alt, score, gene = line.rstrip("\n").split("\t")[:6]
            out[(prefix(chrom), pos, ref, alt)] = (score, gene)
    return out


def load_background(path: Path) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["gene"]] = {k: float(row[k]) for k in ("p50", "p95", "p99") if row.get(k)}
    return out


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    vcf = Path(sys.argv[1])
    bed = Path(sys.argv[2])
    slice_path = Path(sys.argv[3])
    background_path = Path(sys.argv[4])
    dest = Path(sys.argv[5])
    dest.mkdir(parents=True, exist_ok=True)
    windows = load_windows(bed)
    table = load_slice(slice_path)
    background = load_background(background_path)
    per_gene: dict[str, Counter[str]] = {g: Counter() for g in TIER1}
    jury_rows: list[list[object]] = []
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = prefix(fields[0])
            pos, _vid, ref, alt = fields[1:5]
            if set(ref + alt) - set("ACGT") or alt == "*":
                continue
            if not (len(ref) == len(alt) == 1):
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            gene = gene_at(windows, chrom, pos_i)
            if gene is None:
                continue
            bucket = per_gene[gene]
            bucket["snv"] += 1
            hit = table.get((chrom, pos, ref, alt))
            if hit is None:
                bucket["missing"] += 1
                status = "not_in_slice"
                score_s = ""
                score_f: float | None = None
            else:
                score_s, _g = hit
                score_f = float(score_s)
                bucket["scored"] += 1
                status = "exact_allele"
                bg = background.get(gene, {})
                if "p95" in bg and score_f >= bg["p95"]:
                    bucket["ge_p95"] += 1
                if "p99" in bg and score_f >= bg["p99"]:
                    bucket["ge_p99"] += 1
            key = (chrom, pos, ref, alt)
            if key in JURY:
                bg = background.get(gene, {})
                jury_rows.append(
                    [
                        JURY[key],
                        chrom,
                        pos,
                        ref,
                        alt,
                        gene,
                        score_s,
                        status,
                        bg.get("p95", ""),
                        bg.get("p99", ""),
                        (
                            "ge_p99"
                            if score_f is not None and "p99" in bg and score_f >= bg["p99"]
                            else (
                                "ge_p95"
                                if score_f is not None and "p95" in bg and score_f >= bg["p95"]
                                else ("below_p95" if score_f is not None else "unscored")
                            )
                        ),
                        "uncalibrated AlphaGenome splicing; missing is not benign; not a causal claim",
                    ]
                )

    write_tsv(
        dest / "splice_sixgene_counts.tsv",
        ["gene", "proband_snv_in_window", "scored", "missing", "ge_local_p95", "ge_local_p99", "local_p95", "local_p99", "note"],
        [
            [
                gene,
                per_gene[gene]["snv"],
                per_gene[gene]["scored"],
                per_gene[gene]["missing"],
                per_gene[gene]["ge_p95"],
                per_gene[gene]["ge_p99"],
                background.get(gene, {}).get("p95", ""),
                background.get(gene, {}).get("p99", ""),
                "background is all panel SNVs in coding+splice+/-20bp; SpliceAI not used",
            ]
            for gene in TIER1
        ],
    )
    write_tsv(
        dest / "splice_jury.tsv",
        [
            "label",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gene",
            "alphagenome_splicing",
            "lookup_status",
            "local_p95",
            "local_p99",
            "vs_background",
            "note",
        ],
        jury_rows,
    )
    (dest / "caveat.txt").write_text(
        "AlphaGenome splicing is uncalibrated here. A missing row is not evidence of a normal splice site.\n",
        encoding="utf-8",
    )
    print("splice_genes", {g: dict(per_gene[g]) for g in TIER1})


if __name__ == "__main__":
    main()
