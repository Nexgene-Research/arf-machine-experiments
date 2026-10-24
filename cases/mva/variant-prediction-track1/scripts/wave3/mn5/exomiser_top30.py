#!/usr/bin/env python3
"""Take Exomiser .genes.tsv / .variants.tsv and emit a top-30 table plus BUB1B rank."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def main() -> None:
    genes_tsv, out_tsv = Path(sys.argv[1]), Path(sys.argv[2])
    _header, rows = read_tsv(genes_tsv)
    rank_col = None
    gene_col = None
    if not rows:
        out_tsv.write_text(
            "rank\tgene\texomiser_score\texomiser_gene_combined_score\tbub1b_in_top30\tnote\n"
            ".\t.\t.\t.\tno\tExomiser gene table missing; not a negative result\n",
            encoding="utf-8",
        )
        return
    for cand in ("EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE", "combinedScore"):
        if cand in rows[0]:
            rank_col = cand
            break
    for cand in ("GENE_SYMBOL", "gene_symbol", "#GENE_SYMBOL", "gene"):
        if cand in rows[0]:
            gene_col = cand
            break
    if gene_col is None:
        gene_col = next(iter(rows[0]))
    scored = []
    for rec in rows:
        raw = rec.get(rank_col or "", "")
        try:
            score = float(raw)
        except ValueError:
            score = -1.0
        scored.append((score, rec.get(gene_col, "."), rec))
    scored.sort(key=lambda x: x[0], reverse=True)
    top = scored[:30]
    genes30 = {g for _, g, _ in top}
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "rank",
                "gene",
                "exomiser_gene_combined_score",
                "bub1b_in_top30",
                "note",
            ]
        )
        bub = "yes" if "BUB1B" in genes30 else "no"
        for i, (score, gene, _rec) in enumerate(top, start=1):
            w.writerow(
                [
                    i,
                    gene,
                    f"{score:.6g}",
                    bub if i == 1 else "",
                    "phenotype ranker; not missense proof; not a causal claim",
                ]
            )
    bub_rank = next((i for i, (_, g, _) in enumerate(scored, start=1) if g == "BUB1B"), None)
    Path(out_tsv).with_name("bub1b_rank.txt").write_text(
        f"bub1b_rank={bub_rank if bub_rank is not None else 'not_in_table'}\n"
        f"bub1b_in_top30={'yes' if 'BUB1B' in genes30 else 'no'}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
