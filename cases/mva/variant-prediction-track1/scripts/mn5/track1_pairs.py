"""Build Track 1 hypotheses from the scored candidate table.

Het rows that share a gene become compound-het pairs. Homozygous rows are
single-allele hypotheses. A Tier-1 gene with only one passing het allele is
kept as an unresolved second-allele hypothesis. Pairs with alt-AD 0 or GQ
below 20 are dropped.
"""

from __future__ import annotations

import csv
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

TIER1 = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1"}
LOF = (
    "frameshift",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "splice_acceptor",
    "splice_donor",
)


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def alt_reads(ad: str) -> float | None:
    parts = [number(item) for item in (ad or "").split(",")]
    if len(parts) < 2 or parts[1] is None:
        return None
    return parts[1]


def passing(row: dict[str, str]) -> bool:
    gq = number(row.get("gq", ""))
    reads = alt_reads(row.get("ad", ""))
    return gq is not None and gq >= 20 and reads is not None and reads > 0


def allele_score(row: dict[str, str]) -> int:
    consequence = (row.get("consequence") or "").lower()
    clinvar = (row.get("clinvar") or "").lower()
    pathogenic = "pathogenic" in clinvar and "benign" not in clinvar and "conflict" not in clinvar
    lof = row.get("impact") == "HIGH" or any(token in consequence for token in LOF)
    am = number(row.get("am_pathogenicity", ""))
    revel = number(row.get("revel", ""))
    hot = (am is not None and am >= 0.564) or (revel is not None and revel >= 0.5)
    if lof or pathogenic:
        return 100
    if hot and "missense" in consequence:
        return 40
    return 10


def main() -> None:
    src = Path(sys.argv[1]) / "candidates_scored.tsv"
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    hets: dict[str, list[dict[str, str]]] = defaultdict(list)
    homs: list[dict[str, str]] = []
    dropped = {"gq_or_reads": 0}
    for row in rows:
        gt = (row.get("gt") or "").replace("|", "/")
        if not passing(row):
            dropped["gq_or_reads"] += 1
            continue
        if gt in ("0/1", "1/0"):
            hets[row.get("gene") or "."].append(row)
        elif gt == "1/1":
            homs.append(row)

    hypotheses: list[dict[str, str]] = []

    def pack(kind: str, left: dict[str, str], right: dict[str, str] | None) -> dict[str, str]:
        gene = left.get("gene") or "."
        score = allele_score(left) + (allele_score(right) if right else 0)
        if gene in TIER1:
            score += 1000
        item = {
            "hypothesis": kind,
            "gene": gene,
            "tier1": "yes" if gene in TIER1 else "no",
            "prescore": str(score),
            "phase": "unknown",
        }
        for index, allele in ((1, left), (2, right)):
            for field in (
                "chrom",
                "pos",
                "ref",
                "alt",
                "gt",
                "gq",
                "dp",
                "ad",
                "vaf",
                "consequence",
                "impact",
                "clinvar",
                "gnomad_af",
                "af_state",
                "am_pathogenicity",
                "am_class",
                "revel",
            ):
                item[f"{field}_{index}"] = "" if allele is None else allele.get(field, "")
        return item

    pair_genes = 0
    for gene, alleles in hets.items():
        if len(alleles) < 2:
            continue
        pair_genes += 1
        for left, right in itertools.combinations(alleles, 2):
            if (left["chrom"], left["pos"], left["ref"], left["alt"]) == (
                right["chrom"],
                right["pos"],
                right["ref"],
                right["alt"],
            ):
                continue
            hypotheses.append(pack("compound_het", left, right))
    for row in homs:
        hypotheses.append(pack("homozygous", row, None))
    for gene, alleles in hets.items():
        if gene in TIER1 and len(alleles) == 1 and not any(row.get("gene") == gene for row in homs):
            hypotheses.append(pack("tier1_single", alleles[0], None))

    hypotheses.sort(key=lambda item: (-int(item["prescore"]), item["gene"], item["chrom_1"], item["pos_1"]))
    fields = list(hypotheses[0]) if hypotheses else ["hypothesis"]
    with (out / "hypotheses.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(hypotheses)
    top = hypotheses[:60]
    with (out / "top60.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(top)
    summary = {
        "candidates": len(rows),
        "dropped_gq_or_alt_reads": dropped["gq_or_reads"],
        "het_genes_with_pairs": pair_genes,
        "hypotheses": len(hypotheses),
        "compound_het": sum(item["hypothesis"] == "compound_het" for item in hypotheses),
        "homozygous": sum(item["hypothesis"] == "homozygous" for item in hypotheses),
        "tier1_single": sum(item["hypothesis"] == "tier1_single" for item in hypotheses),
        "top60": len(top),
    }
    (out / "hypothesis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
