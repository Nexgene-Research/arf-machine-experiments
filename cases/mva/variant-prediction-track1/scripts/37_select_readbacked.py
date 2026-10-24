"""Independently select legacy read-backed top ten and primary-stage aggregates.

Run only on MN5 against the original full ranked and read-stage tables. Export
only ten selected rows and aggregate counts; no BAM/VCF/reads leave MN5.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def count_table(path: Path, callback) -> None:
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            callback(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked-full", type=Path, required=True)
    parser.add_argument("--read-stage", type=Path, required=True)
    parser.add_argument("--hypotheses", type=Path, required=True)
    parser.add_argument("--scored", type=Path, required=True)
    parser.add_argument("--top60", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    args = parser.parse_args()
    for output in (args.out, args.counts):
        if output.exists():
            raise FileExistsError(output)
    for source in (
        args.ranked_full,
        args.read_stage,
        args.hypotheses,
        args.scored,
        args.top60,
    ):
        if source.stat().st_size < 1:
            raise ValueError(f"empty source: {source}")
    selected = []
    ranked_counts = Counter()
    per_gene_class = Counter()
    hom_per_gene = Counter()
    ranked_genes = set()
    with args.ranked_full.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        needed = {
            "reads_ok",
            "gene",
            "epcr",
            "pile_alt_1",
            "pile_alt_2",
            "chrom_1",
            "pos_1",
            "chrom_2",
            "pos_2",
        }
        if not needed.issubset(reader.fieldnames or ()):
            raise ValueError(
                f"source lacks columns {needed - set(reader.fieldnames or ())}"
            )
        columns = reader.fieldnames
        for row in reader:
            ranked_counts["hypotheses"] += 1
            ranked_counts["capped"] += row["capped"] == "yes"
            ranked_counts["clinvar_joined_rows"] += bool(
                row["clinvar_1"] != "." or row["clinvar_2"] != "."
            )
            ranked_counts["reads_ok_hypotheses"] += row["reads_ok"] == "yes"
            per_gene_class[(row["gene"], row["hypothesis"])] += 1
            if row["hypothesis"] == "homozygous":
                hom_per_gene[row["gene"]] += 1
            if row["hypothesis"] == "compound_het":
                ranked_genes.add(row["gene"])
            if row["reads_ok"] == "yes" and row["hypothesis"] == "compound_het":
                selected.append(row)
    selected.sort(
        key=lambda row: (
            -float(row["epcr"]),
            -(int(row["pile_alt_1"] or 0) + int(row["pile_alt_2"] or 0)),
            row["chrom_1"],
            row["pos_1"],
            row["chrom_2"],
            row["pos_2"],
        )
    )
    selected = selected[:10]
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, delimiter="\t", fieldnames=columns, lineterminator="\r\n"
        )
        writer.writeheader()
        writer.writerows(selected)
    read_counts = Counter()
    read_by_pair = []

    def read_row(row):
        read_counts[row["source"]] += 1
        if row["source"] == "published":
            read_counts["published_sites_pileup_alt"] += int(row["pile_alt_1"] or 0)
        elif row["source"] == "top60":
            read_counts["top60_both_snv"] += all(
                len(row[f"ref_{n}"]) == len(row[f"alt_{n}"]) == 1 for n in ("1", "2")
            )
            if len(read_by_pair) < 2:
                read_by_pair.append(int(row["phase_support"] or 0))

    count_table(args.read_stage, read_row)
    hypothesis_counts = Counter()
    hypothesis_genes = set()

    def hypothesis_row(row):
        hypothesis_counts[row["hypothesis"]] += 1
        if row["hypothesis"] == "compound_het":
            hypothesis_genes.add(row["gene"])

    count_table(args.hypotheses, hypothesis_row)
    scored_counts = Counter()

    def scored_row(row):
        scored_counts["scored_candidates"] += 1
        scored_counts["am_hit"] += bool(row.get("am_pathogenicity"))
        scored_counts["revel_hit"] += bool(row.get("revel"))

    count_table(args.scored, scored_row)
    top60_count = Counter()
    count_table(args.top60, lambda row: top60_count.update(["top60_read_checked"]))
    if (
        ranked_counts["hypotheses"] != sum(hypothesis_counts.values())
        or read_counts["top60"] != top60_count["top60_read_checked"]
        or len(read_by_pair) != 2
    ):
        raise ValueError("primary table count discrepancy")
    counts = {
        "source_sha256": {
            source.name: sha256(source)
            for source in (
                args.ranked_full,
                args.read_stage,
                args.hypotheses,
                args.scored,
                args.top60,
            )
        },
        "output_sha256": sha256(args.out),
        "ranked": dict(ranked_counts),
        "hypothesis_class": dict(hypothesis_counts),
        "selected_rows": len(selected),
        "selected_epcr_cap": max(float(row["epcr"]) for row in selected),
        "read_stage": dict(read_counts),
        "whatshap_spanning_reads_pair0": read_by_pair[0],
        "whatshap_spanning_reads_pair1": read_by_pair[1],
        "scored": dict(scored_counts),
        "top60_read_checked": top60_count["top60_read_checked"],
        "tier1_pair_counts": {
            gene: per_gene_class[(gene, "compound_het")]
            for gene in ("BUB1", "BUB1B", "CENATAC", "CEP57", "MAD1L1", "TRIP13")
        },
        "tier1_homozygous_counts": {
            gene: hom_per_gene[gene]
            for gene in ("BUB1", "BUB1B", "CENATAC", "CEP57", "MAD1L1")
        },
        "pair_genes": len(hypothesis_genes),
        "non_tier1_pair_genes": len(
            hypothesis_genes - {"BUB1", "BUB1B", "CENATAC", "CEP57", "MAD1L1", "TRIP13"}
        ),
    }
    args.counts.write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n")
    print(
        f"selected={len(selected)} from {ranked_counts['hypotheses']} ranked rows; "
        f"hash={counts['output_sha256']}"
    )


if __name__ == "__main__":
    main()
