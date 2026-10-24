#!/usr/bin/env python3
"""Pejaver-style LR for AlphaMissense and REVEL, BUB1B held out.

Streams ClinVar GRCh38 VCF. Joins AlphaMissense and REVEL by exact allele
without loading either table whole: AM is streamed once into a missense-key
set taken from ClinVar (bounded), REVEL zip is scanned for those keys only.
"""

from __future__ import annotations

import csv
import gzip
import math
import sys
import zipfile
from pathlib import Path

TIER1 = {"BUB1B"}
PATH_LABELS = {"pathogenic", "likely_pathogenic"}
BENIGN_LABELS = {"benign", "likely_benign"}


def stars(status: str) -> int:
    text = status.lower()
    if "practice_guideline" in text or "reviewed_by_expert_panel" in text:
        return 3
    if "criteria_provided,_multiple_submitters" in text:
        return 2
    return 0


def clinvar_missense(path: Path) -> list[dict]:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _q, _f, info = fields[:8]
            if len(ref) != 1 or len(alt) != 1:
                continue
            rec = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
            sig = rec.get("CLNSIG", "").lower()
            if stars(rec.get("CLNREVSTAT", "")) < 2:
                continue
            if "conflict" in sig:
                continue
            gene = rec.get("GENEINFO", "").split(":", 1)[0]
            if gene in TIER1:
                continue
            label = None
            if any(name in sig for name in PATH_LABELS) and "benign" not in sig:
                label = "P"
            elif any(name in sig for name in BENIGN_LABELS) and "pathogenic" not in sig:
                label = "B"
            if label is None:
                continue
            chrom_n = chrom if chrom.startswith("chr") else "chr" + chrom
            rows.append(
                {
                    "chrom": chrom_n,
                    "pos": pos,
                    "ref": ref.upper(),
                    "alt": alt.upper(),
                    "gene": gene,
                    "label": label,
                }
            )
    return rows


def chrom_name(value: str) -> str:
    value = value.strip()
    if value.startswith("chr"):
        return value
    if value in {"M", "MT"}:
        return "chrM"
    return "chr" + value


def join_am(path: Path, keys: set[tuple[str, str, str, str]]) -> dict:
    hits = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            key = (chrom_name(parts[0]), parts[1], parts[2].upper(), parts[3].upper())
            if key not in keys:
                continue
            try:
                score = float(parts[8])
            except ValueError:
                continue
            current = hits.get(key)
            if current is None or score > current:
                hits[key] = score
    return hits


def join_revel(path: Path, keys: set[tuple[str, str, str, str]]) -> dict:
    hits = {}
    with zipfile.ZipFile(path) as archive:
        name = "revel_with_transcript_ids"
        if name not in archive.namelist():
            name = next(n for n in archive.namelist() if "revel" in n.lower())
        with archive.open(name) as raw:
            header = True
            for line in raw:
                text = line.decode("utf-8", "replace").rstrip("\n")
                if header:
                    header = False
                    continue
                parts = text.split(",")
                if len(parts) < 8:
                    continue
                key = (chrom_name(parts[0]), parts[2], parts[3].upper(), parts[4].upper())
                if key not in keys:
                    continue
                try:
                    score = float(parts[7])
                except ValueError:
                    continue
                current = hits.get(key)
                if current is None or score > current:
                    hits[key] = score
    return hits


def lr_for_threshold(rows: list[dict], field: str, threshold: float, higher: bool) -> tuple[float, int, int, int, int]:
    tp = fp = fn = tn = 0
    for row in rows:
        score = row.get(field)
        if score is None:
            continue
        called = score >= threshold if higher else score <= threshold
        if row["label"] == "P":
            tp += called
            fn += not called
        else:
            fp += called
            tn += not called
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    # likelihood ratio of a pathogenic call; 0/0 stays undefined
    if spec >= 1 or sens <= 0:
        return float("nan"), tp, fp, fn, tn
    return sens / (1 - spec), tp, fp, fn, tn


def main() -> None:
    clinvar, am, revel, dest = map(Path, sys.argv[1:5])
    alleles = clinvar_missense(clinvar)
    keys = {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in alleles}
    am_hits = join_am(am, keys) if am.is_file() else {}
    revel_hits = join_revel(revel, keys) if revel.is_file() else {}
    for row in alleles:
        key = (row["chrom"], row["pos"], row["ref"], row["alt"])
        row["am"] = am_hits.get(key)
        row["revel"] = revel_hits.get(key)
    # Candidate thresholds around the observed allele (AM 0.9229, REVEL 0.472)
    am_lr, *_ = lr_for_threshold(alleles, "am", 0.9229, True)
    revel_lr, *_ = lr_for_threshold(alleles, "revel", 0.472, True)
    with dest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "predictor",
                "threshold",
                "held_out_genes",
                "n_clinvar_missense",
                "n_with_score",
                "lr_pathogenic_call",
                "candidate_allele",
                "candidate_score",
                "note",
            )
        )
        writer.writerow(
            (
                "AlphaMissense",
                "0.9229",
                "BUB1B",
                len(alleles),
                sum(r["am"] is not None for r in alleles),
                f"{am_lr:.4f}" if am_lr == am_lr else "undefined",
                "chr15-40220612-T-G",
                "0.9229",
                "two-star ClinVar missense; BUB1B excluded; not an ACMG PP3 assignment",
            )
        )
        writer.writerow(
            (
                "REVEL",
                "0.472",
                "BUB1B",
                len(alleles),
                sum(r["revel"] is not None for r in alleles),
                f"{revel_lr:.4f}" if revel_lr == revel_lr else "undefined",
                "chr15-40220612-T-G",
                "0.472",
                "two-star ClinVar missense; BUB1B excluded; not an ACMG PP3 assignment",
            )
        )
    print(f"clinvar={len(alleles)} am_hits={len(am_hits)} revel_hits={len(revel_hits)}")


if __name__ == "__main__":
    main()
