#!/usr/bin/env python3
"""Parse GATK ModelSegments copy-ratio and allele-fraction tables.

A missing segment file is recorded, not invented as zero clones.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_gatk_table(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def finite_float(text: str) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if value != value:  # NaN
        return None
    return value


def main() -> None:
    prefix, out_tsv = Path(sys.argv[1]), Path(sys.argv[2])
    model = read_gatk_table(Path(str(prefix) + ".modelFinal.seg")) or read_gatk_table(
        Path(str(prefix) + ".modelBegin.seg")
    )
    cr = read_gatk_table(Path(str(prefix) + ".cr.seg"))
    af = read_gatk_table(Path(str(prefix) + ".af.seg"))
    rows = []
    source = model or cr or af
    if not source:
        rows.append(
            {
                "chrom": ".",
                "start": ".",
                "end": ".",
                "call": "not_detected_or_unparsed",
                "mean_log2": ".",
                "minor_af": ".",
                "n_af": ".",
                "note": "ModelSegments tables missing or header-only; not a zero-clone claim",
            }
        )
    for rec in source:
        chrom = rec.get("CONTIG") or rec.get("contig") or rec.get("Chromosome") or "."
        start = rec.get("START") or rec.get("start") or rec.get("Start") or "."
        end = rec.get("END") or rec.get("end") or rec.get("End") or "."
        log2_text = (
            rec.get("LOG2_COPY_RATIO_POSTERIOR_50")
            or rec.get("MEAN_LOG2_COPY_RATIO")
            or rec.get("LOG2_COPY_RATIO")
            or "."
        )
        maf_text = (
            rec.get("MINOR_ALLELE_FRACTION_POSTERIOR_50")
            or rec.get("MINOR_ALLELE_FRACTION")
            or rec.get("MEAN_MINOR_ALLELE_FRACTION")
            or rec.get("Segment_Mean")
            or "."
        )
        n_af = rec.get("NUM_POINTS_ALLELE_FRACTION") or rec.get("Num_Probes") or "."
        n_cr = rec.get("NUM_POINTS_COPY_RATIO") or "."
        ratio = None
        log2 = finite_float(log2_text)
        if log2 is not None:
            ratio = 2 ** log2
        maf_f = finite_float(maf_text)
        call = "not_called"
        if ratio is not None:
            call = "copy_ratio_outlier" if ratio < 0.7 or ratio > 1.3 else "copy_ratio_near_diploid"
        elif n_cr in {"0", ".", ""}:
            call = "copy_ratio_not_modelled"
        if maf_f is not None:
            if maf_f < 0.40:
                call = "allelic_imbalance_candidate"
            elif call in {"not_called", "copy_ratio_not_modelled", "copy_ratio_near_diploid"}:
                call = "baf_balanced"
        note = "GATK ModelSegments on BAM allelic counts; blood WGS; not a clone fraction"
        if n_cr in {"0", ".", ""}:
            note += "; copy-ratio points=0 (no CollectReadCounts)"
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "call": call,
                "mean_log2": log2_text if log2 is not None else ".",
                "minor_af": f"{maf_f:.6f}" if maf_f is not None else ".",
                "n_af": n_af,
                "note": note,
            }
        )
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["chrom", "start", "end", "call", "mean_log2", "minor_af", "n_af", "note"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    n_flag = sum(1 for r in rows if "candidate" in r["call"] or "outlier" in r["call"])
    Path(out_tsv).with_name("clone_fraction.txt").write_text(
        "not_estimated\n"
        if n_flag == 0
        else f"segments_flagged={n_flag}; clone_fraction_not_estimated\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
