#!/usr/bin/env python3
"""Summarise bcftools roh RG lines to F_ROH. Does not imply consanguinity or causality."""

from __future__ import annotations

import sys
from pathlib import Path

AUTOSOMES = {f"chr{n}" for n in range(1, 23)}


def prefix(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    if chrom == "MT":
        return "chrM"
    return "chr" + chrom


def autosome_length(fai: Path) -> int:
    total = 0
    for line in fai.read_text(encoding="utf-8").splitlines():
        chrom, length = line.split("\t")[:2]
        if prefix(chrom) in AUTOSOMES:
            total += int(length)
    if total < 2_000_000_000:
        raise SystemExit(f"unexpected autosome length {total}")
    return total


def main() -> None:
    raw = Path(sys.argv[1])
    fai = Path(sys.argv[2])
    dest = Path(sys.argv[3])
    dest.parent.mkdir(parents=True, exist_ok=True)
    genome = autosome_length(fai)
    n_win = 0
    roh_bp = 0
    n_markers = 0
    n_long = 0
    roh_bp_long = 0
    long_min = 1_500_000
    with raw.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("RG"):
                continue
            fields = line.rstrip("\n").split("\t")
            # RG sample chrom start end length nmarkers quality
            chrom = prefix(fields[2])
            if chrom not in AUTOSOMES:
                continue
            length = int(fields[5])
            n_win += 1
            roh_bp += length
            n_markers += int(fields[6])
            if length >= long_min:
                n_long += 1
                roh_bp_long += length
    f_roh = roh_bp / genome if genome else 0.0
    f_roh_long = roh_bp_long / genome if genome else 0.0
    dest.write_text(
        "metric\tvalue\tnote\n"
        f"autosome_bp\t{genome}\tGRCh38 no-alt primary chr1-22\n"
        f"roh_windows\t{n_win}\tbcftools roh RG lines on autosomes\n"
        f"roh_bp\t{roh_bp}\tsum of RG lengths\n"
        f"roh_markers\t{n_markers}\tmarkers inside those windows\n"
        f"F_ROH\t{f_roh:.8g}\tautosomal roh_bp / autosome_bp; not a consanguinity diagnosis\n"
        f"roh_windows_ge_1.5Mb\t{n_long}\tlong windows only\n"
        f"roh_bp_ge_1.5Mb\t{roh_bp_long}\tsum of RG lengths >=1.5 Mb\n"
        f"F_ROH_ge_1.5Mb\t{f_roh_long:.8g}\tlong-window autosomal fraction; still not a pedigree call\n"
        "caveat\tnot_causal\tthis table does not establish causality or a pedigree\n",
        encoding="utf-8",
    )
    print(f"windows={n_win} F_ROH={f_roh:.8g}")


if __name__ == "__main__":
    main()
