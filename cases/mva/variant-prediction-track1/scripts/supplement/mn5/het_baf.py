#!/usr/bin/env python3
"""Stream heterozygous allele balance from an annotated VCF. No full load."""

from __future__ import annotations

import gzip
import statistics
import sys
from pathlib import Path


def main() -> None:
    vcf, dest = Path(sys.argv[1]), Path(sys.argv[2])
    opener = gzip.open if str(vcf).endswith(".gz") else open
    by_chrom: dict[str, list[float]] = {}
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        sample_index = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_index = 9
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if not (chrom.startswith("chr") and (chrom[3:].isdigit() or chrom in {"chrX", "chrY"})):
                continue
            fmt = fields[8].split(":")
            sample = fields[sample_index].split(":")
            rec = dict(zip(fmt, sample))
            gt = rec.get("GT", "./.").replace("|", "/")
            if gt not in {"0/1", "1/0"}:
                continue
            ad = rec.get("AD", "")
            parts = ad.split(",")
            if len(parts) < 2:
                continue
            try:
                ref_n, alt_n = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            depth = ref_n + alt_n
            if depth < 10:
                continue
            by_chrom.setdefault(chrom, []).append(alt_n / depth)
    lines = ["chrom\thets\tmedian_baf\tmean_baf\tnote\n"]
    for chrom, values in sorted(by_chrom.items()):
        lines.append(
            f"{chrom}\t{len(values)}\t{statistics.median(values):.4f}\t"
            f"{statistics.fmean(values):.4f}\tallele_balance_not_clone_fraction\n"
        )
    dest.write_text("".join(lines), encoding="utf-8")
    print(f"chromosomes={len(by_chrom)}")


if __name__ == "__main__":
    main()
