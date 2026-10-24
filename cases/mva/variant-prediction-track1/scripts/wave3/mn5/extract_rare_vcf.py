#!/usr/bin/env python3
"""Subset the annotated VCF to rare HIGH/MODERATE plus spliceai keys for Exomiser.

Prefilter is documented. It is not a diagnosis.
"""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

RARE = 0.01
KEEP_IMPACT = {"HIGH", "MODERATE"}


def info_map(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif item:
            out[item] = ""
    return out


def load_keys(path: Path) -> set[tuple[str, str, str, str]]:
    keys: set[tuple[str, str, str, str]] = set()
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            keys.add((rec["chrom"], rec["pos"], rec["ref"], rec["alt"]))
    return keys


def af_of(info: dict[str, str]) -> float | None:
    for key in ("gnomad_AF", "AF", "gnomAD_AF"):
        raw = info.get(key)
        if raw in (None, "", ".", "nan"):
            continue
        try:
            return float(raw.split(",")[0])
        except ValueError:
            continue
    return None


def impact_of(info: dict[str, str]) -> str:
    if info.get("IMPACT") in KEEP_IMPACT:
        return info["IMPACT"]
    ann = info.get("ANN", "")
    if ann:
        parts = ann.split(",", 1)[0].split("|")
        if len(parts) > 2:
            return parts[2]
    csq = info.get("CSQ", "")
    if csq:
        # VEP: Allele|Consequence|IMPACT|...
        parts = csq.split(",", 1)[0].split("|")
        if len(parts) > 2:
            return parts[2]
    return ""


def main() -> None:
    vcf_in, keys, out_vcf = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    keep = load_keys(keys)
    n_keep = n_rare = n_key = 0
    with gzip.open(vcf_in, "rt", encoding="utf-8") as fin, gzip.open(
        out_vcf, "wt", encoding="utf-8"
    ) as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.split("\t")
            chrom, pos, ref, alt0 = f[0], f[1], f[3], f[4].split(",")[0]
            info = info_map(f[7])
            key = (chrom, pos, ref, alt0)
            take = False
            if key in keep:
                take = True
                n_key += 1
            else:
                impact = impact_of(info)
                af = af_of(info)
                if impact in KEEP_IMPACT and af is not None and af <= RARE:
                    take = True
                    n_rare += 1
            if take:
                fout.write(line)
                n_keep += 1
    Path(str(out_vcf) + ".counts.txt").write_text(
        f"keep={n_keep}\nrare_hm={n_rare}\nkeys={n_key}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
