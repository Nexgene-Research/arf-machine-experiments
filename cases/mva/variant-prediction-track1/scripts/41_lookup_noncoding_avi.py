"""Phase D2: tabix-indexed local AVI lookup of MN5's bounded SNV keys.

No model API or whole-table scan. Reuses the audited BGZF/TBI reader in 28.
"""

from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/track1/rare_noncoding_keys.tsv"
OUT = ROOT / "results/track1/noncoding_avi.tsv"
MODULE = ROOT / "scripts/28_alphagenome_bounded.py"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    spec = importlib.util.spec_from_file_location("avi_bounded", MODULE)
    avi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(avi)
    with INPUT.open(newline="", encoding="utf-8") as fh:
        keys = list(csv.DictReader(fh, delimiter="\t"))
    if not 1 <= len(keys) <= 100:
        raise ValueError("expected bounded, nonempty list of SNV keys")
    index, _ = avi.read_tbi(avi.AVI)
    found = 0
    with avi.AVI.open("rb") as bgzf, OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "gene",
                "grch38_key",
                "avi_raw_score",
                "avi_PHRED",
                "top_percentile",
                "lookup_status",
                "likelihood_ratio",
            )
        )
        for row in keys:
            chrom, pos, ref, alt = row["grch38_key"].split("-")
            if not (len(ref) == len(alt) == 1):
                raise ValueError("D2 requires SNV key")
            record, reason, _ = avi.avi_at(index, bgzf, chrom, int(pos), ref, alt)
            if record:
                phred = float(record[5])
                if not math.isfinite(phred):
                    raise ValueError("invalid AVI PHRED")
                percentile = f"{100 * 10 ** (-phred / 10):.8g}"
                found += 1
            else:
                percentile = ""
            writer.writerow(
                (
                    row["gene"],
                    row["grch38_key"],
                    record[4] if record else "",
                    record[5] if record else "",
                    percentile,
                    "exact_allele" if record else reason,
                    "1",
                )
            )
    print(f"D2 tabix lookups={len(keys)} matched={found}; uncalibrated LR=1")


if __name__ == "__main__":
    main()
