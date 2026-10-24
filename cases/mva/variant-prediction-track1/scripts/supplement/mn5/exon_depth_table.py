#!/usr/bin/env python3
"""Turn a mosdepth regions BED into a small TSV."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path


def main() -> None:
    src, dest = Path(sys.argv[1]), Path(sys.argv[2])
    lines = ["chrom\tstart\tend\tinterval\tdepth\n"]
    with gzip.open(src, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            lines.append("\t".join(cols[:5]) + "\n")
    dest.write_text("".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
