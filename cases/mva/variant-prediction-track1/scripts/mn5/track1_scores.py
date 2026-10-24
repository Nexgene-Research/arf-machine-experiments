"""Join AlphaMissense and REVEL onto the uncapped candidate table.

Streams the full predictor files once and keeps scores only for exact
chrom-pos-ref-alt keys already in candidates_full.tsv. Scores stay
descriptive: the preregistration assigns them LR=1 until calibrated.
"""

from __future__ import annotations

import csv
import gzip
import sys
import zipfile
from collections import Counter
from pathlib import Path


def chrom(value: str) -> str:
    value = value.strip()
    if value.startswith("chr"):
        return "chrM" if value in ("chrMT",) else value
    if value in ("M", "MT"):
        return "chrM"
    return "chr" + value


def key_of(parts: list[str]) -> tuple[str, str, str, str]:
    return (chrom(parts[0]), parts[1], parts[2].upper(), parts[3].upper())


def main() -> None:
    root = Path(sys.argv[1])
    out = Path(sys.argv[2])
    candidates = root / "candidates_full.tsv"
    am_path = Path(sys.argv[3])
    revel_zip = Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)

    with candidates.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    wanted = {key_of([row["chrom"], row["pos"], row["ref"], row["alt"]]) for row in rows}

    am: dict[tuple[str, str, str, str], tuple[float, str, int]] = {}
    with gzip.open(am_path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            key = key_of(parts[:4])
            if key not in wanted:
                continue
            try:
                score = float(parts[8])
            except ValueError:
                continue
            current = am.get(key)
            if current is None or score > current[0]:
                am[key] = (score, parts[9], 1 if current is None else current[2] + 1)
            else:
                am[key] = (current[0], current[1], current[2] + 1)

    revel: dict[tuple[str, str, str, str], tuple[float, int]] = {}
    with zipfile.ZipFile(revel_zip) as archive:
        with archive.open("revel_with_transcript_ids") as raw:
            header = True
            for line in raw:
                text = line.decode("utf-8", "replace").rstrip("\n")
                if header:
                    header = False
                    continue
                parts = text.split(",")
                if len(parts) < 8:
                    continue
                key = (chrom(parts[0]), parts[2], parts[3].upper(), parts[4].upper())
                if key not in wanted:
                    continue
                try:
                    score = float(parts[7])
                except ValueError:
                    continue
                current = revel.get(key)
                if current is None or score > current[0]:
                    revel[key] = (score, 1 if current is None else current[1] + 1)
                else:
                    revel[key] = (current[0], current[1] + 1)

    extra = ["am_pathogenicity", "am_class", "am_transcripts", "revel", "revel_transcripts"]
    coverage: Counter[tuple[str, str]] = Counter()
    with (out / "candidates_scored.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields + extra, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            key = key_of([row["chrom"], row["pos"], row["ref"], row["alt"]])
            hit_am = am.get(key)
            hit_revel = revel.get(key)
            row["am_pathogenicity"] = f"{hit_am[0]:.4f}" if hit_am else ""
            row["am_class"] = hit_am[1] if hit_am else ""
            row["am_transcripts"] = str(hit_am[2]) if hit_am else "0"
            row["revel"] = f"{hit_revel[0]:.3f}" if hit_revel else ""
            row["revel_transcripts"] = str(hit_revel[1]) if hit_revel else "0"
            writer.writerow(row)
            consequence = row.get("consequence") or "."
            coverage[(consequence, "am_present" if hit_am else "am_missing")] += 1
            coverage[(consequence, "revel_present" if hit_revel else "revel_missing")] += 1

    with (out / "score_coverage.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["consequence", "status", "count"])
        for (consequence, status), count in sorted(coverage.items()):
            writer.writerow([consequence, status, count])
    print(
        f"candidates={len(rows)} am_hit={len(am)} revel_hit={len(revel)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
