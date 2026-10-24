"""Write the Track 1 CSV from the MN5 top-10 ranking and run the validator.

Reads only results/track1_top10_ranked.tsv (at most ten rows). Does not load
the full hypothesis table. Never overwrites an existing submission.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "sub" / "track1_submission_template.csv"
RANKED = ROOT / "results" / "track1_top10_ranked.tsv"
OUT = ROOT / "sub" / "Nexgene-Research_wgs-discovery.csv"
REF = ROOT / "ref" / "GRCh38_no_alt.fa"
COLUMNS = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Submission already exists: {OUT}")
    with TEMPLATE.open(newline="", encoding="utf-8-sig") as handle:
        columns = next(csv.reader(handle))
    if columns != COLUMNS:
        raise ValueError(f"Unexpected template columns: {columns}")
    with RANKED.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not 1 <= len(rows) <= 10:
        raise ValueError(f"{len(rows)} ranked rows (need 1-10)")
    ready = []
    prev = float("inf")
    for row in rows:
        epcr = float(row["epcr"])
        if not 0 < epcr <= 1:
            raise ValueError(f"epcr out of range: {row['epcr']!r}")
        if epcr > prev + 1e-12:
            raise ValueError("top10 is not ordered by descending epcr")
        prev = epcr
        chrom2 = (row.get("chrom_2") or "").strip()
        notes = ";".join(
            [
                f"gene={row.get('gene', '')}",
                f"hypothesis={row.get('hypothesis', '')}",
                f"phase={row.get('phase_called', row.get('phase', ''))}",
                f"capped={row.get('capped', '')}",
                f"reads_ok={row.get('reads_ok', '')}",
                "prior=0.1",
                "heuristic_epcr",
            ]
        )
        ready.append(
            {
                "proband_id": "PROBAND01",
                "chrom_1": row["chrom_1"],
                "pos_1": row["pos_1"],
                "ref_1": row["ref_1"],
                "alt_1": row["alt_1"],
                "chrom_2": chrom2,
                "pos_2": row["pos_2"] if chrom2 else "",
                "ref_2": row["ref_2"] if chrom2 else "",
                "alt_2": row["alt_2"] if chrom2 else "",
                "epcr": row["epcr"],
                "finding_type": "primary",
                "notes": notes,
            }
        )
    python = ROOT / "env" / "venv" / "Scripts" / "python.exe"
    validator = ROOT / "scripts" / "validate_submission.py"
    with tempfile.TemporaryDirectory(prefix=".track1_validation_", dir=OUT.parent) as directory:
        trial = Path(directory) / "candidate.csv"
        with trial.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(ready)
        check = subprocess.run(
            [str(python), str(validator), str(trial), str(REF)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        if check.returncode:
            raise ValueError(
                f"validator exit {check.returncode}: {(check.stdout + check.stderr).strip()}"
            )
        os.replace(trial, OUT)
    print(f"OK: {len(ready)} primary rows; {OUT}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, FileExistsError) as error:
        sys.exit(f"ERROR: {error}")
