"""Validate a Track 1 submission CSV against the template and GRCh38.

Usage: env\\venv\\Scripts\\python.exe scripts/validate_submission.py <csv> ref/GRCh38_no_alt.fa
Exit 0 only when every check passes.
"""

import csv
import math
import pathlib
import sys

from pyfaidx import Fasta

ROOT = pathlib.Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "submission" / "track1_submission_template.csv"

if len(sys.argv) != 3:
    sys.exit("usage: validate_submission.py <csv> <GRCh38 fasta>")

csv_path, fasta_path = sys.argv[1], sys.argv[2]
with open(TEMPLATE, newline="", encoding="utf-8-sig") as handle:
    columns = next(csv.reader(handle))
with open(csv_path, newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    header = reader.fieldnames or []
    rows = list(reader)
fasta = Fasta(fasta_path)
errors: list[str] = []

if header != columns:
    errors.append(f"columns differ from template: {header}")
if not 1 <= len(rows) <= 10:
    errors.append(f"{len(rows)} rows (1-10 allowed)")

pairs: set[tuple] = set()
for line, row in enumerate(rows, start=2):
    get = lambda key: (row.get(key) or "").strip()
    if not get("proband_id"):
        errors.append(f"line {line}: empty proband_id")
    alleles = []
    for n in ("1", "2"):
        fields = [get(f"{key}_{n}") for key in ("chrom", "pos", "ref", "alt")]
        if n == "2" and not any(fields):
            continue
        chrom, pos, ref, alt = fields
        if not (chrom.startswith("chr") and pos.isdigit() and int(pos) > 0 and ref and alt and ref != alt):
            errors.append(f"line {line}: malformed variant {n}: {fields}")
            continue
        if set(ref.upper() + alt.upper()) - set("ACGT"):
            errors.append(f"line {line}: non-ACGT allele in variant {n}")
            continue
        if chrom not in fasta:
            errors.append(f"line {line}: contig {chrom} not in reference")
            continue
        start = int(pos) - 1
        genome = fasta[chrom][start : start + len(ref)].seq.upper()
        if genome != ref.upper():
            errors.append(f"line {line}: ref_{n}={ref} but GRCh38 has {genome} at {chrom}:{pos}")
        if len(ref) != len(alt) and ref[0].upper() != alt[0].upper():
            errors.append(f"line {line}: indel in variant {n} lacks the shared VCF padding base")
        alleles.append((chrom, int(pos), ref.upper(), alt.upper()))
    if len(alleles) == 2 and alleles[0] == alleles[1]:
        errors.append(f"line {line}: variant 1 and variant 2 are identical")
    key = tuple(sorted(alleles))
    if key and key in pairs:
        errors.append(f"line {line}: duplicate row")
    pairs.add(key)
    try:
        epcr = float(get("epcr"))
    except ValueError:
        epcr = math.nan
    if not 0 < epcr <= 1:
        errors.append(f"line {line}: epcr {get('epcr')!r} outside (0, 1]")
    if get("finding_type") not in ("primary", "secondary"):
        errors.append(f"line {line}: finding_type {get('finding_type')!r}")

primary = [row for row in rows if (row.get("finding_type") or "").strip() == "primary"]
if rows and not primary:
    errors.append("no primary row")


def score(row: dict) -> float:
    try:
        return float((row.get("epcr") or "").strip())
    except ValueError:
        return 0.0


scores = [score(row) for row in primary]
if scores != sorted(scores, reverse=True):
    errors.append("primary rows are not ordered by epcr, highest first")

print("\n".join(errors) or f"OK: {len(rows)} rows ({len(primary)} primary)")
sys.exit(1 if errors else 0)
