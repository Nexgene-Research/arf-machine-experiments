"""Fail-closed Track 1 CSV writer; requires frozen A and adjudicated ranked pairs.

ranked_pairs.tsv must retain the exact 30_rank_pairs.py A columns and append
epcr,submission_ready,proband_id,acmg_criteria_1,acmg_criteria_2. The latter
criteria may be empty if explicitly reviewed; Stage 32 supplies heuristic
ePCR. Every proposed phenotype explanation is primary. No imputation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import importlib

ranking = importlib.import_module("30_rank_pairs")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", type=Path, default=ROOT / "results/ranked_pairs.tsv")
    ap.add_argument(
        "--candidates", type=Path, default=ROOT / "results/candidates_full.tsv"
    )
    ap.add_argument("--frozen-a", type=Path, default=ROOT / "results/ranked_a.tsv")
    ap.add_argument(
        "--manifest", type=Path, default=ROOT / "results/ranked_a_manifest.json"
    )
    ap.add_argument("--reference", type=Path, default=ROOT / "ref/GRCh38_no_alt.fa")
    ap.add_argument(
        "--out", type=Path, default=ROOT / "sub/Nexgene-Research_wgs-discovery.csv"
    )
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(f"Submission already exists; never overwrite: {args.out}")
    if not args.manifest.is_file() or args.manifest.stat().st_size > 1_000_000:
        raise ValueError("Missing or oversized frozen-A manifest")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if (
        not manifest.get("no_B_C_D_consumed")
        or ranking.digest(args.frozen_a) != manifest["ranked_a_sha256"]
    ):
        raise ValueError(
            "Frozen Angle A hash mismatch or angle independence not documented"
        )
    ranked = ranking.bounded_table(
        args.ranked,
        (
            *ranking.FIELDS,
            "epcr",
            "submission_ready",
            "proband_id",
            "acmg_criteria_1",
            "acmg_criteria_2",
        ),
    )
    # Consume all ranked rows and emit at most ten evidence-supported rows.
    frozen = ranking.bounded_table(args.frozen_a, ranking.FIELDS)
    candidates = ranking.bounded_table(
        args.candidates,
        (
            "gene",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gt",
            "genotype_status",
            "normalization_status",
            "ref_matches_grch38",
        ),
    )
    if ranking.digest(args.candidates) != manifest["inputs_sha256"].get(
        str(args.candidates.resolve())
    ):
        raise ValueError("Candidate input changed since Angle A freeze")
    a_ids = {r["pair_id"]: r for r in frozen}
    if len(a_ids) != len(frozen):
        raise ValueError("Duplicated frozen-A pair IDs")
    c_keys = {ranking.allele_key(r): r for r in candidates}
    if len(c_keys) != len(candidates):
        raise ValueError("Duplicate candidate variant IDs")
    if args.reference.resolve() != Path(manifest["reference_path"]).resolve():
        raise ValueError("Reference changed since frozen Angle A")
    index = Path(str(args.reference) + ".fai")
    if not index.is_file() or (
        manifest["reference_index_sha256"] != "not_available_synthetic_reference"
        and ranking.digest(index) != manifest["reference_index_sha256"]
    ):
        raise ValueError("Reference .fai missing or changed since frozen Angle A")
    fasta = Fasta(str(args.reference), as_raw=True, sequence_always_upper=True)
    if not args.out.parent.is_dir():
        raise FileNotFoundError(f"Submission directory absent: {args.out.parent}")
    seen, ready, prev = set(), [], float("inf")
    for r in ranked:
        pair_id = r["pair_id"]
        if pair_id in seen or pair_id not in a_ids:
            raise ValueError(f"Unknown or duplicate ranked pair ID: {pair_id}")
        seen.add(pair_id)
        source = a_ids[pair_id]
        if any(
            r[k] != source[k]
            for k in ranking.FIELDS
            if k not in ("p_path_1", "p_path_2")
        ):
            raise ValueError(f"Stage 32 changed frozen-A evidence: {pair_id}")
        if r["submission_ready"] not in ("yes", "no"):
            raise ValueError(
                f"Stage 32 must explicitly adjudicate submission_ready=yes/no: {pair_id}"
            )
        if not r["epcr"].strip():
            raise ValueError(f"Missing adjudicated epcr: {pair_id}")
        value = ranking.probability(r["epcr"], "epcr")
        if value > prev:
            raise ValueError("ranked_pairs.tsv is not ordered by descending epcr")
        prev = value
        if not r["proband_id"].strip():
            raise ValueError(f"Missing non-synthetic proband_id: {pair_id}")
        for n in (1, 2):
            key = r[f"variant_{n}"]
            if n == 2 and not key:
                continue  # homozygous/hemizygous/dominant hypothesis has ONE allele
            if (
                key not in c_keys
                or r[f"normalization_status_{n}"] != "bcftools_norm_GRCh38"
            ):
                raise ValueError(f"Unverified candidate allele: {key}")
            ranking.verify_allele(key, fasta)
        if r["submission_ready"] == "yes" and r.get("phase") != "cis":
            ready.append(r)
    if not ready:
        raise ValueError(
            "No ranked hypothesis with normalized GRCh38 alleles and positive ePCR"
        )
    ready = ready[:10]
    template = ROOT / "sub/track1_submission_template.csv"
    with template.open(newline="", encoding="utf-8-sig") as f:
        columns = next(csv.reader(f))
    expected = [
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
    if columns != expected:
        raise ValueError("Unexpected submission template schema")
    # Write in system scratch; only publish output after project's real validator exits 0.
    with tempfile.TemporaryDirectory(
        prefix=".track1_validation_", dir=args.out.parent
    ) as directory:
        trial = Path(directory) / "candidate.csv"
        with trial.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            for r in ready:
                row = {column: r[column] for column in columns if column in r}
                row["finding_type"] = (
                    "secondary"
                    if r.get("hypothesis_track") == "incidental"
                    else "primary"
                )
                row["notes"] = (
                    f"pair_id={r['pair_id']};phase={r['phase']};epcr_status={r.get('epcr_status', 'provisional')};Angle A frozen"
                )
                writer.writerow(row)
        command = [
            str(ROOT / "env/venv/Scripts/python.exe"),
            str(ROOT / "scripts/validate_submission.py"),
            str(trial),
            str(args.reference),
        ]
        check = subprocess.run(
            command, cwd=str(ROOT), capture_output=True, text=True, check=False
        )
        if check.returncode:
            raise ValueError(
                f"Project validator exit {check.returncode}: {(check.stdout + check.stderr).strip()}"
            )
        os.replace(trial, args.out)
    print(
        f"Validated {len(ready)} genuine rows; project validator exit=0; submission={args.out}"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as error:
        sys.exit(f"ERROR: {error}")
