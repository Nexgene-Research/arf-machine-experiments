#!/usr/bin/env python3
"""Build a tiny VCF of jury + rare-noncoding keys and score with SpliceAI.

Missing score is not benign. This table does not establish causality.
"""

from __future__ import annotations

import os

import csv
import gzip
import subprocess
import sys
from pathlib import Path


def contig_headers(ref: Path) -> list[str]:
    fai = Path(str(ref) + ".fai")
    lines = []
    for line in fai.read_text(encoding="utf-8").splitlines():
        chrom, length, *_ = line.split("\t")
        lines.append(f"##contig=<ID={chrom},length={length}>")
    if not lines:
        raise FileNotFoundError(f"empty or missing fasta index: {fai}")
    return lines


def write_vcf(keys: Path, ref: Path, out_vcf: Path) -> list[dict[str, str]]:
    records = []
    with keys.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            records.append(rec)
    lines = [
        "##fileformat=VCFv4.2",
        f"##reference={ref}",
        *contig_headers(ref),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for rec in records:
        lines.append(
            f"{rec['chrom']}\t{rec['pos']}\t{rec['label']}\t{rec['ref']}\t{rec['alt']}\t.\t.\t."
        )
    out_vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return records


def parse_spliceai(info: str) -> dict[str, str]:
    for item in info.split(";"):
        if item.startswith("SpliceAI="):
            # ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            parts = item.split("=", 1)[1].split(",")[0].split("|")
            if len(parts) >= 6:
                return {
                    "allele": parts[0],
                    "symbol": parts[1],
                    "DS_AG": parts[2],
                    "DS_AL": parts[3],
                    "DS_DG": parts[4],
                    "DS_DL": parts[5],
                }
    return {}


def main() -> None:
    keys, ref, out_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf = out_dir / "keys.vcf"
    records = write_vcf(keys, ref, vcf)
    scored = out_dir / "keys.spliceai.vcf"
    shared = Path(os.environ["MVA_BASE"]) / "jobs/track1_wave3/spliceai_venv"
    venv = shared if (shared / "bin" / "python").is_file() else (out_dir / "venv")
    py = venv / "bin" / "python"
    wheels = Path(os.environ["MVA_BASE"]) / "jobs/track1_wave3/wheels"
    if not py.is_file():
        subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
        pip = venv / "bin" / "pip"
        pip_cmd = [str(pip), "install", "--quiet"]
        if wheels.is_dir():
            pip_cmd += ["--no-index", "--find-links", str(wheels)]
        pip_cmd += ["spliceai", "tensorflow-cpu"]
        proc = subprocess.run(
            pip_cmd,
            capture_output=True,
            timeout=3600,
        )
        (out_dir / "pip.log").write_bytes(proc.stdout + proc.stderr)
        if proc.returncode:
            tsv = out_dir / "spliceai_keys.tsv"
            with tsv.open("w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(
                    fh,
                    fieldnames=[
                        "label",
                        "chrom",
                        "pos",
                        "ref",
                        "alt",
                        "gene",
                        "DS_AG",
                        "DS_AL",
                        "DS_DG",
                        "DS_DL",
                        "max_delta",
                        "status",
                        "note",
                    ],
                    delimiter="\t",
                )
                w.writeheader()
                for rec in records:
                    w.writerow(
                        {
                            **{k: rec[k] for k in ("label", "chrom", "pos", "ref", "alt", "gene")},
                            "DS_AG": ".",
                            "DS_AL": ".",
                            "DS_DG": ".",
                            "DS_DL": ".",
                            "max_delta": ".",
                            "status": "blocked_pip",
                            "note": "spliceai/tensorflow-cpu install failed; missing is not benign",
                        }
                    )
            (out_dir / "caveat.txt").write_text(
                "SpliceAI pip install failed. See pip.log. Scores not invented.\n",
                encoding="utf-8",
            )
            return
    cmd = [
        str(py),
        "-m",
        "spliceai",
        "-I",
        str(vcf),
        "-O",
        str(scored),
        "-R",
        str(ref),
        "-A",
        "grch38",
    ]
    proc = subprocess.run(cmd, capture_output=True, timeout=3600)
    (out_dir / "spliceai.log").write_bytes(proc.stdout + proc.stderr)
    scored_map: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if scored.is_file():
        text = scored.read_text(encoding="utf-8")
        for line in text.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            parsed = parse_spliceai(f[7] if len(f) > 7 else "")
            scored_map[(f[0], f[1], f[3], f[4])] = parsed
    tsv = out_dir / "spliceai_keys.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "label",
                "chrom",
                "pos",
                "ref",
                "alt",
                "gene",
                "DS_AG",
                "DS_AL",
                "DS_DG",
                "DS_DL",
                "max_delta",
                "status",
                "note",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for rec in records:
            key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
            parsed = scored_map.get(key, {})
            deltas = []
            for name in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
                try:
                    value = float(parsed.get(name, "nan"))
                except ValueError:
                    continue
                if value == value:
                    deltas.append(value)
            max_d = max(deltas) if deltas else None
            w.writerow(
                {
                    **{k: rec[k] for k in ("label", "chrom", "pos", "ref", "alt", "gene")},
                    "DS_AG": parsed.get("DS_AG", "."),
                    "DS_AL": parsed.get("DS_AL", "."),
                    "DS_DG": parsed.get("DS_DG", "."),
                    "DS_DL": parsed.get("DS_DL", "."),
                    "max_delta": f"{max_d:.4f}" if max_d is not None else ".",
                    "status": "scored" if parsed else "missing",
                    "note": "SpliceAI; missing is not benign; not a causal claim",
                }
            )
    (out_dir / "caveat.txt").write_text(
        "SpliceAI on a bounded key set only. Missing row is not a normal splice site.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
