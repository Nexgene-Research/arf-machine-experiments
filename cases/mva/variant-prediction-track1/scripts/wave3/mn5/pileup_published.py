#!/usr/bin/env python3
"""Genotype newly mapped published MVA alleles from this BAM.

Does not re-query the four SNVs already pileup'd alt=0.
Missing call is not absence. This table does not establish causality.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    sites, bam, ref, bcf_sif, out_dir = (Path(a) for a in sys.argv[1:6])
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf_out = out_dir / "published_called.vcf.gz"
    cmd = [
        "singularity",
        "exec",
        "--cleanenv",
        str(bcf_sif),
        "bcftools",
        "mpileup",
        "-f",
        str(ref),
        "-T",
        str(sites),
        "-a",
        "FORMAT/AD,FORMAT/DP",
        "-Ou",
        str(bam),
    ]
    call = [
        "singularity",
        "exec",
        "--cleanenv",
        str(bcf_sif),
        "bcftools",
        "call",
        "-m",
        "-Oz",
        "-o",
        str(vcf_out),
    ]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.run(call, stdin=p1.stdout, capture_output=True)
    if p1.stdout:
        p1.stdout.close()
    p1.wait()
    (out_dir / "pileup.stderr.txt").write_bytes(p2.stderr)
    if p2.returncode:
        raise SystemExit(f"bcftools call failed: {p2.stderr[-800:]!r}")
    subprocess.run(
        [
            "singularity",
            "exec",
            "--cleanenv",
            str(bcf_sif),
            "bcftools",
            "index",
            "-t",
            str(vcf_out),
        ],
        check=True,
    )
    wanted = {}
    with sites.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, ident, ref, alt, *_ = line.rstrip("\n").split("\t")
            wanted[(chrom, pos, ref, alt)] = ident
    wanted_pos = {(c, p): ident for (c, p, _r, _a), ident in wanted.items()}
    rows = []
    opener = gzip.open if vcf_out.suffixes[-1] == ".gz" else open
    with opener(vcf_out, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ident, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            fmt = fields[8].split(":")
            smp = fields[9].split(":")
            rec = dict(zip(fmt, smp))
            key = (chrom, pos, ref, alt)
            pos_key = (chrom, pos)
            label = wanted.get(key) or wanted_pos.get(pos_key, ident)
            rows.append(
                {
                    "label": label,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gt": rec.get("GT", "."),
                    "dp": rec.get("DP", "."),
                    "ad": rec.get("AD", "."),
                    "filter": fields[6],
                    "note": "native genotype; not BAMSurgeon; not a causal claim",
                }
            )
    seen_pos = {(r["chrom"], r["pos"]) for r in rows}
    for key, ident in wanted.items():
        chrom, pos, ref, alt = key
        if (chrom, pos) in seen_pos:
            continue
        rows.append(
                {
                    "label": ident,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gt": "no_call",
                    "dp": ".",
                    "ad": ".",
                    "filter": ".",
                    "note": "not emitted by bcftools call; missing is not absence",
                }
            )
    tsv = out_dir / "published_site_gt_new.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["label", "chrom", "pos", "ref", "alt", "gt", "dp", "ad", "filter", "note"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    (out_dir / "inputs.sha256").write_text(
        f"{sha256(sites)}  {sites.name}\n",
        encoding="utf-8",
    )
    (out_dir / "caveat.txt").write_text(
        "Native pileup of newly mapped published alleles only. "
        "Does not re-query the four SNVs already alt=0. Not a diagnosis.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
