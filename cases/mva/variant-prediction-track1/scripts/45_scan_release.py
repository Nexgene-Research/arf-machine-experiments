"""Scan a staged public release for controlled files and credential patterns.

Scans bounded text; refuses binary genomes and large files without reading them.
"""

from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RELEASE = ROOT if (ROOT / "submission").is_dir() else ROOT / "release/Nexgene-AI-track1"
DENY_EXT = (".fastq", ".fq", ".bam", ".cram", ".vcf", ".bcf", ".docx", ".tbi", ".csi")
PATTERNS = (
    re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(rb"(?i)(?:hf_|ghp_|github_pat_)[a-zA-Z0-9_]{24,}"),
    re.compile(
        rb"(?i)(?:api[_-]?key|access[_-]?token|client[_-]?secret)\s*[:=]\s*['\"]?[A-Za-z0-9._-]{20,}"
    ),
    re.compile(rb"(?i)(?:transfer1|alogin1|MN5-ACC)\.bsc\.es"),
    re.compile(rb"(?i)(?:@(?:transfer1|alogin1)\.bsc\.es|BEGIN PGP " rb"PRIVATE)"),
    re.compile(rb"/gpfs/scratch/ehpc926/nexi928147"),
    re.compile(rb"(?i)\b(?:Astra|GPT-6|Azure|OpenScience)\b"),
)


def main() -> int:
    issues = []
    if not RELEASE.is_dir():
        print("missing release")
        return 1
    for path in RELEASE.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        rel = path.relative_to(RELEASE).as_posix()
        if path.stat().st_size > 50_000_000 or any(
            part in path.name.lower() for part in DENY_EXT
        ):
            issues.append(f"prohibited filename/size: {rel}")
            continue
        if path.suffix in {".xlsx", ".png", ".pdf"}:
            continue
        if path.name == "45_scan_release.py":
            continue
        content = path.read_bytes()
        for rx in PATTERNS:
            if rx.search(content):
                issues.append(f"credential/host pattern: {rel}")
                break
        if rel.startswith("results/") and (
            b"WGS_EX2312012" in content or b"@HD\tVN" in content
        ):
            issues.append(f"restricted patient source metadata: {rel}")
        if rel.endswith("_track1_report.md") and b"v1 ranked" in content.lower():
            issues.append(f"process wording in report: {rel}")
    print(
        "OK: controlled-data and credential scan clean"
        if not issues
        else "\n".join(issues)
    )
    return int(bool(issues))


if __name__ == "__main__":
    sys.exit(main())
