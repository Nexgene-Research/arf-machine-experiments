"""Run all Track 1 freeze checks and commit a scanned local-only release.

Does not push or mutate legacy v1. Writes a new final receipt only on success.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release/Nexgene-AI-track1"
RECEIPT = ROOT / "results/track1/final_release_receipt.json"
PYTHON = ROOT / "env/venv/Scripts/python.exe"
CSV = ROOT / "sub/Nexgene-AI_mn5-evidence-v2.csv"
REPORT = ROOT / "results/Nexgene-AI_track1_report.md"
BOOK = ROOT / "sub/Nexgene-AI_methods_description.xlsx"
V1 = ROOT / "sub/Nexgene-Research_wgs-discovery.csv"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked(cmd: list[str], cwd: Path = ROOT) -> str:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    output = (result.stdout + result.stderr).strip()
    if result.returncode:
        raise RuntimeError(
            f"exit={result.returncode} command={cmd[0]} {cmd[1:]}: {output[:2000]}"
        )
    print(f"exit=0 {cmd[0]} {cmd[1:]}: {output[:250]}")
    return output


def main() -> None:
    if RECEIPT.exists():
        raise FileExistsError(RECEIPT)
    if sha(V1) != "41bd90a7fbc9738ad670943a41bc5ee7437d6902c5de5b3dde1ecc844b9bcecb":
        raise ValueError("fallback v1 changed")
    validated = checked(
        [
            str(PYTHON),
            "scripts/validate_submission.py",
            str(CSV),
            "ref/GRCh38_no_alt.fa",
        ]
    )
    release_checked = checked([str(PYTHON), "scripts/36_release_checks.py"])
    scan = checked([str(PYTHON), "scripts/45_scan_release.py"])
    with (RELEASE / "MANIFEST.sha256").open(encoding="utf-8") as fh:
        manifest = [line.rstrip("\n").split("  ", 1) for line in fh]
    for expected, name in manifest:
        if expected != sha(RELEASE / name):
            raise ValueError(f"manifest mismatch: {name}")
    spec = importlib.util.spec_from_file_location(
        "release_tests", RELEASE / "tests/test_synthetic.py"
    )
    if spec is None or spec.loader is None:
        raise ImportError("release synthetic tests cannot be imported")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in sorted(vars(module)):
        if name.startswith("test_"):
            getattr(module, name)()
    tests = [name for name in vars(module) if name.startswith("test_")]
    with CSV.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError("no rows")
    review = (ROOT / "results/redteam/internal_review_v2.md").read_text(encoding="utf-8")
    response = (ROOT / "results/redteam/responses_v2.md").read_text(encoding="utf-8")
    blocking = "**BLOCKING" in review and "Scientific freeze is a no-go" in response
    if not blocking:
        raise ValueError(
            "scientific disposition not documented; do not certify submission"
        )
    if (RELEASE / ".git").exists():
        raise FileExistsError(
            "release git already initialized; inspect it rather than overwrite"
        )
    checked(["git", "init"], RELEASE)
    checked(["git", "symbolic-ref", "HEAD", "refs/heads/release-v2"], RELEASE)
    intended = [name for _, name in manifest] + ["MANIFEST.sha256"]
    checked(["git", "add", "--", *intended], RELEASE)
    status = checked(["git", "status", "--short"], RELEASE)
    checked(["git", "diff", "--cached", "--check"], RELEASE)
    checked(
        ["git", "commit", "-m", "Track 1 audit archive - no scientific freeze"], RELEASE
    )
    commit = checked(["git", "rev-parse", "HEAD"], RELEASE)
    files = {
        path.name: sha(path)
        for path in (CSV, REPORT, BOOK, RELEASE / "MANIFEST.sha256")
    }
    info = {
        "filenames_sha256": files,
        "row_count": len(rows),
        "top_gene": rows[0]["notes"].split(";")[0].removeprefix("gene="),
        "phase_status": next(
            (
                v.removeprefix("phase=")
                for v in rows[0]["notes"].split(";")
                if v.startswith("phase=")
            ),
            "unknown",
        ),
        "validator_output": validated,
        "release_check_output": release_checked,
        "tests": {name: "passed" for name in tests},
        "controlled_data_scan": scan,
        "red_team_status": "one fresh-context internal review completed early; second-allele eligibility and score objections recorded",
        "scientific_freeze_status": "NO-GO: v2 CSV format-valid, not frozen-eligible",
        "files_cleared_for_upload": [],
        "operational_fallback": "sub/Nexgene-Research_wgs-discovery.csv preserved and validator-accepted; old intronic causal claims not endorsed",
        "session_spend_usd": "unavailable: harness did not expose a mission cost meter",
        "v1_sha256_unchanged": sha(V1),
        "commit_hash": commit,
        "git_status_pre_commit": status,
    }
    RECEIPT.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(
        f"archive checks passed, SCIENTIFIC FREEZE NO-GO; receipt {RECEIPT}; commit {commit}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        sys.exit(f"FREEZE FAILED: {err}")
