"""Stage, submit, poll and collect MN5 rigor-pass jobs.

Uses Git OpenSSH. Does not write into release/Nexgene-AI-track1.
"""

from __future__ import annotations

import os

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import subprocess

DEST = pathlib.Path(__file__).resolve().parents[1]
SSH = r"C:\Program Files\Git\usr\bin\ssh.exe"
SFTP = r"C:\Program Files\Git\usr\bin\sftp.exe"
OPTS = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=30"]
REMOTE_JOBS = str(pathlib.Path(os.environ["MVA_BASE"]) / "jobs/track1_rigor")
REMOTE_OUT = str(pathlib.Path(os.environ["MVA_BASE"]) / "out/track1_rigor")
HOST_LOGIN = "<user>@<mn5-login>"
HOST_XFER = "<user>@<mn5-transfer>"

UPLOAD = [
    "scripts/mn5/common.sh",
    "scripts/mn5/rigor_scan.py",
    "scripts/mn5/rigor_pid.py",
    "scripts/mn5/rigor_roh.py",
    "scripts/mn5/rigor_splice.py",
    "scripts/mn5/rigor_clinvar.py",
    "scripts/mn5/rigor_scan.sbatch",
    "scripts/mn5/rigor_pid.sbatch",
    "scripts/mn5/rigor_roh.sbatch",
    "scripts/mn5/rigor_roh_summarize.sbatch",
    "scripts/mn5/rigor_splice.sbatch",
    "scripts/mn5/rigor_clinvar.sbatch",
]

REFS = [
    "refs/ACMG_SF_v3.2.genes.txt",
    "refs/tier1_sixgene_coding_splice20.bed",
    "refs/alphagenome_splicing_sixgene_slice.tsv",
    "refs/alphagenome_splicing_sixgene_background.tsv",
]

COLLECT = (
    "vaf_counts.tsv",
    "vaf_by_gene_six.tsv",
    "vaf_jury_alleles.tsv",
    "clinvar_plp_summary.tsv",
    "clinvar_plp_genes.tsv",
    "clinvar_plp_jury.tsv",
    "rare_counts.tsv",
    "rare_sixgene.tsv",
    "rare_acmg_sf.tsv",
    "pid_pair_counts.tsv",
    "pid_sixgene.tsv",
    "pid_jury.tsv",
    "roh_summary.tsv",
    "splice_sixgene_counts.tsv",
    "splice_jury.tsv",
    "scan_counts.json",
    "python.version",
    "bcftools.version",
    "inputs.sha256",
    "caveat.txt",
    "exit_code.txt",
    "run.log",
    "success",
)


def ssh(host: str, command: str, data: bytes | None = None, timeout: int = 180) -> str:
    proc = subprocess.run(
        [SSH, *OPTS, host, command],
        input=data,
        capture_output=True,
        timeout=timeout,
    )
    if proc.returncode:
        raise RuntimeError(
            f"SSH {proc.returncode}: {proc.stderr.decode(errors='replace')[-800:]}"
        )
    return proc.stdout.decode(errors="replace")


def lf_bytes(path: pathlib.Path) -> bytes:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    if b"\r" in data:
        raise ValueError(f"CR remains: {path}")
    return data


def job_id(text: str) -> str:
    match = re.search(r"Submitted batch job (\d+)", text)
    if not match:
        raise RuntimeError(f"no job id in: {text!r}")
    return match.group(1)


def put(src: pathlib.Path, remote_name: str, evidence: dict[str, str]) -> None:
    data = lf_bytes(src)
    digest = hashlib.sha256(data).hexdigest()
    remote = f"{REMOTE_JOBS}/{remote_name}"
    ssh(HOST_XFER, f"cat > '{remote}'", data)
    observed = ssh(HOST_XFER, f"sha256sum '{remote}'").split()[0]
    if observed != digest:
        raise RuntimeError(f"hash mismatch {remote}")
    evidence[remote_name] = digest


def stage() -> None:
    ssh(HOST_LOGIN, f"mkdir -p '{REMOTE_JOBS}' '{REMOTE_OUT}'")
    evidence: dict[str, str] = {}
    for rel in UPLOAD:
        put(DEST / rel, pathlib.Path(rel).name, evidence)
    for rel in REFS:
        put(DEST / rel, pathlib.Path(rel).name, evidence)
    receipt = DEST / "logs" / "mn5_stage.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"staged": list(evidence)}, indent=2))


def submit_more() -> None:
    ledger = DEST / "logs" / "mn5_submit.tsv"
    rows = ledger.read_text(encoding="utf-8") if ledger.is_file() else "job\tslurm_id\tutc\n"
    for name in ("rigor_clinvar.sbatch", "rigor_roh_summarize.sbatch"):
        text = ssh(HOST_LOGIN, f"sbatch --chdir='{REMOTE_JOBS}' '{REMOTE_JOBS}/{name}'")
        ident = job_id(text)
        print(name, ident, text.strip())
        rows += f"{name}\t{ident}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n"
    ledger.write_text(rows, encoding="utf-8")


def submit() -> None:
    rows = ["job\tslurm_id\tutc\n"]
    scan_text = ssh(
        HOST_LOGIN,
        f"sbatch --chdir='{REMOTE_JOBS}' '{REMOTE_JOBS}/rigor_scan.sbatch'",
    )
    scan_id = job_id(scan_text)
    print("rigor_scan.sbatch", scan_id, scan_text.strip())
    rows.append(
        f"rigor_scan.sbatch\t{scan_id}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n"
    )
    independent = ["rigor_roh.sbatch", "rigor_splice.sbatch"]
    for name in independent:
        text = ssh(HOST_LOGIN, f"sbatch --chdir='{REMOTE_JOBS}' '{REMOTE_JOBS}/{name}'")
        ident = job_id(text)
        print(name, ident, text.strip())
        rows.append(f"{name}\t{ident}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n")
    pid_text = ssh(
        HOST_LOGIN,
        f"sbatch --dependency=afterok:{scan_id} --chdir='{REMOTE_JOBS}' "
        f"'{REMOTE_JOBS}/rigor_pid.sbatch'",
    )
    pid_id = job_id(pid_text)
    print("rigor_pid.sbatch", pid_id, pid_text.strip())
    rows.append(
        f"rigor_pid.sbatch\t{pid_id}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n"
    )
    (DEST / "logs" / "mn5_submit.tsv").write_text("".join(rows), encoding="utf-8")


def status() -> None:
    print(ssh(HOST_LOGIN, "squeue -u \"$USER\" -o '%.18i %.12j %.8T %.10M %.6D %R'"))
    print(ssh(HOST_LOGIN, f"ls -ld {REMOTE_OUT}/* 2>/dev/null | tail"))


def collect() -> None:
    names = ",".join(repr(n) for n in COLLECT)
    script = (
        "python3 - <<'PY'\n"
        "import pathlib\n"
        f"root=pathlib.Path('{REMOTE_OUT}')\n"
        f"names=({names})\n"
        "for p in root.rglob('*'):\n"
        "    if p.is_file() and p.name in names and p.stat().st_size < 50_000_000:\n"
        "        if p.name == 'keys_impact.tsv':\n"
        "            continue\n"
        "        print(p)\n"
        "PY"
    )
    paths = [line.strip() for line in ssh(HOST_LOGIN, script).splitlines() if line.strip()]
    dest_root = DEST / "results"
    dest_root.mkdir(parents=True, exist_ok=True)
    batch = ["lcd " + str(dest_root).replace("\\", "/")]
    for remote in paths:
        parent = pathlib.Path(remote).parent.name
        batch.append(f"get {remote} {parent}__{pathlib.Path(remote).name}")
    proc = subprocess.run(
        [SFTP, *OPTS, "-b", "-", HOST_XFER],
        input="\n".join(batch) + "\n",
        capture_output=True,
        text=True,
        timeout=300,
    )
    (DEST / "logs" / "mn5_collect.log").write_text(
        proc.stdout + "\n" + proc.stderr, encoding="utf-8"
    )
    if proc.returncode:
        raise RuntimeError(proc.stderr[-800:])
    print(f"collected={len(paths)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode", choices=("stage", "submit", "submit_more", "status", "collect")
    )
    args = parser.parse_args()
    {
        "stage": stage,
        "submit": submit,
        "submit_more": submit_more,
        "status": status,
        "collect": collect,
    }[args.mode]()


if __name__ == "__main__":
    main()
