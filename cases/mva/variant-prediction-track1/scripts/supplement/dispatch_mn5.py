"""Stage, submit, poll and collect MN5 supplement jobs.

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
import time

DEST = pathlib.Path(__file__).resolve().parents[1]
PROJECT = DEST.parent.parent
SSH = r"C:\Program Files\Git\usr\bin\ssh.exe"
SFTP = r"C:\Program Files\Git\usr\bin\sftp.exe"
OPTS = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=30"]
REMOTE_JOBS = str(pathlib.Path(os.environ["MVA_BASE"]) / "jobs/track1_supplement")
REMOTE_OUT = str(pathlib.Path(os.environ["MVA_BASE"]) / "out/track1_supplement")
HOST_LOGIN = "<user>@<mn5-login>"
HOST_XFER = "<user>@<mn5-transfer>"

UPLOAD = [
    "scripts/mn5/common.sh",
    "scripts/mn5/l1_consensus.sbatch",
    "scripts/mn5/l1_extract.sbatch",
    "scripts/mn5/l1_extract_windows.py",
    "scripts/mn5/l3_manta.sbatch",
    "scripts/mn5/l3_delly.sbatch",
    "scripts/mn5/l3_depth.sbatch",
    "scripts/mn5/filter_sv_windows.py",
    "scripts/mn5/exon_depth_table.py",
    "scripts/mn5/het_baf.py",
    "scripts/mn5/het_baf.sbatch",
    "scripts/mn5/l6_missense.py",
    "scripts/mn5/l6_missense.sbatch",
    "scripts/mn5/calibrate_predictors.py",
    "scripts/mn5/calibrate.sbatch",
    "scripts/mn5/alphagenome_status.py",
    "scripts/mn5/alphagenome_infer.py",
    "scripts/mn5/alphagenome_pilot.sbatch",
    "scripts/mn5/hg002_null.py",
    "scripts/mn5/hg002_null.sbatch",
    "scripts/mn5/bamsurgeon.sbatch",
    "scripts/mn5/l2_locus.sbatch",
]


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


def stage() -> None:
    ssh(HOST_LOGIN, f"mkdir -p '{REMOTE_JOBS}' '{REMOTE_OUT}'")
    extra = [
        (DEST / "scripts" / "mn5" / "l1_consensus_merge.py", "l1_consensus_merge.py"),
        (DEST / "computed" / "l3_sv" / "tier1_sixgene_exon_utr.bed", "tier1_sixgene_exon_utr.bed"),
        (DEST / "compiled" / "alleles" / "alphagenome_queue.tsv", "alphagenome_queue.tsv"),
    ]
    files = [(DEST / rel, pathlib.Path(rel).name) for rel in UPLOAD]
    files.extend(extra)
    evidence = {}
    for src, name in files:
        data = lf_bytes(src)
        digest = hashlib.sha256(data).hexdigest()
        remote = f"{REMOTE_JOBS}/{name}"
        ssh(HOST_XFER, f"cat > '{remote}'", data)
        observed = ssh(HOST_XFER, f"sha256sum '{remote}'").split()[0]
        if observed != digest:
            raise RuntimeError(f"hash mismatch {remote}")
        evidence[name] = digest
    receipt = DEST / "logs" / "mn5_stage.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"staged": list(evidence)}, indent=2))


def submit() -> None:
    cpu = [
        "l1_consensus.sbatch",
        "l3_manta.sbatch",
        "l3_delly.sbatch",
        "l3_depth.sbatch",
        "het_baf.sbatch",
        "calibrate.sbatch",
        "hg002_null.sbatch",
        "bamsurgeon.sbatch",
        "l2_locus.sbatch",
    ]
    gpu = [
        "l6_missense.sbatch",
        "alphagenome_pilot.sbatch",
    ]
    rows = ["job\tslurm_id\tutc\n"]
    ledger = DEST / "logs" / "mn5_submit.tsv"
    l1_id = None
    for name in cpu + gpu:
        extra = ""
        if name == "l1_extract.sbatch" and l1_id:
            extra = f" --dependency=afterok:{l1_id}"
        text = ssh(
            HOST_LOGIN,
            f"sbatch{extra} --chdir='{REMOTE_JOBS}' '{REMOTE_JOBS}/{name}'",
        )
        ident = job_id(text)
        print(name, ident, text.strip())
        rows.append(f"{name}\t{ident}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n")
        ledger.write_text("".join(rows), encoding="utf-8")
        if name == "l1_consensus.sbatch":
            l1_id = ident
        time.sleep(1)
    if l1_id:
        text = ssh(
            HOST_LOGIN,
            f"sbatch --dependency=afterok:{l1_id} --chdir='{REMOTE_JOBS}' "
            f"'{REMOTE_JOBS}/l1_extract.sbatch'",
        )
        ident = job_id(text)
        print("l1_extract.sbatch", ident, text.strip())
        rows.append(
            f"l1_extract.sbatch\t{ident}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n"
        )
    (DEST / "logs" / "mn5_submit.tsv").write_text("".join(rows), encoding="utf-8")


def status() -> None:
    print(ssh(HOST_LOGIN, "squeue -u \"$USER\" -o '%.18i %.12j %.8T %.10M %.6D %R'"))


def collect() -> None:
    script = (
        "python3 - <<'PY'\n"
        "import pathlib\n"
        f"root=pathlib.Path('{REMOTE_OUT}')\n"
        "names=('caller_discordance.tsv','het_baf.tsv','exon_depth.tsv',"
        "'manta_tier1_sv.tsv','delly_tier1_sv.tsv','asn1002lys.tsv',"
        "'indel_inference.tsv','predictor_lr.tsv','positive_controls.tsv',"
        "'null_hg002.tsv','l2_locus.tsv','exit_code.txt')\n"
        "for p in root.rglob('*'):\n"
        "    if p.is_file() and p.name in names and p.stat().st_size < 50_000_000:\n"
        "        print(p)\n"
        "PY"
    )
    paths = [line.strip() for line in ssh(HOST_LOGIN, script).splitlines() if line.strip()]
    dest_root = DEST / "computed" / "mn5_collect"
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
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("stage", "submit", "status", "collect"))
    args = ap.parse_args()
    {"stage": stage, "submit": submit, "status": status, "collect": collect}[args.mode]()


if __name__ == "__main__":
    main()
