"""Stage, submit, poll and collect MN5 wave-3 jobs.

Uses Git OpenSSH. Does not write into the three frozen Track 1 packs.
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
REMOTE_JOBS = str(pathlib.Path(os.environ["MVA_BASE"]) / "jobs/track1_wave3")
REMOTE_OUT = str(pathlib.Path(os.environ["MVA_BASE"]) / "out/track1_wave3")
HOST_LOGIN = "<user>@<mn5-login>"
HOST_XFER = "<user>@<mn5-transfer>"

UPLOAD = [
    "scripts/mn5/common.sh",
    "scripts/mn5/pileup_published.py",
    "scripts/mn5/sixgene_cnv.py",
    "scripts/mn5/mosaic_table.py",
    "scripts/mn5/spliceai_keys.py",
    "scripts/mn5/extract_rare_vcf.py",
    "scripts/mn5/exomiser_top30.py",
    "scripts/mn5/wave3_pileup.sbatch",
    "scripts/mn5/wave3_cnv.sbatch",
    "scripts/mn5/wave3_verifybamid.sbatch",
    "scripts/mn5/wave3_mosaic.sbatch",
    "scripts/mn5/wave3_spliceai.sbatch",
    "scripts/mn5/wave3_exomiser.sbatch",
]

REFS = [
    "refs/published_sites.vcf",
    "refs/published_alleles.tsv",
    "refs/hpo.txt",
    "refs/gene_body_2kb.bed",
    "refs/spliceai_keys.tsv",
    "refs/tier1_sixgene_exon_utr.bed",
    "refs/verifybamid/1000g.phase3.10k.b38.vcf.gz.dat.UD",
    "refs/verifybamid/1000g.phase3.10k.b38.vcf.gz.dat.mu",
    "refs/verifybamid/1000g.phase3.10k.b38.vcf.gz.dat.bed",
]

COLLECT = (
    "published_site_gt_new.tsv",
    "sixgene_cnv.tsv",
    "mosaic_segments.tsv",
    "clone_fraction.txt",
    "spliceai_keys.tsv",
    "exomiser_top30.tsv",
    "bub1b_rank.txt",
    "verify.selfSM",
    "caveat.txt",
    "exit_code.txt",
    "run.log",
    "inputs.sha256",
    "success",
    "pip.log",
    "spliceai.log",
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


def submit() -> None:
    rows = ["job\tslurm_id\tutc\n"]
    for name in (
        "wave3_pileup.sbatch",
        "wave3_cnv.sbatch",
        "wave3_verifybamid.sbatch",
        "wave3_mosaic.sbatch",
        "wave3_spliceai.sbatch",
        "wave3_exomiser.sbatch",
    ):
        text = ssh(
            HOST_LOGIN,
            f"sbatch --chdir='{REMOTE_JOBS}' '{REMOTE_JOBS}/{name}'",
        )
        ident = job_id(text)
        print(name, ident, text.strip())
        rows.append(
            f"{name}\t{ident}\t{dt.datetime.now(dt.timezone.utc).isoformat()}\n"
        )
    (DEST / "logs" / "mn5_submit.tsv").write_text("".join(rows), encoding="utf-8")


def status() -> None:
    print(ssh(HOST_LOGIN, "squeue -u \"$USER\" -o '%.18i %.12j %.8T %.10M %.6D %R'"))
    print(
        ssh(
            HOST_LOGIN,
            f"ls -ld {REMOTE_OUT}/* 2>/dev/null | tail",
        )
    )


def collect() -> None:
    names = ",".join(repr(n) for n in COLLECT)
    script = (
        "python3 - <<'PY'\n"
        "import pathlib\n"
        f"root=pathlib.Path('{REMOTE_OUT}')\n"
        f"names=({names})\n"
        "for p in root.rglob('*'):\n"
        "    if p.is_file() and p.name in names and p.stat().st_size < 50_000_000:\n"
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
    parser.add_argument("mode", choices=("stage", "submit", "status", "collect"))
    args = parser.parse_args()
    {
        "stage": stage,
        "submit": submit,
        "status": status,
        "collect": collect,
    }[args.mode]()


if __name__ == "__main__":
    main()
