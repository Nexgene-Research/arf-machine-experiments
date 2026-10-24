"""Build auditable Track 1 lineage from recorded hashes and current scripts.

Never hashes patient genomes on the laptop; the MN5 job receipts supply hashes.
Unknown hashes are explicitly reported as not independently recorded.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/lineage.md"
STEPS = [
    (
        "L1 alignment",
        "47349918_1-4",
        "scripts/mn5/l1_20261015_align.sbatch",
        "FASTQ lane hashes: results/l1_20261015/provenance_20261015T211105Z.tsv",
        "BAM hash: remote L1 receipt",
        "bwa-mem2: remote receipt",
    ),
    (
        "L1 merge / markdup",
        "47353088",
        "scripts/mn5/l1_20261015_merge.sbatch",
        "aligned lanes: remote L1 receipt",
        "marked.bam: remote L1 receipt",
        "Picard: remote receipt",
    ),
    (
        "gnomAD fastpath",
        "47392756",
        "scripts/mn5/fastpath_run_v2.sbatch",
        "source VCF: results/fastpath_v2_full_collect.json",
        "4aa21a3477584295a08fdfb22a577b60612d1d7b6cfc88a367d737a89b4eee5c",
        "72b84cccdaf50c705530aa532dcb7f1cabc4b803da5a00bc48d153c65f27f889",
    ),
    (
        "critical annotation",
        "47397018",
        "scripts/mn5/critical_annotation_20261016.sbatch",
        "4aa21a3477584295a08fdfb22a577b60612d1d7b6cfc88a367d737a89b4eee5c",
        "6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d",
        "72b84cccdaf50c705530aa532dcb7f1cabc4b803da5a00bc48d153c65f27f889",
    ),
    (
        "candidate scores",
        "47397018 downstream",
        "scripts/mn5/track1_scores.py",
        "6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d",
        "7d0f2efd998fa92e44a98d998a5a480e460417e8c43712994684f2359519cd61",
        "MN5 provenance not independently recorded",
    ),
    (
        "candidate pairs",
        "47397018 downstream",
        "scripts/mn5/track1_pairs.py",
        "0d115c7145d96c9318dc7cf3d0cdec077e5d339c09c969bd95079bf8f4516e1c",
        "372c4dffe80861c92fa1a149572c6212661dd29c852c6e68e7ab6b7814593841",
        "MN5 provenance not independently recorded",
    ),
    (
        "read screen",
        "47397018 downstream",
        "scripts/mn5/track1_reads.py",
        "372c4dffe80861c92fa1a149572c6212661dd29c852c6e68e7ab6b7814593841",
        "read-stage output hash not independently recorded",
        "WhatsHap image: results/inputs.sha256 if used",
    ),
    (
        "ePCR ranking",
        "47397018 downstream",
        "scripts/mn5/track1_epcr.py",
        "read-stage output hash not independently recorded",
        "0c19ae9b2256e18186d651efd2e957f1f010db6dbba1d6ff6719c2b831fdff6a",
        "system python3",
    ),
    (
        "legacy CSV",
        "local",
        "scripts/35_write_mn5_csv.py",
        "abcd7194e1ac61d25bd9def66ab96d93df45ae824ebd71f8399f19eaa480f981",
        "41bd90a7fbc9738ad670943a41bc5ee7437d6902c5de5b3dde1ecc844b9bcecb",
        "project Python venv",
    ),
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    output = Path(sys.argv[1]) if len(sys.argv) == 2 else OUT
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for label, job, script, before, after, container in STEPS:
        path = ROOT / script
        rows.append(
            f"| {label} | {job} | `{script}` | `{sha(path) if path.is_file() else 'not archived locally'}` | {before} | {after} | {container} |"
        )
    receipt = json.loads(
        (ROOT / "results/track1/bcd_47403447_receipt.json").read_text()
    )
    outputs = receipt["outputs"]
    bcd_script = "scripts/mn5/track1_v2_bcd.py"
    bcd_sbatch = "scripts/mn5/track1_v2_bcd.sbatch"
    rows.append(
        f"| v2 bounded allele/phase/second hit | 47403447 | `{bcd_sbatch}`; `{bcd_script}` | `{sha(ROOT / bcd_sbatch)}`; `{sha(ROOT / bcd_script)}` | annotated VCF `6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d`; BAM hash not independently recorded | allele validation `{outputs['allele_validation.tsv']['sha256']}`; phase `{outputs['phase_pair_evidence.tsv']['sha256']}`; other outputs in `results/track1/bcd_47403447_receipt.json` | bcftools `72b84cccdaf50c705530aa532dcb7f1cabc4b803da5a00bc48d153c65f27f889`; WhatsHap image digest not independently recorded |"
    )
    rows.append(
        f"| v2 eligibility and CSV | local | `scripts/41_write_v2.py` | `{sha(ROOT / 'scripts/41_write_v2.py')}` | v1 top10 `abcd7194e1ac61d25bd9def66ab96d93df45ae824ebd71f8399f19eaa480f981`; allele `{outputs['allele_validation.tsv']['sha256']}`; phase `{outputs['phase_pair_evidence.tsv']['sha256']}` | `{sha(ROOT / 'sub/Nexgene-AI_mn5-evidence-v2.csv')}` | project Python venv |"
    )
    output.write_text(
        "# Track 1 lineage (audited 2026-10-16)\n\n"
        "Source job identifiers, FASTA/VCF hash receipts: `results/track1/funnel.json`, `results/inputs.sha256`, `results/outputs.sha256`. Textual 'not independently recorded' is an evidence gap, not a checksum. Script SHA-256 is measured from the project file by `scripts/40_lineage.py`. Read FASTQ/BAM/VCF stays on MN5. **There are two upstream branches, not one fabricated FASTQ-to-VCF chain:** supplied Sentieon variant calls feed the gnomAD fastpath and critical annotation; FASTQs independently feed the bwa-mem2/markdup BAM used for v2 allele recall and phase. Both converge only at v2 adjudication on the *same individual*. Thus the listed L1 jobs did not create the supplied source VCF.\n\n"
        "| Step | Slurm job ID | Script path | Script SHA-256 | Input SHA-256 / receipt | Output SHA-256 / receipt | Container SHA-256 / receipt |\n|---|---|---|---|---|---|---|\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {output}: {len(rows)} lineage rows")


if __name__ == "__main__":
    main()
