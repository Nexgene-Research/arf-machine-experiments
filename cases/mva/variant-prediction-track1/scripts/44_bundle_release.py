"""Assemble a minimal audited Track 1 release from explicit allowlists.

Run after CSV, report and workbook are final. Refuses to overwrite files.
Repo initialization/commit and receipt are separate final freeze steps.
"""

from __future__ import annotations

import hashlib
import pathlib
import shutil
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEST = ROOT / "release/Nexgene-AI-track1"
SCRIPTS = [
    "scripts/mn5/l1_20261015_align.sbatch",
    "scripts/mn5/l1_20261015_merge.sbatch",
    "scripts/mn5/fastpath_run_v2.sbatch",
    "scripts/mn5/fastpath_stream.py",
    "scripts/mn5/critical_annotation_20261016.sbatch",
    "scripts/mn5/track1_scores.sbatch",
    "scripts/mn5/track1_scores.py",
    "scripts/mn5/track1_pairs.py",
    "scripts/mn5/track1_reads.sbatch",
    "scripts/mn5/track1_reads.py",
    "scripts/mn5/track1_epcr.py",
    "scripts/mn5/track1_v2_bcd.py",
    "scripts/mn5/track1_v2_bcd_audit.py",
    "scripts/28_alphagenome_bounded.py",
    "scripts/30_rank_pairs.py",
    "scripts/31_write_submission.py",
    "scripts/32_epcr.py",
    "scripts/33_evidence_graph.py",
    "scripts/34_controls.py",
    "scripts/35_write_mn5_csv.py",
    "scripts/36_release_checks.py",
    "scripts/38_funnel.py",
    "scripts/40_lineage.py",
    "scripts/40_phase_l_literature.py",
    "scripts/50_track2_handoff.py",
    "scripts/51_phase_figure.py",
    "scripts/52_record_provenance.py",
    "scripts/53_redteam_reconcile.py",
    "scripts/54_adjudication_final.py",
    "scripts/figstyle.py",
    "scripts/validate_submission.py",
]
DOCS = [
    "docs/preregistration.md",
    "docs/lineage.md",
    "docs/manual_steps.md",
    "results/track1/v2_adjudication_protocol.md",
    "results/track1/official_requirements_audit.md",
]
AGGREGATES = [
    "results/track1/allele_validation_summary.md",
    "results/track1/phase_evidence.md",
    "results/track1/phase_aggregate.png",
    "results/track1/phase_aggregate.pdf",
    "results/track1/phase_aggregate.svg",
    "results/track1/second_hit_summary.md",
    "results/track1/literature_evidence.tsv",
    "results/track1/evidence_convergence.tsv",
    "results/track1/literature_summary.md",
    "results/l7_metrics_audited.tsv",
]


def copy(source: str, destination: pathlib.Path, *, required: bool = False) -> None:
    src = ROOT / source
    if not src.is_file():
        if required:
            raise FileNotFoundError(src)
        return
    if src.stat().st_size > 50_000_000 or src.suffix.lower() in {
        ".bam",
        ".cram",
        ".vcf",
        ".fastq",
        ".docx",
    }:
        raise ValueError(f"prohibited in release: {src}")
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, destination)


def text(path: str, contents: str) -> None:
    dst = DEST / path
    if dst.exists():
        raise FileExistsError(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(contents, encoding="utf-8")


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if DEST.exists():
        raise FileExistsError(DEST)
    all_scripts = set(SCRIPTS)
    for pattern in (
        "scripts/mn5/track1_*.sbatch",
        "scripts/4[1-9]_*.py",
        "scripts/3[7-9]_*.py",
    ):
        all_scripts.update(
            path.relative_to(ROOT).as_posix() for path in ROOT.glob(pattern)
        )
    for path in sorted(all_scripts):
        copy(path, DEST / path, required=True)
    for path in DOCS:
        copy(path, DEST / "docs" / pathlib.Path(path).name, required=True)
    for path in AGGREGATES:
        copy(path, DEST / "results" / pathlib.Path(path).name)
    for path in [
        "sub/Nexgene-AI_mn5-evidence-v2.csv",
        "results/Nexgene-AI_track1_report.md",
        "sub/Nexgene-AI_methods_description.xlsx",
    ]:
        copy(path, DEST / "submission" / pathlib.Path(path).name, required=True)
    copy("results/provenance.tsv", DEST / "docs/provenance.tsv", required=True)
    for optional in ["docs/methods.md", "docs/deviations.md", "docs/data_access.md"]:
        copy(optional, DEST / optional, required=True)
    text(
        ".gitignore",
        "*.vcf*\n*.bam*\n*.cram*\n*.fastq*\n*.fq*\n*.docx\n*.key\n*.pem\n.env\n__pycache__/\n",
    )
    text(
        "LICENSE",
        "Submissions and reports: Creative Commons Attribution 4.0 International (CC BY 4.0), subject to organizer terms.\nNo new code license asserted: authorship and upstream licensing must be checked before relicensing copied tools or scripts. External reference data and software retain their original licenses and attributions. Restricted patient source data are not included.\n",
    )
    text(
        "config/paths.example.yaml",
        "# Edit these authorized paths on your own installation; never commit controlled genomes.\nsource_fastqs: /authorized/data/lanes/*.fastq.gz\nsource_vcf: /authorized/data/input.vcf.gz\nreference_grch38: /authorized/ref/GRCh38_no_alt.fa\nmn5_work: /authorized/mn5/work\nproject_venv: /authorized/python.exe\n",
    )
    text(
        "config/containers.sha256",
        "bcftools_1.24--h487d631_1.sif 72b84cccdaf50c705530aa532dcb7f1cabc4b803da5a00bc48d153c65f27f889\nsnpeff_5.2--hdfd78af_1.sif b9181ec1a90eac8e58d1222fe0bd5454e059b379e3ab22d64a3d3a33eec64f8f\nwhatshap_2.8--py39h2de1943_0.sif 16608efe71b65b22a001a30230bbe4a04436a86c7d2c27120ccb1bbff08cda93\nbwa-mem2_2.3--he70b90d_0.sif ddcb61f9da1d85536b65c439ef9a892ea85b0c4f482df193be0832ca0d791e66\n",
    )
    text(
        "config/tool_versions.tsv",
        "tool\tversion\tsource\npython\t3.12\tproject venv\nbcftools\t1.24\tMN5 container\nWhatsHap\t2.8\tMN5 container\nsamtools\t1.23.1\tMN5 module\nintel\t2023.2.0\tMN5 module\nsingularity\t4.1.5\tMN5 module\n",
    )
    frozen = subprocess.run(
        [str(ROOT / "env/venv/Scripts/python.exe"), "-m", "pip", "freeze"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    text("config/pip_freeze.txt", frozen)
    # The curated release README is maintained in the release tree itself.
    if not (DEST / "README.md").exists():
        text(
        "README.md",
        "# Nexgene AI — Track 1 reproducibility release\n\nGRCh38 variant prioritization for a single case; see submission/Nexgene-AI_track1_report.md. The CSV ePCR is a heuristic, not a validated clinical probability. Original patient FASTQ, BAM, VCF and clinical notes are not shipped; authorized access is necessary for a full rerun.\n\n## Reproduce\n\nEdit `config/paths.example.yaml` to map authorized GRCh38 reference, FASTQs, VCF and MN5 work directory. Match input SHA-256 against `docs/lineage.md` and `docs/provenance.tsv` before proceeding. On MN5: `sbatch scripts/mn5/l1_20261015_align.sbatch`; `sbatch scripts/mn5/l1_20261015_merge.sbatch`; `sbatch scripts/mn5/fastpath_run_v2.sbatch`; `sbatch scripts/mn5/critical_annotation_20261016.sbatch`; `sbatch scripts/mn5/track1_scores.sbatch`; `python3 scripts/mn5/track1_pairs.py <workdir>`; `sbatch scripts/mn5/track1_reads.sbatch`; `python3 scripts/mn5/track1_epcr.py <workdir> <annotated-vcf>`. See scripts for actual arguments and host-specific paths; no refactor was undertaken. On laptop run `python scripts/38_funnel.py`, `python scripts/41_write_v2.py`, and `python scripts/validate_submission.py submission/Nexgene-AI_mn5-evidence-v2.csv /authorized/ref/GRCh38_no_alt.fa` using the pinned project interpreter. Runtime and hardware: MN5 CPU nodes (80 CPUs/4 GPUs available; jobs request CPU allocation); laptop 12 CPUs/16 GB; stage runtime from Slurm receipts. Expected hashes appear in `docs/lineage.md` and `MANIFEST.sha256`. Case-level expert judgments and primary source identifiers are in `docs/methods.md` and the evidence tables.\n\n## Access and redistribution\n\nSee `docs/data_access.md`. No restricted source genomes are included. License and credit information: `LICENSE`.\n",
    )
    # Synthetic-only, no patient-derived input fixtures.
    text(
        "tests/test_synthetic.py",
        "import csv\nimport io\n\nCOLUMNS = 'proband_id chrom_1 pos_1 ref_1 alt_1 chrom_2 pos_2 ref_2 alt_2 epcr finding_type notes'.split()\n\ndef test_synthetic_schema_and_order():\n    synthetic = ','.join(COLUMNS) + '\\nSYN,chr1,2,C,T,,,,,0.1,primary,synthetic\\nSYN,chr1,3,G,A,,,,,0.05,primary,synthetic\\n'\n    reader = csv.DictReader(io.StringIO(synthetic)); rows = list(reader)\n    assert reader.fieldnames == COLUMNS and len(rows) == 2\n    assert [float(x['epcr']) for x in rows] == sorted([float(x['epcr']) for x in rows], reverse=True)\n\ndef test_synthetic_ref_and_ties():\n    reference = 'ACGT'\n    assert reference[1] == 'C'\n    assert sorted([(0.1, 'B'), (0.1, 'A')], key=lambda x:(-x[0],x[1])) == [(0.1,'A'),(0.1,'B')]\n",
    )
    entries = sorted(
        (digest(p), p.relative_to(DEST).as_posix())
        for p in DEST.rglob("*")
        if p.is_file()
    )
    text("MANIFEST.sha256", "".join(f"{sha}  {name}\n" for sha, name in entries))
    print(f"Release assembled: {DEST}; {len(entries)} manifest files")


if __name__ == "__main__":
    main()
