"""Arrange collected MN5 tables, write provenance and MANIFEST.sha256."""

from __future__ import annotations

import os

import hashlib
import pathlib
import shutil

DEST = pathlib.Path(__file__).resolve().parents[1]
RESULTS = DEST / "results"

MAP = {
    "vaf_counts.tsv": "vaf",
    "vaf_by_gene_six.tsv": "vaf",
    "vaf_jury_alleles.tsv": "vaf",
    "clinvar_plp_summary.tsv": "clinvar",
    "clinvar_plp_genes.tsv": "clinvar",
    "clinvar_plp_jury.tsv": "clinvar",
    "rare_counts.tsv": "rare",
    "rare_sixgene.tsv": "rare",
    "rare_acmg_sf.tsv": "rare",
    "pid_pair_counts.tsv": "pid",
    "pid_sixgene.tsv": "pid",
    "pid_jury.tsv": "pid",
    "roh_summary.tsv": "roh",
    "splice_sixgene_counts.tsv": "splice",
    "splice_jury.tsv": "splice",
    "scan_counts.json": "vaf",
}

# Scan wrote empty CLNSIG tables; prefer the dedicated ClinVar join and long-window ROH.
PREFIX_RANK = {
    "clinvar_": 30,
    "roh_sum_": 20,
    "pid_": 10,
    "splice_": 10,
    "scan_": 5,
    "roh_": 1,
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    placed: list[pathlib.Path] = []
    chosen: dict[str, tuple[int, pathlib.Path]] = {}
    for path in RESULTS.glob("*__*"):
        if not path.is_file():
            continue
        prefix, name = path.name.split("__", 1)
        folder = MAP.get(name)
        if folder is None:
            continue
        rank = 0
        for key, value in PREFIX_RANK.items():
            if prefix.startswith(key.rstrip("_")) or prefix.startswith(key):
                rank = max(rank, value)
        prev = chosen.get(name)
        if prev is None or rank >= prev[0]:
            chosen[name] = (rank, path)
    for name, (_rank, path) in sorted(chosen.items()):
        dest = RESULTS / MAP[name] / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
        placed.append(dest)

    logs = DEST / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    empty = RESULTS / "scan_47437284__clinvar_plp_summary.tsv"
    if empty.is_file():
        shutil.copy2(empty, logs / "original_annotate_clnsig_empty.tsv")
    absent = RESULTS / "clinvar" / "clinvar_plp_jury_absent.tsv"
    absent.parent.mkdir(parents=True, exist_ok=True)
    absent.write_text(
        "label\tchrom\tpos\tref\talt\tgene\tclinvar_class\tnote\n"
        "BUB1B_p.Asn1002Lys\tchr15\t40220612\tT\tG\tBUB1B\tnot_in_clinvar_exact_allele\t"
        "no CHROM+POS+REF+ALT hit after chr-prefixing ClinVar; not proof of absence from every ClinVar file; not a causal claim\n",
        encoding="utf-8",
    )
    placed.append(absent)

    for pattern in ("*__run.log", "*__exit_code.txt", "*__success", "*__python.version", "*__bcftools.version"):
        for path in RESULTS.glob(pattern):
            dest = logs / path.name
            shutil.copy2(path, dest)

    provenance = DEST / "provenance.tsv"
    rows = ["item\tversion\tsource\tsha256\tnotes"]
    source = DEST / "refs" / "SOURCE.tsv"
    if source.is_file():
        for line in source.read_text(encoding="utf-8").splitlines()[1:]:
            if line.strip():
                parts = line.split("\t")
                rows.append(
                    f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[-1]}"
                )
    submit = DEST / "logs" / "mn5_submit.tsv"
    if submit.is_file():
        rows.append(
            f"mn5_submit.tsv\tacc_ehpc\t{submit.as_posix()}\t{sha256(submit)}\t"
            "scan 47437284; roh 47437286; splice 47437289; pid 47437290; clinvar 47437514; roh_sum 47437515"
        )
    stage = DEST / "logs" / "mn5_stage.json"
    if stage.is_file():
        rows.append(
            f"mn5_stage.json\tstaged scripts/refs\t{stage.as_posix()}\t{sha256(stage)}\t"
            "SHA-256 verified on <mn5-transfer>"
        )
    rows.append(
        "annotated_vcf\tcritical_annotation_47397018\t"
        + str(pathlib.Path(os.environ["MVA_BASE"]) / "out/critical_annotation_47397018/with_clinvar.vcf.gz")
        + "\t"
        "6f7c4255f931ff4e6d1f403ebe696cda73a4158c2890a49efc21b69de3e5bd6d\t"
        "not copied to the laptop; this table does not establish causality"
    )
    for path in placed:
        rows.append(
            f"{path.relative_to(DEST).as_posix()}\tMN5 collect\t{path.as_posix()}\t"
            f"{sha256(path)}\tthis table does not establish causality"
        )
    provenance.write_text("\n".join(rows) + "\n", encoding="utf-8")

    manifest = DEST / "MANIFEST.sha256"
    hashed = []
    for path in sorted(DEST.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(DEST).as_posix()
        if rel == "MANIFEST.sha256":
            continue
        if path.stat().st_size >= 50_000_000:
            raise SystemExit(f"oversized laptop file {rel}")
        hashed.append(f"{sha256(path)}  {rel}")
    manifest.write_text("\n".join(hashed) + "\n", encoding="utf-8")
    print(f"placed={len(placed)} manifest={len(hashed)}")


if __name__ == "__main__":
    main()
