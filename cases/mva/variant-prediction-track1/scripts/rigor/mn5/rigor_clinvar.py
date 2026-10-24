#!/usr/bin/env python3
"""Intersect the proband VCF with ClinVar after prefixing ClinVar contigs.

The annotated with_clinvar.vcf.gz has empty CLNSIG because ClinVar uses
unprefixed names (15) while the proband uses chr15. This job joins explicitly.
Missing CLNSIG is not benign. This table does not establish causality.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

TIER1 = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1")
JURY = {
    ("chr15", "40209701", "T", "G"): "BUB1B_p.Leu737Ter",
    ("chr15", "40220612", "T", "G"): "BUB1B_p.Asn1002Lys",
}


def load_genes(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def prefix(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"MT", "M"}:
        return "chrM"
    return "chr" + chrom


def chrom_ok(chrom: str) -> bool:
    rest = chrom[3:] if chrom.startswith("chr") else chrom
    if rest in {"X", "Y"}:
        return True
    return rest.isdigit() and 1 <= int(rest) <= 22


def info_map(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
        else:
            out[item] = ""
    return out


def gt_class(gt: str) -> str:
    norm = gt.replace("|", "/")
    if norm in {"0/1", "1/0"}:
        return "het"
    if norm in {"1/1"}:
        return "hom"
    if norm in {"1", "1/."}:
        return "hemi"
    return "other"


def ann_gene(info: dict[str, str]) -> str:
    raw = info.get("ANN", "")
    first = raw.split(",", 1)[0] if raw else ""
    parts = first.split("|")
    return parts[3] if len(parts) > 3 else "."


def clinvar_gene(geneinfo: str) -> str:
    if not geneinfo or geneinfo == ".":
        return "."
    return geneinfo.split("|", 1)[0].split(":", 1)[0]


def plp_class(clnsig: str) -> str:
    if clnsig in ("", ".", None):
        return "none"
    cleaned = clnsig.replace(" ", "_").lower()
    for sep in ("|", "&", ",", "/"):
        cleaned = cleaned.replace(sep, " ")
    tokens = [tok.strip("_") for tok in cleaned.split() if tok]
    pathogenic = False
    benign = False
    conflicting = False
    for tok in tokens:
        if tok in {"pathogenic", "likely_pathogenic"}:
            pathogenic = True
        elif tok in {"benign", "likely_benign"}:
            benign = True
        elif "conflicting" in tok:
            conflicting = True
    if conflicting or (pathogenic and benign):
        return "conflicting"
    if pathogenic:
        return "P/LP"
    return "other"


def stars(revstat: str) -> int:
    s = (revstat or "").lower()
    if "practice_guideline" in s:
        return 4
    if "reviewed_by_expert_panel" in s:
        return 3
    if "multiple_submitters" in s and "no_conflicts" in s:
        return 2
    if "criteria_provided" in s:
        return 1
    return 0


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def load_clinvar(path: Path) -> tuple[dict[tuple[str, str, str, str], tuple[str, str, str]], Counter[str]]:
    table: dict[tuple[str, str, str, str], tuple[str, str, str]] = {}
    stats: Counter[str] = Counter()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            stats["clinvar_records"] += 1
            fields = line.rstrip("\n").split("\t")
            chrom = prefix(fields[0])
            if not chrom_ok(chrom):
                stats["clinvar_dropped_contig"] += 1
                continue
            pos, _vid, ref, alt = fields[1:5]
            if alt == "*" or "," in alt:
                stats["clinvar_multiallelic_or_symbolic"] += 1
                continue
            inf = info_map(fields[7])
            clnsig = inf.get("CLNSIG", ".")
            if clnsig in {"", "."}:
                continue
            table[(chrom, pos, ref, alt)] = (
                clnsig,
                inf.get("CLNREVSTAT", "."),
                clinvar_gene(inf.get("GENEINFO", ".")),
            )
            stats["clinvar_with_clnsig"] += 1
    return table, stats


def main() -> None:
    vcf = Path(sys.argv[1])
    clin = Path(sys.argv[2])
    acmg_path = Path(sys.argv[3])
    dest = Path(sys.argv[4])
    dest.mkdir(parents=True, exist_ok=True)
    acmg = load_genes(acmg_path)
    table, stats = load_clinvar(clin)
    plp_genes: dict[str, Counter[str]] = defaultdict(Counter)
    plp_star: dict[str, int] = defaultdict(int)
    jury_rows: dict[str, dict[str, str]] = {}
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            stats["proband_records"] += 1
            fields = line.rstrip("\n").split("\t")
            chrom = prefix(fields[0])
            if not chrom_ok(chrom):
                continue
            pos, _vid, ref, alt = fields[1:5]
            hit = table.get((chrom, pos, ref, alt))
            if hit is None:
                continue
            stats["overlap"] += 1
            clnsig, revstat, cv_gene = hit
            klass_plp = plp_class(clnsig)
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}
            gclass = gt_class(fmt.get("GT", "."))
            inf = info_map(fields[7])
            gene = ann_gene(inf)
            if gene in {".", ""}:
                gene = cv_gene
            key = (chrom, pos, ref, alt)
            if klass_plp == "P/LP":
                stats["clinvar_plp"] += 1
                plp_genes[gene][gclass] += 1
                plp_genes[gene]["n"] += 1
                plp_star[gene] = max(plp_star[gene], stars(revstat))
                if gclass == "hom":
                    stats["clinvar_plp_hom"] += 1
                if gene in TIER1:
                    stats["clinvar_plp_tier1"] += 1
                if gene in acmg:
                    stats["clinvar_plp_acmg_sf"] += 1
            elif klass_plp == "conflicting":
                stats["clinvar_conflicting"] += 1
            if key in JURY:
                jury_rows[JURY[key]] = {
                    "label": JURY[key],
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene,
                    "gt": fmt.get("GT", "."),
                    "filter": fields[6],
                    "clinvar": clnsig,
                    "clinvar_class": klass_plp,
                    "clinvar_stars": str(stars(revstat)),
                    "clinvar_gene": cv_gene,
                    "note": "joined on chr-prefixed CHROM+POS+REF+ALT; not a causal claim",
                }
    write_tsv(
        dest / "clinvar_plp_summary.tsv",
        ["metric", "count", "note"],
        [
            ["clinvar_records_with_clnsig", stats["clinvar_with_clnsig"], "ClinVar GRCh38 after chr prefix; alt contigs dropped"],
            ["proband_exact_overlap", stats["overlap"], "same chrom,pos,ref,alt"],
            ["clinvar_plp_alleles", stats["clinvar_plp"], "CLNSIG Pathogenic or Likely_pathogenic, not conflicting"],
            ["clinvar_plp_hom", stats["clinvar_plp_hom"], "P/LP with GT 1/1"],
            ["clinvar_plp_genes", len(plp_genes), "distinct gene symbols"],
            ["clinvar_plp_tier1_alleles", stats["clinvar_plp_tier1"], "P/LP inside six MVA genes"],
            ["clinvar_plp_acmg_sf_alleles", stats["clinvar_plp_acmg_sf"], "P/LP whose gene is on ACMG SF v3.2"],
            ["clinvar_conflicting", stats["clinvar_conflicting"], "not counted as P/LP"],
            [
                "original_with_clinvar_clnsig_nonempty",
                0,
                "bcftools annotate did not join: ClinVar CHROM=15 vs proband CHROM=chr15",
            ],
        ],
    )
    plp_rows = []
    for gene, bucket in sorted(plp_genes.items(), key=lambda kv: (-kv[1]["n"], kv[0])):
        plp_rows.append(
            [
                gene,
                bucket["n"],
                bucket["het"],
                bucket["hom"],
                bucket["hemi"],
                plp_star[gene],
                "ACMG_SF_v3.2" if gene in acmg else "not_SF",
                "tier1" if gene in TIER1 else "other",
                "gene-level only; this is not a secondary-finding report",
            ]
        )
    write_tsv(
        dest / "clinvar_plp_genes.tsv",
        ["gene", "n_plp", "het", "hom", "hemi", "max_stars", "acmg_sf", "panel", "note"],
        plp_rows,
    )
    jury_header = [
        "label",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "gt",
        "filter",
        "clinvar",
        "clinvar_class",
        "clinvar_stars",
        "clinvar_gene",
        "note",
    ]
    write_tsv(
        dest / "clinvar_plp_jury.tsv",
        jury_header,
        [[row[h] for h in jury_header] for _label, row in sorted(jury_rows.items())],
    )
    (dest / "inputs.sha256").write_text(
        f"{hashlib.sha256(acmg_path.read_bytes()).hexdigest()}  {acmg_path}\n"
        f"vcf_path  {vcf}\n"
        f"clinvar_path  {clin}\n",
        encoding="utf-8",
    )
    (dest / "caveat.txt").write_text(
        "ClinVar join after chr-prefixing. Empty original CLNSIG is a contig mismatch, not evidence of absence. Not causal.\n",
        encoding="utf-8",
    )
    print(
        f"clinvar_clnsig={stats['clinvar_with_clnsig']} overlap={stats['overlap']} "
        f"plp={stats['clinvar_plp']} genes={len(plp_genes)}"
    )


if __name__ == "__main__":
    main()
