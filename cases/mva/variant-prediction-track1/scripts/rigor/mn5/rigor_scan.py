#!/usr/bin/env python3
"""One pass over the annotated VCF: VAF window, ClinVar P/LP, rare counts.

Does not dump genome-wide coordinates. Jury-report alleles may be named.
Missing AF is unqueried, never rare. This table does not establish causality.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

TIER1 = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1")
JURY = {
    ("chr15", "40209701", "T", "G"): "BUB1B_p.Leu737Ter",
    ("chr15", "40220612", "T", "G"): "BUB1B_p.Asn1002Lys",
}
HET_LO, HET_HI, HOM_MIN = 0.25, 0.75, 0.85
RARE_AF = 0.01
HIGH_MOD = {"HIGH", "MODERATE"}


def load_genes(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


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


def num(value: str | None) -> float | None:
    if value in (None, "", ".", "nan", "NaN"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def vaf_of(ad: str) -> float | None:
    parts = ad.split(",")
    if len(parts) < 2:
        return None
    try:
        ref_n, alt_n = int(parts[0]), int(parts[1])
    except ValueError:
        return None
    total = ref_n + alt_n
    if total <= 0:
        return None
    return alt_n / total


def chrom_ok(chrom: str) -> bool:
    if chrom.startswith("chr"):
        rest = chrom[3:]
    else:
        rest = chrom
        chrom = "chr" + chrom
    if rest in {"X", "Y", "M", "MT"}:
        return rest != "MT"
    return rest.isdigit() and 1 <= int(rest) <= 22


def prefix(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else "chr" + chrom


def gt_class(gt: str) -> str:
    norm = gt.replace("|", "/")
    if norm in {"0/1", "1/0"}:
        return "het"
    if norm in {"1/1"}:
        return "hom"
    if norm in {"1", "1/."}:
        return "hemi"
    return "other"


def ann_best(info: dict[str, str]) -> tuple[str, str, str, str]:
    raw = info.get("ANN", "")
    first = raw.split(",", 1)[0] if raw else ""
    parts = first.split("|")
    gene = parts[3] if len(parts) > 3 else "."
    effect = parts[1] if len(parts) > 1 else "."
    impact = parts[2] if len(parts) > 2 else "."
    hgvs = parts[9] if len(parts) > 9 else "."
    return gene, effect, impact, hgvs


def plp_class(clnsig: str) -> str:
    if clnsig in ("", ".", None):
        return "none"
    tokens = [tok.strip("_") for tok in clnsig.replace(" ", "_").lower().replace("/", ",").split(",") if tok]
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


def main() -> None:
    vcf = Path(sys.argv[1])
    acmg_path = Path(sys.argv[2])
    dest = Path(sys.argv[3])
    dest.mkdir(parents=True, exist_ok=True)
    acmg = load_genes(acmg_path)
    opener = gzip.open if str(vcf).endswith(".gz") else open
    counts: Counter[str] = Counter()
    vaf_six: dict[str, Counter[str]] = {g: Counter() for g in TIER1}
    rare_six: dict[str, Counter[str]] = {g: Counter() for g in TIER1}
    rare_acmg: dict[str, Counter[str]] = defaultdict(Counter)
    plp_genes: dict[str, Counter[str]] = defaultdict(Counter)
    plp_star: dict[str, int] = defaultdict(int)
    rare_het_genes: Counter[str] = Counter()
    rare_hom_genes: Counter[str] = Counter()
    unq_het_genes: Counter[str] = Counter()
    jury_rows: dict[str, dict[str, str]] = {}
    keys_path = dest / "keys_impact.tsv"
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle, keys_path.open(
        "w", encoding="utf-8", newline=""
    ) as keys_fh:
        keys_w = csv.writer(keys_fh, delimiter="\t", lineterminator="\n")
        keys_w.writerow(["chrom", "pos", "ref", "alt", "gene", "impact", "effect", "gt", "filter"])
        sample_index = 9
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, _qual, flt = fields[:7]
            counts["records"] += 1
            chrom = prefix(chrom)
            if not chrom_ok(chrom):
                counts["dropped_contig"] += 1
                continue
            if alt == "*" or set(ref + alt) - set("ACGT"):
                counts["symbolic"] += 1
                continue
            inf = info_map(fields[7])
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}
            gt = fmt.get("GT", ".")
            klass = gt_class(gt)
            ad = fmt.get("AD", ".")
            vaf = vaf_of(ad)
            gene, effect, impact, hgvs = ann_best(inf)
            af = num(inf.get("gnomad_AF"))
            faf = num(inf.get("gnomad_faf95max"))
            observed = [x for x in (af, faf) if x is not None]
            maxaf = max(observed) if observed else None
            af_state = "unqueried" if maxaf is None else "observed"
            rare = maxaf is not None and maxaf <= RARE_AF
            coding = impact in HIGH_MOD
            passed = flt == "PASS"
            clnsig = inf.get("CLNSIG", ".")
            revstat = inf.get("CLNREVSTAT", ".")
            plp = plp_class(clnsig)
            key = (chrom, pos, ref, alt)
            if passed:
                counts["pass"] += 1
            if klass == "het":
                counts["het"] += 1
                if vaf is None:
                    counts["het_vaf_missing"] += 1
                elif vaf < HET_LO or vaf > HET_HI:
                    counts["het_vaf_out"] += 1
                    if passed:
                        counts["het_vaf_out_pass"] += 1
                else:
                    counts["het_vaf_ok"] += 1
            elif klass == "hom":
                counts["hom"] += 1
                if vaf is None:
                    counts["hom_vaf_missing"] += 1
                elif vaf < HOM_MIN:
                    counts["hom_vaf_low"] += 1
                    if passed:
                        counts["hom_vaf_low_pass"] += 1
                else:
                    counts["hom_vaf_ok"] += 1
            if gene in vaf_six and klass in {"het", "hom"}:
                bucket = vaf_six[gene]
                bucket[klass] += 1
                if klass == "het" and vaf is not None and (vaf < HET_LO or vaf > HET_HI):
                    bucket["het_vaf_out"] += 1
                if klass == "hom" and vaf is not None and vaf < HOM_MIN:
                    bucket["hom_vaf_low"] += 1
            if plp == "P/LP":
                counts["clinvar_plp"] += 1
                plp_genes[gene][klass] += 1
                plp_genes[gene]["n"] += 1
                plp_star[gene] = max(plp_star[gene], stars(revstat))
                if klass == "hom":
                    counts["clinvar_plp_hom"] += 1
                if gene in TIER1:
                    counts["clinvar_plp_tier1"] += 1
            elif plp == "conflicting":
                counts["clinvar_conflicting"] += 1
            if coding and klass == "het":
                if rare:
                    rare_het_genes[gene] += 1
                    counts["rare_coding_het"] += 1
                elif af_state == "unqueried":
                    unq_het_genes[gene] += 1
                    counts["unqueried_coding_het"] += 1
            if coding and klass == "hom":
                if rare:
                    rare_hom_genes[gene] += 1
                    counts["rare_coding_hom"] += 1
                elif af_state == "unqueried":
                    counts["unqueried_coding_hom"] += 1
            if gene in rare_six and coding:
                if rare and klass == "het":
                    rare_six[gene]["rare_het"] += 1
                if rare and klass == "hom":
                    rare_six[gene]["rare_hom"] += 1
                if af_state == "unqueried" and klass == "het":
                    rare_six[gene]["unqueried_het"] += 1
                if af_state == "unqueried" and klass == "hom":
                    rare_six[gene]["unqueried_hom"] += 1
            if gene in acmg and coding:
                if rare and klass == "het":
                    rare_acmg[gene]["rare_het"] += 1
                if rare and klass == "hom":
                    rare_acmg[gene]["rare_hom"] += 1
                if af_state == "unqueried" and klass == "het":
                    rare_acmg[gene]["unqueried_het"] += 1
            if coding and klass == "het":
                keys_w.writerow([chrom, pos, ref, alt, gene, impact, effect, gt, flt])
                counts["keys_impact_het"] += 1
            if key in JURY:
                jury_rows[JURY[key]] = {
                    "label": JURY[key],
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene,
                    "gt": gt,
                    "filter": flt,
                    "vaf": "" if vaf is None else f"{vaf:.6g}",
                    "vaf_window": (
                        "missing"
                        if vaf is None
                        else (
                            "het_ok"
                            if klass == "het" and HET_LO <= vaf <= HET_HI
                            else (
                                "het_out"
                                if klass == "het"
                                else ("hom_ok" if klass == "hom" and vaf >= HOM_MIN else "hom_low")
                            )
                        )
                    ),
                    "impact": impact,
                    "effect": effect,
                    "hgvs": hgvs,
                    "af_state": af_state,
                    "gnomad_af": inf.get("gnomad_AF", "."),
                    "clinvar": clnsig,
                    "clinvar_class": plp,
                    "clinvar_stars": str(stars(revstat)),
                    "note": "jury-report allele; this row is not a causal claim",
                }

    two_plus = sum(1 for _gene, n in rare_het_genes.items() if n >= 2)
    two_plus_unq = sum(1 for _gene, n in unq_het_genes.items() if n >= 2)
    hom_genes = sum(1 for n in rare_hom_genes.values() if n)
    counts["genes_two_plus_rare_coding_het"] = two_plus
    counts["genes_two_plus_unqueried_coding_het"] = two_plus_unq
    counts["genes_rare_coding_hom"] = hom_genes
    tier_two = [g for g in TIER1 if rare_six[g]["rare_het"] >= 2]
    counts["tier1_genes_two_plus_rare_het"] = len(tier_two)

    write_tsv(
        dest / "vaf_counts.tsv",
        ["metric", "count", "note"],
        [
            ["records", counts["records"], "primary contigs after dropping symbolic alleles"],
            ["het", counts["het"], "GT 0/1 or 1/0"],
            ["hom", counts["hom"], "GT 1/1"],
            ["het_vaf_ok", counts["het_vaf_ok"], "het VAF inside 0.25-0.75"],
            ["het_vaf_out", counts["het_vaf_out"], "het VAF <0.25 or >0.75"],
            ["het_vaf_out_pass", counts["het_vaf_out_pass"], "same, FILTER=PASS"],
            ["het_vaf_missing", counts["het_vaf_missing"], "het without usable AD"],
            ["hom_vaf_ok", counts["hom_vaf_ok"], "hom VAF >=0.85"],
            ["hom_vaf_low", counts["hom_vaf_low"], "hom VAF <0.85"],
            ["hom_vaf_low_pass", counts["hom_vaf_low_pass"], "same, FILTER=PASS"],
            ["hom_vaf_missing", counts["hom_vaf_missing"], "hom without usable AD"],
        ],
    )
    write_tsv(
        dest / "vaf_by_gene_six.tsv",
        ["gene", "het", "hom", "het_vaf_out", "hom_vaf_low", "note"],
        [
            [
                gene,
                vaf_six[gene]["het"],
                vaf_six[gene]["hom"],
                vaf_six[gene]["het_vaf_out"],
                vaf_six[gene]["hom_vaf_low"],
                "counts only; coordinates omitted except jury alleles",
            ]
            for gene in TIER1
        ],
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
        "vaf",
        "vaf_window",
        "impact",
        "effect",
        "hgvs",
        "af_state",
        "gnomad_af",
        "clinvar",
        "clinvar_class",
        "clinvar_stars",
        "note",
    ]
    write_tsv(
        dest / "vaf_jury_alleles.tsv",
        jury_header,
        [[row[h] for h in jury_header] for _label, row in sorted(jury_rows.items())],
    )
    write_tsv(
        dest / "clinvar_plp_summary.tsv",
        ["metric", "count", "note"],
        [
            ["clinvar_plp_alleles", counts["clinvar_plp"], "CLNSIG Pathogenic or Likely_pathogenic, not conflicting"],
            ["clinvar_plp_hom", counts["clinvar_plp_hom"], "P/LP with GT 1/1"],
            ["clinvar_plp_genes", len(plp_genes), "distinct ANN gene symbols"],
            ["clinvar_plp_tier1_alleles", counts["clinvar_plp_tier1"], "P/LP inside six MVA genes"],
            ["clinvar_conflicting", counts["clinvar_conflicting"], "not counted as P/LP"],
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
    write_tsv(
        dest / "clinvar_plp_jury.tsv",
        jury_header,
        [[row[h] for h in jury_header] for _label, row in sorted(jury_rows.items())],
    )
    write_tsv(
        dest / "rare_counts.tsv",
        ["metric", "count", "note"],
        [
            ["rare_af_le_0.01_coding_het", counts["rare_coding_het"], "HIGH/MODERATE; AF observed and <=0.01"],
            ["rare_af_le_0.01_coding_hom", counts["rare_coding_hom"], "HIGH/MODERATE; AF observed and <=0.01"],
            ["unqueried_coding_het", counts["unqueried_coding_het"], "HIGH/MODERATE het with no AF tag; not rare"],
            ["unqueried_coding_hom", counts["unqueried_coding_hom"], "HIGH/MODERATE hom with no AF tag; not rare"],
            ["genes_two_plus_rare_coding_het", two_plus, "count of genes, not a coordinate dump"],
            ["genes_two_plus_unqueried_coding_het", two_plus_unq, "do not treat as rare compound-het genes"],
            ["genes_rare_coding_hom", hom_genes, "genes with at least one rare HIGH/MODERATE hom"],
            ["tier1_genes_two_plus_rare_het", len(tier_two), ",".join(tier_two) or "none"],
        ],
    )
    write_tsv(
        dest / "rare_sixgene.tsv",
        ["gene", "rare_het", "rare_hom", "unqueried_het", "unqueried_hom", "two_plus_rare_het", "note"],
        [
            [
                gene,
                rare_six[gene]["rare_het"],
                rare_six[gene]["rare_hom"],
                rare_six[gene]["unqueried_het"],
                rare_six[gene]["unqueried_hom"],
                int(rare_six[gene]["rare_het"] >= 2),
                "unqueried is not rare; this table does not establish biallelic disease",
            ]
            for gene in TIER1
        ],
    )
    acmg_rows = []
    for gene, bucket in sorted(rare_acmg.items()):
        if not any(bucket[k] for k in ("rare_het", "rare_hom", "unqueried_het")):
            continue
        acmg_rows.append(
            [
                gene,
                bucket["rare_het"],
                bucket["rare_hom"],
                bucket["unqueried_het"],
                "P/LP gene-level is in clinvar_plp_genes.tsv; this is not a returnable SF call",
            ]
        )
    write_tsv(
        dest / "rare_acmg_sf.tsv",
        ["gene", "rare_coding_het", "rare_coding_hom", "unqueried_coding_het", "note"],
        acmg_rows,
    )
    (dest / "scan_counts.json").write_text(json.dumps(dict(counts), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (dest / "inputs.sha256").write_text(
        f"{hashlib.sha256(acmg_path.read_bytes()).hexdigest()}  {acmg_path}\n"
        f"vcf_path  {vcf}\n",
        encoding="utf-8",
    )
    (dest / "caveat.txt").write_text(
        "These tables do not establish causality, phase, or a reportable secondary finding.\n",
        encoding="utf-8",
    )
    print(f"records={counts['records']} plp_genes={len(plp_genes)} keys={counts['keys_impact_het']}")


if __name__ == "__main__":
    main()
