#!/usr/bin/env python3
"""Controlled-data MN5-only VCF preparation and blinded, bounded triage.

The full original and normalized VCF must remain on MN5. No gene/phenotype names
enter candidate selection. All VCF rows (including missing GT and non-PASS)
contribute to per-contig QC; caps are explicitly counted, not hidden exclusions.
"""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import heapq
import json
from pathlib import Path
import re
import sys

CANONICAL = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
HEAD = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "orig_chrom",
    "orig_pos",
    "orig_ref",
    "orig_alts",
    "orig_alt_index",
    "orig_gt",
    "orig_ad_ref",
    "orig_ad_alt",
    "filter",
    "qual",
    "gt",
    "gq",
    "dp",
    "ad_ref",
    "ad_alt",
    "vaf_ad",
    "gnomad_af",
    "gnomad_ac",
    "ref_check",
    "reason",
]


def textopen(path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.endswith(".gz")
        else open(path, encoding="utf-8")
    )


def prepare(args):
    n = 0
    allowed = {
        line.split("\t", 1)[0]
        for line in Path(args.reference_fai).read_text().splitlines()
    }
    unsupported = defaultdict(lambda: {"count": 0, "examples": []})
    with textopen(args.input) as src:
        for line in src:
            if line.startswith("##contig=<ID="):
                line = re.sub(
                    r"^(##contig=<ID=)([^,>]+)",
                    lambda m: m[1] + (m[2] if m[2].startswith("chr") else "chr" + m[2]),
                    line,
                )
            elif line.startswith("#CHROM"):
                sys.stdout.write(
                    '##INFO=<ID=FP_ORIG,Number=1,Type=String,Description="Original VCF chrom:pos:ref:all-ALT order before split/left alignment">\n'
                )
                sys.stdout.write(
                    '##INFO=<ID=FP_ALTIDX,Number=A,Type=Integer,Description="One-based original ALT index (pre-split)">\n'
                )
                sys.stdout.write(
                    '##INFO=<ID=FP_OLDGT,Number=1,Type=String,Description="Original unsplit GT">\n'
                )
                sys.stdout.write(
                    '##INFO=<ID=FP_OLDADREF,Number=1,Type=Integer,Description="Original unsplit reference allele AD">\n'
                )
                sys.stdout.write(
                    '##INFO=<ID=FP_OLDADALT,Number=A,Type=Integer,Description="Original unsplit per-ALT allele AD">\n'
                )
            elif not line.startswith("#"):
                col = line.rstrip("\n").split("\t")
                if len(col) < 10:
                    raise ValueError(
                        "Expected one-sample VCF: malformed or missing sample column"
                    )
                chrom, pos, _, ref, alts = col[:5]
                if not re.fullmatch(r"[A-Za-z0-9_.-]+", chrom):
                    raise ValueError("Unsafe original contig name")
                if any(x in alts for x in ";| "):
                    raise ValueError("Unsupported ALT string for provenance INFO")
                col[0] = chrom if chrom.startswith("chr") else "chr" + chrom
                if col[0] not in allowed:
                    record = unsupported[chrom]
                    record["count"] += 1
                    if len(record["examples"]) < 3:
                        record["examples"].append(
                            {
                                "chrom": chrom,
                                "pos": pos,
                                "ref": ref,
                                "alt": alts,
                                "filter": col[6],
                                "format": col[8],
                                "sample": col[9],
                            }
                        )
                    n += 1
                    if args.limit and n >= args.limit:
                        break
                    continue
                vals = dict(zip(col[8].split(":"), col[9].split(":")))
                oldad = vals.get("AD", ".").split(",")
                if len(oldad) != len(alts.split(",")) + 1:
                    oldad = ["."] * (len(alts.split(",")) + 1)
                prefix = (
                    "FP_ORIG="
                    + chrom
                    + ":"
                    + pos
                    + ":"
                    + ref
                    + ":"
                    + alts.replace(",", "~")
                )
                prefix += ";FP_ALTIDX=" + ",".join(
                    str(i) for i in range(1, len(alts.split(",")) + 1)
                )
                prefix += ";FP_OLDGT=" + vals.get("GT", ".")
                if oldad[0] != ".":
                    prefix += (
                        ";FP_OLDADREF="
                        + oldad[0]
                        + ";FP_OLDADALT="
                        + ",".join(oldad[1:])
                    )
                col[7] = prefix if col[7] == "." else col[7] + ";" + prefix
                line = "\t".join(col) + "\n"
                n += 1
            sys.stdout.write(line)
            if args.limit and n >= args.limit:
                break
    print("prepared_original_records=" + str(n), file=sys.stderr)
    Path(args.unsupported).write_text(
        json.dumps(unsupported, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


def info_parse(field):
    return dict(x.split("=", 1) if "=" in x else (x, "") for x in field.split(";"))


def numeric(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def sort_score(row):
    """Phenotype-blind order; failing FILTER never omitted, penalized only in ranking."""
    filt = int(row["filter"] in ("PASS", "."))
    gq = numeric(row["gq"]) or 0
    dp = numeric(row["dp"]) or 0
    va = numeric(row["vaf_ad"]) or 0
    qual = numeric(row["qual"]) or 0
    return (
        int(row["gnomad_af"] in (".", "0", "0.0")),
        filt,
        min(va, 1 - va) if va else 0,
        min(gq, 99),
        min(dp, 120),
        min(qual, 1000),
    )


def summarize(args):
    destination = Path(args.output)
    if destination.exists():
        raise FileExistsError(destination)
    destination.mkdir(parents=True)
    qc, drops = defaultdict(Counter), Counter()
    kept = defaultdict(list)
    serial = 0
    with textopen(args.input) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 10:
                raise ValueError("Expected one-sample VCF")
            chrom, pos, _, ref, alt, qual, filt, infotxt, fmt, sample = col[:10]
            info = info_parse(infotxt)
            if "FP_ORIG" not in info:
                raise ValueError("Untraceable input: FP_ORIG missing")
            orig = info["FP_ORIG"].split(":", 3)
            if len(orig) != 4:
                raise ValueError("Malformed FP_ORIG")
            old_alts = orig[3].split("~")
            idx = int(info["FP_ALTIDX"])
            if idx > len(old_alts):
                raise ValueError("Original ALT index out of range")
            values = dict(zip(fmt.split(":"), sample.split(":")))
            gt = values.get("GT", ".")
            dp, gq = values.get("DP", "."), values.get("GQ", ".")
            ad = values.get("AD", ".").split(",")
            adref = ad[0] if len(ad) > 0 else "."
            adalt = ad[1] if len(ad) > 1 else "."
            denominator = (numeric(adref) or 0) + (numeric(adalt) or 0)
            vaf = (numeric(adalt) or 0) / denominator if denominator else None
            af = info.get("gnomad_AF", info.get("GNOMAD_AF", "."))
            ac = info.get("gnomad_AC", ".")
            stats = qc[chrom]
            stats["normalized_alleles"] += 1
            stats["pass" if filt in ("PASS", ".") else "nonpass"] += 1
            stats["nocall" if "." in re.split(r"[/|]", gt) else "called"] += 1
            stats[
                "multiallelic_origin" if len(old_alts) > 1 else "single_alt_origin"
            ] += 1
            stats["snv" if len(ref) == len(alt) == 1 else "indel_or_complex"] += 1
            if numeric(dp) is not None:
                stats["dp_sum"] += int(float(dp))
                stats["dp_observed"] += 1
                if float(dp) < 10:
                    stats["dp_lt10"] += 1
            if vaf is not None:
                stats[
                    "vaf_"
                    + (
                        "0_05"
                        if vaf < 0.05
                        else "05_20"
                        if vaf < 0.2
                        else "20_40"
                        if vaf < 0.4
                        else "40_60"
                        if vaf < 0.6
                        else "60_80"
                        if vaf < 0.8
                        else "80_100"
                    )
                ] += 1
            if (
                re.fullmatch(r"[0-9][0-9]?", chrom.removeprefix("chr"))
                and chrom.removeprefix("chr") not in CANONICAL
            ):
                drops["noncanonical_contig"] += 1
            allele_present = (numeric(adalt) or 0) > 0 or bool(
                re.search(r"(^|[/|])1($|[/|])", gt)
            )
            if not allele_present:
                drops["no_alt_support_in_vcf_fields"] += 1
                continue
            if numeric(af) is not None and float(af) > 0.01:
                drops["population_af_gt_0.01"] += 1
                continue
            reason = (
                "rare_af" if numeric(af) is not None else "af_missing_not_assumed_rare"
            )
            if filt not in ("PASS", "."):
                reason += ";nonpass"
            if vaf is not None and vaf < 0.2:
                reason += ";low_vaf"
            if "." in re.split(r"[/|]", gt):
                reason += ";no_call"
            if numeric(dp) is not None and float(dp) < 10:
                reason += ";low_dp"
            row = dict(
                zip(
                    HEAD,
                    (
                        chrom,
                        pos,
                        ref,
                        alt,
                        *orig,
                        str(idx),
                        info.get("FP_OLDGT", "."),
                        info.get("FP_OLDADREF", "."),
                        info.get("FP_OLDADALT", "."),
                        filt,
                        qual,
                        gt,
                        gq,
                        dp,
                        adref,
                        adalt,
                        f"{vaf:.5f}" if vaf is not None else ".",
                        af,
                        ac,
                        "bcftools_norm_c_e",
                        reason,
                    ),
                )
            )
            serial += 1
            item = (sort_score(row), serial, row)
            if len(kept[chrom]) < args.per_contig:
                heapq.heappush(kept[chrom], item)
            elif item > kept[chrom][0]:
                heapq.heapreplace(kept[chrom], item)
                drops["cap_excluded_" + chrom] += 1
            else:
                drops["cap_excluded_" + chrom] += 1
    selected = []
    for chrom, rows in sorted(kept.items()):
        rows.sort(reverse=True)
        # Report truncations explicitly. Never interpret absent allele as reference.
        selected.extend(item[2] for item in rows)
        qc[chrom]["shortlisted"] = min(len(rows), args.per_contig)
    with (destination / "blind_candidates.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=HEAD)
        writer.writeheader()
        writer.writerows(selected)
    with (destination / "per_chrom_qc.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        fields = ["chrom", *sorted(set(k for q in qc.values() for k in q))]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({"chrom": c, **qc[c]} for c in sorted(qc))
    (destination / "exclusion_counts.json").write_text(
        json.dumps(drops, sort_keys=True, indent=2) + "\n"
    )
    (destination / "manifest.json").write_text(
        json.dumps(
            {
                "source": args.input,
                "blind_selection": True,
                "per_contig_cap": args.per_contig,
                "candidates": len(selected),
                "normalized_records": sum(x["normalized_alleles"] for x in qc.values()),
                "annotation_missing_warning": "Population AF missing != rare; no consequence annotation in blind stage",
            },
            indent=2,
        )
        + "\n"
    )
    assert (destination / "blind_candidates.tsv").stat().st_size < 45_000_000, (
        "Candidate table exceeds transfer cap; retain on MN5 and adjust per-contig cap in a new run"
    )


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--input", required=True)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--reference-fai", required=True)
    p.add_argument("--unsupported", required=True)
    p.set_defaults(func=prepare)
    p = sub.add_parser("summarize")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--per-contig", type=int, default=500)
    p.set_defaults(func=summarize)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
