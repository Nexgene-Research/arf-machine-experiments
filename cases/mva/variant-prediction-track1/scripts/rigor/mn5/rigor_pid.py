#!/usr/bin/env python3
"""GATK PID/PGT labels for HIGH/MODERATE heterozygous keys.

Caller-phased PID is not WhatsHap and does not prove physical haplotypy.
This table does not establish causality.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

TIER1 = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1")
JURY = {
    ("chr15", "40209701", "T", "G"): "BUB1B_p.Leu737Ter",
    ("chr15", "40220612", "T", "G"): "BUB1B_p.Asn1002Lys",
}


def prefix(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else "chr" + chrom


def hap_of(pgt: str) -> str | None:
    if pgt == "0|1":
        return "1"
    if pgt == "1|0":
        return "0"
    if pgt == "1|1":
        return "both"
    return None


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def load_keys(path: Path) -> dict[tuple[str, str, str, str], dict[str, str]]:
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[(row["chrom"], row["pos"], row["ref"], row["alt"])] = row
    return out


def pair_call(a_hap: str | None, b_hap: str | None) -> str:
    if a_hap in (None, "both") or b_hap in (None, "both"):
        return "unknown"
    if a_hap == b_hap:
        return "cis"
    return "trans"


def main() -> None:
    keys_path = Path(sys.argv[1])
    gatk_root = Path(sys.argv[2])
    dest = Path(sys.argv[3])
    dest.mkdir(parents=True, exist_ok=True)
    keys = load_keys(keys_path)
    counts: Counter[str] = Counter()
    pid_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    jury_hit: dict[str, dict[str, str]] = {}
    chroms = [f"chr{n}" for n in range(1, 23)] + ["chrX"]
    for chrom in chroms:
        vcf = gatk_root / chrom / "normalized.vcf.gz"
        if not vcf.is_file():
            counts[f"missing_{chrom}"] += 1
            continue
        opener = gzip.open
        with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                rec_chrom = prefix(fields[0])
                pos, _vid, ref, alt = fields[1:5]
                key = (rec_chrom, pos, ref, alt)
                meta = keys.get(key)
                if meta is None:
                    continue
                fmt = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}
                pid = fmt.get("PID", ".")
                pgt = fmt.get("PGT", ".")
                gt = fmt.get("GT", ".")
                counts["matched_keys"] += 1
                rec = {
                    "chrom": rec_chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": meta["gene"],
                    "impact": meta["impact"],
                    "pid": pid,
                    "pgt": pgt,
                    "gt": gt,
                    "hap": hap_of(pgt) or "",
                }
                if key in JURY:
                    jury_hit[JURY[key]] = rec
                if pid not in {"", ".", None}:
                    counts["with_pid"] += 1
                    pid_groups[(rec_chrom, pid)].append(rec)
                else:
                    counts["without_pid"] += 1

    pair_counts: Counter[str] = Counter()
    six: dict[str, Counter[str]] = {g: Counter() for g in TIER1}
    for (_chrom, pid), members in pid_groups.items():
        by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
        for rec in members:
            by_gene[rec["gene"]].append(rec)
        for gene, recs in by_gene.items():
            if len(recs) < 2:
                continue
            pair_counts["gene_pid_blocks_ge2"] += 1
            if gene in six:
                six[gene]["pid_blocks_ge2"] += 1
            for i, left in enumerate(recs):
                for right in recs[i + 1 :]:
                    call = pair_call(left["hap"] or None, right["hap"] or None)
                    pair_counts[call] += 1
                    if gene in six:
                        six[gene][call] += 1

    write_tsv(
        dest / "pid_pair_counts.tsv",
        ["metric", "count", "note"],
        [
            ["matched_high_moderate_hets", counts["matched_keys"], "intersection of GATK GT with scan keys"],
            ["with_pid", counts["with_pid"], "FORMAT/PID present"],
            ["without_pid", counts["without_pid"], "unphased by the caller"],
            ["gene_pid_blocks_ge2", pair_counts["gene_pid_blocks_ge2"], "same PID, same gene, >=2 hets"],
            ["cis_pairs", pair_counts["cis"], "same PGT haplotype; caller label only"],
            ["trans_pairs", pair_counts["trans"], "opposite PGT haplotype; caller label only"],
            ["unknown_pairs", pair_counts["unknown"], "PID without interpretable PGT"],
        ],
    )
    write_tsv(
        dest / "pid_sixgene.tsv",
        ["gene", "pid_blocks_ge2", "cis_pairs", "trans_pairs", "unknown_pairs", "note"],
        [
            [
                gene,
                six[gene]["pid_blocks_ge2"],
                six[gene]["cis"],
                six[gene]["trans"],
                six[gene]["unknown"],
                "does not replace WhatsHap unknown; not physical phase",
            ]
            for gene in TIER1
        ],
    )
    stop = jury_hit.get("BUB1B_p.Leu737Ter")
    miss = jury_hit.get("BUB1B_p.Asn1002Lys")
    same = (
        "not_both_in_gatk"
        if not (stop and miss)
        else (
            "same_pid"
            if stop["pid"] not in {".", ""} and stop["pid"] == miss["pid"]
            else "different_or_missing_pid"
        )
    )
    call = "unknown"
    if stop and miss and same == "same_pid":
        call = pair_call(stop["hap"] or None, miss["hap"] or None)
    write_tsv(
        dest / "pid_jury.tsv",
        [
            "allele_a",
            "allele_b",
            "a_pid",
            "b_pid",
            "a_pgt",
            "b_pgt",
            "pid_status",
            "caller_phase",
            "note",
        ],
        [
            [
                "BUB1B_p.Leu737Ter",
                "BUB1B_p.Asn1002Lys",
                stop["pid"] if stop else ".",
                miss["pid"] if miss else ".",
                stop["pgt"] if stop else ".",
                miss["pgt"] if miss else ".",
                same,
                call,
                "GATK PID/PGT only; WhatsHap remains unknown; not a causal claim",
            ]
        ],
    )
    (dest / "caveat.txt").write_text(
        "GATK PID/PGT is the caller's local phase set, not physical haplotypy and not causality.\n",
        encoding="utf-8",
    )
    print(f"matched={counts['matched_keys']} cis={pair_counts['cis']} trans={pair_counts['trans']}")


if __name__ == "__main__":
    main()
