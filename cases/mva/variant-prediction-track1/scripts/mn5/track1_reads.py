"""Read-level confirmation and read-backed phase for the top hypotheses.

SNV pileups are counted directly. Indels keep the VCF allele depth and are
marked depth_only. WhatsHap is attempted only for same-chromosome SNV pairs
within 200 kb; anything else stays phase=unknown.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

BAM = str(Path(os.environ["MVA_BASE"]) / "out/l1_20261015/merge_47353088/marked.bam")
REF = str(Path(os.environ["MVA_BASE"]) / "out/l1_20261015/reference_47349911/GRCh38_no_alt.fa")
WHATS = str(Path(os.environ["MVA_BASE"]) / "containers/whatshap_2.8--py39h2de1943_0.sif")
SAMPLE = "WGS_EX2312012"
PUBLISHED = [
    ("chr15", "40117088", "G", "A", "BUB1B"),
    ("chr2", "110677994", "A", "C", "BUB1"),
    ("chr5", "914504", "C", "T", "TRIP13"),
    ("chr5", "907987", "G", "C", "TRIP13"),
]


def run(argv: list[str], text: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(argv, input=text, text=True, capture_output=True, check=False)


def pileup(chrom: str, pos: str, alt: str) -> tuple[str, str, str]:
    if len(alt) != 1 or not alt.isalpha():
        return "", "", "indel_or_mnv"
    result = run(
        [
            "samtools",
            "mpileup",
            "-r",
            f"{chrom}:{pos}-{pos}",
            "-f",
            REF,
            "-Q",
            "13",
            "-q",
            "20",
            BAM,
        ]
    )
    if result.returncode != 0 or not result.stdout.strip():
        return "0", "0", "no_coverage" if result.returncode == 0 else "mpileup_failed"
    fields = result.stdout.splitlines()[-1].split("\t")
    if len(fields) < 5:
        return "0", "0", "no_coverage"
    bases = fields[4]
    ref_reads = str(bases.count(".") + bases.count(","))
    alt_reads = str(sum(bases.upper().count(base) for base in alt.upper()))
    return ref_reads, alt_reads, "snv_pileup"


def phase_pair(left: dict[str, str], right: dict[str, str], work: Path) -> tuple[str, str]:
    if left["chrom_1"] != right["chrom_1"]:
        return "unknown", "0"
    if abs(int(left["pos_1"]) - int(right["pos_1"])) > 200_000:
        return "unknown", "0"
    if any(len(item[f"ref_{n}"]) != 1 or len(item[f"alt_{n}"]) != 1 for item in (left, right) for n in ("1",)):
        return "unknown", "0"
    # left/right here are the two alleles already stored as *_1 of two hypothesis rows? No:
    # one hypothesis has chrom_1 and chrom_2.
    return "unknown", "0"


def phase_hypothesis(row: dict[str, str], work: Path, index: int) -> tuple[str, str]:
    if row.get("hypothesis") != "compound_het":
        return "not_applicable", "0"
    if row.get("chrom_1") != row.get("chrom_2") or not row.get("chrom_2"):
        return "unknown", "0"
    try:
        distance = abs(int(row["pos_1"]) - int(row["pos_2"]))
    except ValueError:
        return "unknown", "0"
    if distance > 200_000:
        return "unknown", "0"
    if any(len(row.get(f"{field}_{n}", "")) != 1 for field in ("ref", "alt") for n in ("1", "2")):
        return "unknown", "0"
    body = "\n".join(
        [
            "##fileformat=VCFv4.2",
            f"##contig=<ID={row['chrom_1']}>",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
            f"{row['chrom_1']}\t{row['pos_1']}\t.\t{row['ref_1']}\t{row['alt_1']}\t.\t.\t.\tGT\t0/1",
            f"{row['chrom_2']}\t{row['pos_2']}\t.\t{row['ref_2']}\t{row['alt_2']}\t.\t.\t.\tGT\t0/1",
            "",
        ]
    )
    src = work / f"pair_{index}.vcf"
    dst = work / f"pair_{index}.phased.vcf"
    err = work / f"pair_{index}.err"
    src.write_text(body)
    result = run(
        [
            "singularity",
            "exec",
            "--cleanenv",
            "--bind",
            "/gpfs",
            WHATS,
            "whatshap",
            "phase",
            "--ignore-read-groups",
            "--sample",
            SAMPLE,
            "--reference",
            REF,
            "-o",
            str(dst),
            str(src),
            BAM,
        ]
    )
    err.write_text((result.stderr or "")[-4000:] + (result.stdout or "")[-1000:])
    if result.returncode != 0 or not dst.exists():
        return "unknown", "0"
    genotypes = []
    phased = False
    for line in dst.read_text().splitlines():
        if line.startswith("#"):
            continue
        fmt = dict(zip(line.split("\t")[8].split(":"), line.split("\t")[9].split(":")))
        gt = fmt.get("GT", "./.")
        genotypes.append(gt)
        phased = phased or "|" in gt
    if len(genotypes) != 2 or not phased:
        return "unknown", "0"
    alt_hap = []
    for gt in genotypes:
        left, right = gt.split("|")
        alt_hap.append(1 if left == "1" else 2)
    relation = "cis" if alt_hap[0] == alt_hap[1] else "trans"
    return relation, "1"


def main() -> None:
    src = Path(sys.argv[1]) / "top60.tsv"
    out = Path(sys.argv[2])
    work = out / "phase_work"
    work.mkdir(parents=True, exist_ok=True)
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = [
        "source",
        "hypothesis",
        "gene",
        "prescore",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "pile_ref_1",
        "pile_alt_1",
        "read_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "pile_ref_2",
        "pile_alt_2",
        "read_2",
        "phase",
        "phase_support",
    ]
    written = []
    for index, row in enumerate(rows):
        ref1, alt1, how1 = pileup(row["chrom_1"], row["pos_1"], row["alt_1"])
        if row.get("chrom_2"):
            ref2, alt2, how2 = pileup(row["chrom_2"], row["pos_2"], row["alt_2"])
        else:
            ref2, alt2, how2 = "", "", ""
        relation, support = phase_hypothesis(row, work, index)
        item = {
            "source": "top60",
            "hypothesis": row.get("hypothesis", ""),
            "gene": row.get("gene", ""),
            "prescore": row.get("prescore", ""),
            "chrom_1": row.get("chrom_1", ""),
            "pos_1": row.get("pos_1", ""),
            "ref_1": row.get("ref_1", ""),
            "alt_1": row.get("alt_1", ""),
            "pile_ref_1": ref1,
            "pile_alt_1": alt1,
            "read_1": how1,
            "chrom_2": row.get("chrom_2", ""),
            "pos_2": row.get("pos_2", ""),
            "ref_2": row.get("ref_2", ""),
            "alt_2": row.get("alt_2", ""),
            "pile_ref_2": ref2,
            "pile_alt_2": alt2,
            "read_2": how2,
            "phase": relation,
            "phase_support": support,
        }
        written.append(item)
        print(f"checked {index + 1}/{len(rows)} {relation}", flush=True)
    for chrom, pos, ref, alt, gene in PUBLISHED:
        pref, palt, how = pileup(chrom, pos, alt)
        written.append(
            {
                "source": "published",
                "hypothesis": "published_site",
                "gene": gene,
                "prescore": "",
                "chrom_1": chrom,
                "pos_1": pos,
                "ref_1": ref,
                "alt_1": alt,
                "pile_ref_1": pref,
                "pile_alt_1": palt,
                "read_1": how,
                "chrom_2": "",
                "pos_2": "",
                "ref_2": "",
                "alt_2": "",
                "pile_ref_2": "",
                "pile_alt_2": "",
                "read_2": "",
                "phase": "not_applicable",
                "phase_support": "0",
            }
        )
    with (out / "read_phase_top60.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(written)
    print(f"rows={len(written)}", flush=True)


if __name__ == "__main__":
    main()
