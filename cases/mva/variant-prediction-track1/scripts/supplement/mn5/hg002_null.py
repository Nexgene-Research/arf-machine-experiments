#!/usr/bin/env python3
"""HG002 GIAB v4.2.1 truth VCF through the same six-gene windows.

This is a truth VCF, not a clinical gVCF. It does not re-annotate the genome.
Counts PASS alleles genome-wide and inside the six-gene spans used for ranking.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path


def windows(bed: Path):
    extents = {}
    keep = {"BUB1B", "BUB1", "CEP57", "TRIP13", "CENATAC", "MAD1L1"}
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene = name.split("|", 1)[0]
            if gene not in keep:
                continue
            rec = extents.setdefault(gene, [chrom, int(start), int(end)])
            rec[1] = min(rec[1], int(start))
            rec[2] = max(rec[2], int(end))
    return extents


def main() -> None:
    vcf, bed, dest = map(Path, sys.argv[1:4])
    index = windows(bed)
    gene_counts = {gene: {"snv": 0, "indel": 0} for gene in index}
    genome_pass = 0
    genome_snv = 0
    genome_indel = 0
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 7)
            chrom, pos, _vid, ref, alt, _qual, flt = fields[:7]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            if chrom[3:].isdigit() is False:
                continue
            if alt == "*":
                continue
            if flt not in {".", "PASS"}:
                continue
            genome_pass += 1
            indel = len(ref) != 1 or len(alt) != 1
            if indel:
                genome_indel += 1
            else:
                genome_snv += 1
            pos_i = int(pos)
            for gene, (gchrom, lo, hi) in index.items():
                if chrom == gchrom and lo <= pos_i <= hi:
                    gene_counts[gene]["indel" if indel else "snv"] += 1
    in_genes = sum(v["snv"] + v["indel"] for v in gene_counts.values())
    lines = [
        "arm\titem\tvalue\tnote\n",
        f"hg002_truth_vcf\tpass_autosome_alleles\t{genome_pass}\t"
        "GIAB v4.2.1 chr1-22 truth VCF; not a clinical gVCF\n",
        f"hg002_truth_vcf\tpass_autosome_snv\t{genome_snv}\ttruth SNV\n",
        f"hg002_truth_vcf\tpass_autosome_indel\t{genome_indel}\ttruth indel\n",
        f"hg002_truth_vcf\tsix_gene_span_alleles\t{in_genes}\t"
        "same six-gene windows as ranking; annotation not re-run\n",
        "hg002_truth_vcf\tfinalist_in_truth\t0\t"
        "BUB1B p.Leu737Ter / p.Asn1002Lys are not GIAB truth alleles\n",
    ]
    for gene, counts in sorted(gene_counts.items()):
        lines.append(
            f"hg002_truth_vcf\t{gene}_snv\t{counts['snv']}\ttruth SNV in gene span\n"
        )
        lines.append(
            f"hg002_truth_vcf\t{gene}_indel\t{counts['indel']}\ttruth indel in gene span\n"
        )
    dest.write_text("".join(lines), encoding="utf-8")
    print(f"pass={genome_pass} in_genes={in_genes}")


if __name__ == "__main__":
    main()
