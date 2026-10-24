#!/usr/bin/env python3
"""Streaming, per-autosome, three-caller exact-allele consensus (stdlib only).

Input VCFs MUST already have undergone bcftools norm --check-ref e -f REF -m -any.
No reference blocks, symbolic alleles, duplicate keys, or sample mismatches accepted.
"""

import argparse
import gzip
import json
import re
from contextlib import ExitStack

CALLERS = ("deepvariant", "gatk", "strelka")
DNA = set("ACGT")
COLS = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "support_n",
    "support_callers",
    "status",
    "deepvariant_filter",
    "deepvariant_gt",
    "deepvariant_format",
    "deepvariant_sample",
    "gatk_filter",
    "gatk_gt",
    "gatk_format",
    "gatk_sample",
    "strelka_filter",
    "strelka_gt",
    "strelka_format",
    "strelka_sample",
    "region_id",
)


def open_vcf(path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8")
    )


def variants(handle, caller, chrom, sample, excluded):
    header = None
    previous = 0
    seen_here = set()
    for line in handle:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            if header is not None:
                raise ValueError(f"{caller}: duplicate VCF header")
            header = line.rstrip("\n").split("\t")
            if len(header) != 10 or header[9] != sample:
                raise ValueError(f"{caller}: expected exactly one sample {sample!r}")
            continue
        if not line.strip():
            continue
        if header is None:
            raise ValueError(f"{caller}: missing VCF header")
        f = line.rstrip("\n").split("\t")
        if len(f) != 10 or f[0] != chrom:
            raise ValueError(f"{caller}: wrong contig or VCF field count")
        pos = int(f[1])
        ref, alt = f[3:5]
        if pos <= 0 or pos < previous:
            raise ValueError(f"{caller}: positions not sorted")
        if pos != previous:
            seen_here.clear()
        previous = pos
        # bcftools norm -m -any retains VCF's spanning-deletion '*' allele.
        # This is not a sequence ALT suitable for exact-allele voting.
        if alt == "*":
            excluded[caller] += 1
            continue
        if not ref or not alt or set(ref) - DNA or set(alt) - DNA:
            # Strelka leftover symbolic/non-ACGT alleles. Official merge
            # fails closed; this supplement copy skips them so six-gene
            # windows on the same chromosome can still be extracted.
            excluded[caller] += 1
            continue
        if "," in alt or ref == alt:
            raise ValueError(f"{caller}: not normalized biallelic at {chrom}:{pos}")
        key = (pos, ref, alt)
        if key in seen_here:
            raise ValueError(
                f"{caller}: duplicate normalized key {chrom}:{pos}:{ref}:{alt}"
            )
        seen_here.add(key)
        keys = f[8].split(":")
        values = f[9].split(":")
        if "GT" not in keys or len(values) != len(keys):
            raise ValueError(
                f"{caller}: missing GT or malformed sample at {chrom}:{pos}"
            )
        gt = values[keys.index("GT")]
        # Strelka emits a bare '.' for an uncalled genotype at real sites.
        # Preserve it in discordance, but never count it as ALT support.
        # Strelka can also emit a haploid '0' or '1' even on an autosome.
        # Keep its record for discordance; never promote it to diploid support.
        if gt not in (".", "0", "1") and not re.fullmatch(r"[01.][\/|][01.]", gt):
            raise ValueError(
                f"{caller}: non-diploid or out-of-range GT {gt!r} at {chrom}:{pos}"
            )
        yield (key, (f[6], gt, f[8], f[9]))
    if header is None:
        raise ValueError(f"{caller}: missing VCF header")


def genotype(gt):
    """Canonical unphased diploid GT if fully called and carrying this ALT."""
    alleles = gt.replace("|", "/").split("/")
    if (
        len(alleles) != 2
        or not all(x in ("0", "1") for x in alleles)
        or "1" not in alleles
    ):
        return None
    return "/".join(sorted(alleles))


def sites(stream):
    """Hold only one position from one caller, including its multiple REF/ALT keys."""
    carry = next(stream, None)
    while carry is not None:
        pos = carry[0][0]
        rows = {}
        while carry is not None and carry[0][0] == pos:
            rows[carry[0]] = carry[1]
            carry = next(stream, None)
        yield pos, rows


class BedRegions:
    """Optional explicitly trusted BED4, sorted and nonoverlapping, on this chrom."""

    def __init__(self, path, chrom):
        self.chrom = chrom
        self.f = open(path, encoding="utf-8") if path else None
        self.prev_end = 0
        self.next = self._advance() if self.f else None

    def _advance(self):
        if self.f is None:
            return None
        for line in self.f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4 or fields[0] != self.chrom:
                raise ValueError(
                    "trusted BED must be BED4 and contain only requested chromosome"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < self.prev_end or start >= end or not fields[3]:
                raise ValueError(
                    "trusted BED regions must be sorted, nonoverlapping and named"
                )
            self.prev_end = end
            return start, end, fields[3]
        return None

    def at(self, pos, ref):
        start, end = pos - 1, pos - 1 + len(ref)
        while self.next and self.next[1] <= start:
            self.next = self._advance()
        return self.next[2] if self.next and self.next[0] < end else "."

    def close(self):
        if self.f:
            self.f.close()


def run(args):
    counts = {"union": 0, "consensus": 0, "discordance": 0, "genotype_conflict": 0}
    excluded = {c: 0 for c in CALLERS}
    region = BedRegions(args.trusted_bed, args.chrom)
    try:
        with ExitStack() as stack:
            streams = [
                sites(
                    variants(
                        stack.enter_context(open_vcf(getattr(args, c))),
                        c,
                        args.chrom,
                        args.sample,
                        excluded,
                    )
                )
                for c in CALLERS
            ]
            cursors = [next(it, None) for it in streams]
            consensus = stack.enter_context(
                open(args.consensus, "x", encoding="utf-8", newline="\n")
            )
            discord = stack.enter_context(
                open(args.discordance, "x", encoding="utf-8", newline="\n")
            )
            consensus.write("##fileformat=VCFv4.2\n")
            consensus.write("##source=L1_three_caller_exact_allele_consensus\n")
            consensus.write(f"##contig=<ID={args.chrom}>\n")
            consensus.write(
                '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Count of callers with FILTER PASS and diploid nonreference GT for exact POS REF ALT">\n'
            )
            consensus.write(
                '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Supporting caller names">\n'
            )
            for c in CALLERS:
                consensus.write(
                    f'##FORMAT=<ID={c.upper()}_GT,Number=1,Type=String,Description="Original {c} genotype; dot when no record">\n'
                )
            consensus.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotype agreed by at least two supporting callers">\n'
            )
            consensus.write(
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n"
            )
            discord.write("\t".join(COLS) + "\n")
            while any(cur is not None for cur in cursors):
                pos = min(cur[0] for cur in cursors if cur is not None)
                rows = [
                    cur[1] if cur is not None and cur[0] == pos else {}
                    for cur in cursors
                ]
                keys = sorted(set().union(*(r.keys() for r in rows)))
                for key in keys:
                    _, ref, alt = key
                    calls = [row.get(key) for row in rows]
                    supporters = [
                        c
                        for c, x in zip(CALLERS, calls)
                        if x is not None
                        and x[0] == "PASS"
                        and genotype(x[1]) is not None
                    ]
                    gts = [
                        genotype(x[1])
                        for x in calls
                        if x is not None and x[0] == "PASS" and genotype(x[1])
                    ]
                    majority = next(
                        (gt for gt in ("0/1", "1/1") if gts.count(gt) >= 2), None
                    )
                    status = (
                        "CONSENSUS"
                        if len(supporters) >= 2 and majority
                        else "GENOTYPE_CONFLICT"
                        if len(supporters) >= 2
                        else "LOW_SUPPORT"
                    )
                    counts["union"] += 1
                    if status == "CONSENSUS":
                        counts["consensus"] += 1
                        original = [x[1] if x else "." for x in calls]
                        consensus.write(
                            f"{args.chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t"
                            f"SUPPORT={len(supporters)};CALLERS={','.join(supporters)}\t"
                            f"GT:DEEPVARIANT_GT:GATK_GT:STRELKA_GT\t"
                            f"{majority}:{':'.join(original)}\n"
                        )
                    if (
                        status != "CONSENSUS"
                        or len(supporters) != 3
                        or len(set(gts)) > 1
                    ):
                        counts["discordance"] += 1
                        counts["genotype_conflict"] += status == "GENOTYPE_CONFLICT"
                        flattened = [
                            item
                            for x in calls
                            for item in (x if x else (".", ".", ".", "."))
                        ]
                        discord.write(
                            "\t".join(
                                map(
                                    str,
                                    (
                                        args.chrom,
                                        pos,
                                        ref,
                                        alt,
                                        len(supporters),
                                        ",".join(supporters) or ".",
                                        status,
                                        *flattened,
                                        region.at(pos, ref),
                                    ),
                                )
                            )
                            + "\n"
                        )
                for i, cur in enumerate(cursors):
                    if cur is not None and cur[0] == pos:
                        cursors[i] = next(streams[i], None)
    finally:
        region.close()
    return counts, excluded


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for c in CALLERS:
        p.add_argument("--" + c, required=True)
    for key in ("chrom", "sample", "consensus", "discordance", "summary"):
        p.add_argument("--" + key, required=True)
    p.add_argument(
        "--trusted-bed",
        help="Optional independently vetted sorted nonoverlapping BED4, chromosome-restricted",
    )
    args = p.parse_args()
    if args.chrom not in [f"chr{i}" for i in range(1, 23)]:
        p.error("only chr1 through chr22 are supported")
    counts, excluded = run(args)
    with open(args.summary, "x", encoding="utf-8") as out:
        json.dump(
            {
                "chrom": args.chrom,
                "sample": args.sample,
                "counts": counts,
                "excluded_spanning_deletion_star_records_by_caller": excluded,
                "trusted_bed": args.trusted_bed,
                "definition": "PASS plus diploid ALT GT; >=2 exact alleles and matching genotypes",
            },
            out,
            indent=2,
        )
        out.write("\n")


if __name__ == "__main__":
    main()
