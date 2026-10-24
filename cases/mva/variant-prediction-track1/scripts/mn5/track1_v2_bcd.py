#!/usr/bin/env python3
"""Bounded, aggregate-only Tier-1 B/C/D screen; raw intermediate VCF/BAM stay MN5.

Inputs are a locked six-gene GRCh38 MANE model, exon BED, candidate rows,
and the indexed annotated VCF. All generated aggregate outputs contain neither
read identifiers nor sequences. Read-level work is restricted to six loci.
"""

import csv
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys

B = Path(os.environ["MVA_BASE"])
BAM = B / "out/l1_20261015/merge_47353088/marked.bam"
FA = B / "out/l1_20261015/reference_47349911/GRCh38_no_alt.fa"
VCF = B / "out/critical_annotation_47397018/with_clinvar.vcf.gz"
BC = B / "containers/bcftools_1.24--h487d631_1.sif"
WH = B / "containers/whatshap_2.8--py39h2de1943_0.sif"
SAMPLE = "WGS_EX2312012"
GENES = ("BUB1B", "BUB1", "CEP57", "TRIP13", "CENATAC", "MAD1L1")
ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])
RAW = OUT / "mn5_only_raw"
RAW.mkdir(parents=True, exist_ok=False)


def run(args, output=None):
    print("RUN", " ".join(map(str, args[:9])), flush=True)
    if output is None:
        p = subprocess.run(list(map(str, args)), text=True, capture_output=True)
    else:
        with open(output, "w") as fp:
            p = subprocess.run(
                list(map(str, args)), stdout=fp, stderr=subprocess.PIPE, text=True
            )
    if p.returncode:
        raise RuntimeError(f"{args[0]} exit={p.returncode}: {p.stderr[-1800:]}")
    return p.stdout if output is None else ""


def bc(*args, output=None):
    return run(
        ("singularity", "exec", "--cleanenv", "--bind", "/gpfs", BC, "bcftools", *args),
        output,
    )


def wh(*args):
    return run(
        ("singularity", "exec", "--cleanenv", "--bind", "/gpfs", WH, "whatshap", *args)
    )


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def table(name, data, columns):
    target = OUT / name
    with target.open("x", newline="") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(data)


def info_map(s):
    return {p.split("=", 1)[0]: p.split("=", 1)[1] for p in s.split(";") if "=" in p}


def variant(line):
    p = line.rstrip("\n").split("\t")
    i = info_map(p[7])
    fmt = dict(zip(p[8].split(":"), p[9].split(":")))
    return {
        "chrom": p[0],
        "pos": int(p[1]),
        "ref": p[3],
        "alt": p[4],
        "filter": p[6],
        "qual": p[5],
        "info": i,
        "fmt": fmt,
        "key": "-".join((p[0], p[1], p[3], p[4])),
    }


def numeric(x, fallback=-1):
    try:
        return float(x)
    except (ValueError, TypeError):
        return fallback


def annot(v, gene, mane):
    entries = [x.split("|") for x in v["info"].get("ANN", "").split(",")]
    hits = [p for p in entries if len(p) > 6 and p[3] == gene and p[6] == mane]
    p = hits[0] if hits else []
    nmd = (
        "unknown_no_MANE_ANN"
        if not p
        else (
            "predicted_NMD_review_transcript"
            if ("stop_gained" in p[1] or "frameshift" in p[1])
            else "not_inferred"
        )
    )
    return (
        p[1] if p else "unknown_MANE_effect",
        p[9] if len(p) > 9 else "",
        p[10] if len(p) > 10 else "",
        nmd,
    )


def read_events(fields, sites):
    """Count anchor and simple SNV/indel events without exporting SAM records."""
    cigar = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
    seq = fields[9]
    bq = fields[10]
    refpos = int(fields[3])
    qpos = 0
    bases = {}
    insertions = {}
    deletions = {}
    aligned = 0
    soft = 0
    for n0, op in cigar:
        n = int(n0)
        if op in "M=X":
            for pos in sites:
                if refpos <= pos < refpos + n:
                    idx = qpos + pos - refpos
                    bases[pos] = (seq[idx], ord(bq[idx]) - 33)
            refpos += n
            qpos += n
            aligned += n
        elif op == "I":
            insertions[refpos - 1] = seq[qpos : qpos + n]
            qpos += n
        elif op in "DN":
            if op == "D":
                deletions[refpos - 1] = n
            refpos += n
        elif op == "S":
            qpos += n
            soft += n
    return bases, insertions, deletions, soft / max(1, aligned + soft)


def classify(v, bases, ins, dels):
    pos = v["pos"]
    ref = v["ref"]
    alt = v["alt"]
    base = bases.get(pos)
    if not base or base[1] < 20:
        return None
    if len(ref) == len(alt) == 1:
        return (
            "alt"
            if base[0].upper() == alt
            else ("ref" if base[0].upper() == ref else None)
        )
    if len(alt) > len(ref) and alt.startswith(ref):
        return (
            "alt"
            if ins.get(pos + len(ref) - 1) == alt[len(ref) :]
            else (
                "ref"
                if base[0].upper() == ref[0] and pos + len(ref) - 1 not in ins
                else None
            )
        )
    if len(ref) > len(alt) and ref.startswith(alt):
        return (
            "alt"
            if dels.get(pos + len(alt) - 1) == len(ref) - len(alt)
            else (
                "ref"
                if base[0].upper() == ref[0] and pos + len(alt) - 1 not in dels
                else None
            )
        )
    return None


def read_metrics(gene, chrom, start, end, target):
    counts = {
        v["key"]: {
            "ref": 0,
            "alt": 0,
            "forward": 0,
            "reverse": 0,
            "soft": 0,
            "sa": 0,
            "xa": 0,
        }
        for v in target
    }
    sites = {v["pos"] for v in target}
    signal = {
        "reads": 0,
        "discordant_pairs": 0,
        "supplementary": 0,
        "softclip_20bp": 0,
        "SA_tag": 0,
    }
    # Stream SAM one record at a time; never persist read names or sequences.
    p = subprocess.Popen(
        ["samtools", "view", "-F", "0x504", str(BAM), f"{chrom}:{start}-{end}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        signal["reads"] += 1
        if flag & 0x800:
            signal["supplementary"] += 1
        if flag & 0x2 == 0 and flag & 0x1:
            signal["discordant_pairs"] += 1
        cigar = f[5]
        clip = sum(int(n) for n in re.findall(r"(\d+)S", cigar))
        if clip >= 20:
            signal["softclip_20bp"] += 1
        if any(t.startswith("SA:Z:") for t in f[11:]):
            signal["SA_tag"] += 1
        if int(f[4]) < 60 or (flag & 0x400):
            continue
        b, i, d, soft = read_events(f, sites)
        for v in target:
            kind = classify(v, b, i, d)
            if kind is None:
                continue
            c = counts[v["key"]]
            c[kind] += 1
            if kind == "alt":
                c["reverse" if flag & 0x10 else "forward"] += 1
                c["soft"] += int(soft >= 0.15)
                c["sa"] += int(any(t.startswith("SA:Z:") for t in f[11:]))
                c["xa"] += int(any(t.startswith("XA:Z:") for t in f[11:]))
    stderr = p.stderr.read()
    rc = p.wait()
    if rc:
        raise RuntimeError(f"samtools {gene} exit={rc} {stderr[-1500:]}")
    return counts, signal


def parse_phase(path, pairs, gene, chrom):
    phased = {}
    blocks = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        v = variant(line)
        gt = v["fmt"].get("GT", "./.")
        ps = v["fmt"].get("PS", ".")
        if "|" in gt and ps not in (".", ""):
            phased[v["key"]] = (ps, gt)
            blocks.setdefault(ps, []).append(v["pos"])
    records = []
    for a, b in pairs:
        x = phased.get(a)
        y = phased.get(b)
        same = bool(x and y and x[0] == y[0])
        p = blocks[x[0]] if same else []
        relation = "unknown"
        if same:
            relation = "cis" if x[1].index("1") == y[1].index("1") else "trans"
        records.append(
            {
                "gene": gene,
                "allele_1": a,
                "allele_2": b,
                "same_block": str(same).lower(),
                "phase": relation,
                "block_start": min(p) if p else "",
                "block_end": max(p) if p else "",
                "block_span_bp": max(p) - min(p) + 1 if p else "",
                "informative_fragments": "not_reported_by_whatshap_VCF",
                "link_type": "chained_or_direct_unresolved" if same else "none",
            }
        )
    return records, blocks


def main():
    model = json.loads((ROOT / "gene_models.json").read_text())
    assert model["assembly"] == "GRCh38"
    models = {g["gene"]: g for g in model["gene_models"]}
    assert set(GENES) <= set(models)
    target = rows(ROOT / "tier1_bounded_calls.tsv")
    top = rows(ROOT / "track1_top10_ranked.tsv")
    assert all(t["region_gene"] in GENES for t in target)
    assert bc("query", "-l", VCF).strip() == SAMPLE
    exon = [
        line.rstrip("\n").split("\t")
        for line in open(ROOT / "exon_utr.bed")
        if line.startswith("chr")
    ]
    validation = []
    concordance = []
    screen = []
    matrix = []
    keys = []
    blocksall = []
    phaseall = []
    for gene in GENES:
        print("STAGE GENE", gene, flush=True)
        m = models[gene]
        chrom = m["chrom"]
        lo = int(m["gene_start_1based"])
        hi = int(m["gene_end_1based_inclusive"])
        start = max(1, lo - 100000)
        end = hi + 100000
        span = f"{chrom}:{start}-{end}"
        # Normalization checks the FASTA REF, left-aligns indels and splits multiallelics.
        source = RAW / f"{gene}.normalized.vcf"
        cmd = f"singularity exec --cleanenv --bind /gpfs {BC} bcftools view -r {span} -Ou {VCF} | singularity exec --cleanenv --bind /gpfs {BC} bcftools norm -f {FA} -m -any -c e -Ov -o {source}"
        run(["bash", "-o", "pipefail", "-c", cmd])
        var = {}
        for line in open(source):
            if line.startswith("#"):
                continue
            v = variant(line)
            var[v["key"]] = v
        q = [
            v
            for v in var.values()
            if v["fmt"].get("GT") in ("0/1", "1/0")
            and v["filter"] in ("PASS", ".")
            and numeric(v["qual"]) >= 20
            and numeric(v["fmt"].get("GQ")) >= 20
            and numeric(v["fmt"].get("DP")) >= 10
        ]
        # ALL passing heterozygous sites, regardless of rarity, AF or consequence.
        wanted = {v["key"] for v in q}
        regional = RAW / f"{gene}.all_hets.vcf"
        with regional.open("x") as fp:
            for line in open(source):
                if line.startswith("#") or variant(line)["key"] in wanted:
                    fp.write(line)
        phased = RAW / f"{gene}.phased.vcf"
        blocklist = RAW / f"{gene}.blocks.tsv"
        if q:
            wh(
                "phase",
                "--ignore-read-groups",
                "--sample",
                SAMPLE,
                "--reference",
                FA,
                "-o",
                phased,
                regional,
                BAM,
            )
            wh("stats", "--sample", SAMPLE, "--block-list", blocklist, phased)
        else:
            phased.write_text(regional.read_text())
            blocklist.write_text("")
        genepairs = []
        for r in top:
            if r["gene"] != gene or r["hypothesis"] != "compound_het":
                continue
            a = "-".join((r["chrom_1"], r["pos_1"], r["ref_1"], r["alt_1"]))
            b = "-".join((r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"]))
            genepairs.append((a, b))
        p, block = parse_phase(phased, genepairs, gene, chrom)
        phaseall += p
        for ps, positions in block.items():
            blocksall.append(
                {
                    "gene": gene,
                    "chrom": chrom,
                    "block_id_local": ps,
                    "start_1based": min(positions),
                    "end_1based": max(positions),
                    "span_bp": max(positions) - min(positions) + 1,
                    "phased_hets": len(positions),
                    "total_high_quality_hets_gene_window": len(q),
                }
            )
        # Recall from the separately aligned BAM, on the same individual.
        recall = RAW / f"{gene}.recall.vcf"
        c = f"singularity exec --cleanenv --bind /gpfs {BC} bcftools mpileup -f {FA} -r {span} -q 20 -Q 20 -Ou {BAM} | singularity exec --cleanenv --bind /gpfs {BC} bcftools call -m -v -Ov -o {recall}"
        run(["bash", "-o", "pipefail", "-c", c])
        recalled = {}
        for line in open(recall):
            if line.startswith("#"):
                continue
            v = variant(line)
            recalled[v["key"]] = v
        gene_target = []
        for row in target:
            if row["region_gene"] != gene:
                continue
            k = "-".join((row["chr"], row["pos"], row["ref"], row["alt"]))
            gene_target.append(
                var.get(
                    k,
                    {
                        "key": k,
                        "pos": int(row["pos"]),
                        "ref": row["ref"],
                        "alt": row["alt"],
                    },
                )
            )
        # Add all rare gene-body calls, including low-VAF and non-PASS, capped only in exported table.
        rare = [
            v
            for v in var.values()
            if lo <= v["pos"] <= hi
            and (numeric(v["info"].get("gnomad_AF"), 0) <= 0.01)
            and v["fmt"].get("GT") not in ("0/0", "./.", ".")
        ]
        extra = {v["key"]: v for v in gene_target}
        for v in rare:
            extra[v["key"]] = v
        # Separate alignment metrics evaluated on all rare gene-body alleles; no genome-wide pass.
        metrics, signal = read_metrics(
            gene, chrom, max(1, lo - 50000), hi + 50000, list(extra.values())
        )
        for v in extra.values():
            k = v["key"]
            c = metrics[k]
            original = var.get(k)
            re_call = recalled.get(k)
            meta = original or {
                "fmt": {},
                "info": {},
                "filter": "missing",
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
            }
            effect, hgvsc, hgvsp, nmd = annot(meta, gene, m["transcript"])
            allele = {
                "gene": gene,
                "grch38_key": k,
                "ref_check": "pass" if original else "not_exact_after_normalization",
                "source_gt": meta["fmt"].get("GT", "missing"),
                "source_filter": meta["filter"],
                "source_DP": meta["fmt"].get("DP", ""),
                "source_AD": meta["fmt"].get("AD", ""),
                "read_DP_MQ60_BQ20": c["ref"] + c["alt"],
                "read_ref_MQ60_BQ20": c["ref"],
                "read_alt_MQ60_BQ20": c["alt"],
                "alt_forward": c["forward"],
                "alt_reverse": c["reverse"],
                "alt_softclip_fraction": round(c["soft"] / c["alt"], 4)
                if c["alt"]
                else "",
                "alt_SA_count": c["sa"],
                "alt_XA_count": c["xa"],
                "recall_gt": re_call["fmt"].get("GT", "not_called")
                if re_call
                else "not_called",
                "recall_DP": re_call["fmt"].get("DP", "") if re_call else "",
                "caller_relation": (
                    "concordant"
                    if re_call and re_call["fmt"].get("GT") == meta["fmt"].get("GT")
                    else "discordant_or_not_called"
                ),
                "gnomad_AF": meta["info"].get("gnomad_AF", ""),
                "gnomad_HomAlt": meta["info"].get("gnomad_HomAlt", ""),
                "ClinVar_CLNSIG": meta["info"].get("CLNSIG", ""),
                "ClinVar_CLNREVSTAT": meta["info"].get("CLNREVSTAT", ""),
                "MANE_transcript": m["transcript"],
                "MANE_consequence": effect,
                "MANE_HGVSc": hgvsc,
                "MANE_HGVSp": hgvsp,
                "NMD_expectation": nmd,
                "AM_REVEL": "not_available_in_this_VCF_uncalibrated_LR_1",
            }
            validation.append(allele)
            if v in gene_target or k in {x["key"] for x in gene_target}:
                concordance.append(
                    {
                        key: allele[key]
                        for key in (
                            "gene",
                            "grch38_key",
                            "source_gt",
                            "recall_gt",
                            "caller_relation",
                            "source_DP",
                            "recall_DP",
                            "read_alt_MQ60_BQ20",
                        )
                    }
                )
            if (
                numeric(meta["info"].get("gnomad_AF"), 0) <= 0.01
                and effect
                in (
                    "intron_variant",
                    "upstream_gene_variant",
                    "downstream_gene_variant",
                    "5_prime_UTR_variant",
                    "3_prime_UTR_variant",
                    "intergenic_region",
                )
                and len(v["ref"]) == len(v["alt"]) == 1
            ):
                keys.append(
                    {
                        "gene": gene,
                        "grch38_key": k,
                        "gnomad_AF": meta["info"].get("gnomad_AF", ""),
                        "MANE_consequence": effect,
                        "read_alt_MQ60_BQ20": c["alt"],
                    }
                )
        # Depth: independent exon-level samtools depth, compare gene-body median to +/-50kb flanks.
        exonrows = [x for x in exon if x[3].split("|")[0] == gene]
        exdepth = {x[3]: [] for x in exonrows}
        flanks = []
        dp = subprocess.Popen(
            [
                "samtools",
                "depth",
                "-aa",
                "-r",
                f"{chrom}:{max(1, lo - 50000)}-{hi + 50000}",
                str(BAM),
            ],
            stdout=subprocess.PIPE,
            text=True,
        )
        for line in dp.stdout:
            _, pos0, dep = line.split()[:3]
            pos = int(pos0)
            d = int(dep)
            if pos < lo or pos > hi:
                flanks.append(d)
            for x in exonrows:
                if int(x[1]) < pos <= int(x[2]):
                    exdepth[x[3]].append(d)
        if dp.wait():
            raise RuntimeError("samtools depth failed " + gene)
        flank = statistics.median(flanks) if flanks else 0
        for x in exonrows:
            values = exdepth[x[3]]
            med = statistics.median(values) if values else 0
            screen.append(
                {
                    "gene": gene,
                    "segment": x[3],
                    "chrom": chrom,
                    "start0": x[1],
                    "end0": x[2],
                    "median_depth": med,
                    "flank_median_depth": flank,
                    "depth_ratio": round(med / flank, 4) if flank else "",
                    "discordant_pair_reads_gene_50kb": signal["discordant_pairs"],
                    "supplementary_reads_gene_50kb": signal["supplementary"],
                    "SA_tag_reads_gene_50kb": signal["SA_tag"],
                    "softclip_20bp_reads_gene_50kb": signal["softclip_20bp"],
                    "total_reads_gene_50kb": signal["reads"],
                    "interpretation": "screen_only_not_CNV_or_SV_confirmation",
                }
            )
        gvals = [a for a in validation if a["gene"] == gene]
        matrix.append(
            {
                "gene": gene,
                "rare_gene_body_alleles_including_nonPASS": len(rare),
                "bounded_candidate_alleles": len(gene_target),
                "all_high_quality_hets_gene_100kb": len(q),
                "phase_blocks": len(block),
                "candidate_pairs": len(genepairs),
                "candidate_pairs_trans": sum(x["phase"] == "trans" for x in p),
                "candidate_pairs_cis": sum(x["phase"] == "cis" for x in p),
                "candidate_pairs_unknown": sum(x["phase"] == "unknown" for x in p),
                "bounded_alleles_alt_read_supported": sum(
                    int(a["read_alt_MQ60_BQ20"]) > 0 for a in gvals
                ),
                "lowest_exon_flank_depth_ratio": min(
                    (
                        float(x["depth_ratio"])
                        for x in screen
                        if x["gene"] == gene and x["depth_ratio"] != ""
                    ),
                    default="",
                ),
                "interpretation": "screen_only; trans only if in same phased block",
            }
        )
        print(
            "STAGE DONE",
            gene,
            "hets",
            len(q),
            "rare",
            len(rare),
            "blocks",
            len(block),
            flush=True,
        )
    valfields = list(validation[0])
    table("allele_validation.tsv", validation[:300], valfields)
    table("caller_concordance.tsv", concordance, list(concordance[0]))
    table(
        "phase_blocks_summary.tsv",
        blocksall,
        [
            "gene",
            "chrom",
            "block_id_local",
            "start_1based",
            "end_1based",
            "span_bp",
            "phased_hets",
            "total_high_quality_hets_gene_window",
        ],
    )
    table(
        "phase_pair_evidence.tsv",
        phaseall,
        [
            "gene",
            "allele_1",
            "allele_2",
            "same_block",
            "phase",
            "block_start",
            "block_end",
            "block_span_bp",
            "informative_fragments",
            "link_type",
        ],
    )
    table("second_hit_screen.tsv", screen, list(screen[0]))
    table("tier1_evidence_matrix.tsv", matrix, list(matrix[0]))
    table(
        "rare_noncoding_keys.tsv",
        keys[:200],
        ["gene", "grch38_key", "gnomad_AF", "MANE_consequence", "read_alt_MQ60_BQ20"],
    )
    (OUT / "allele_validation_summary.md").write_text(
        "# Tier-1 allele authenticity\n\n"
        + f"Authenticated against GRCh38; independent bwa-mem2 alignment with bcftools caller, same individual. {len(validation)} alleles inspected ({min(300, len(validation))} emitted), {len(concordance)} prespecified bounded alleles. MAPQ>=60 BQ>=20; duplicate and secondary alignments excluded. Missing call is not absence.\n\n"
        + "Read flags summarize SA/XA and soft clips; one caller is not independent biological replication. See tables for per-allele values.\n"
    )
    (OUT / "phase_evidence.md").write_text(
        "# Read-backed regional phase\n\n"
        + f"WhatsHap 2.8 with ALL high-quality heterozygotes within each of six genes ±100 kb; {len(blocksall)} phased blocks, {len(phaseall)} prespecified pairs assessed.\n\n"
        + "\n".join(
            f"- {x['gene']} {x['allele_1']} / {x['allele_2']}: {x['phase']}; same block {x['same_block']}; span {x['block_span_bp'] or 'unknown'} bp; direct vs chained and informative fragment counts unresolved (Whatshap VCF alone cannot distinguish)."
            for x in phaseall
        )
        + "\n"
    )
    (OUT / "second_hit_summary.md").write_text(
        "# Six-gene second-hit screen\n\n"
        + f"{len(screen)} exon/UTR intervals; {len(keys[:200])} rare SNV noncoding keys exported. Depth ratio is median segment depth / median ±50 kb non-gene flank depth.\n\n"
        + "Discordant-pair, supplementary, SA and long-soft-clip counts are nonspecific gene-window signals, **not CNV or SV confirmations**. NMD label is tentative: transcript-specific exon-junction position review is needed. AM and REVEL are unavailable in annotated VCF, uncalibrated LR=1.\n"
    )
    print("AGGREGATE_COMPLETE", len(validation), len(screen), len(phaseall), flush=True)


if __name__ == "__main__":
    main()
