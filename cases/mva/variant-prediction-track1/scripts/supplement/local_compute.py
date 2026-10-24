"""Local, bounded computations for the Track 1 evidence supplement.

Writes only under release/track1-evidence-supplement/computed.
Never opens a patient VCF, FASTQ or BAM. Does not write the frozen release.
"""

from __future__ import annotations

import os

import csv
import gzip
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
PROJECT = DEST.parent.parent
AVI_MOD = PROJECT / "scripts" / "28_alphagenome_bounded.py"
GRAPH_MOD = PROJECT / "scripts" / "33_evidence_graph.py"
FASTA = PROJECT / "ref" / "GRCh38_no_alt.fa"

# GRCh38 primary-assembly centromere spans (UCSC cytoBand, public).
# p arm is [0, centro_start); q arm is [centro_end, chrom_end).
CENTROMERE = {
    "chr1": (121700000, 125100000),
    "chr2": (91800000, 96000000),
    "chr3": (87800000, 94000000),
    "chr4": (48200000, 51800000),
    "chr5": (46100000, 51400000),
    "chr6": (58500000, 62600000),
    "chr7": (58100000, 62100000),
    "chr8": (43200000, 47200000),
    "chr9": (42200000, 45500000),
    "chr10": (38000000, 41600000),
    "chr11": (51000000, 55800000),
    "chr12": (33200000, 37800000),
    "chr13": (16500000, 18900000),
    "chr14": (16100000, 18200000),
    "chr15": (17500000, 20500000),
    "chr16": (35300000, 38400000),
    "chr17": (22700000, 27400000),
    "chr18": (15400000, 21500000),
    "chr19": (24200000, 28100000),
    "chr20": (25700000, 30400000),
    "chr21": (10900000, 13000000),
    "chr22": (13700000, 17400000),
    "chrX": (58100000, 63800000),
    "chrY": (10300000, 13200000),
}

STOP = ("chr15", 40209701, "T", "G")
MISSENSE = ("chr15", 40220612, "T", "G")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def phase_span() -> None:
    start, end = STOP[1], MISSENSE[1]
    distance = end - start
    phase_path = DEST / "compiled" / "phase" / "phase_pair_evidence.tsv"
    whatshap = "file_absent"
    if phase_path.is_file():
        with phase_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        matching = [
            row
            for row in rows
            if {row.get("pos_1"), row.get("pos_2")} == {str(start), str(end)}
            or {row.get("position_1"), row.get("position_2")} == {str(start), str(end)}
        ]
        if not matching:
            # pair table uses chrom:pos keys in some revisions
            matching = [
                row
                for row in rows
                if "40209701" in str(row) and "40220612" in str(row)
            ]
        statuses = sorted({row.get("phase") or row.get("phase_call") or "" for row in matching})
        whatshap = ",".join(s for s in statuses if s) or "no_matching_pair_row"
    write_tsv(
        DEST / "computed" / "phase_span" / "span.tsv",
        [
            "chrom",
            "pos_stop",
            "pos_missense",
            "distance_bp",
            "short_read_spanning_possible",
            "whatshap_status",
            "note",
        ],
        [
            {
                "chrom": "chr15",
                "pos_stop": start,
                "pos_missense": end,
                "distance_bp": distance,
                "short_read_spanning_possible": "no",
                "whatshap_status": whatshap,
                "note": "10911 bp exceeds paired-end insert; no parent or long-read data",
            }
        ],
    )
    print(f"phase_span distance_bp={distance} whatshap={whatshap}")


def coding_avi_splice() -> None:
    avi = load_module(AVI_MOD, "avi_bounded")
    if not avi.AVI.is_file() or not avi.SPLICE.is_file():
        raise FileNotFoundError("AlphaGenome tables missing")
    from pyfaidx import Fasta

    fasta = Fasta(str(FASTA), as_raw=True, sequence_always_upper=True)
    index, _ = avi.read_tbi(avi.AVI)
    size = avi.SPLICE.stat().st_size
    rows = []
    with avi.AVI.open("rb") as bgzf, avi.SPLICE.open("rb") as splice:
        for chrom, pos, ref, alt, label in (
            (*STOP, "BUB1B_p.Leu737Ter"),
            (*MISSENSE, "BUB1B_p.Asn1002Lys"),
        ):
            observed = fasta[chrom][pos - 1 : pos - 1 + len(ref)]
            if observed != ref:
                raise ValueError(f"REF mismatch {chrom}:{pos} {ref} vs {observed}")
            hit, avi_reason, _ = avi.avi_at(index, bgzf, chrom, pos, ref, alt)
            splice_value, splice_reason, _ = avi.splice_at(
                splice, size, chrom, pos, ref, alt
            )
            percentile = ""
            if hit:
                phred = float(hit[5])
                percentile = f"{100 * 10 ** (-phred / 10):.8g}"
            rows.append(
                {
                    "label": label,
                    "grch38_key": f"{chrom}-{pos}-{ref}-{alt}",
                    "avi_raw_score": hit[4] if hit else "",
                    "avi_PHRED": hit[5] if hit else "",
                    "top_percentile": percentile,
                    "avi_status": "exact_allele" if hit else avi_reason,
                    "splice_value": splice_value or "",
                    "splice_status": "exact_allele" if splice_value else splice_reason,
                    "likelihood_ratio": "1",
                    "note": "uncalibrated table lookup; missing is not benign",
                }
            )
    write_tsv(
        DEST / "computed" / "l5_alphagenome" / "coding_snv_scores.tsv",
        list(rows[0]),
        rows,
    )
    print(f"coding_avi_splice rows={len(rows)}")


def ranked_pairs() -> Path:
    src = DEST / "compiled" / "ranking" / "top10_reconstructed.tsv"
    out = DEST / "computed" / "graph" / "ranked_pairs.tsv"
    with src.open(encoding="utf-8", newline="") as handle:
        incoming = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for index, row in enumerate(incoming, start=1):
        v1 = f"{row['chrom_1']}-{row['pos_1']}-{row['ref_1']}-{row['alt_1']}"
        v2 = ""
        if row.get("chrom_2") and row.get("pos_2"):
            v2 = f"{row['chrom_2']}-{row['pos_2']}-{row['ref_2']}-{row['alt_2']}"
        rows.append(
            {
                "pair_id": f"{row['gene']}:{v1};{v2}" if v2 else f"{row['gene']}:{v1}",
                "gene": row["gene"],
                "variant_1": v1,
                "variant_2": v2,
                "epcr": row.get("epcr", ""),
                "epcr_status": "provisional_v1_reconstructed",
                "rank": str(index),
            }
        )
    write_tsv(out, list(rows[0]), rows)
    return out


def evidence_graph() -> None:
    ranked = ranked_pairs()
    graph = load_module(GRAPH_MOD, "evidence_graph")
    sys.argv = [
        str(GRAPH_MOD),
        "--ranked",
        str(ranked),
        "--b",
        str(DEST / "compiled" / "literature" / "angle_b_phenotype.tsv"),
        "--c",
        str(DEST / "compiled" / "literature" / "angle_c_mechanism.tsv"),
        "--d",
        str(DEST / "compiled" / "literature" / "literature_cases.tsv"),
        "--out-dir",
        str(DEST / "computed" / "graph"),
    ]
    graph.main()
    print("evidence_graph written")


def arm_depth() -> None:
    candidates = [
        PROJECT / "mn5_out" / "l4_full_mosdepth_47355682" / "genome_100kb.regions.bed.gz",
        DEST / "computed" / "dosage_baf" / "genome_100kb.regions.bed.gz",
    ]
    bed = next((path for path in candidates if path.is_file()), None)
    note_path = DEST / "computed" / "dosage_baf" / "arm_depth.STATUS.txt"
    if bed is None:
        note_path.write_text(
            "100 kb mosdepth bed is not on this laptop. Fetch the 207 kB "
            "genome_100kb.regions.bed.gz from MN5 "
            + str(Path(os.environ["MVA_BASE"]) / "out/l4_full_qc/mosdepth_47355682/")
            + " into computed/dosage_baf/ and rerun local_compute.py. "
            "Do not fetch the BAM.\n",
            encoding="utf-8",
        )
        print("arm_depth deferred: mosdepth bed absent")
        return
    totals: dict[tuple[str, str], list[float]] = {}
    opener = gzip.open if bed.suffix == ".gz" or bed.name.endswith(".gz") else open
    with opener(bed, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start_s, end_s, *rest = line.split()
            if chrom not in CENTROMERE:
                continue
            start, end = int(start_s), int(end_s)
            depth = float(rest[0]) if rest else 0.0
            centro = CENTROMERE[chrom]
            mid = (start + end) / 2
            arm = "p" if mid < centro[0] else ("q" if mid >= centro[1] else "centro")
            totals.setdefault((chrom, arm), []).append(depth)
    autosomal = [
        sorted(values)[len(values) // 2]
        for (chrom, arm), values in totals.items()
        if chrom.startswith("chr") and chrom[3:].isdigit() and arm in {"p", "q"}
    ]
    auto_median = sorted(autosomal)[len(autosomal) // 2] if autosomal else float("nan")
    rows = []
    for (chrom, arm), values in sorted(totals.items()):
        ordered = sorted(values)
        median = ordered[len(ordered) // 2]
        rows.append(
            {
                "chrom": chrom,
                "arm": arm,
                "bins": len(values),
                "median_depth": f"{median:.4f}",
                "ratio_to_autosome_arm_median": (
                    f"{median / auto_median:.4f}" if auto_median else ""
                ),
                "source": str(bed),
                "note": "depth track; not a mosaic call or clone fraction",
            }
        )
    write_tsv(
        DEST / "computed" / "dosage_baf" / "arm_depth.tsv",
        list(rows[0]),
        rows,
    )
    print(f"arm_depth arms={len(rows)} auto_median={auto_median}")


def filter_exon_bed() -> None:
    src = PROJECT / "results" / "l3_priority_exon_utr.bed"
    dest = DEST / "computed" / "l3_sv" / "tier1_sixgene_exon_utr.bed"
    keep = {"BUB1B", "BUB1", "CEP57", "TRIP13", "CENATAC", "MAD1L1"}
    rows = []
    with src.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            gene = cols[3].split("|", 1)[0]
            if gene in keep:
                rows.append(line if line.endswith("\n") else line + "\n")
    dest.write_text("".join(rows), encoding="utf-8")
    print(f"exon_bed rows={len(rows)}")


def main() -> None:
    phase_span()
    coding_avi_splice()
    evidence_graph()
    arm_depth()
    filter_exon_bed()


if __name__ == "__main__":
    main()
