"""Bounded, read-only GRCh38 AlphaGenome SNV lookups; never open patient VCF.

Usage: env/venv/Scripts/python.exe scripts/28_alphagenome_bounded.py
       env/venv/Scripts/python.exe scripts/28_alphagenome_bounded.py --self-test
Only the small tier1_bounded_calls.tsv is an input. No panel is used as a
surrogate for the full table, and the plain splicing TSV's unrelated .tbi is
never opened. Missing predictions are not benign evidence.
"""

import argparse
import csv
import gzip
import hashlib
import io
import struct
import tempfile
import zlib
from pathlib import Path

from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
CALLS = ROOT / "results/tier1_bounded_calls.tsv"
AVI = Path(
    r"F:\data\avi_scores_snvs_tabix\alphagenome_variant_impact_score_snvs.tsv.gz"
)
SPLICE = Path(
    r"F:\data\combined_splicing_snvs_tabix\combined_alphagenome_splicing_snvs.tsv"
)
OUT = ROOT / "results/alphagenome_scores.tsv"
NOTES = ROOT / "results/alphagenome_lookup_notes.md"
MAX_CALLS = 100
MAX_LINE = 1024
MAX_SEEKS = 50


def order_key(chrom, pos):
    c = chrom.removeprefix("chr")
    rank = (
        int(c)
        if c.isdecimal() and 1 <= int(c) <= 22
        else {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c)
    )
    if rank is None:
        raise ValueError(f"Unrecognized GRCh38 contig {chrom!r}")
    return (rank, int(pos))


def read_tbi(path):
    # Bound index inflation independent of the size of the BGZF data file.
    with gzip.open(str(path) + ".tbi", "rb") as handle:
        raw = handle.read(32 * 1024 * 1024 + 1)
    if len(raw) > 32 * 1024 * 1024 or raw[:4] != b"TBI\x01":
        raise ValueError("Invalid or oversized AVI TBI")
    buf = io.BytesIO(raw)
    buf.seek(4)
    nref, fmt, col_seq, col_beg, col_end, meta, skip, nbytes = struct.unpack(
        "<8i", buf.read(32)
    )
    if not (
        0 < nref < 1000
        and 0 < nbytes < 100000
        and fmt & 0xFFFF in (0, 2)  # generic or VCF preset for TSV-like fields
        and col_seq == 1
        and col_beg == 2
    ):
        raise ValueError("Unexpected AVI tabix columns / index")
    names = buf.read(nbytes).rstrip(b"\0").decode("ascii").split("\0")
    if len(names) != nref:
        raise ValueError("TBI contig count mismatch")
    indexes = {}
    for name in names:
        (nbin,) = struct.unpack("<i", buf.read(4))
        bins = {}
        if not 0 <= nbin <= 100000:
            raise ValueError("Invalid TBI bin count")
        for _ in range(nbin):
            bin_id, nchunk = struct.unpack("<Ii", buf.read(8))
            if not 0 <= nchunk <= 100000:
                raise ValueError("Invalid TBI chunk count")
            chunks = [struct.unpack("<QQ", buf.read(16)) for _ in range(nchunk)]
            if bin_id != 37450:  # index statistics pseudo-bin is not genomic data
                bins[bin_id] = chunks
        (ninterval,) = struct.unpack("<i", buf.read(4))
        if not 0 <= ninterval <= 2000000:
            raise ValueError("Invalid TBI interval count")
        buf.seek(8 * ninterval, 1)
        indexes[name] = bins
    return indexes, (fmt, col_seq, col_beg, col_end, meta, skip)


def bgzf_block(handle, addr):
    handle.seek(addr)
    head = handle.read(18)
    if len(head) != 18 or head[:4] != b"\x1f\x8b\x08\x04" or head[12:14] != b"BC":
        raise ValueError("Invalid BGZF header")
    size = struct.unpack_from("<H", head, 16)[0] + 1
    if size < 26 or size > 65536:
        raise ValueError("Invalid BGZF block size")
    packed = head + handle.read(size - 18)
    if len(packed) != size:
        raise ValueError("Truncated BGZF block")
    data = zlib.decompress(packed, wbits=31)
    if len(data) > 65536 or len(data) != struct.unpack_from("<I", packed, size - 4)[0]:
        raise ValueError("BGZF uncompressed size mismatch")
    return data, size


def avi_at(index, handle, chrom, pos, ref, alt):
    name = chrom if chrom in index else chrom.removeprefix("chr")
    if name not in index:
        return None, "contig_not_indexed", 0
    p = pos - 1
    bins = [
        0,
        1 + (p >> 26),
        9 + (p >> 23),
        73 + (p >> 20),
        585 + (p >> 17),
        4681 + (p >> 14),
    ]
    chunks = [chunk for bin_id in bins for chunk in index[name].get(bin_id, ())]
    if not chunks:
        return None, "position_not_indexed", 0
    lo, hi = min(a for a, b in chunks), max(b for a, b in chunks)
    # Linear-offset trimming for coarse parent bins, if needed, can be added
    # later; hard read cap prevents a stray wide bin from scanning the table.
    addr, pending, count = lo >> 16, b"", 0
    while addr <= (hi >> 16) + 1:
        data, size = bgzf_block(handle, addr)
        count += 1
        if count > 96:
            raise ValueError("AVI region exceeds 96 BGZF blocks; refusing broad scan")
        if count == 1:
            within = lo & 0xFFFF
            if within > len(data) or (within and data[within - 1 : within] != b"\n"):
                raise ValueError("TBI virtual offset not on BGZF line boundary")
            data = data[within:]
        parts = (pending + data).split(b"\n")
        pending = parts.pop()
        for line in parts:
            if not line or line.startswith(b"#"):
                continue
            fields = line.rstrip(b"\r").decode("ascii").split("\t")
            if len(fields) != 6:
                raise ValueError("Invalid AVI record")
            row_key = (fields[0], int(fields[1]))
            if row_key == (name, pos) and fields[2:4] == [ref, alt]:
                float(fields[4])
                float(fields[5])
                return fields, "", count
            if fields[0] == name and int(fields[1]) > pos:
                return None, "allele_not_in_avi", count
        addr += size
    return None, "allele_not_in_avi", count


def aligned_line(handle, offset, size, stats):
    """Seek to a byte midpoint; return the first *complete* following row."""
    stats["seeks"] += 1
    if stats["seeks"] > MAX_SEEKS:
        raise ValueError("Plain-TSV per-position seek budget exceeded")
    handle.seek(offset)
    if offset:
        partial = handle.readline(MAX_LINE + 1)
        if len(partial) > MAX_LINE:
            raise ValueError("Oversized partial TSV line")
    start = handle.tell()
    if start >= size:
        return size, size, None
    line = handle.readline(MAX_LINE + 1)
    if len(line) > MAX_LINE or not line.endswith(b"\n"):
        raise ValueError("Oversized or incomplete TSV row")
    stats["bytes"] += len(line) + (start - offset)
    fields = line.rstrip(b"\r\n").decode("ascii").split("\t")
    if len(fields) != 5:
        raise ValueError("Invalid splicing TSV row")
    return start, handle.tell(), fields


def splice_at(handle, size, chrom, pos, ref, alt):
    """Lower bound by physical byte offsets, then inspect *all* same-position alleles."""
    target = order_key(chrom, pos)
    lo, hi = 0, size
    stats = {"seeks": 0, "bytes": 0, "same_position": 0}
    while hi - lo > MAX_LINE:
        mid = (lo + hi) // 2
        start, end, fields = aligned_line(handle, mid, size, stats)
        if (
            fields is None
            or fields[0].startswith("#")
            or order_key(fields[0], fields[1]) >= target
        ):
            hi = mid if fields is None else start
        else:
            lo = end
        if lo > hi:
            raise ValueError("Splicing table not sorted or binary interval inverted")
    stats["seeks"] += 1
    handle.seek(lo)
    # lo always denotes a physical line start (including zero). At most one
    # short final interval is read, plus the 4 alleles at a SNV position.
    while True:
        line = handle.readline(MAX_LINE + 1)
        stats["bytes"] += len(line)
        if not line:
            return None, "position_not_in_splicing", stats
        if len(line) > MAX_LINE:
            raise ValueError("Oversized splicing TSV row")
        fields = line.rstrip(b"\r\n").decode("ascii").split("\t")
        if fields[0].startswith("#"):
            continue
        row_key = order_key(fields[0], fields[1])
        if row_key > target:
            reason = (
                "allele_not_in_splicing"
                if stats["same_position"]
                else "position_not_in_splicing"
            )
            return None, reason, stats
        if row_key == target:
            stats["same_position"] += 1
            if fields[2:4] == [ref, alt]:
                if len(fields) != 5:
                    raise ValueError("Malformed splicing score")
                float(fields[4])
                return fields[4], "", stats
        if stats["bytes"] > 128 * 1024 or stats["same_position"] > 16:
            raise ValueError("Splicing lookup exceeded bounded local scan")


def test_search():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "mini.tsv"
        lines = [b"#CHROM\tPOS\tREF\tALT\talphagenome_splicing\n"]
        rows = []
        for chrom in ("chr1", "chr2", "chr10", "chr15", "chrX"):
            for pos in range(1, 400):
                for ref, alt in (("A", "C"), ("A", "G"), ("A", "T")):
                    value = f"{pos / 1000:.4f}"
                    lines.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{value}\n".encode())
                    rows.append((chrom, pos, ref, alt, value))
        path.write_bytes(b"".join(lines))
        with path.open("rb") as f:
            for chrom, pos, ref, alt, value in (
                rows[0],
                rows[-1],
                rows[1700],
                rows[2401],
                rows[4700],
            ):
                found, reason, stats = splice_at(
                    f, path.stat().st_size, chrom, pos, ref, alt
                )
                assert (found, reason) == (value, ""), (chrom, pos, found, reason)
                assert stats["seeks"] <= MAX_SEEKS
            assert (
                splice_at(f, path.stat().st_size, "chr15", 200, "A", "G")[0] == "0.2000"
            )
            assert (
                splice_at(f, path.stat().st_size, "chr2", 240, "C", "G")[1]
                == "allele_not_in_splicing"
            )
            assert (
                splice_at(f, path.stat().st_size, "chr15", 40168264, "G", "A")[1]
                == "position_not_in_splicing"
            )

    # Independently construct two small BGZF members; a variant line starts
    # in one member and terminates in the next. The query must join both.
    def pack_bgzf(plain):
        compressor = zlib.compressobj(level=6, wbits=-15)
        payload = compressor.compress(plain) + compressor.flush()
        size = 18 + len(payload) + 8
        return (
            b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
            + struct.pack("<H", size - 1)
            + payload
            + struct.pack("<II", zlib.crc32(plain), len(plain))
        )

    first = b"chr15\t40168263\tA\tC\t0.1\t1\nchr15\t40168264\tG"
    second = b"\tA\t0.06781\t6.91273\nchr15\t40168264\tG\tC\t0.1\t1\nchr15\t40168265\tG\tC\t0.2\t2\n"
    packed_first = pack_bgzf(first)
    key_bin = 4681 + ((40168264 - 1) >> 14)
    synthetic_index = {
        "chr15": {key_bin: [(0, (len(packed_first) << 16) | len(second))]}
    }
    with io.BytesIO(packed_first + pack_bgzf(second)) as handle:
        hit, reason, blocks = avi_at(
            synthetic_index, handle, "chr15", 40168264, "G", "A"
        )
        assert hit is not None and (hit[4:], reason, blocks) == (
            ["0.06781", "6.91273"],
            "",
            2,
        )
        hit, reason, blocks = avi_at(
            synthetic_index, handle, "chr15", 40168264, "G", "T"
        )
        assert hit is None and reason == "allele_not_in_avi"
    return "synthetic checks: TSV first/last/middle, natural contig order, duplicate alleles, missing allele/position; BGZF cross-block line and absent allele PASS"


def run():
    test_result = test_search()
    with CALLS.open(encoding="utf-8", newline="") as f:
        calls = list(csv.DictReader(f, delimiter="\t"))
    if not 0 < len(calls) <= MAX_CALLS:
        raise ValueError("Unexpected bounded input size")
    if len({(r["chr"], r["pos"], r["ref"], r["alt"]) for r in calls}) != len(calls):
        raise ValueError("Repeated input variant")
    if not AVI.is_file() or not SPLICE.is_file():
        raise FileNotFoundError("Required AlphaGenome data missing")
    index, tabix_cols = read_tbi(AVI)
    if not ({"chr15", "15"} & set(index)):
        raise ValueError("AVI lacks chr15 indexed contig")
    fasta = Fasta(
        str(ROOT / "ref/GRCh38_no_alt.fa"), as_raw=True, sequence_always_upper=True
    )
    # Verify natural chromosome ordering at fixed small offsets, not a scan.
    sample_keys = []
    size = SPLICE.stat().st_size
    with SPLICE.open("rb") as splice:
        header = splice.readline(MAX_LINE + 1).rstrip(b"\r\n")
        if header != b"#CHROM\tPOS\tREF\tALT\talphagenome_splicing":
            raise ValueError("Unexpected plain TSV header")
        for i in range(16):
            start, end, fields = aligned_line(
                splice, size * i // 16, size, {"seeks": 0, "bytes": 0}
            )
            if fields is not None and not fields[0].startswith("#"):
                sample_keys.append(order_key(fields[0], fields[1]))
        if sample_keys != sorted(sample_keys):
            raise ValueError(
                "Fixed-offset probes contradict natural chromosome-position order"
            )
        rows, max_seeks, splice_bytes, avi_blocks = [], 0, 0, 0
        with AVI.open("rb") as avi:
            for call in calls:
                chrom, pos, ref, alt = (
                    call["chr"],
                    int(call["pos"]),
                    call["ref"],
                    call["alt"],
                )
                if (
                    not chrom.startswith("chr")
                    or pos < 1
                    or fasta[chrom][pos - 1 : pos - 1 + len(ref)] != ref
                ):
                    raise ValueError(
                        f"Bad contig/position/GRCh38 REF for {chrom}:{pos}:{ref}:{alt}"
                    )
                snv = (
                    len(ref) == len(alt) == 1
                    and ref != alt
                    and ref in "ACGT"
                    and alt in "ACGT"
                )
                row = dict(
                    variant=f"{chrom}-{pos}-{ref}-{alt}",
                    source="AVI_SNV_BGZF_TBI;splicing_plain_TSV_binary_seek",
                    raw_score="",
                    phred="",
                    top_percentile="",
                    splice_value="",
                    avi_missing_reason="",
                    splice_missing_reason="",
                    provenance="GRCh38;exact_REF_ALT;uncalibrated_descriptive_LR=1",
                )
                if snv:
                    hit, reason, blocks = avi_at(index, avi, chrom, pos, ref, alt)
                    avi_blocks += blocks
                    row["avi_missing_reason"] = reason
                    if hit:
                        raw, phred = float(hit[4]), float(hit[5])
                        row.update(
                            raw_score=hit[4],
                            phred=hit[5],
                            top_percentile=f"{100 * 10 ** (-phred / 10):.9g}",
                        )
                    splice_value, reason, stats = splice_at(
                        splice, size, chrom, pos, ref, alt
                    )
                    max_seeks = max(max_seeks, stats["seeks"])
                    splice_bytes += stats["bytes"]
                    row["splice_missing_reason"] = reason
                    if splice_value is not None:
                        row["splice_value"] = splice_value
                else:
                    row["avi_missing_reason"] = "not_SNV_or_monomorphic;AVI_SNV_only"
                    row["splice_missing_reason"] = (
                        "not_SNV_or_monomorphic;splicing_SNV_table_only"
                    )
                rows.append(row)
    # Write output only after all checks and bounded lookups complete.
    with OUT.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    digest = hashlib.sha256(OUT.read_bytes()).hexdigest()
    avi_hits = sum(bool(r["raw_score"]) for r in rows)
    splice_hits = sum(bool(r["splice_value"]) for r in rows)
    NOTES.write_text(
        "# Bounded AlphaGenome lookups (2026-10-16)\n\n"
        f"Run: `env\\venv\\Scripts\\python.exe scripts/28_alphagenome_bounded.py`. Input: `results/tier1_bounded_calls.tsv` ({len(rows)} exact GRCh38 alleles); no patient VCF or genome-wide table scan. "
        f"Only SNVs queried: {sum(len(r['ref']) == len(r['alt']) == 1 for r in calls)}; indels await MN5 bcftools normalization/inference.\n\n"
        f"AVI exact-allele hits: {avi_hits}; splicing exact-allele hits: {splice_hits}; maximal per-SNV plain-TSV seeks: {max_seeks}/{MAX_SEEKS}; "
        f"total plain-TSV bytes read in queries: {splice_bytes}; BGZF blocks: {avi_blocks}; fixed-offset order probes: {len(sample_keys)}. "
        f"Tabix index columns (format,seq,beg,end,meta,skip): {tabix_cols}. "
        "Splicing `.tsv.gz.tbi` was never opened because it does not index the plain `.tsv`; panel subset was not used. "
        "Natural chr1..chr22,X,Y order established only at fixed offset probes, not proved globally.\n\n"
        f"Verification: {test_result}. Every REF matched local GRCh38 FASTA. AVI top percentile = 100 × 10^(-PHRED/10); "
        "it is **not** the splicing percentile. Scores are uncalibrated, descriptive only, LR=1; missing values cannot count as benign. "
        "No phenotype, phase, or causality is inferred from these values.\n\n"
        "Queue suggestion (do not edit shared `results/alphagenome_queue.tsv`): consider exact SNV records with absent AVI/splicing separately for local AlphaGenome inference only after confirming the model outputs match these score definitions; "
        "all five existing indels remain queued, with repeat alleles normalized/deduplicated on MN5 first.\n\n"
        f"Sources: `{AVI}` (BGZF indexed, file size {AVI.stat().st_size} bytes); `{AVI}.tbi` (index, size {Path(str(AVI) + '.tbi').stat().st_size} bytes); "
        f"`{SPLICE}` (plain ordered TSV, size {size} bytes); local FASTA `ref/GRCh38_no_alt.fa`. "
        f"SHA256 of final scores TSV: `{digest}`. Full source-table hashes not computed (bounded-I/O requirement).\n",
        encoding="utf-8",
    )
    print(
        f"calls={len(rows)} avi_hits={avi_hits} splice_hits={splice_hits} max_seeks={max_seeks} splice_bytes={splice_bytes} avi_blocks={avi_blocks} scores_sha256={digest}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Only exercise synthetic bounded plain-TSV search",
    )
    args = parser.parse_args()
    if args.self_test:
        print(test_search())
    else:
        run()
