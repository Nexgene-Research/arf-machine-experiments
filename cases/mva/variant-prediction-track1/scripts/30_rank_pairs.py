"""Freeze genotype/rarity/effect-only Angle A on bounded, normalized GRCh38 tables.

Required TSV schemas (additional columns allowed):
 candidates: gene,chrom,pos,ref,alt,gt,genotype_status,normalization_status,
             ref_matches_grch38
 hypotheses: pair_id,gene,variant_1,variant_2,phase,phase_evidence,
             hypothesis_track (optional: compound_het/homozygous/hemizygous/
             dominant_mosaic/blended)
 read-check: variant,alt_reads,genotype_status (missing rows remain provisional)
 scores: variant,raw_score,phred,splice_value (descriptive only, never LR).
  variant and variant_1/2 are exact chrN-POS-REF-ALT keys. No input posterior:
  Stage 32 computes P_path and ePCR after this independent A freeze.
  score_a is a descriptive read-support ordering, not a probability.
 "--dry-run" still checks evidence but permits writing only to a synthetic_
 file or session-scratch path, and never creates a submission.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys

from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 50_000_000
MAX_ROWS = 100_000
MAX_NEW_SNVS = 100
KEY = re.compile(
    r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))-([1-9][0-9]*)-([ACGT]+)-([ACGT]+)$"
)
BASE = (
    "pair_id",
    "gene",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "gt_1",
    "gt_2",
    "alt_reads_1",
    "alt_reads_2",
    "phase",
    "score_a",
    "angle_a_status",
)
EXTRA = (
    "variant_1",
    "variant_2",
    "p_path_1",
    "p_path_2",
    "read_status_1",
    "read_status_2",
    "normalization_status_1",
    "normalization_status_2",
    "phase_evidence",
    "avi_phred_1",
    "avi_phred_2",
    "splice_value_1",
    "splice_value_2",
    "hypothesis_track",
    "proband_id",
    "evidence_status",
)
FIELDS = BASE + EXTRA
CRITERIA = "exact GRCh38 REF, normalized allele, nonfailing genotype; missing/zero read support retained as provisional; trans validated only with phase evidence; cis downranked; descriptive A score not P_path; no B/C/D; tie by read evidence and lexicographic key"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bounded_table(path, fields):
    if not path.is_file():
        raise FileNotFoundError(f"Required bounded input absent: {path}")
    if path.stat().st_size > MAX_BYTES:
        raise ValueError(f"Genome-scale or oversized input refused (>50 MB): {path}")
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", strict=True)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or any(not s for s in header):
            raise ValueError(f"Duplicate or empty TSV header: {path}")
        missing = sorted(set(fields) - set(header))
        if missing:
            raise ValueError(f"Schema error {path}: missing columns {missing}")
        rows = []
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Schema error {path}: malformed TSV row")
            rows.append(row)
            if len(rows) > MAX_ROWS:
                raise ValueError(f"Genome-scale row count refused: {path}")
    if not rows:
        raise ValueError(f"Empty evidence table: {path}")
    return rows


def key_parts(value):
    m = KEY.fullmatch(value)
    if not m or m.group(3) == m.group(4):
        raise ValueError(f"Invalid exact normalized VCF key: {value!r}")
    return m.group(1), int(m.group(2)), m.group(3), m.group(4)


def allele_key(row):
    return "-".join(row[c] for c in ("chrom", "pos", "ref", "alt"))


def verify_allele(key, fasta):
    chrom, pos, ref, alt = key_parts(key)
    if chrom not in fasta or pos + len(ref) - 1 > len(fasta[chrom]):
        raise ValueError(f"GRCh38 reference interval absent for {key}")
    actual = str(fasta[chrom][pos - 1 : pos - 1 + len(ref)]).upper()
    if actual != ref:
        raise ValueError(f"GRCh38 REF mismatch: {key}; reference={actual}")
    if len(ref) != len(alt):
        if ref[0] != alt[0]:
            raise ValueError(f"VCF padding base missing: {key}")
        # A normalized VCF indel is both minimally represented and left aligned.
        if len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            raise ValueError(f"Nonminimal indel representation: {key}")
        if pos > 1:
            prior = str(fasta[chrom][pos - 2 : pos - 1]).upper()
            if (
                ref[-1] == alt[-1]
                or (len(ref) == 1 and alt[-1] == prior)
                or (len(alt) == 1 and ref[-1] == prior)
            ):
                raise ValueError(f"Indel not left aligned: {key}")
    elif len(ref) > 1:
        raise ValueError(
            f"MNV/complex allele unsupported without decomposed normalized SNV: {key}"
        )
    return chrom, pos, ref, alt


def probability(value, label):
    try:
        n = float(value)
    except ValueError as e:
        raise ValueError(f"Invalid {label}: {value!r}") from e
    if not math.isfinite(n) or not 0 < n <= 1:
        raise ValueError(f"Invalid {label} outside (0,1]: {value!r}")
    return n


def read_score_table(path):
    rows = bounded_table(path, ("variant", "raw_score", "phred", "splice_value"))
    result = {}
    for row in rows:
        key_parts(row["variant"])
        if row["variant"] in result:
            raise ValueError(f"Duplicate score key: {row['variant']}")
        result[row["variant"]] = row
    return result


def incremental_scores(keys, reference, dry_run):
    if not keys:
        return []
    if len(keys) > MAX_NEW_SNVS:
        raise ValueError("New SNV count exceeds bounded 100-allele lookup budget")
    if dry_run:
        raise ValueError(
            "Dry run with new SNVs needs a synthetic complete --scores table; external score files are not queried"
        )
    spec = importlib.util.spec_from_file_location(
        "alphagenome_bounded_28", ROOT / "scripts/28_alphagenome_bounded.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not module.AVI.is_file() or not module.SPLICE.is_file():
        raise FileNotFoundError("New SNVs require both local AlphaGenome score tables")
    index, _ = module.read_tbi(module.AVI)
    size = module.SPLICE.stat().st_size
    output = []
    with module.AVI.open("rb") as avi, module.SPLICE.open("rb") as splice:
        header = splice.readline(module.MAX_LINE + 1).rstrip(b"\r\n")
        if header != b"#CHROM\tPOS\tREF\tALT\talphagenome_splicing":
            raise ValueError("Unexpected AlphaGenome splicing TSV schema")
        for key in sorted(keys):
            chrom, pos, ref, alt = key_parts(key)
            row, reason, _ = module.avi_at(index, avi, chrom, pos, ref, alt)
            splice_value, splice_reason, _ = module.splice_at(
                splice, size, chrom, pos, ref, alt
            )
            output.append(
                dict(
                    variant=key,
                    raw_score=row[4] if row else "",
                    phred=row[5] if row else "",
                    splice_value=splice_value or "",
                    avi_missing_reason=reason,
                    splice_missing_reason=splice_reason,
                    source="28_alphagenome_bounded_exact_lookup",
                    provenance="GRCh38;uncalibrated_descriptive_LR=1",
                )
            )
    return output


def render(rows, fields):
    import io

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--candidates", type=Path, default=ROOT / "results/candidates_full.tsv"
    )
    p.add_argument(
        "--pairs", type=Path, default=ROOT / "results/compound_het_candidates.tsv"
    )
    p.add_argument("--read-check", type=Path, default=ROOT / "results/read_check.tsv")
    p.add_argument(
        "--scores", type=Path, default=ROOT / "results/alphagenome_scores.tsv"
    )
    p.add_argument("--reference", type=Path, default=ROOT / "ref/GRCh38_no_alt.fa")
    p.add_argument("--out", type=Path, default=ROOT / "results/ranked_a.tsv")
    p.add_argument(
        "--manifest", type=Path, default=ROOT / "results/ranked_a_manifest.json"
    )
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    incremental = a.out.with_name("alphagenome_scores_incremental.tsv")
    destinations = [a.out, a.manifest, incremental]
    if len(set(destinations)) != len(destinations) or any(
        path.exists() for path in destinations
    ):
        raise FileExistsError(
            "Frozen Angle A destination already exists; never overwrite"
        )
    if a.dry_run and any(
        "synthetic_" not in path.name and "ses_" not in str(path)
        for path in (a.out, a.manifest)
    ):
        raise ValueError(
            "Dry run requires synthetic_ output names or session scratch destinations"
        )
    candidates = bounded_table(
        a.candidates,
        (
            "gene",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gt",
            "genotype_status",
            "normalization_status",
            "ref_matches_grch38",
        ),
    )
    pairs = bounded_table(
        a.pairs,
        ("pair_id", "gene", "variant_1", "variant_2", "phase", "phase_evidence"),
    )
    reads = bounded_table(a.read_check, ("variant", "alt_reads", "genotype_status"))
    score = read_score_table(a.scores)
    if not a.dry_run and not Path(str(a.reference) + ".fai").is_file():
        raise FileNotFoundError(
            "Reference .fai missing; refusing unbounded whole-FASTA index build on laptop"
        )
    fasta = Fasta(str(a.reference), as_raw=True, sequence_always_upper=True)
    by_key = {}
    for c in candidates:
        key = allele_key(c)
        verify_allele(key, fasta)
        if key in by_key:
            raise ValueError(f"Duplicate candidate: {key}")
        if (
            c["normalization_status"] != "bcftools_norm_GRCh38"
            or c["ref_matches_grch38"] != "true"
        ):
            raise ValueError(
                f"Missing exact bcftools GRCh38 normalization/REF proof: {key}"
            )
        by_key[key] = c
    read_by_key = {}
    for r in reads:
        key_parts(r["variant"])
        if r["variant"] in read_by_key:
            raise ValueError(f"Duplicate read-check key: {r['variant']}")
        read_by_key[r["variant"]] = r
    seen, ranked, new_snvs = set(), [], set()
    for pair in pairs:
        pair_id = pair["pair_id"].strip()
        if not pair_id or pair_id in seen:
            raise ValueError(f"Empty/duplicate pair_id: {pair_id!r}")
        seen.add(pair_id)
        k1, k2 = pair["variant_1"], pair["variant_2"]
        if k1 not in by_key or (k2 and k2 not in by_key) or (k2 and k1 == k2):
            raise ValueError(
                f"Hypothesis references absent/identical candidate alleles: {pair_id}"
            )
        c1, c2 = by_key[k1], by_key[k2] if k2 else None
        track = pair.get("hypothesis_track", "compound_het" if k2 else "homozygous")
        if track not in (
            "compound_het",
            "homozygous",
            "hemizygous",
            "dominant_mosaic",
            "blended",
        ):
            raise ValueError(f"Invalid hypothesis track {track}: {pair_id}")
        if (track == "compound_het") != bool(k2) and track != "blended":
            if not (
                not k2 and track in ("homozygous", "hemizygous", "dominant_mosaic")
            ):
                raise ValueError(f"Track/allele count conflict: {pair_id}")
        if (
            c1["gene"] != pair["gene"]
            or (c2 is not None and c2["gene"] != pair["gene"] and track != "blended")
            or not pair["gene"]
        ):
            raise ValueError(f"Pair gene mismatch: {pair_id}")
        if pair["phase"] not in ("cis", "trans", "unknown") or (
            pair["phase"] in ("cis", "trans") and not pair["phase_evidence"].strip()
        ):
            raise ValueError(f"Missing explicit phase status/evidence: {pair_id}")
        evidence = []
        for key, c in ((k1, c1), *([(k2, c2)] if c2 is not None else [])):
            r = read_by_key.get(key)
            if r is None:
                evidence.append(({"genotype_status": "missing"}, 0))
                continue
            try:
                alt_reads = int(r["alt_reads"]) if r["alt_reads"].strip() else 0
            except ValueError as e:
                raise ValueError(f"Invalid ALT reads: {key}") from e
            if alt_reads < 0:
                raise ValueError(f"Negative ALT reads: {key}")
            evidence.append((r, alt_reads))
        if any(
            c["genotype_status"].lower()
            in ("fail", "failed", "invalid", "ref_mismatch")
            or r["genotype_status"].lower() in ("fail", "failed", "invalid")
            for c, (r, n) in zip((c1, *([c2] if c2 else [])), evidence)
        ):
            continue  # explicit failing genotype, not merely missing read evidence
        if track == "compound_het" and any(
            c["gt"].replace("|", "/") not in ("0/1", "1/0") for c in (c1, c2)
        ):
            continue
        if track == "homozygous" and c1["gt"].replace("|", "/") != "1/1":
            continue
        for key in (k1, *([k2] if k2 else [])):
            chrom, pos, ref, alt = key_parts(key)
            if len(ref) == len(alt) == 1 and key not in score:
                new_snvs.add(key)
        confident = all(
            c["genotype_status"] == "validated"
            and r["genotype_status"] == "validated"
            and n > 0
            for c, (r, n) in zip((c1, *([c2] if c2 else [])), evidence)
        )
        value = (1 if confident else 0.1) * {"trans": 1, "unknown": 0.5, "cis": 0.01}[
            pair["phase"]
        ]
        row = {
            "pair_id": pair_id,
            "gene": pair["gene"],
            "phase": pair["phase"],
            "score_a": format(value, ".12g"),
            "angle_a_status": "verified"
            if confident
            else "provisional_missing_read_or_genotype",
            "hypothesis_track": track,
            "proband_id": pair.get("proband_id", ""),
            "evidence_status": "read_supported"
            if confident
            else "unverified_critical_read_link",
        }
        for i, (key, c, read) in enumerate(
            ((k1, c1, evidence[0]), *([(k2, c2, evidence[1])] if c2 else [])), 1
        ):
            chrom, pos, ref, alt = key_parts(key)
            row.update(
                {
                    f"{n}_{i}": v
                    for n, v in (
                        ("chrom", chrom),
                        ("pos", str(pos)),
                        ("ref", ref),
                        ("alt", alt),
                        ("gt", c["gt"]),
                        ("alt_reads", str(read[1])),
                        ("variant", key),
                        ("p_path", ""),
                        ("read_status", read[0]["genotype_status"]),
                        ("normalization_status", c["normalization_status"]),
                    )
                }
            )
        row["phase_evidence"] = pair["phase_evidence"]
        if not c2:
            for col in (
                "chrom",
                "pos",
                "ref",
                "alt",
                "gt",
                "alt_reads",
                "variant",
                "p_path",
                "read_status",
                "normalization_status",
            ):
                row[f"{col}_2"] = ""
        ranked.append(row)
    if not ranked:
        raise ValueError(
            "No eligible pair with two validated genotypes and positive ALT reads; no frozen A"
        )
    lookup = incremental_scores(new_snvs, a.reference, a.dry_run)
    score.update({r["variant"]: r for r in lookup})
    for row in ranked:
        for i in (1, 2):
            entry = score.get(row[f"variant_{i}"], {})
            row[f"avi_phred_{i}"] = entry.get("phred", "")
            row[f"splice_value_{i}"] = entry.get("splice_value", "")
    ranked.sort(
        key=lambda r: (
            -float(r["score_a"]),
            -(int(r["alt_reads_1"]) + int(r["alt_reads_2"])),
            min(r["variant_1"], r["variant_2"]),
            max(r["variant_1"], r["variant_2"]),
            r["pair_id"],
        )
    )
    data = render(ranked, FIELDS)
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "criteria": CRITERIA,
        "inputs_sha256": {
            str(path.resolve()): digest(path)
            for path in (a.candidates, a.pairs, a.read_check, a.scores)
        },
        "reference_path": str(a.reference.resolve()),
        "reference_index_sha256": digest(Path(str(a.reference) + ".fai"))
        if Path(str(a.reference) + ".fai").is_file()
        else "not_available_synthetic_reference",
        "ranked_a_sha256": hashlib.sha256(data).hexdigest(),
        "ranked_a_path": str(a.out.resolve()),
        "incremental_score_keys": sorted(new_snvs),
        "eligible_pairs": len(ranked),
        "excluded_pairs": len(pairs) - len(ranked),
        "no_B_C_D_consumed": True,
    }
    for path in (a.out, a.manifest):
        if not path.parent.is_dir():
            raise FileNotFoundError(f"Output directory absent: {path.parent}")
    if lookup:
        incremental.write_bytes(
            render(
                lookup,
                (
                    "variant",
                    "raw_score",
                    "phred",
                    "splice_value",
                    "avi_missing_reason",
                    "splice_missing_reason",
                    "source",
                    "provenance",
                ),
            )
        )
    a.out.write_bytes(data)
    a.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Frozen A: {len(ranked)} eligible / {len(pairs) - len(ranked)} excluded; sha256={manifest['ranked_a_sha256']}; incremental_SNV={len(lookup)}"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"ERROR: {error}")
