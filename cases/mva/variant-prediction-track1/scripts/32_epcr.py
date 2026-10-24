"""Preregistered variant evidence arithmetic and *gated* case-level ePCR.

Input contract: ranked A TSV has pair_id,gene,chrom_1,pos_1,ref_1,alt_1,
chrom_2,pos_2,ref_2,alt_2,gt_1,gt_2,alt_reads_1,alt_reads_2,
phase,score_a,angle_a_status. Empty allele 2 means a single allele.
Never interpret score_a or B/C/D descriptive points as likelihood ratios.

Criterion TSV: variant_key,criterion,strength,source_id,audit_status,
computational_family (optional). Unique (variant_key,criterion) is mandatory.
Allowed strength labels are very_strong, strong, moderate, supporting,
benign_strong, benign_supporting. BA1 is not assigned a generic LR;
PP3/BP4 require two-star gene-heldout predictor-specific calibration.
Case links TSV: pair_id,link,q,source_id,calibration_status (heldout_real or
synthetic), exactly phenotype,mechanism,read,inheritance per pair. These q's
MUST originate in independently measured case-level calibration; never enter
descriptive scores here. Synthetic mode never writes to clinical results.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = (
    "pair_id gene chrom_1 pos_1 ref_1 alt_1 chrom_2 pos_2 ref_2 alt_2 "
    "gt_1 gt_2 alt_reads_1 alt_reads_2 phase score_a angle_a_status"
).split()
LINKS = ("phenotype", "mechanism", "read", "inheritance")
LR = {
    "very_strong": 350.0,
    "strong": 18.7,
    "moderate": 4.3,
    "supporting": 2.08,
    "benign_strong": 1 / 18.7,
    "benign_supporting": 1 / 2.08,
}
TIER1 = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1"}


def table(path, required=(), maximum=100000):
    if Path(path).stat().st_size > 50_000_000:
        raise ValueError(f"{path}: bounded table exceeds 50 MB; stage on MN5")
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(
                f"{path}: missing columns {sorted(set(required) - set(reader.fieldnames or []))}"
            )
        out = []
        for row in reader:
            if len(out) >= maximum:
                raise ValueError(
                    f"{path}: table exceeds bounded {maximum} rows; run large passes on MN5"
                )
            out.append(row)
    return out, reader.fieldnames


def write(path, columns, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def allele_key(row, idx):
    c, p, ref, alt = (
        row[f"{part}_{idx}"].strip() for part in ("chrom", "pos", "ref", "alt")
    )
    if not any((c, p, ref, alt)):
        return ""
    if (
        not all((c, p, ref, alt))
        or not c.startswith("chr")
        or not p.isdigit()
        or int(p) < 1
    ):
        raise ValueError(f"{row['pair_id']}: incomplete/non-GRCh38-style allele {idx}")
    if any(x not in "ACGT" for x in ref + alt) or ref == alt:
        raise ValueError(f"{row['pair_id']}: invalid REF/ALT at allele {idx}")
    return f"{c}-{int(p)}-{ref}-{alt}"


def posterior(evidence, prior, multiplier=1.0, source_counts=None):
    odds = prior / (1 - prior)
    for item in evidence:
        odds *= float(item.get("empirical_lr") or LR[item["strength"]]) ** (
            multiplier
            * (
                source_counts.get(item["source_id"], 0)
                if source_counts is not None
                else 1
            )
        )
    return odds / (1 + odds)


def phase_probability(row, unknown=0.5):
    phase = row["phase"].strip().lower()
    if phase in {"trans", "cis"}:
        proof = row.get("phase_evidence", "").lower()
        if not any(
            x in proof
            for x in (
                "physical",
                "familial",
                "family",
                "verified",
                "validated",
                "read-backed",
                "readbacked",
            )
        ):
            raise ValueError(
                f"{row['pair_id']}: phase={phase} needs documented physical/familial validation; no unproved trans inference"
            )
        return 1.0 if phase == "trans" else 0.0
    if phase in {"trans_verified", "trans_physical", "trans_familial"}:
        return 1.0
    if phase in {"cis_verified", "cis_physical", "cis_familial"}:
        return 0.0
    if phase in {"unknown", "unphased", ""}:
        return unknown
    raise ValueError(
        f"{row['pair_id']}: unrecognized phase {phase!r}; use verified trans/cis or unknown"
    )


def pick_angles(row, b, c, d):
    # Independent A plus nonempty independent gene-level B/C/D; not their scores.
    gene = row["gene"]
    count = int(
        row["angle_a_status"].strip().lower()
        in {
            "verified",
            "supported",
            "frozen_a_validated_genotypes_descriptive_predictions_only",
        }
    )
    count += int(gene in b) + int(gene in c) + int(gene in d)
    return count


def read_angles(args):
    result = []
    for path in (args.angle_b, args.angle_c, args.angle_d):
        if path is None or not Path(path).exists():
            result.append(set())
            continue
        rows, _ = table(path, ("gene",))
        result.append({r["gene"].strip() for r in rows if r["gene"].strip()})
    return result


def empirical_computational_calibration(path, keys):
    """Learn data thresholds on calibration genes; assess LR on held-out genes.

    No Pejaver-published numeric threshold is presumed or substituted. A
    predictor/score direction must be supplied and constant per stratum.
    This is an optional empirical mapping for a two-star nonconflicting
    ClinVar P/LP vs B/LB reference, NEVER a patient-specific calibration.
    """
    if path is None or not Path(path).exists():
        return {}
    fields = (
        "variant_key",
        "gene",
        "label",
        "stars",
        "conflict",
        "split",
        "predictor",
        "score",
        "orientation",
        "consequence",
        "af_bin",
    )
    rows, _ = table(path, fields)
    per_group = defaultdict(lambda: defaultdict(list))
    seen = set()
    for r in rows:
        if (
            r["variant_key"] in keys
            or r["variant_key"] in seen
            or r["gene"] in TIER1
            or r["conflict"] != "no"
            or r["label"] not in ("P", "LP", "B", "LB")
            or int(r["stars"]) < 2
        ):
            if r["variant_key"] in seen:
                raise ValueError(
                    "duplicate ClinVar allele in computational calibration"
                )
            continue
        seen.add(r["variant_key"])
        if r["split"] not in ("calibration", "heldout") or r["orientation"] not in (
            "high",
            "low",
        ):
            raise ValueError(
                "calibration needs explicit split and validated score orientation"
            )
        value = float(r["score"])
        if not math.isfinite(value) or not all(
            r[k] for k in ("gene", "predictor", "consequence", "af_bin")
        ):
            raise ValueError(
                "missing calibrated predictor score, consequence or AF bin"
            )
        group = (r["predictor"], r["consequence"], r["af_bin"], r["orientation"])
        per_group[group][r["split"]].append(
            (r["gene"], r["label"] in ("P", "LP"), value)
        )
    result = {}
    for group, splits in per_group.items():
        train, held = splits["calibration"], splits["heldout"]
        train_genes = {x[0] for x in train}
        held_genes = {x[0] for x in held}
        if train_genes & held_genes:
            raise ValueError(
                f"gene leakage between calibration and heldout for {group}"
            )
        if any(
            sum(x[1] == cls for x in train) < 10
            or sum(x[1] == cls for x in held) < 5
            or len({x[0] for x in held if x[1] == cls}) < 5
            for cls in (True, False)
        ):
            continue
        direction = 1 if group[3] == "high" else -1
        thresholds = sorted({direction * x[2] for x in train})
        if len(thresholds) > 256:
            thresholds = [
                thresholds[i * (len(thresholds) - 1) // 255] for i in range(256)
            ]
        # Threshold selection uses calibration genes only; deterministic ties.
        n_p = sum(x[1] for x in train)
        n_b = len(train) - n_p
        threshold = max(
            thresholds,
            key=lambda t: (
                sum(x[1] and direction * x[2] >= t for x in train) / n_p
                - sum(not x[1] and direction * x[2] >= t for x in train) / n_b,
                -t,
            ),
        )
        counts = {
            (label, hit): sum(
                x[1] == label and (direction * x[2] >= threshold) == hit for x in held
            )
            for label in (True, False)
            for hit in (True, False)
        }
        n_hp = sum(x[1] for x in held)
        n_hb = len(held) - n_hp
        lr_high = ((counts[True, True] + 1) / (n_hp + 2)) / (
            (counts[False, True] + 1) / (n_hb + 2)
        )
        lr_low = ((counts[True, False] + 1) / (n_hp + 2)) / (
            (counts[False, False] + 1) / (n_hb + 2)
        )
        result[group] = (threshold, lr_high, lr_low, len(train), len(held))
    return result


def load_evidence(path, keys, calibration=None):
    evidence = defaultdict(list)
    if path is None or not Path(path).exists():
        return evidence
    rows, _ = table(
        path, ("variant_key", "criterion", "strength", "source_id", "audit_status")
    )
    seen = set()
    calibration = calibration or {}
    used_sources = defaultdict(set)
    for row in rows:
        key, criterion = row["variant_key"], row["criterion"].upper()
        if key not in keys:
            continue
        if (key, criterion) in seen:
            raise ValueError(
                f"duplicate criterion {key}/{criterion} (do not double count)"
            )
        seen.add((key, criterion))
        strength = row["strength"]
        if strength not in LR or not row["source_id"].strip():
            raise ValueError(f"invalid criterion strength/source for {key}/{criterion}")
        if row["audit_status"] != "verified":
            row["included"] = "no_unverified"
        elif criterion in {"PP3", "BP4"} or row.get("computational_family", "").strip():
            if criterion not in {"PP3", "BP4"}:
                row["included"] = "no_uncalibrated_computational"
                continue
            if strength != (
                "supporting" if criterion == "PP3" else "benign_supporting"
            ):
                raise ValueError(
                    f"{key}/{criterion}: computational evidence cannot exceed supporting"
                )
            group = tuple(
                row.get(name, "")
                for name in ("predictor", "consequence", "af_bin", "orientation")
            )
            model = calibration.get(group)
            if model is None or not row.get("score"):
                row["included"] = "no_uncalibrated_computational"
                continue
            threshold, high, low, n_train, n_held = model
            case_score = float(row["score"])
            if not math.isfinite(case_score):
                raise ValueError(f"{key}/{criterion}: predictor score not finite")
            hit = (1 if group[3] == "high" else -1) * case_score >= threshold
            ratio = high if hit else low
            if not (ratio > 1 if criterion == "PP3" else ratio < 1):
                row["included"] = "no_predictor_direction_disagrees_criterion"
                continue
            row["empirical_lr"] = fmt(
                min(ratio, 2.08) if ratio > 1 else max(ratio, 1 / 2.08)
            )
            row["calibration_threshold"] = fmt(threshold)
            row["n_calibration"] = n_train
            row["n_heldout"] = n_held
            row["included"] = "yes_empirical_gene_heldout"
            evidence[key].append(row)
            used_sources[key, row["source_id"]].add(criterion)
        elif criterion == "BA1":
            row["included"] = "no_BA1_standalone_LR_mapping"
        elif not (criterion.startswith(("PVS", "PS", "PM", "PP", "BS", "BP"))):
            row["included"] = "no_unknown_criterion"
        elif criterion.startswith(("BA", "BS", "BP")) != strength.startswith("benign"):
            raise ValueError(f"criterion/strength direction conflict {key}/{criterion}")
        else:
            row["included"] = "yes"
            evidence[key].append(row)
            used_sources[key, row["source_id"]].add(criterion)
    for key in keys:
        comp = [r for r in evidence[key] if r["criterion"].upper() in ("PP3", "BP4")]
        if len(comp) > 1:
            chosen = max(
                comp,
                key=lambda r: (abs(math.log(float(r["empirical_lr"]))), r["source_id"]),
            )
            for r in comp:
                if r is not chosen:
                    evidence[key].remove(r)
                    r["included"] = "no_correlated_computational_family"
                    used_sources[key, r["source_id"]].discard(r["criterion"].upper())
    for (key, source), criteria in used_sources.items():
        if len(criteria) > 1:
            raise ValueError(
                f"{key}: source {source} assigned multiple independent criteria: {sorted(criteria)}"
            )
    return evidence, rows


def verify_case_calibration(path, pair_ids):
    """Hard gate: independent labeled, gene-blocked held-out case controls.

    This does not estimate a q or validate its reliability: a separate audited
    calibrator must have produced the supplied q values and their source ID.
    """
    if path is None or not Path(path).exists():
        return set(), "independent labeled case calibration table absent"
    rows, _ = table(
        path, ("source_id", "link", "gene", "label", "partition", "case_id")
    )
    by_source = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for r in rows:
        if r["case_id"] in pair_ids:
            raise ValueError(
                "case-calibration contains target case identifier (label leakage)"
            )
        if (
            r["link"] not in LINKS
            or r["label"] not in {"positive", "negative"}
            or r["partition"] not in {"train", "heldout"}
        ):
            raise ValueError(
                "case-calibration requires recognized link, binary label and train/heldout partition"
            )
        if (
            not r["source_id"]
            or not r["gene"]
            or not r["case_id"]
            or r["gene"] in TIER1
        ):
            raise ValueError(
                "case-calibration requires source/case/gene and excludes Tier 1 genes"
            )
        by_source[r["source_id"]][r["link"]][r["partition"], r["label"]].add(r["gene"])
    eligible = set()
    for source, links in by_source.items():
        if all(
            len(links[link].get((partition, label), set())) >= 5
            for link in LINKS
            for partition in ("train", "heldout")
            for label in ("positive", "negative")
        ):
            for link in LINKS:
                train = set().union(
                    *(
                        links[link].get(("train", label), set())
                        for label in ("positive", "negative")
                    )
                )
                held = set().union(
                    *(
                        links[link].get(("heldout", label), set())
                        for label in ("positive", "negative")
                    )
                )
                if train & held:
                    raise ValueError(
                        "case-calibration gene leakage across train/heldout"
                    )
            eligible.add(source)
    return (
        eligible,
        ""
        if eligible
        else "insufficient gene-blocked train/heldout labeled real cases per link and class",
    )


def case_links(path, pair_ids, synthetic, eligible_sources):
    if path is None or not Path(path).exists():
        return {}, "case-level calibrated q links absent"
    rows, _ = table(path, ("pair_id", "link", "q", "source_id", "calibration_status"))
    by_pair = defaultdict(dict)
    for r in rows:
        if r["pair_id"] not in pair_ids:
            continue
        if r["link"] not in LINKS or r["link"] in by_pair[r["pair_id"]]:
            raise ValueError(f"invalid or duplicated link {r['pair_id']}/{r['link']}")
        if not r["source_id"].strip() or r["calibration_status"] != (
            "synthetic" if synthetic else "heldout_real"
        ):
            raise ValueError(
                f"unverified calibration provenance {r['pair_id']}/{r['link']}"
            )
        if not synthetic and r["source_id"] not in eligible_sources:
            return (
                {},
                "case q source lacks independent gene-blocked held-out labeled calibration",
            )
        value = float(r["q"])
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"invalid calibrated q {r['pair_id']}/{r['link']}")
        by_pair[r["pair_id"]][r["link"]] = (value, r["source_id"])
    return by_pair, ""


def fmt(x):
    return format(x, ".9g")


def score(
    row, evidence, links, prior=0.1, unknown=0.5, multiplier=1.0, source_counts=None
):
    keys = [k for k in (allele_key(row, 1), allele_key(row, 2)) if k]
    probabilities = [
        posterior(
            evidence.get(key, []),
            prior,
            multiplier,
            source_counts,
        )
        for key in keys
    ]
    q_variant = (
        min(probabilities) * phase_probability(row, unknown)
        if len(keys) == 2
        else probabilities[0]
    )
    if not links or set(links) != set(LINKS):
        # Preregistered unverified-link cap: descriptive provisional ePCR,
        # NEVER an estimate of a clinically validated probability.
        return q_variant, min(0.1, q_variant)
    # q_variant is a variant/pathogenicity heuristic, not a causal probability.
    q_case = min([q_variant] + [links[link][0] for link in LINKS])
    return q_variant, min(1.0, max(sys.float_info.min, q_case))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked-a", type=Path, default=ROOT / "results/ranked_a.tsv")
    ap.add_argument(
        "--angle-b", type=Path, default=ROOT / "results/angle_b_phenotype.tsv"
    )
    ap.add_argument(
        "--angle-c", type=Path, default=ROOT / "results/angle_c_mechanism.tsv"
    )
    ap.add_argument(
        "--angle-d", type=Path, default=ROOT / "results/literature_cases.tsv"
    )
    ap.add_argument(
        "--criteria", type=Path, default=ROOT / "results/allele_criteria.tsv"
    )
    ap.add_argument(
        "--clinvar-calibration",
        type=Path,
        default=ROOT / "results/clinvar_computational_calibration.tsv",
        help="optional two-star nonconflicting ClinVar rows, with gene-disjoint calibration/heldout splits",
    )
    ap.add_argument(
        "--case-links", type=Path, default=ROOT / "results/case_links_calibrated.tsv"
    )
    ap.add_argument(
        "--case-calibration",
        type=Path,
        default=ROOT / "results/case_calibration_heldout.tsv",
    )
    ap.add_argument("--output-dir", type=Path, default=ROOT / "results")
    ap.add_argument(
        "--synthetic",
        action="store_true",
        help="synthetic evidence only, output outside clinical results",
    )
    ap.add_argument("--bootstrap-reps", type=int, default=400)
    ap.add_argument("--seed", type=int, default=32026)
    args = ap.parse_args(argv)
    if args.synthetic and args.output_dir.resolve() == (ROOT / "results").resolve():
        ap.error("synthetic mode requires --output-dir separate from clinical results")
    if args.bootstrap_reps < 100:
        ap.error("bootstrap-reps must be at least 100")
    pairs, cols = table(args.ranked_a, BASE)
    if not args.synthetic:
        manifest = args.ranked_a.with_name("ranked_a_manifest.json")
        if not manifest.is_file():
            raise ValueError("Frozen A manifest absent")
        locked = json.loads(manifest.read_text(encoding="utf-8"))
        if hashlib.sha256(args.ranked_a.read_bytes()).hexdigest() != locked[
            "ranked_a_sha256"
        ] or not locked.get("no_B_C_D_consumed"):
            raise ValueError("Frozen A checksum/independence mismatch")
    if not pairs:
        raise ValueError(
            "ranked A contains no verified allele pairs; refuse empty outputs"
        )
    destinations = tuple(
        args.output_dir / filename
        for filename in (
            "ranked_pairs.tsv",
            "epcr.tsv",
            "acmg_evidence.tsv",
            "epcr_sensitivity.tsv",
            "epcr_bootstrap.tsv",
        )
    )
    if any(path.exists() for path in destinations):
        raise FileExistsError(
            "32 refuses overwrite of ranked/evidence/sensitivity/bootstrap outputs; use new --output-dir"
        )
    if len({r["pair_id"] for r in pairs}) != len(pairs) or any(
        not r["pair_id"] for r in pairs
    ):
        raise ValueError("pair_id must be nonempty and unique")
    keys = {allele_key(r, i) for r in pairs for i in (1, 2)} - {""}
    calibration = empirical_computational_calibration(args.clinvar_calibration, keys)
    loaded = load_evidence(args.criteria, keys, calibration)
    evidence, audit = loaded if isinstance(loaded, tuple) else (loaded, [])
    eligible, calibration_gap = (
        (set(), "")
        if args.synthetic
        else verify_case_calibration(
            args.case_calibration, {r["pair_id"] for r in pairs}
        )
    )
    by_pair, gap = case_links(
        args.case_links, {r["pair_id"] for r in pairs}, args.synthetic, eligible
    )
    if not args.synthetic and calibration_gap:
        by_pair = {}
        gap = calibration_gap
    angles = read_angles(args)
    ranked, sensitivity, bootstrap = [], [], []
    rng = random.Random(args.seed)
    for r in pairs:
        links = by_pair.get(r["pair_id"], {})
        n = pick_angles(r, *angles)
        # Exactly preregistered: for n>=1 A=min(1,1+0.1(n-1)) = 1.
        angle_multiplier = min(1.0, 1 + 0.1 * (n - 1)) if n else 0.0
        qv, epcr = score(r, evidence, links)
        if n == 0:
            n = 1  # the frozen A row is itself an angle, even if read support is unverified
            angle_multiplier = 1.0
        epcr *= angle_multiplier
        if (
            r.get("evidence_status") != "read_supported"
            or r.get("angle_a_status") != "verified"
        ):
            epcr = min(epcr, 0.1)
        if r.get("phase") == "cis":
            epcr = min(
                epcr, 1e-12
            )  # a cis pair cannot explain a two-hit recessive mechanism
        status = (
            "synthetic_scored"
            if args.synthetic
            else (
                "heuristic_case_links_calibrated"
                if set(links) == set(LINKS)
                else "provisional_unverified_link_cap_0.1"
            )
        )
        keys_in_pair = [allele_key(r, i) for i in (1, 2) if allele_key(r, i)]
        ps = [posterior(evidence.get(key, []), 0.1) for key in keys_in_pair]
        ranked.append(
            {
                **r,
                "epcr": fmt(max(sys.float_info.min, epcr)),
                "p_path_1": fmt(ps[0]),
                "p_path_2": fmt(ps[1]) if len(ps) > 1 else "",
                "p_pair": fmt(qv),
                "variant_posterior": fmt(qv),
                "independent_angle_count": n,
                "epcr_status": status,
                "submission_ready": "yes"
                if epcr > 0 and r.get("phase") != "cis"
                else "no",
                "proband_id": r.get("proband_id", ""),
                "acmg_criteria_1": ";".join(
                    sorted(x["criterion"] for x in evidence.get(keys_in_pair[0], []))
                ),
                "acmg_criteria_2": ";".join(
                    sorted(x["criterion"] for x in evidence.get(keys_in_pair[1], []))
                )
                if len(keys_in_pair) > 1
                else "",
            }
        )
        for prior in (0.01, 0.1, 0.2):
            for factor in (0.8, 1.0, 1.2):
                for phase_p in (0.25, 0.5, 0.75):
                    sv, se = score(r, evidence, links, prior, phase_p, factor)
                    sensitivity.append(
                        {
                            "pair_id": r["pair_id"],
                            "prior": fmt(prior),
                            "strength_loglr_scale": fmt(factor),
                            "unknown_phase_probability": fmt(phase_p),
                            "variant_posterior": fmt(sv),
                            "epcr": fmt(
                                max(
                                    sys.float_info.min,
                                    min(
                                        se * angle_multiplier,
                                        0.1
                                        if status
                                        == "provisional_unverified_link_cap_0.1"
                                        or r.get("evidence_status") != "read_supported"
                                        else 1.0,
                                    ),
                                )
                            )
                            if r.get("phase") != "cis"
                            else "1e-12",
                            "epcr_status": status,
                        }
                    )
        sources = sorted(
            {
                e["source_id"]
                for k in (allele_key(r, 1), allele_key(r, 2))
                for e in evidence.get(k, [])
            }
        )
        if len(sources) >= 2:
            # Cluster resampling of independent criterion sources; q links remain fixed,
            # so this is conditional uncertainty, not a calibrated case confidence interval.
            for rep in range(args.bootstrap_reps):
                sample = [rng.choice(sources) for _ in sources]
                counts = {src: sample.count(src) for src in sources}
                bv, be = score(r, evidence, links, source_counts=counts)
                assert be is not None
                bootstrap.append(
                    {
                        "pair_id": r["pair_id"],
                        "replicate": rep,
                        "resampled_sources": ";".join(sample),
                        "variant_posterior": fmt(bv),
                        "epcr": fmt(
                            min(
                                be * angle_multiplier,
                                0.1
                                if status == "provisional_unverified_link_cap_0.1"
                                or r.get("evidence_status") != "read_supported"
                                else 1.0,
                            )
                        ),
                        "bootstrap_scope": "criterion_sources_only_links_fixed",
                    }
                )
    ranked.sort(
        key=lambda r: (
            -float(r["epcr"]) if r["epcr"] else math.inf,
            -sum(
                int(r[f"alt_reads_{i}"]) if r[f"alt_reads_{i}"].isdigit() else 0
                for i in (1, 2)
            ),
            allele_key(r, 1),
            allele_key(r, 2),
            r["pair_id"],
        )
    )
    # Provisional scores are bounded heuristic prioritization values, not
    # clinical probabilities. Critical links need observed independent evidence.
    out = args.output_dir
    write(
        out / "ranked_pairs.tsv",
        list(
            dict.fromkeys(
                list(cols)
                + [
                    "epcr",
                    "p_pair",
                    "submission_ready",
                    "proband_id",
                    "acmg_criteria_1",
                    "acmg_criteria_2",
                    "variant_posterior",
                    "independent_angle_count",
                    "epcr_status",
                ]
            )
        ),
        ranked,
    )
    write(
        out / "epcr.tsv",
        [
            "pair_id",
            "gene",
            "p_path_1",
            "p_path_2",
            "p_pair",
            "epcr",
            "independent_angle_count",
            "epcr_status",
        ],
        ranked,
    )
    write(
        out / "acmg_evidence.tsv",
        [
            "variant_key",
            "criterion",
            "strength",
            "source_id",
            "audit_status",
            "computational_family",
            "included",
            "predictor",
            "score",
            "consequence",
            "af_bin",
            "orientation",
            "empirical_lr",
            "calibration_threshold",
            "n_calibration",
            "n_heldout",
        ],
        audit,
    )
    write(
        out / "epcr_sensitivity.tsv",
        [
            "pair_id",
            "prior",
            "strength_loglr_scale",
            "unknown_phase_probability",
            "variant_posterior",
            "epcr",
            "epcr_status",
        ],
        sensitivity,
    )
    write(
        out / "epcr_bootstrap.tsv",
        [
            "pair_id",
            "replicate",
            "resampled_sources",
            "variant_posterior",
            "epcr",
            "bootstrap_scope",
        ],
        bootstrap,
    )
    print(
        f"pairs={len(pairs)} scored={sum(bool(r['epcr']) for r in ranked)} bootstrap_replicates={len(bootstrap)}; {gap or 'links supplied'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
