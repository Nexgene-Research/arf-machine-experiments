"""Frozen-pipeline controls; fail closed when real control inputs are absent.

Manifest JSON schema: {
  "hg002_hpo_terms": ["HP:...", ...],
  "original_hpo_terms": ["HP:...", ...],
  "hg002": {"ranked_a": "...", "scored_dir": "...", "commands": [[PYTHON, ABSOLUTE_30, ...], [PYTHON, ABSOLUTE_32, ...]]},
  "positive": {...same...}, "case": {...same...},
  "permutations": [{"terms": [...], "ranked_a": "...", "scored_dir": "...",
                     "commands": [[PYTHON, ABSOLUTE_10, ...], [PYTHON, ABSOLUTE_30, ...], [PYTHON, ABSOLUTE_32, ...]]}, ...],
  "hpo_depth_universe": "TSV with term,depth_bin",
  "published_pairs": "TSV with case_id,gene,assembly,chrom_1,pos_1,ref_1,alt_1,
     chrom_2,pos_2,ref_2,alt_2,gt_1,gt_2,phase,normalized,ref_checked,source_id",
  "positive_candidate_table": "30-format candidate TSV physically containing BOTH
     exact normalized published-case allele rows (case_id additional metadata)"
}
All command arrays are explicit, no shell, and must invoke the frozen scripts.
Manifest is control metadata only, never a claim of real genotype or causal rank.
The HG002 and positive runs must pass through the actual script 30 then script 32
exactly as the case run does; each HPO permutation reruns 10 then 30 then 32.
No attempt to reconstruct genotype tables from an incomplete case catalogue.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = {
    "10": ROOT / "scripts/10_phenotype.py",
    "30": ROOT / "scripts/30_rank_pairs.py",
    "32": ROOT / "scripts/32_epcr.py",
}
ORIGINAL = (
    "HP:0002859",
    "HP:0000121",
    "HP:0004322",
    "HP:0001508",
    "HP:0003202",
    "HP:0001622",
    "HP:0001518",
)
SEED = 34026


def read_tsv(path, fields):
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not set(fields) <= set(reader.fieldnames):
            raise ValueError(
                f"{path}: missing required columns {sorted(set(fields) - set(reader.fieldnames or []))}"
            )
        out = []
        for row in reader:
            if len(out) >= 100000:
                raise ValueError("controls table exceeds laptop bound; stage on MN5")
            out.append(row)
    return out


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ensure_commands(item, stages, frozen):
    commands = item.get("commands", [])
    if len(commands) != len(stages):
        raise ValueError(
            "control must run each frozen stage once and in order: " + ",".join(stages)
        )
    for cmd, stage in zip(commands, stages):
        if (
            not isinstance(cmd, list)
            or len(cmd) < 2
            or Path(cmd[1]).resolve() != SCRIPT[stage].resolve()
        ):
            raise ValueError(
                f"control stage {stage} must invoke exact script {SCRIPT[stage]}"
            )
        if Path(cmd[0]).resolve() != (ROOT / "env/venv/Scripts/python.exe").resolve():
            raise ValueError(
                "all controls must use the same pinned project interpreter"
            )
        if stage == "30":
            for arg in (
                "--candidates",
                "--pairs",
                "--read-check",
                "--scores",
                "--out",
                "--manifest",
            ):
                if arg not in cmd or cmd.index(arg) + 1 == len(cmd):
                    raise ValueError(
                        f"30 controls require explicit {arg} to prevent patient input/output reuse"
                    )
            if (
                Path(cmd[cmd.index("--out") + 1]).resolve()
                != Path(item["ranked_a"]).resolve()
            ):
                raise ValueError("30 --out must equal manifest ranked_a")
        if stage == "32" and (
            "--synthetic" in cmd or "--output-dir" not in cmd or "--ranked-a" not in cmd
        ):
            raise ValueError(
                "real control 32 invocation requires --ranked-a and --output-dir without --synthetic"
            )
        if stage == "32":
            for arg in (
                "--case-calibration",
                "--case-links",
                "--criteria",
                "--angle-b",
            ):
                if (
                    arg not in cmd
                    or cmd.index(arg) + 1 == len(cmd)
                    or not Path(cmd[cmd.index(arg) + 1]).exists()
                ):
                    raise ValueError(
                        f"32 real controls need supplied existing {arg}, not default/missing calibration"
                    )
            if (
                Path(cmd[cmd.index("--output-dir") + 1]).resolve()
                != Path(item["scored_dir"]).resolve()
            ):
                raise ValueError(
                    "32 --output-dir must equal unique manifest scored_dir"
                )
            if (
                Path(cmd[cmd.index("--ranked-a") + 1]).resolve()
                != Path(item["ranked_a"]).resolve()
            ):
                raise ValueError("32 must consume 30 ranked A for this control")
        if stage == "10":
            if "--hpo-terms-file" not in cmd or "--output" not in cmd:
                raise ValueError(
                    "permutation stage 10 must accept explicit HPO terms and separate output; otherwise it repeats proband HPO"
                )
            terms_path = Path(cmd[cmd.index("--hpo-terms-file") + 1])
            if sorted(terms_path.read_text(encoding="utf-8").split()) != sorted(
                item["terms"]
            ):
                raise ValueError(
                    "permutation 10 HPO input differs from fixed-seed depth-matched terms"
                )
            if "--hpo-terms-file" not in SCRIPT["10"].read_text(encoding="utf-8"):
                raise ValueError(
                    "frozen 10_phenotype.py lacks --hpo-terms-file; cannot produce genuine HPO permutation control"
                )
        if digest(SCRIPT[stage]) != frozen[stage]:
            raise ValueError(f"frozen pipeline script {stage} changed")


def check_published(path, candidate_table):
    fields = (
        "case_id",
        "gene",
        "assembly",
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
        "phase",
        "normalized",
        "ref_checked",
        "source_id",
    )
    rows = read_tsv(path, fields)
    if not rows:
        raise ValueError(
            "no published pairs with both exact resolved alleles and case-specific genotypes"
        )
    for r in rows:
        if (
            r["assembly"] != "GRCh38"
            or r["normalized"] != "yes"
            or r["ref_checked"] != "yes"
            or r["phase"]
            not in {
                "trans_verified",
                "trans_familial",
                "trans_physical",
            }
            or any(r[f"gt_{i}"].replace("|", "/") not in {"0/1", "1/0"} for i in (1, 2))
            or (r["chrom_1"], r["pos_1"], r["ref_1"], r["alt_1"])
            == (r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"])
            or not all(r[k] for k in fields)
            or any(
                not r[f"chrom_{i}"].startswith("chr")
                or not r[f"pos_{i}"].isdigit()
                or not r[f"gt_{i}"]
                for i in (1, 2)
            )
        ):
            raise ValueError(
                f"published pair {r['case_id']} unresolved; cannot inject guessed normalized allele or genotype"
            )
    injected = read_tsv(
        candidate_table,
        (
            "case_id",
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
    for r in rows:
        for i in (1, 2):
            matching = [
                x
                for x in injected
                if x["case_id"] == r["case_id"]
                and x["gene"] == r["gene"]
                and all(
                    x[k] == r[f"{k}_{i}"] for k in ("chrom", "pos", "ref", "alt", "gt")
                )
            ]
            if (
                len(matching) != 1
                or matching[0]["genotype_status"] != "validated"
                or matching[0]["normalization_status"] != "bcftools_norm_GRCh38"
                or matching[0]["ref_matches_grch38"] != "true"
            ):
                raise ValueError(
                    f"published pair {r['case_id']} allele {i}: absent or unvalidated in 30 candidate table"
                )
    return rows


def check_permutations(manifest):
    original = manifest.get("original_hpo_terms")
    if not isinstance(original, list) or sorted(original) != sorted(ORIGINAL):
        raise ValueError(
            "original phenotype must equal locked 7 proband HPO terms; parental history excluded"
        )
    if manifest.get("hg002_hpo_terms") != original:
        raise ValueError("HG002 must use exact proband HPO terms unchanged")
    universe = read_tsv(manifest["hpo_depth_universe"], ("term", "depth_bin"))
    depth = {}
    for r in universe:
        if r["term"] in depth and depth[r["term"]] != r["depth_bin"]:
            raise ValueError("conflicting ontology depth assignments")
        depth[r["term"]] = r["depth_bin"]
    if any(t not in depth for t in original):
        raise ValueError("original terms missing from ontology depth universe")
    bins = Counter(depth[t] for t in original)
    rng = random.Random(SEED)
    by_bin = {b: sorted(t for t, d in depth.items() if d == b) for b in bins}
    if any(len(by_bin[b]) < count for b, count in bins.items()):
        raise ValueError("not enough ontology terms within original broad-depth bins")
    permutations = manifest.get("permutations", [])
    if len(permutations) < 100:
        raise ValueError(
            "at least 100 fixed-seed depth-matched HPO permutations are required"
        )
    for idx, item in enumerate(permutations):
        expected = sorted(
            t for b, n in sorted(bins.items()) for t in rng.sample(by_bin[b], n)
        )
        if sorted(item.get("terms", [])) != expected:
            raise ValueError(
                f"permutation {idx} differs from fixed seed {SEED} depth-matched draw"
            )
    return len(permutations)


def scored(path):
    rows = read_tsv(Path(path) / "ranked_pairs.tsv", ("pair_id", "epcr", "epcr_status"))
    if not rows or any(
        not r["epcr"] or r["epcr_status"] != "heuristic_case_links_calibrated"
        for r in rows
    ):
        raise ValueError(
            "unscored ePCR; numeric null and ranks prohibited until real case calibration exists"
        )
    vals = [float(r["epcr"]) for r in rows]
    if any(not 0 < v <= 1 for v in vals):
        raise ValueError("invalid ePCR value")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--manifest", type=Path, help="real-label control manifest with command arrays"
    )
    ap.add_argument("--output-dir", type=Path, default=ROOT / "results/controls")
    ap.add_argument(
        "--preflight",
        action="store_true",
        help="check inputs without launching pipeline",
    )
    args = ap.parse_args(argv)
    if not args.manifest or not args.manifest.exists():
        ap.error(
            "control manifest absent: need HG002 genotype/candidate table, published exact normalized case-pair genotypes inserted into 30 candidate input, labeled calibration, and >=100 ontology-depth-matched HPO permutations; no numeric null/rank"
        )
    m = json.loads(args.manifest.read_text(encoding="utf-8"))
    if not all(p.exists() for p in SCRIPT.values()):
        raise ValueError(
            "frozen 10/30/32 scripts not all present; cannot use proxy control"
        )
    frozen = {k: digest(v) for k, v in SCRIPT.items()}
    n = check_permutations(m)
    pairs = check_published(m["published_pairs"], m["positive_candidate_table"])
    items = [m[k] for k in ("case", "hg002", "positive")]
    items += m["permutations"]
    directories = [Path(item["scored_dir"]).resolve() for item in items]
    if len(set(directories)) != len(directories):
        raise ValueError(
            "case/control/permutation outputs must use distinct directories"
        )
    for i, item in enumerate(items):
        ensure_commands(item, ("10", "30", "32") if i >= 3 else ("30", "32"), frozen)
        if (
            Path(item["ranked_a"]).exists()
            or (Path(item["scored_dir"]) / "ranked_pairs.tsv").exists()
        ):
            raise FileExistsError(
                "control stage outputs already exist; use new per-run paths to avoid overwriting earlier results"
            )
        if i >= 3:
            phenotype_command, _, epcr_command = item["commands"]
            if (
                Path(
                    phenotype_command[phenotype_command.index("--output") + 1]
                ).resolve()
                != Path(epcr_command[epcr_command.index("--angle-b") + 1]).resolve()
            ):
                raise ValueError(
                    "HPO permutation must feed its own newly scored B table to 32"
                )

    def angle_path(item):
        command = item["commands"][-1]
        return Path(command[command.index("--angle-b") + 1]).resolve()

    if len({angle_path(item) for item in items[:3]}) != 1:
        raise ValueError(
            "HG002, positive controls and case must use the same frozen proband phenotype B table"
        )
    cmd_30_positive = items[2]["commands"][0]
    if str(Path(m["positive_candidate_table"]).resolve()) not in [
        str(Path(x).resolve()) for x in cmd_30_positive[2:]
    ]:
        raise ValueError(
            "positive candidate table must be passed directly to script 30"
        )
    submitted = read_tsv(
        Path(cmd_30_positive[cmd_30_positive.index("--pairs") + 1]),
        ("pair_id", "gene", "variant_1", "variant_2"),
    )
    for pair in pairs:
        first = "-".join(pair[k] for k in ("chrom_1", "pos_1", "ref_1", "alt_1"))
        second = "-".join(pair[k] for k in ("chrom_2", "pos_2", "ref_2", "alt_2"))
        if not any(
            row["pair_id"] == pair["case_id"]
            and row["gene"] == pair["gene"]
            and {row["variant_1"], row["variant_2"]} == {first, second}
            for row in submitted
        ):
            raise ValueError(
                f"published pair {pair['case_id']} absent as exact 30 pair candidate"
            )
    read_checks = read_tsv(
        Path(cmd_30_positive[cmd_30_positive.index("--read-check") + 1]),
        ("variant", "alt_reads", "genotype_status", "read_provenance"),
    )
    for pair in pairs:
        for i in (1, 2):
            key = "-".join(
                pair[k] for k in (f"chrom_{i}", f"pos_{i}", f"ref_{i}", f"alt_{i}")
            )
            matching = [r for r in read_checks if r["variant"] == key]
            if (
                len(matching) != 1
                or matching[0]["genotype_status"] != "validated"
                or matching[0]["read_provenance"] != "observed_real_molecules"
                or not matching[0]["alt_reads"].isdigit()
                or int(matching[0]["alt_reads"]) < 1
            ):
                raise ValueError(
                    f"{key}: 30 requires independently validated real ALT reads; cannot fabricate VCF-only positive-control read counts"
                )
    if args.preflight:
        print(
            f"preflight OK: {n} depth-matched permutations, {len(pairs)} fully specified published pairs, frozen 10/30/32 scripts; no numeric control claimed"
        )
        return 0
    if (args.output_dir / "control_null.tsv").exists():
        raise FileExistsError(
            "control_null.tsv already exists; rerun to a new output directory"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for item in items:
        for cmd in item["commands"]:
            if (
                digest(Path(cmd[1]))
                != frozen[
                    next(
                        k
                        for k, v in SCRIPT.items()
                        if Path(cmd[1]).resolve() == v.resolve()
                    )
                ]
            ):
                raise ValueError("frozen script changed during control execution")
            subprocess.run(cmd, cwd=ROOT, check=True)
        scored(item["scored_dir"])  # fail before consuming an unscored output
    case = scored(items[0]["scored_dir"])
    hg002 = scored(items[1]["scored_dir"])
    positive = scored(items[2]["scored_dir"])
    nulls = [scored(item["scored_dir"]) for item in items[3:]]
    positive_ranks = {r["pair_id"]: index for index, r in enumerate(positive, 1)}
    if not {p["case_id"] for p in pairs} <= set(positive_ranks):
        raise ValueError(
            "positive control pair IDs not retained by 30/32; rank not reportable"
        )
    positive_summary = ";".join(
        f"{p['case_id']}:{positive_ranks[p['case_id']]}" for p in pairs
    )
    thresholds = sorted({float(r["epcr"]) for r in case}, reverse=True)
    path = args.output_dir / "control_null.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            (
                "null_regime",
                "threshold",
                "empirical_tail_p",
                "false_candidates_per_genome",
                "n_null_genomes_permutations",
                "positive_pair_ranks",
            )
        )
        for regime, samples in (
            ("hg002_fixed_hpo", [hg002]),
            ("proband_depth_matched_hpo_permutations", nulls),
        ):
            null_max = [max(float(r["epcr"]) for r in rows) for rows in samples]
            for t in thresholds:
                w.writerow(
                    (
                        regime,
                        format(t, ".9g"),
                        format(
                            (1 + sum(v >= t for v in null_max)) / (1 + len(null_max)),
                            ".9g",
                        ),
                        format(
                            sum(
                                sum(float(r["epcr"]) >= t for r in rows)
                                for rows in samples
                            )
                            / len(samples),
                            ".9g",
                        ),
                        len(samples),
                        positive_summary,
                    )
                )
    print(
        f"{path}: {n} permutations + HG002; empirical p and false candidates calculated only from actual 30/32 outputs"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
