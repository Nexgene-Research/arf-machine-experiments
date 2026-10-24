"""Preregistered heuristic ePCR for Track 1 hypotheses.

Automatic criteria only. AlphaMissense and REVEL stay descriptive (LR=1).
An unverified read or phase link caps case-level ePCR at 0.1. This is not a
validated clinical probability.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import sys
from pathlib import Path

TIER1 = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1"}
LOF_GENES = {"BUB1B", "TRIP13"}
LOF = ("frameshift", "stop_gained", "stop_lost", "start_lost", "splice_acceptor", "splice_donor")
STRENGTH = {
    "very_strong": 350.0,
    "strong": 18.7,
    "moderate": 4.3,
    "supporting": 2.08,
    "benign_strong": 1 / 18.7,
    "benign_supporting": 1 / 2.08,
}
PRIORS = (0.01, 0.1, 0.2)


def logit(p: float) -> float:
    p = min(max(p, 1e-9), 1 - 1e-9)
    return math.log(p / (1 - p))


def logistic(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1 / (1 + z)
    z = math.exp(x)
    return z / (1 + z)


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def chrom_of(value: str) -> str:
    if value.startswith("chr"):
        return value
    return "chrM" if value in ("M", "MT") else "chr" + value


def stars(review: str) -> int:
    text = (review or "").lower()
    if "expert_panel" in text or "practice_guideline" in text:
        return 3
    if "multiple_submitters" in text and "conflict" not in text:
        return 2
    if "single_submitter" in text:
        return 1
    return 0


def criteria(row: dict[str, str], index: int, clin: dict[tuple[str, str, str, str], tuple[str, str]]) -> list[tuple[str, str]]:
    gene = row.get("gene", "")
    consequence = (row.get(f"consequence_{index}") or "").lower()
    impact = row.get(f"impact_{index}") or ""
    af = number(row.get(f"gnomad_af_{index}", ""))
    state = row.get(f"af_state_{index}") or ""
    key = (
        chrom_of(row.get(f"chrom_{index}", "")),
        row.get(f"pos_{index}", ""),
        (row.get(f"ref_{index}") or "").upper(),
        (row.get(f"alt_{index}") or "").upper(),
    )
    sig, review = clin.get(key, ((row.get(f"clinvar_{index}") or ""), ""))
    found: list[tuple[str, str]] = []
    lof = impact == "HIGH" or any(token in consequence for token in LOF)
    if lof and gene in LOF_GENES:
        found.append(("PVS1", "very_strong"))
    if state == "missing" or (af is not None and af <= 1e-4):
        found.append(("PM2", "moderate"))
    if stars(review) >= 2 and "conflict" not in sig.lower():
        low = sig.lower()
        if "likely_pathogenic" in low:
            found.append(("PS1", "supporting"))
        elif "pathogenic" in low and "benign" not in low:
            found.append(("PS1", "strong"))
        elif "likely_benign" in low:
            found.append(("BP6", "benign_supporting"))
        elif low.startswith("benign") or ";benign" in low or low == "benign":
            found.append(("BA1", "benign_strong"))
    if gene in TIER1:
        found.append(("PP4", "supporting"))
    return found


def path_prob(bits: list[tuple[str, str]], prior: float) -> float:
    total = logit(prior)
    seen: set[str] = set()
    for name, strength in bits:
        if name in seen:
            continue
        seen.add(name)
        total += math.log(STRENGTH[strength])
    return logistic(total)


def main() -> None:
    root = Path(sys.argv[1])
    vcf = Path(sys.argv[2])
    with (root / "hypotheses.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    phase: dict[tuple[str, str, str, str], dict[str, str]] = {}
    phase_path = root / "read_phase_top60.tsv"
    if phase_path.exists():
        with phase_path.open(newline="") as handle:
            for item in csv.DictReader(handle, delimiter="\t"):
                if item.get("source") != "top60":
                    continue
                phase[(item["chrom_1"], item["pos_1"], item["chrom_2"], item["pos_2"])] = item
    wanted = set()
    for row in rows:
        for index in ("1", "2"):
            if row.get(f"chrom_{index}"):
                wanted.add(
                    (
                        chrom_of(row[f"chrom_{index}"]),
                        row[f"pos_{index}"],
                        row[f"ref_{index}"].upper(),
                        row[f"alt_{index}"].upper(),
                    )
                )
    clin: dict[tuple[str, str, str, str], tuple[str, str]] = {}
    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            key = (chrom_of(parts[0]), parts[1], parts[3].upper(), parts[4].upper())
            if key not in wanted:
                continue
            info = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    name, value = item.split("=", 1)
                    info[name] = value
            clin[key] = (info.get("CLNSIG", ""), info.get("CLNREVSTAT", ""))
            if len(clin) == len(wanted):
                break

    ranked = []
    for row in rows:
        evidence = phase.get((row["chrom_1"], row["pos_1"], row.get("chrom_2", ""), row.get("pos_2", "")), {})
        relation = evidence.get("phase") or row.get("phase") or "unknown"
        if relation == "not_applicable":
            relation = "homozygous" if row["hypothesis"] == "homozygous" else "unknown"
        bits1 = criteria(row, 1, clin)
        bits2 = criteria(row, 2, clin) if row.get("chrom_2") else []
        if relation == "trans":
            bits1.append(("PM3", "strong"))
            if bits2:
                bits2.append(("PM3", "strong"))
        p1 = {prior: path_prob(bits1, prior) for prior in PRIORS}
        p2 = {prior: path_prob(bits2, prior) for prior in PRIORS} if bits2 else None
        p_trans = {"trans": 1.0, "unknown": 0.5, "cis": 0.05}.get(relation, 1.0)
        verified = bool(evidence) and evidence.get("read_1") == "snv_pileup" and (
            row["hypothesis"] != "compound_het" or evidence.get("read_2") == "snv_pileup"
        )
        reads_ok = verified and all(
            (evidence.get(f"pile_alt_{n}") or "0") not in ("", "0")
            for n in ("1", "2")
            if row.get(f"chrom_{n}")
        )
        phase_ok = row["hypothesis"] != "compound_het" or relation in ("trans", "cis")
        for prior in PRIORS:
            if p2 is None:
                q = p1[prior]
            else:
                q = min(p1[prior], p2[prior]) * p_trans
            epcr = min(max(q, 1e-6), 1.0)
            if not reads_ok or not phase_ok:
                epcr = min(epcr, 0.1)
            if prior == 0.1:
                row_out = dict(row)
                row_out.update(
                    {
                        "phase_called": relation,
                        "p_path_1": f"{p1[prior]:.6g}",
                        "p_path_2": "" if p2 is None else f"{p2[prior]:.6g}",
                        "p_trans": str(p_trans),
                        "q_variant": f"{q:.6g}",
                        "epcr": f"{epcr:.6g}",
                        "capped": "yes" if epcr < q - 1e-12 else "no",
                        "reads_ok": "yes" if reads_ok else "no",
                        "criteria_1": ",".join(f"{name}:{strength}" for name, strength in bits1),
                        "criteria_2": ",".join(f"{name}:{strength}" for name, strength in bits2),
                        "pile_alt_1": evidence.get("pile_alt_1", ""),
                        "pile_alt_2": evidence.get("pile_alt_2", ""),
                    }
                )
                ranked.append(row_out)
    ranked.sort(
        key=lambda item: (
            -float(item["epcr"]),
            -(int(item.get("pile_alt_1") or 0) + int(item.get("pile_alt_2") or 0)),
            item.get("chrom_1", ""),
            item.get("pos_1", ""),
            item.get("chrom_2", ""),
            item.get("pos_2", ""),
        )
    )
    fields = list(ranked[0]) if ranked else ["epcr"]
    with (root / "ranked_pairs_real.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(ranked)
    with (root / "top10_ranked.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(ranked[:10])
    summary = {
        "hypotheses": len(ranked),
        "capped": sum(item["capped"] == "yes" for item in ranked),
        "reads_checked_in_top": sum(1 for item in ranked[:60] if item["reads_ok"] == "yes"),
        "clinvar_joined": len(clin),
        "top_epcr": [float(item["epcr"]) for item in ranked[:10]],
        "top_genes": [item["gene"] for item in ranked[:10]],
        "top_hypothesis": [item["hypothesis"] for item in ranked[:10]],
    }
    (root / "epcr_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
