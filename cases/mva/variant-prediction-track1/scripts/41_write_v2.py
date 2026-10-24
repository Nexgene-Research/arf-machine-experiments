"""Apply frozen eligibility and produce one strictly provisional Track 1 CSV.

Inputs are bounded old ranking plus independently read-validated and phased
candidate endpoints. Preserve v1; validate trial against the real GRCh38 FASTA
before publishing. Cannot promote scores or infer trans from two heterozygotes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/track1_top10_ranked.tsv"
ALLELES = ROOT / "results/track1/allele_validation.tsv"
PHASE = ROOT / "results/track1/phase_pair_evidence.tsv"
OUT = ROOT / "sub/Nexgene-AI_mn5-evidence-v2.csv"
CARDS = ROOT / "results/track1/v2_evidence_cards.md"
AUDIT = ROOT / "results/track1/v2_adjudication.tsv"
REF = ROOT / "ref/GRCh38_no_alt.fa"
PROTOCOL = ROOT / "results/track1/v2_adjudication_protocol.md"
LOCKED_HASH = "259fc684603c1b52882e1c749ee3d9d2fa792e2e2fc82abc851bf354dd6d22eb"
HEADER = "proband_id chrom_1 pos_1 ref_1 alt_1 chrom_2 pos_2 ref_2 alt_2 epcr finding_type notes".split()
ESTABLISHED_AR = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1"}


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(row: dict[str, str], n: int) -> str:
    return "-".join(row[f"{part}_{n}"] for part in ("chrom", "pos", "ref", "alt"))


def main() -> None:
    if any(p.exists() for p in (OUT, CARDS, AUDIT)):
        raise FileExistsError("v2 outputs already exist; preserve existing version")
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != LOCKED_HASH:
        raise ValueError("adjudication protocol hash changed")
    original = load(BASE)
    calls = {a["grch38_key"]: a for a in load(ALLELES)}
    phases = {frozenset((r["allele_1"], r["allele_2"])): r for r in load(PHASE)}
    if len(original) != 10:
        raise ValueError("expected bounded v1 evidence table")
    decisions, admitted = [], []
    for row in original:
        gene = row["gene"]
        keys = (key(row, 1), key(row, 2))
        mapped = [calls.get(k) for k in keys]
        reasons = []
        if gene not in ESTABLISHED_AR:
            reasons.append(
                "gene-disease inheritance fit not sufficiently established for classical MVA"
            )
        for n, allele in enumerate(mapped, 1):
            if allele is None:
                reasons.append(f"allele {n} not authenticated in B+C+D output")
                continue
            if allele["ref_check"] != "pass" or allele["source_gt"] != "0/1":
                reasons.append(f"allele {n} REF/heterozygosity not verified")
            if (
                allele["caller_relation"] != "concordant"
                or allele["recall_gt"] != "0/1"
            ):
                reasons.append(
                    f"allele {n} not concordant on separate alignment/caller"
                )
            if (
                min(
                    int(allele["read_alt_MQ60_BQ20"]),
                    int(allele["alt_forward"]),
                    int(allele["alt_reverse"]),
                )
                < 1
            ):
                reasons.append(f"allele {n} lacks stringent bidirectional alt reads")
            af = allele["gnomad_AF"]
            if not af or af == ".":
                reasons.append(f"allele {n} population AF unresolved")
            elif float(af) > 0.01:
                reasons.append(f"allele {n} AF>0.01; no founder exception")
            if int(allele["alt_SA_count"]) or int(allele["alt_XA_count"]):
                reasons.append(f"allele {n} supplementary/XA paralog flag unresolved")
        limiting = row["consequence_2"].lower()
        # For these ten rows, LoF is allele 1 only in BUB1B; a second
        # intron/upstream has no RNA/functional validation and fails eligibility.
        if gene == "BUB1B" and "stop_gained" in row["consequence_1"].lower():
            if "missense" not in limiting:
                reasons.append(
                    "second allele has no established damaging molecular effect"
                )
        else:
            reasons.append("no supported LoF plus plausible coding second-allele pair")
        phase = phases.get(frozenset(keys))
        relation = phase["phase"] if phase else "unknown"
        if relation == "cis":
            reasons.append("cis-proven pair cannot be compound heterozygote")
        if phase is None:
            reasons.append("no regional phase entry")
        if reasons:
            decisions.append((";".join(keys), gene, "excluded", "; ".join(reasons), ""))
            continue
        if relation not in ("unknown", "trans"):
            raise ValueError(f"unrecognized relation {relation}")
        original_epcr = float(row["epcr"])
        if original_epcr <= 0 or original_epcr > 0.1:
            raise ValueError("unexpected legacy capped heuristic")
        # No new independent calibrated LR, and unknown phase/second-allele
        # function remains critical, hence an unmodified 0.1 cap.
        score = min(original_epcr, 0.1)
        notes = (
            f"gene={gene};mechanism=LoF_plus_putative_kinase_missense;"
            f"phase={relation};tie=caller_gene_molecular;"
            "limitation=missense_effect_and_trans_unverified;epcr=heuristic"
        )
        admitted.append(
            {
                "proband_id": "PROBAND01",
                **{
                    f"{part}_{n}": row[f"{part}_{n}"]
                    for n in ("1", "2")
                    for part in ("chrom", "pos", "ref", "alt")
                },
                "epcr": f"{score:.6g}",
                "finding_type": "primary",
                "notes": notes,
            }
        )
        decisions.append(
            (
                ";".join(keys),
                gene,
                "admitted_provisionally",
                "two independently recalled bidirectional alleles; putative kinase missense not experimentally validated; phase unknown",
                f"{score:.6g}",
            )
        )
    if not 1 <= len(admitted) <= 10:
        raise ValueError(
            "no eligible provisional v2 causal scenario; preserve valid v1 fallback"
        )
    admitted.sort(
        key=lambda r: (
            -float(r["epcr"]),
            r["chrom_1"],
            int(r["pos_1"]),
            r["chrom_2"],
            int(r["pos_2"]),
        )
    )
    with tempfile.TemporaryDirectory(
        prefix="track1_v2_trial_", dir=OUT.parent
    ) as dirname:
        trial = Path(dirname) / "trial.csv"
        with trial.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=HEADER, lineterminator="\n")
            writer.writeheader()
            writer.writerows(admitted)
        check = subprocess.run(
            [
                str(ROOT / "env/venv/Scripts/python.exe"),
                str(ROOT / "scripts/validate_submission.py"),
                str(trial),
                str(REF),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        if check.returncode:
            raise RuntimeError(
                f"trial validator exit={check.returncode}: {check.stdout} {check.stderr}"
            )
        if OUT.exists():
            raise FileExistsError(OUT)
        os.replace(trial, OUT)
    with AUDIT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("grch38_pair", "gene", "decision", "reason", "epcr"))
        writer.writerows(decisions)
    winning = admitted[0]
    CARDS.write_text(
        "# v2 evidence cards (one provisionally admitted hypothesis)\n\n"
        f"**BUB1B GRCh38 chr15:{winning['pos_1']}:{winning['ref_1']}>{winning['alt_1']} and chr15:{winning['pos_2']}:{winning['ref_2']}>{winning['alt_2']}**. First allele MANE NM_001211.6:c.2210T>G, p.Leu737Ter; second allele MANE NM_001211.6:c.3006T>G, p.Asn1002Lys, within annotated kinase-domain interval, effect unverified. Both 0/1 in the supplied callset and independently recalled from a separately aligned BAM of the **same individual**, with bidirectional MQ≥60/BQ≥20 alt reads (`results/track1/allele_validation.tsv`); no SA/XA alt flags in these two records. gnomAD AF values as archived in that table, with a ClinVar warning about site-level frequency quality for the stop. Physical phase is **unknown**, distinct regional blocks (`results/track1/phase_pair_evidence.tsv`); no trans claim. Molecular mechanism hypothesis: BUB1B loss plus a potential hypomorph, not proven; no allele-specific RNA or protein assay. Gene–disease support: Hanks et al. DOI 10.1038/ng1449 and Ochiai et al. DOI 10.1073/pnas.1317008111 (latter abstract in this review). Phenotype-only BUB1B rank and gaps: `results/angle_b_phenotype.tsv`; no absent feature is asserted. Exact stop is ClinVar VCV000533901.9 aggregate P/LP two stars, but disease-specific P one-star and no independently justified PS1 or PM3 (`results/track1/clinvar_exact_audit.tsv`). Published cases use other alleles; a MAD1L1-driven tumor mechanism is a competing explanation if a valid MAD1L1 genotype is found. Strongest objection: p.Asn1002Lys may be neutral and the two BUB1B alleles may be cis. Confirm via trio/long reads, RNA/protein, independent functional assays and orthogonal SV/CNV analysis. ePCR remains capped at 0.1 heuristic.\n\n"
        "All excluded v1 hypotheses and allele-specific reasons: `results/track1/v2_adjudication.tsv`; Phase L matrix: `results/track1/evidence_convergence.tsv`. Gene-level literature does not establish allele-level function.\n",
        encoding="utf-8",
    )
    print(
        f"OK: admitted={len(admitted)} excluded={len(decisions) - len(admitted)}; trial validator {check.stdout.strip()}; final {OUT}"
    )


if __name__ == "__main__":
    main()
