"""Map collected MN5 TSVs into computed/ destinations and merge null rows."""

from __future__ import annotations

import csv
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
SRC = DEST / "computed" / "mn5_collect"
MAP = {
    "caller_discordance.tsv": DEST / "computed" / "l1_consensus" / "caller_discordance.tsv",
    "het_baf.tsv": DEST / "computed" / "dosage_baf" / "het_baf.tsv",
    "exon_depth.tsv": DEST / "computed" / "l3_sv" / "exon_depth.tsv",
    "manta_tier1_sv.tsv": DEST / "computed" / "l3_sv" / "manta_tier1_sv.tsv",
    "delly_tier1_sv.tsv": DEST / "computed" / "l3_sv" / "delly_tier1_sv.tsv",
    "asn1002lys.tsv": DEST / "computed" / "l6_missense" / "asn1002lys.tsv",
    "indel_inference.tsv": DEST / "computed" / "l5_alphagenome" / "indel_inference.tsv",
    "predictor_lr.tsv": DEST / "computed" / "calibration" / "predictor_lr.tsv",
    "positive_controls.tsv": DEST / "computed" / "positive_controls" / "positive_controls.tsv",
    "l2_locus.tsv": DEST / "computed" / "l1_consensus" / "l2_locus.tsv",
}


def main() -> None:
    if not SRC.is_dir():
        print("no collect directory")
        return
    for name, dest in MAP.items():
        src = SRC / name
        if not src.is_file():
            continue
        if name == "asn1002lys.tsv" and dest.is_file() and dest.stat().st_size > src.stat().st_size:
            print(f"kept local {dest.relative_to(DEST)}")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(src.read_bytes())
        print(f"placed {name} -> {dest.relative_to(DEST)}")
    sv_dest = DEST / "computed" / "l3_sv" / "tier1_sv.tsv"
    parts = []
    for name in ("manta_tier1_sv.tsv", "delly_tier1_sv.tsv"):
        path = DEST / "computed" / "l3_sv" / name
        if path.is_file():
            parts.append(path)
    if parts:
        header = None
        rows = []
        for path in parts:
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                header = reader.fieldnames
                rows.extend(reader)
        with sv_dest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    hg = SRC / "null_hg002.tsv"
    null = DEST / "computed" / "null" / "null_runs.tsv"
    if hg.is_file() and null.is_file():
        existing = null.read_text(encoding="utf-8")
        extra = hg.read_text(encoding="utf-8")
        extra_lines = [line for line in extra.splitlines() if line and not line.startswith("arm\t")]
        if extra_lines and "hg002_truth_vcf" not in existing:
            null.write_text(existing.rstrip() + "\n" + "\n".join(extra_lines) + "\n", encoding="utf-8")
            print("merged hg002 null into null_runs.tsv")


if __name__ == "__main__":
    main()
