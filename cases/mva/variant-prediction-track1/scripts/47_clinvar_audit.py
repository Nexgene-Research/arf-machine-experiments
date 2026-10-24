"""Save direct ClinVar exact-allele review without inventing independent ACMG LR.

Verified on ClinVar VCV000533901.9, accessed 2026-10-16; queried with only
public allele notation/Variation ID. Supplementary MyVariant mirror is older.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/track1/clinvar_exact_audit.tsv"
RECORDS = [
    (
        "chr15:40209701:T:G",
        "VCV000533901.9",
        "P/LP",
        "2",
        "RCV000641226.9 MVA1 P one-star; RCV003324778.1 condition not provided LP one-star",
        "2024-10-09",
        "Labcorp/Invitae;GeneDx",
        "none asserted: published PMID citations relate to BUB1B loss of function, not necessarily this exact allele; quality warning for gnomAD site",
        "https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/",
    ),
    (
        "chr15:40220612:T:G",
        "not found in public MyVariant ClinVar mirror",
        "unresolved",
        "unresolved",
        "no exact submission verified by this targeted mirror query; not proof of absence",
        "",
        "",
        "do not apply PS1/PM3",
        "https://myvariant.info/v1/variant/chr15:g.40220612T%3EG?assembly=hg38&fields=clinvar",
    ),
]


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "grch38_allele",
                "public_record",
                "aggregate_clinsig",
                "aggregate_stars",
                "condition_specific_status",
                "last_evaluated",
                "submitters",
                "criterion_limit",
                "source_url",
            )
        )
        writer.writerows(RECORDS)
    print(f"Audited {len(RECORDS)} public allele queries; PS1/PM3 not added")


if __name__ == "__main__":
    main()
