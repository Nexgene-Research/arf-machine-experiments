"""Archive publicly verified MANE crosswalks for two GRCh38 alleles.

Manual source verification: Clinical Genome Allele Registry per-allele JSON,
retrieved 2026-10-16. No case context was sent to that public service.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/track1/variant_crosswalk.tsv"
ROWS = [
    (
        "chr15:40209701:T:G",
        "NM_001211.6:c.2210T>G",
        "NP_001202.5:p.Leu737Ter",
        "MANE Select 1.3",
        "https://reg.genome.network/allele?hgvs=NC_000015.10:g.40209701T%3EG",
        "stop gained; exon/NMD not established",
    ),
    (
        "chr15:40220612:T:G",
        "NM_001211.6:c.3006T>G",
        "NP_001202.5:p.Asn1002Lys",
        "MANE Select 1.3",
        "https://reg.genome.network/allele?hgvs=NC_000015.10:g.40220612T%3EG",
        "within UniProt O60566 canonical kinase aa766–1050; interface/effect not established",
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
                "mane_hgvs",
                "mane_protein",
                "mane_version",
                "public_source",
                "interpretation_limit",
            )
        )
        writer.writerows(ROWS)
    print(f"Archived {len(ROWS)} independently looked-up MANE crosswalks")


if __name__ == "__main__":
    main()
