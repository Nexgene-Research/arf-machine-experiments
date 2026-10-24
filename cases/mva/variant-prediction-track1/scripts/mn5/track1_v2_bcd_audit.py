#!/usr/bin/env python3
"""End-to-end audit of the immutable Phase B+C+D bounded aggregate tables."""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results/track1"
RECEIPT = R / "bcd_47403447_receipt.json"


def read(name):
    with (R / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for name, meta in receipt["outputs"].items():
        f = Path(meta["path"])
        assert f.is_file() and f.stat().st_size == meta["bytes"]
        assert hashlib.sha256(f.read_bytes()).hexdigest() == meta["sha256"]
    calls = read("allele_validation.tsv")
    pairs = read("phase_pair_evidence.tsv")
    recall = read("caller_concordance.tsv")
    matrix = read("tier1_evidence_matrix.tsv")
    exons = read("second_hit_screen.tsv")
    blocks = read("phase_blocks_summary.tsv")
    keys = read("rare_noncoding_keys.tsv")
    bykey = {row["grch38_key"]: row for row in calls}
    ends = sorted({x[field] for x in pairs for field in ("allele_1", "allele_2")})
    missing = [x for x in ends if x not in bykey]
    concordant = [x for x in recall if x["caller_relation"] == "concordant"]
    q = []
    for key in (
        "chr15-40209701-T-G",
        "chr15-40216470-A-G",
        "chr15-40220612-T-G",
        "chr2-110668376-C-T",
        "chr2-110671177-G-A",
    ):
        x = bykey[key]
        q.append(
            f"- `{key}`: GT {x['source_gt']} / recall {x['recall_gt']}; MQ60/BQ20 alt {x['read_alt_MQ60_BQ20']} (forward {x['alt_forward']}, reverse {x['alt_reverse']}); AF {x['gnomad_AF']}; {x['MANE_consequence']}.\n"
        )
    path = R / "bcd_47403447_audit.md"
    assert not path.exists(), f"refuse overwriting {path}"
    path.write_text(
        "# Independent B+C+D aggregate integrity and scope audit\n\n"
        f"Slurm {receipt['state'].strip()}; frozen protocol SHA-256 `{receipt['protocol_sha256']}`. "
        "All ten retrieved files match remote size and SHA-256 in the receipt.\n\n"
        f"- Allele rows {len(calls)}; prespecified bounded candidates {len(recall)}, "
        f"concordant independent BAM bcftools calls {len(concordant)}, "
        f"discordant or not called {len(recall) - len(concordant)}. This is the same sample, not replication.\n"
        f"- WhatsHap blocks {len(blocks)} from six Tier-1 loci; candidate pairs {len(pairs)}: "
        f"trans {sum(x['phase'] == 'trans' for x in pairs)}, cis {sum(x['phase'] == 'cis' for x in pairs)}, "
        f"unknown {sum(x['phase'] == 'unknown' for x in pairs)}. No phase may be inferred from separate 0/1 calls.\n"
        f"- Gene rows {len(matrix)}; exon/UTR depth rows {len(exons)}; noncoding SNV keys {len(keys)}. "
        "Depth and anomalous read counts are screens, not confirmed CNV/SV.\n"
        f"- Unique pair alleles {len(ends)}; alleles not authenticated by this output "
        f"{len(missing)}: {', '.join('`' + k + '`' for k in missing) if missing else 'none'}. "
        "Do not claim authentication of these missing endpoints.\n\n"
        "## Principal observations (reported directly from the pipeline table)\n\n"
        + "".join(q)
        + "\n"
        "## Scope cautions\n\n"
        "The rare-gene count also includes alleles with missing gnomAD AF (missing is not zero frequency); "
        "re-check AF before assigning rarity. Short-read blocks do not resolve any tested pair. "
        "WhatsHap VCF does not expose informative-fragment counts or direct versus chained support; "
        "these fields remain unknown. Read flags and exon depth need orthogonal follow-up; "
        "genotype discovery outside the six bounded gene windows was not performed. "
        "The BUB1B upstream pair endpoint listed above was not in the bounded allele-read table.\n",
        encoding="utf-8",
    )
    print(path)
    print(
        "caller_concordant",
        len(concordant),
        "of",
        len(recall),
        "pairs_unknown",
        sum(x["phase"] == "unknown" for x in pairs),
        "missing_pair_alleles",
        missing,
    )


if __name__ == "__main__":
    main()
