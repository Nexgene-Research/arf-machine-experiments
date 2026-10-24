"""Verify bounded Phase-I reconstruction against untouched legacy artifacts."""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "results"
TRACK = ROOT / "track1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    legacy = ROOT / "track1_top10_ranked.tsv"
    rebuilt = TRACK / "top10_reconstructed.tsv"
    primary = json.loads((TRACK / "primary_counts_reconstructed.json").read_text())
    funnel = json.loads((TRACK / "funnel_primary_reconstructed.json").read_text())
    exclusions = tsv(ROOT / "exclusions.tsv")
    new_exclusions = tsv(TRACK / "exclusions_primary_reconstructed_v2.tsv")
    old_keys = {(r["gene"], r["hypothesis_class"]): r["n"] for r in exclusions}
    new_keys = {(r["gene"], r["hypothesis_class"]): r["n"] for r in new_exclusions}
    assert legacy.read_bytes() == rebuilt.read_bytes()
    assert len(tsv(rebuilt)) == primary["selected_rows"] == 10
    assert primary["output_sha256"] == digest(rebuilt)
    assert not funnel["legacy_mismatches"]
    assert primary["pair_genes"] == primary["non_tier1_pair_genes"] + 5
    result = {
        "legacy_sha256": digest(legacy),
        "reconstructed_sha256": digest(rebuilt),
        "byte_identical": True,
        "reconstructed_rows": len(tsv(rebuilt)),
        "primary_counts_sha256": digest(TRACK / "primary_counts_reconstructed.json"),
        "primary_script_sha256": digest(
            ROOT.parent / "scripts/37_select_readbacked.py"
        ),
        "job_log_sha256": digest(
            ROOT.parent / "logs/track1_reconstruction_47406904.log"
        ),
        "funnel_rebuilt_fields": len(funnel["reconstructed"]),
        "funnel_field_mismatches": funnel["legacy_mismatches"],
        "funnel_not_rebuilt": funnel["unreconstructed_fields"],
        "legacy_exclusion_rows": len(exclusions),
        "aggregate_exclusion_rows": len(new_exclusions),
        "shared_exclusion_count_differences": {
            f"{gene}/{kind}": {"legacy": old_keys[(gene, kind)], "primary": n}
            for (gene, kind), n in new_keys.items()
            if (gene, kind) in old_keys and n != old_keys[(gene, kind)]
        },
        "legacy_exclusion_keys_missing_from_bounded_reconstruction": [
            f"{gene}/{kind}" for gene, kind in old_keys if (gene, kind) not in new_keys
        ],
        "upstream_job_id": "47397018",
        "reconstruction_job_id": "47406904",
    }
    target = TRACK / "legacy_reconstruction_verification_v2.json"
    if target.exists():
        raise FileExistsError(target)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
