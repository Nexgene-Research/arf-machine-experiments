#!/usr/bin/env python3
"""Local AlphaGenome probe: import, optional SNV parity, then queued indels.

Never invents a score. Empty indel rows stay empty when the package or
parity check fails.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


PARITY_SNV = {
    "chrom": "chr15",
    "pos": "40209701",
    "ref": "T",
    "alt": "G",
    "tabulated_raw": 1.873,
}


def main() -> None:
    dest, status_path = Path(sys.argv[1]), Path(sys.argv[2])
    queue = Path(sys.argv[3]) if len(sys.argv) > 3 else None
    status = {
        "importable": False,
        "parity": None,
        "indels": "not_run",
        "note": "",
    }
    rows = [
        [
            "variant",
            "chrom",
            "pos",
            "ref",
            "alt",
            "source",
            "raw_score",
            "note",
        ]
    ]
    try:
        import alphagenome  # type: ignore

        status["importable"] = True
        status["version"] = getattr(alphagenome, "__version__", "unknown")
    except Exception as exc:
        status["error"] = f"{type(exc).__name__}: {str(exc)[:300]}"
        status["note"] = "research package not importable; weights may still be staged"
        rows.append(
            [
                "queued_BUB1B_indels",
                ".",
                ".",
                ".",
                ".",
                "local_alphagenome",
                "",
                "package not importable; no score invented",
            ]
        )
        dest.write_text("\t".join(rows[0]) + "\n" + "\n".join("\t".join(r) for r in rows[1:]) + "\n")
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        print(status)
        return

    status["parity"] = "not_run_no_local_scorer"
    status["note"] = (
        "alphagenome imports; this environment has no documented local "
        "variant-scoring entrypoint. Queued indels stay empty."
    )
    status["indels"] = "held"
    if queue and queue.is_file():
        with queue.open(encoding="utf-8", newline="") as handle:
            for rec in csv.DictReader(handle, delimiter="\t"):
                rows.append(
                    [
                        rec.get("gene", "BUB1B") + "_indel",
                        rec.get("chr", rec.get("chrom", ".")),
                        rec.get("pos", "."),
                        rec.get("ref", "."),
                        rec.get("alt", "."),
                        "local_alphagenome",
                        "",
                        "held until a local scorer and SNV parity exist",
                    ]
                )
    dest.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8"
    )
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(status)


if __name__ == "__main__":
    main()
