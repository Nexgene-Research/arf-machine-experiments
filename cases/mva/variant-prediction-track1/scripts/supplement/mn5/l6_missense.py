#!/usr/bin/env python3
"""Score BUB1B p.Asn1002Lys only if a staged structure contains that residue.

Fails closed: writes a coverage row, never invents a stability or ESM number.
"""

from __future__ import annotations

import os

import csv
import sys
from pathlib import Path

RESIDUE = 1002
WT, MUT = "N", "K"
PDB_DIR = Path(os.environ["MVA_BASE"]) / "stage/T6/pdb"
if not PDB_DIR.is_dir():
    PDB_DIR = Path(os.environ["MVA_BASE"]) / "models/pdb"


def chain_sequence(path: Path) -> dict[str, dict[int, str]]:
    aa = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
    }
    chains: dict[str, dict[int, str]] = {}
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            res = line[17:20]
            chain = line[21]
            try:
                seq = int(line[22:26])
            except ValueError:
                continue
            if res in aa:
                chains.setdefault(chain, {})[seq] = aa[res]
    return chains


def main() -> None:
    dest = Path(sys.argv[1])
    pdb_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else PDB_DIR
    rows = []
    pdbs = sorted(pdb_dir.glob("*.pdb")) if pdb_dir.is_dir() else []
    for pdb in pdbs:
        chains = chain_sequence(pdb)
        for chain, seq in chains.items():
            residue = seq.get(RESIDUE)
            rows.append(
                {
                    "structure": pdb.name,
                    "chain": chain,
                    "residue_number": RESIDUE,
                    "residue_in_structure": residue or "",
                    "expected_wt": WT,
                    "mutation": f"{WT}{RESIDUE}{MUT}",
                    "covers_asn1002": "yes" if residue == WT else "no",
                    "thermoMPNN": "",
                    "esm1v": "",
                    "boltz2": "",
                    "note": (
                        "residue present; GPU scorer not invoked in this script"
                        if residue == WT
                        else "residue absent or not asparagine; no score invented"
                    ),
                }
            )
    if not rows:
        rows.append(
            {
                "structure": "",
                "chain": "",
                "residue_number": RESIDUE,
                "residue_in_structure": "",
                "expected_wt": WT,
                "mutation": f"{WT}{RESIDUE}{MUT}",
                "covers_asn1002": "no",
                "thermoMPNN": "",
                "esm1v": "",
                "boltz2": "",
                "note": f"no PDB files under {pdb_dir}",
            }
        )
    with dest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    covered = sum(row["covers_asn1002"] == "yes" for row in rows)
    print(f"structures={len(pdbs)} covering={covered}")


if __name__ == "__main__":
    main()
