"""Hash computed outputs and append SOURCE records. Does not invent numbers."""

from __future__ import annotations

import hashlib
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
CLAIMS = {
    "computed/phase_span/span.tsv": "10911 bp between the two BUB1B SNVs; WhatsHap unknown",
    "computed/l5_alphagenome/coding_snv_scores.tsv": "AVI and splicing table hits for stop and missense",
    "computed/graph/evidence_graph.json": "Evidence graph from compiled B/C/D tables and reconstructed pairs",
    "computed/dosage_baf/arm_depth.tsv": "Arm depth from 100 kb mosdepth; not a clone fraction",
    "computed/primekg/kg_rediscovery.tsv": "PrimeKG leave-one-out of five known MVA gene-disease links",
    "computed/null/null_runs.tsv": "HPO-weight permutation null plus HG002 truth counts when present",
    "computed/not_available.md": "Claims that cannot be computed from present data",
    "computed/l3_sv/tier1_sixgene_exon_utr.bed": "Six-gene exon/UTR intervals used for MN5 SV and depth filters",
    "computed/l6_missense/asn1002lys.tsv": "PDB coverage of BUBR1 Asn1002; scores blank when uncovered or unrun",
    "computed/l1_consensus/caller_discordance.tsv": "Tier-1 window discordance from three finished callers",
    "computed/l3_sv/tier1_sv.tsv": "Manta/Delly records overlapping the six genes",
    "computed/l3_sv/exon_depth.tsv": "mosdepth exon means for the six-gene BED",
    "computed/dosage_baf/het_baf.tsv": "Per-chromosome heterozygous allele balance",
    "computed/calibration/predictor_lr.tsv": "ClinVar two-star LR for AM/REVEL with BUB1B held out",
    "computed/l5_alphagenome/indel_inference.tsv": "Local AlphaGenome indel attempt after import/parity probe",
    "computed/positive_controls/positive_controls.tsv": "BAMSurgeon recovery of GRCh38-mapped published alleles",
    "computed/l3_sv/exon_depth_by_gene.tsv": "Per-gene min/median mosdepth over the six-gene exon BED",
    "computed/l1_consensus/target_alleles.tsv": "Three-caller support for the two BUB1B coding SNVs",
    "computed/l3_sv/manta_tier1_sv.tsv": "Manta records whose span overlaps a six-gene exon",
    "computed/l3_sv/delly_tier1_sv.tsv": "Delly records whose span overlaps a six-gene exon",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = DEST / "MANIFEST.sha256"
    seen = set()
    if manifest.is_file():
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if line.strip():
                seen.add(line.split()[-1])
    extra = []
    for rel, claim in CLAIMS.items():
        path = DEST / rel
        if not path.is_file():
            continue
        sha = digest(path)
        note = path.with_name(path.name + ".SOURCE.txt")
        note.write_text(
            f"source\t{path.resolve()}\n"
            f"sha256\t{sha}\n"
            f"bytes\t{path.stat().st_size}\n"
            f"claim\t{claim}\n",
            encoding="utf-8",
        )
        if rel not in seen:
            extra.append(f"{sha}  {rel}")
    if extra:
        with manifest.open("a", encoding="utf-8") as handle:
            handle.write("\n".join(extra) + "\n")
    print(f"ledger_extra={len(extra)}")


if __name__ == "__main__":
    main()
