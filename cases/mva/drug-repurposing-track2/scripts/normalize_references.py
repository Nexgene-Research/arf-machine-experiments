"""Recompute conditional reference-allele context; never normalize a patient genotype."""

import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import DATA, REFERENCE


def main():
    base = REFERENCE
    gb = (base / "NM_001211.6.gb").read_text(encoding="utf-8")
    seq = "".join(re.findall("[acgt]+", gb.split("ORIGIN")[1].split("//")[0])).upper()
    cds_match = re.search(r"\n\s+CDS\s+(\d+)\.\.(\d+)", gb)
    assert cds_match is not None, "Expected single contiguous RefSeq CDS annotation"
    cds_start, cds_end = map(int, cds_match.groups())
    cds = seq[cds_start - 1 : cds_end]
    tx = json.loads((base / "bub1b_transcript.json").read_text())
    exons = sorted(tx["Exon"], key=lambda e: e["start"])
    assert tx["strand"] == 1
    offset = 0
    intervals = []
    for i, exon in enumerate(exons, 1):
        length = exon["end"] - exon["start"] + 1
        intervals.append((i, offset + 1, offset + length, exon["start"], exon["end"]))
        offset += length
    assert offset == len(seq), (
        "RefSeq and Ensembl transcript lengths differ; inspect alignment before interpreting"
    )
    variants = []
    # Exact public ClinVar records, not reverse translation from patient protein labels.
    for label, c, ref, alt, genomic, vcv in [
        ("p.Leu737Ter", 2210, "T", "G", 40209701, "VCV000533901"),
        ("p.Asn1002Lys", 3006, "T", "A", 40220612, "VCV004600147"),
    ]:
        assert cds[c - 1] == ref
        pos = cds_start + c - 1
        exon = next(e for e in intervals if e[1] <= pos <= e[2])
        assert exon[3] + pos - exon[1] == genomic
        codon_start = ((c - 1) // 3) * 3
        wt = cds[codon_start : codon_start + 3]
        mutant = list(wt)
        mutant[(c - 1) % 3] = alt
        mutant = "".join(mutant)
        flanking = (
            cds[max(0, codon_start - 12) : codon_start]
            + "["
            + mutant
            + "]"
            + cds[codon_start + 3 : codon_start + 15]
        )
        variants.append(
            dict(
                public_record=vcv,
                protein=label,
                transcript="NM_001211.6",
                protein_accession="NP_001202.5",
                cdna=f"c.{c}{ref}>{alt}",
                genomic=f"NC_000015.10:g.{genomic}{ref}>{alt}",
                exon_number_in_mane=exon[0],
                transcript_position_1based=pos,
                distance_to_exon_end_nt=exon[2] - pos,
                distance_to_final_exon_junction_nt=intervals[-2][2] - pos,
                reference_codon=wt,
                public_allele_codon=mutant,
                context_5prime_to_3prime_DNA=flanking,
                patient_identity="Unverified: requires original laboratory report and phase",
            )
        )
    result = dict(
        reference_only=True,
        ensembl=f"{tx['id']}.{tx['version']}",
        protein_length_aa=tx["Translation"]["length"],
        exon_count=len(exons),
        cds_start_1based=cds_start,
        cds_end_1based=cds_end,
        variants=variants,
        input_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in base.iterdir()
            if p.is_file()
        },
    )
    DATA.mkdir(exist_ok=True)
    (DATA / "reference_context.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
