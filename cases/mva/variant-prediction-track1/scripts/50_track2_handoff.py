"""Produce bounded research-only Track 2 mechanism handoff from final CSV.

No treatment recommendation. Interpretations remain hypotheses until causal
alleles, phase, functional impact and pediatric safety are independently tested.
"""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "sub/Nexgene-AI_mn5-evidence-v2.csv"
OUT = ROOT / "results"
SOURCES = {
    "BUB1B": "10.1038/ng1449;10.1073/pnas.1317008111;ClinVar:VCV000533901.9;UniProt:O60566",
    "BUB1": "10.1126/sciadv.abk0114;UniProt:O43683",
    "TRIP13": "10.1038/ng.3883",
    "CEP57": "10.1038/ng.822",
    "CENATAC": "10.15252/embj.2020106536",
    "MAD1L1": "10.1126/sciadv.abq5914",
}
FUNCTION = {
    "BUB1B": (
        "candidate biallelic loss plus possible missense hypomorph; not proven",
        "spindle assembly checkpoint and mitotic checkpoint complex",
        "premature chromosome segregation / mosaic whole-chromosome imbalance if functional loss verified",
        "BUBR1 protein abundance, checkpoint arrest after spindle perturbation, live-cell chromosome-segregation fidelity",
        "MAD2–CDC20 checkpoint signalling; BUB1/MAD1L1/TRIP13 network",
    ),
    "BUB1": (
        "candidate biallelic splicing/regulatory loss; no functional evidence for listed introns",
        "kinetochore spindle signalling and cohesion",
        "chromosome segregation / cohesion changes only if allele effect verified",
        "BUB1 RNA splice assay, BUB1 protein, chromosome segregation and checkpoint assays",
        "BUB1B/MAD1L1 downstream checkpoint network",
    ),
    "TRIP13": (
        "candidate biallelic loss; allele mechanism must be verified",
        "p31comet–TRIP13 ATPase checkpoint recycling",
        "spindle-checkpoint release defect or missegregation",
        "TRIP13 ATPase/protein expression and chromosome-segregation assays",
        "p31comet–MAD2; MCC disassembly",
    ),
    "CEP57": (
        "candidate biallelic loss; allele mechanism must be verified",
        "centrosomal microtubule anchoring",
        "microtubule attachment and chromosome imbalance",
        "centrosome localization and mitotic karyotyping",
        "centrosomal microtubule attachment nodes",
    ),
    "CENATAC": (
        "candidate biallelic minor-spliceosome impairment; RNA confirmation needed",
        "U12 minor-intron spliceosome",
        "minor-intron retention and congression defects",
        "RNA minor-intron retention panel and wild-type rescue",
        "minor-intron-dependent mitotic transcripts",
    ),
    "MAD1L1": (
        "candidate biallelic loss; allele mechanism must be verified",
        "kinetochore MAD1–MAD2 checkpoint",
        "checkpoint escape and mosaic chromosomal imbalance",
        "MAD1/MAD2 protein recruitment and checkpoint arrest",
        "MAD2–CDC20 checkpoint signaling",
    ),
}


def save(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    if path.exists():
        raise FileExistsError(path)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def main() -> None:
    if not (OUT / "track1/final_release_receipt.json").is_file():
        raise ValueError("Track 1 freeze receipt required before Track 2 handoff")
    import json

    receipt = json.loads(
        (OUT / "track1/final_release_receipt.json").read_text(encoding="utf-8")
    )
    if not str(receipt.get("scientific_freeze_status", "")).startswith("NO-GO"):
        raise ValueError("this research-only handoff is scoped to the recorded no-go")
    with CSV.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    graph, targets, gaps, lines = [], [], [], []
    for ix, row in enumerate(rows, 1):
        notes = dict(
            token.split("=", 1) for token in row["notes"].split(";") if "=" in token
        )
        gene = notes.get("gene", "")
        if gene not in FUNCTION:
            raise ValueError(f"no primary mechanism summary for {gene}")
        variants = ";".join(
            ":".join(row[f"{part}_{n}"] for part in ("chrom", "pos", "ref", "alt"))
            for n in ("1", "2")
            if row[f"chrom_{n}"]
        )
        mechanism, pathway, cellular, readouts, network = FUNCTION[gene]
        phase = notes.get("phase", "unknown")
        confidence = "NOT frozen-eligible; exploratory only; phase=" + phase
        domain = (
            "MANE BUB1B p.Asn1002Lys in kinase-domain interval; effect unknown"
            if "40220612" in variants
            else "site-specific domain or effect unresolved"
        )
        direction = "restore/stabilize/compensate impaired function only if loss confirmed; avoid inhibition of a deficient checkpoint"
        risk = "pediatric developmental/toxicology evaluation required; mitotic checkpoint perturbation may increase aneuploidy or tumor risk; no medicine proposed"
        graph.append(
            [
                str(ix),
                gene,
                variants,
                confidence,
                mechanism,
                domain,
                pathway,
                cellular,
                network,
                SOURCES[gene],
                "gene mechanism established; proposed allelic effect exploratory",
            ]
        )
        targets.append(
            [
                str(ix),
                gene,
                network,
                direction,
                readouts,
                risk,
                SOURCES[gene],
                "exploratory research nodes, not medications",
            ]
        )
        gaps.append(
            [
                str(ix),
                gene,
                variants,
                "independent genotype and phase; RNA/protein effect of second allele; cellular rescue; tumor versus constitutional phenotype; pediatric safety",
                "long-read/trio phasing; RNA assays; protein/localization; isogenic rescue; orthogonal karyotype",
                "no treatment recommendation",
            ]
        )
        lines.append(
            f"{ix}. **{gene}** — {variants} (GRCh38); {confidence}. Predicted mechanism: {mechanism}. Domain: {domain}. Pathway: {pathway}. Cellular readouts: {readouts}. Intervention *direction to investigate* only after validation: {direction}. Safety: {risk}. Sources: {SOURCES[gene]}. Open causal links: {gaps[-1][3]}."
        )
    save(
        OUT / "track2_mechanism_graph.tsv",
        [
            "scenario",
            "gene",
            "normalized_grch38_variants",
            "confidence",
            "putative_mechanism",
            "domain",
            "pathway",
            "cellular_phenotype",
            "network_nodes",
            "sources",
            "evidence_status",
        ],
        graph,
    )
    save(
        OUT / "track2_targets.tsv",
        [
            "scenario",
            "gene",
            "research_nodes",
            "hypothesis_direction",
            "rescue_readouts",
            "contraindications_pediatric_safety",
            "sources",
            "status",
        ],
        targets,
    )
    save(
        OUT / "track2_evidence_gaps.tsv",
        [
            "scenario",
            "gene",
            "normalized_grch38_variants",
            "open_links",
            "discriminating_assays",
            "caution",
        ],
        gaps,
    )
    report = OUT / "track2_drug_repurposing_handoff.md"
    if report.exists():
        raise FileExistsError(report)
    report.write_text(
        "# Track 2 mechanism handoff — research only\n\n**There are no surviving frozen-eligible causal scenarios.** The one row below is a rejected, exploratory allele-pair hypothesis preserved solely to define discriminating experiments, not a validated diagnosis or clinical-confidence value. No treatment is recommended. DrugCentral, ChEMBL, Open Targets, ClinicalTrials.gov and DGIdb are **future** research sources, not queried to support a medicine here.\n\n"
        + "\n\n".join(lines)
        + "\n\nSee `results/track2_mechanism_graph.tsv`, `track2_targets.tsv`, and `track2_evidence_gaps.tsv`. Established gene-pathway findings are separated from exploratory case-specific mechanism hypotheses in the tables.\n",
        encoding="utf-8",
    )
    print(
        f"Track 2 handoff: zero eligible scenarios, {len(rows)} exploratory audit rows, no drugs"
    )


if __name__ == "__main__":
    main()
