"""Export auditable evidence links from frozen A and Stage 32.

No diagram renderer is available; figure placeholder is a caption, not a
fabricated illustration. Edges remain provisional until real read validation.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path, required):
    if path.stat().st_size > 50_000_000:
        raise ValueError("Bounded graph input exceeds laptop file limit")
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(
                f"{path}: missing {set(required) - set(reader.fieldnames or [])}"
            )
        return list(reader)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", type=Path, default=ROOT / "results/ranked_pairs.tsv")
    ap.add_argument("--b", type=Path, default=ROOT / "results/angle_b_phenotype.tsv")
    ap.add_argument("--c", type=Path, default=ROOT / "results/angle_c_mechanism.tsv")
    ap.add_argument("--d", type=Path, default=ROOT / "results/literature_cases.tsv")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results")
    args = ap.parse_args()
    ranked = read(args.ranked, ("pair_id", "gene", "variant_1", "epcr"))
    b = {
        r["gene"]: r for r in read(args.b, ("gene", "hypothesis", "evidence_reference"))
    }
    c = {
        r["gene"]: r
        for r in read(args.c, ("gene", "pathway_node", "primary_source_ids"))
    }
    d = {}
    for r in read(args.d, ("gene", "source_doi_pmid", "case_id")):
        d.setdefault(r["gene"], []).append(r)
    edges, nodes = [], {}

    def node(kind, name):
        nid = f"{kind}:{name}"
        nodes[nid] = {"id": nid, "kind": kind, "label": name}
        return nid

    def edge(a, z, relation, source, strength, status):
        edges.append(
            dict(
                source_node=a,
                target_node=z,
                relation=relation,
                source=source,
                access_date="2026-10-16",
                strength=strength,
                status=status,
            )
        )

    for row in ranked:
        gene = node("gene", row["gene"])
        hypothesis = node("hypothesis", row["pair_id"])
        edge(
            hypothesis,
            gene,
            "proposed_gene",
            str(args.ranked),
            "unvalidated_case_link",
            row.get("epcr_status", "provisional"),
        )
        for i in (1, 2):
            key = row.get(f"variant_{i}", "")
            if not key:
                continue
            allele = node("variant", key)
            effect = node("molecular_effect", f"{key}:unverified")
            edge(
                allele,
                effect,
                "effect_requires_annotation",
                str(args.ranked),
                "unknown",
                "unverified",
            )
            edge(effect, gene, "gene_link", str(args.ranked), "unknown", "unverified")
            edge(
                allele,
                hypothesis,
                "proposed_allele",
                str(args.ranked),
                "computational",
                row.get("evidence_status", "provisional"),
            )
        if row["gene"] in c:
            mechanism = c[row["gene"]]
            pathway = node("pathway", mechanism["pathway_node"])
            edge(
                gene,
                pathway,
                "mechanism",
                mechanism["primary_source_ids"],
                mechanism.get("link_strength", "qualitative"),
                "literature",
            )
        if row["gene"] in b:
            phenotype = b[row["gene"]]
            clinical = node("clinical_feature", phenotype["hypothesis"])
            edge(
                gene,
                clinical,
                "phenotype_overlap",
                phenotype["evidence_reference"],
                "qualitative",
                "literature",
            )
        for lit in d.get(row["gene"], []):
            paper = node("paper", lit["source_doi_pmid"])
            edge(
                paper,
                gene,
                "published_case",
                lit["source_doi_pmid"],
                "primary_source",
                "case_context_not_patient_genotype",
            )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "evidence_graph.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        w = csv.DictWriter(
            stream,
            fieldnames=[
                "source_node",
                "target_node",
                "relation",
                "source",
                "access_date",
                "strength",
                "status",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(edges)
    manifest = {
        "nodes": sorted(nodes.values(), key=lambda n: n["id"]),
        "edges": edges,
        "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.ranked, args.b, args.c, args.d)
        },
        "figure": "figure_unavailable; requested caption in report; no diagram generated",
    }
    (args.out_dir / "evidence_graph.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"graph nodes={len(nodes)} edges={len(edges)}")


if __name__ == "__main__":
    main()
