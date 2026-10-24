"""Leave-one-out audit of known MVA gene-disease edges in PrimeKG.

Streams edges.csv. Loads nodes.tab only. After hiding a known gene's disease
edges, ranks remaining genes by how many of those diseases they still touch.
The held-out gene cannot score on a hidden direct edge; rediscovery is the
rank it would take among other genes that retain a path of length 1 to the
same diseases. No novel-gene list is written.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
PROJECT = DEST.parent.parent
NODES = PROJECT / "stage" / "T6" / "kg" / "nodes.tab"
EDGES = PROJECT / "stage" / "T6" / "kg" / "edges.csv"
KNOWN = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1")


def load_nodes(path: Path):
    by_index = {}
    gene_index = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            idx = int(row["node_index"].strip('"'))
            kind = row["node_type"].strip('"')
            name = row["node_name"].strip('"')
            by_index[idx] = (kind, name)
            if kind.startswith("gene") and name in KNOWN:
                gene_index[name] = idx
    return by_index, gene_index


def neighbors_of(path: Path, seeds: set[int]) -> dict[int, set[int]]:
    neigh: dict[int, set[int]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            x, y = int(row["x_index"]), int(row["y_index"])
            if x in seeds or y in seeds:
                neigh[x].add(y)
                neigh[y].add(x)
    return neigh


def main() -> None:
    if not NODES.is_file() or not EDGES.is_file():
        raise FileNotFoundError("PrimeKG nodes/edges missing under stage/T6/kg")
    nodes, gene_index = load_nodes(NODES)
    missing = [name for name in KNOWN if name not in gene_index]
    first = neighbors_of(EDGES, set(gene_index.values()))
    disease_ids = set()
    for idx in gene_index.values():
        for nb in first.get(idx, ()):
            kind, _name = nodes.get(nb, ("", ""))
            if "disease" in kind.lower():
                disease_ids.add(nb)
    second = neighbors_of(EDGES, disease_ids)
    dest = DEST / "computed" / "primekg" / "kg_rediscovery.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "gene",
        "gene_index",
        "n_disease_neighbors",
        "hidden_disease_names",
        "other_genes_sharing_any_hidden_disease",
        "best_competitor",
        "best_competitor_shared_diseases",
        "unique_to_those_diseases",
        "note",
    ]
    rows = [header]
    for gene, idx in gene_index.items():
        hidden = [nb for nb in first.get(idx, ()) if nb in disease_ids]
        names = [nodes[i][1] for i in hidden]
        scores: dict[int, int] = {}
        for disease in hidden:
            for other in second.get(disease, ()):
                kind, name = nodes.get(other, ("", ""))
                if not kind.startswith("gene") or other == idx:
                    continue
                scores[other] = scores.get(other, 0) + 1
        competitor = ""
        shared = 0
        if scores:
            other, shared = max(scores.items(), key=lambda item: item[1])
            competitor = nodes[other][1]
        rows.append(
            [
                gene,
                str(idx),
                str(len(hidden)),
                ";".join(names),
                str(len(scores)),
                competitor,
                str(shared),
                "yes" if hidden and not scores else "no",
                "direct gene-disease edges only; embeddings not fit; no novel gene emitted"
                + (f"; missing={missing}" if missing else ""),
            ]
        )
    dest.write_text("\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")
    print(
        f"primekg known={len(gene_index)} missing={missing} diseases={len(disease_ids)}"
    )


if __name__ == "__main__":
    main()
