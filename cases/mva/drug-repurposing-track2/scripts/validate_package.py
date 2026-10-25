"""Independent checks on final static exports, not validation of clinical efficacy."""

import csv
import hashlib
import json
import platform
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import DATA, DOCS, ROOT

EXCLUDED_DIRS = {
    "__pycache__",
    ".ruff_cache",
    ".git",
    "assets",
    "pdf",
}


def read_csv(rel):
    path = ROOT / rel
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = list(reader)
    assert columns and rows, rel
    assert len(set(columns)) == len(columns), rel
    for row in rows:
        assert None not in row and all(
            v is not None and v.strip() for v in row.values()
        ), (rel, row)
        assert not any(
            re.fullmatch(
                r"(?i)(todo|tbd|placeholder|lorem ipsum|dummy|nan|null)", v.strip()
            )
            for v in row.values()
        ), rel
    return columns, rows


def main():
    required = [
        "docs/report.md",
        "docs/appendices.md",
        "docs/claim_audit.md",
        "docs/DELIVERABLES.md",
        "data/nodes.csv",
        "data/edges.csv",
        "data/evidence.csv",
        "data/paths.csv",
        "data/graph.md",
        "data/candidates.csv",
        "data/rejections.csv",
        "data/search_log.csv",
        "data/source_coverage.csv",
        "data/datasets.csv",
        "data/research_tools.csv",
        "data/phenotype_terms.csv",
        "data/reference_context.json",
        "data/population_annotations.json",
        "data/package_summary.json",
        "curation/evidence_overrides.json",
        "curation/curation.json",
        "curation/bridge_updates.json",
        "curation/graph_curation.json",
        "curation/candidate_curation.json",
        "curation/search_curation.json",
        "README.md",
        "environment.txt",
        "scripts/normalize_references.py",
        "scripts/build_package.py",
        "scripts/validate_package.py",
        "scripts/package_results.py",
        "scripts/paths.py",
        "scripts/render_pdf.py",
    ]
    for name in required:
        assert (ROOT / name).is_file() and (ROOT / name).stat().st_size > 0, name
    tables = {
        Path(name).name: read_csv(name) for name in required if name.endswith(".csv")
    }
    nodes = tables["nodes.csv"][1]
    edges = tables["edges.csv"][1]
    evidence = tables["evidence.csv"][1]
    paths = tables["paths.csv"][1]
    candidates = tables["candidates.csv"][1]
    ids = {}
    for name, rows, key in [
        ("nodes", nodes, "node_id"),
        ("edges", edges, "edge_id"),
        ("evidence", evidence, "evidence_id"),
        ("paths", paths, "path_id"),
    ]:
        seq = [r[key] for r in rows]
        assert len(seq) == len(set(seq)), name
        ids[name] = set(seq)
    node_types = {r["node_type"] for r in nodes}
    assert {
        "variant",
        "allele",
        "transcript",
        "protein_state_or_function",
        "complex_or_interaction",
        "biological_process",
        "cell_phenotype",
        "tissue_or_clinical_phenotype",
        "drug",
        "metabolite",
        "target",
        "safety_outcome",
        "experimental_model",
        "evidence_source",
    } <= node_types
    assert {"ABUND", "FUNC", "NORMALMODEL", "CANCERMODEL"} <= ids["nodes"]
    for edge in edges:
        assert edge["source"] in ids["nodes"] and edge["target"] in ids["nodes"], edge
        assert set(edge["evidence_ids"].split(";")) <= ids["evidence"], edge
        assert edge["evidence_state"] in {
            "measured",
            "inferred",
            "hypothesized",
            "contradicted",
        }, edge
    for row in evidence:
        assert row["url"].startswith("https://"), row["evidence_id"]
        for field in [
            "rigor",
            "causal_directness",
            "variant_model_relevance",
            "exposure_feasibility",
            "reproducibility",
            "clinical_context_relevance",
        ]:
            assert row[field] in {
                "high",
                "moderate",
                "low",
                "unknown",
                "not_applicable",
                "conditional",
            }, (row["evidence_id"], field)
    em = {r["edge_id"]: r for r in edges}
    for path in paths:
        es = [em[x] for x in path["edge_ids"].split(";")]
        assert all(a["target"] == b["source"] for a, b in zip(es, es[1:])), path[
            "path_id"
        ]
        assert int(path["hop_count"]) == len(es)
        assert int(path["unproven_links"]) == sum(
            e["evidence_state"] != "measured" for e in es
        )
        transfers = set(
            t for e in es for t in e["transfer_assumptions"].split(";") if ":" in t
        )
        assert int(path["context_transfer_count"]) == len(transfers)
        for prefix, field in [
            ("species:", "species_transfers"),
            ("tissue:", "tissue_transfers"),
            ("genotype:", "genotype_transfers"),
        ]:
            assert int(path[field]) == sum(t.startswith(prefix) for t in transfers)
    needed = [
        "Candidate",
        "Approval status/date",
        "Objective",
        "Axis",
        "Mechanistic path",
        "Variant condition",
        "Evidence anchors",
        "Counterevidence",
        "Weakest link",
        "Exposure",
        "Key risks",
        "Tier",
        "Decisive test",
        "Kill criterion",
    ]
    assert set(needed) <= set(tables["candidates.csv"][0])
    for c in candidates:
        assert c["Tier"] in {"1", "2", "3", "X"} and c["Objective"] in {
            "A",
            "B",
            "C",
            "D",
            "E",
        }
        assert set(c["Evidence anchors"].split(";")) <= ids["evidence"]
    assert not any(c["Tier"] == "1" for c in candidates)
    assert all(c["Objective"] in {"C", "E"} for c in candidates if c["Tier"] == "2")
    report = (DOCS / "report.md").read_text(encoding="utf-8")
    appendix = (DOCS / "appendices.md").read_text(encoding="utf-8")
    for letter in "ABCDEFGHI":
        assert "### Axis " + letter + ":" in appendix, letter
    for n in range(1, 14):
        assert f"| {n} |" in appendix, n
    for eid in set(re.findall(r"\bE\d{2}\b", report + appendix)):
        assert eid in ids["evidence"], eid
    markdown_files = [
        DOCS / "report.md",
        DOCS / "appendices.md",
        ROOT / "README.md",
    ]
    local_link_targets = []
    for md_path in markdown_files:
        text = md_path.read_text(encoding="utf-8")
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if re.match(r"https?://", target):
                continue
            resolved = (md_path.parent / target).resolve()
            local_link_targets.append(resolved)
            rel = str(resolved.relative_to(ROOT.resolve())).replace("\\", "/")
            if rel != "data/validation.json":
                assert resolved.is_file(), f"{md_path.name} -> {target}"
    r = json.loads((DATA / "reference_context.json").read_text())
    assert r["reference_only"] is True
    assert (
        r["variants"][0]["reference_codon"] == "TTA"
        and r["variants"][0]["public_allele_codon"] == "TGA"
    )
    assert (
        r["variants"][1]["cdna"] == "c.3006T>A"
        and "Unverified" in r["variants"][1]["patient_identity"]
    )
    summary = json.loads((DATA / "package_summary.json").read_text())
    assert (
        summary["evidence_claims"] == len(evidence)
        and summary["nodes"] == len(nodes)
        and summary["edges"] == len(edges)
    )
    assert summary["paths"] == len(paths) and summary["candidates"] == len(candidates)
    assert summary["tiers"] == dict(Counter(c["Tier"] for c in candidates))
    files = [
        p
        for p in ROOT.rglob("*")
        if p.is_file()
        and p.name not in {"validation.json", "research_package.zip"}
        and not p.name.endswith("_err.txt")
        and not EXCLUDED_DIRS & set(p.relative_to(ROOT).parts)
    ]
    result = {
        "status": "passed",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "validation_scope": "Static public-evidence package only; no patient biology, clinical benefit or efficacy validated",
        "summary": summary,
        "checks": [
            "required files and nonempty CSV schemas",
            "no empty/dummy-value cells in CSV exports",
            "unique identifiers and referential integrity",
            "graph states/types and continuous paths",
            "independent link/transfer counts",
            "candidate evidence and tier rules",
            "all axes A-I and 13 decision rows",
            "report evidence IDs and local links",
            "conditional reference codon sanity check",
            "final file hashes",
        ],
        "files": {
            str(p.relative_to(ROOT)).replace("\\", "/"): {
                "bytes": p.stat().st_size,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
            for p in sorted(files)
        },
    }
    (DATA / "validation.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    for target in local_link_targets:
        assert target.is_file(), target
    print(
        json.dumps(
            {
                "status": result["status"],
                "summary": summary,
                "files_checked": len(files),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
