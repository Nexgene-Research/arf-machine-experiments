"""Deterministic export of manually sourced claims; no efficacy calculations or web calls."""

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import CURATION, DATA, REFERENCE


def load(name):
    return json.loads((CURATION / name).read_text(encoding="utf-8"))


def write_csv(name, rows, columns=None):
    if columns is None:
        columns = list(rows[0])
    DATA.mkdir(exist_ok=True)
    with (DATA / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def mapped(columns, rows):
    for row in rows:
        assert len(row) == len(columns), (len(row), len(columns), row[0])
    return [dict(zip(columns, row)) for row in rows]


def main():
    c = load("curation.json")
    update = load("bridge_updates.json")
    evidence = []
    for raw in c["evidence"] + c["additional_evidence"] + update["evidence"]:
        row = dict(c["evidence_defaults"])
        row.update(raw)
        evidence.append(row)
    for eid, drug, setid, indication, revision, risk in c["regulatory"]:
        row = dict(c["evidence_defaults"])
        row.update(
            evidence_id=eid,
            claim=f"{drug}: {indication}. {risk}",
            identifier=f"DailyMed setid:{setid}; source date:{revision}",
            url=f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={setid}",
            source_type="regulatory",
            access="full_text_sections",
            locator="Label indications, warnings, selected pediatric/interactions sections; date is revision/web update, not original approval date",
            model_species="Human clinical label",
            genotype="Not BUB1B-selected",
            intervention_comparator=f"{drug} approved-use/safety record; no MVA comparator",
            exposure="Approved indication context; exact MVA rescue exposure not established",
            endpoint="Approval, formulation, indication and safety",
            effect_magnitude=risk,
            clinical_context_relevance="moderate",
            study_family=f"label_{drug}",
        )
        evidence.append(row)
    overrides = load("evidence_overrides.json")
    for row in evidence:
        row.update(overrides.get(row["evidence_id"], {}))
    evidence.sort(key=lambda r: int(r["evidence_id"][1:]))
    evidence_columns = [
        "evidence_id",
        "claim",
        "identifier",
        "url",
        "source_type",
        "access",
        "read_scope",
        "locator",
        "model_species",
        "tissue_cell_state",
        "genotype",
        "intervention_comparator",
        "exposure",
        "duration",
        "endpoint",
        "effect_magnitude",
        "replication_statistics",
        "target_specificity",
        "orthogonal_validation",
        "negative_conflicting",
        "case_relevance",
        "transfer_assumptions",
        "funding_conflicts",
        "rigor",
        "causal_directness",
        "variant_model_relevance",
        "exposure_feasibility",
        "reproducibility",
        "clinical_context_relevance",
        "study_family",
        "access_date",
    ]
    write_csv("evidence.csv", evidence, evidence_columns)
    g = load("graph_curation.json")
    node_columns = ["node_id", "node_type", "label", "context", "description"]
    nodes = mapped(node_columns, g["nodes"] + update["nodes"])
    for row in evidence:
        nodes.append(
            dict(
                node_id="SRC_" + row["evidence_id"],
                node_type="evidence_source",
                label=row["evidence_id"] + " " + row["identifier"],
                context=row["model_species"],
                description=row["url"],
            )
        )
    edge_columns = [
        "edge_id",
        "source",
        "target",
        "relationship",
        "sign",
        "evidence_ids",
        "evidence_state",
        "species",
        "tissue",
        "genotype",
        "cell_state",
        "exposure",
        "limitations",
        "transfer_assumptions",
    ]
    edges = mapped(edge_columns, g["edges"] + update["edges"])
    write_csv("nodes.csv", nodes, node_columns)
    write_csv("edges.csv", edges, edge_columns)
    edge_map = {r["edge_id"]: r for r in edges}
    node_map = {r["node_id"]: r for r in nodes}
    path_columns = [
        "path_id",
        "category",
        "objective",
        "edge_ids",
        "hypothesis",
        "weakest_link",
        "decisive_test",
        "feedback_timing_competing_effects",
    ]
    paths = mapped(path_columns, g["paths"] + update["paths"])
    for p in paths:
        selected = [edge_map[x] for x in p["edge_ids"].split(";")]
        assert all(
            a["target"] == b["source"] for a, b in zip(selected, selected[1:])
        ), p["path_id"]
        transfers = sorted(
            {
                t
                for e in selected
                for t in e["transfer_assumptions"].split(";")
                if ":" in t
            }
        )
        p["hop_count"] = len(selected)
        p["unproven_links"] = sum(e["evidence_state"] != "measured" for e in selected)
        p["hypothesized_links"] = sum(
            e["evidence_state"] == "hypothesized" for e in selected
        )
        p["inferred_links"] = sum(e["evidence_state"] == "inferred" for e in selected)
        p["contradicted_links"] = sum(
            e["evidence_state"] == "contradicted" for e in selected
        )
        p["context_transfer_count"] = len(transfers)
        p["species_transfers"] = sum(t.startswith("species:") for t in transfers)
        p["tissue_transfers"] = sum(t.startswith("tissue:") for t in transfers)
        p["genotype_transfers"] = sum(t.startswith("genotype:") for t in transfers)
        p["transfer_events"] = (
            ";".join(transfers) or "No transfer specified; inspect edge context"
        )
        p["direction_check"] = (
            "Directed topology checked; biological net sign is conditional, not multiplied"
        )
    write_csv("paths.csv", paths)
    cc = load("candidate_curation.json")
    candidates = mapped(cc["columns"], cc["rows"])
    write_csv("candidates.csv", candidates, cc["columns"])
    rejections = []
    for r in candidates:
        if r["Tier"] == "X":
            rejections.append(
                dict(
                    candidate=r["Candidate"],
                    objective=r["Objective"],
                    reason=r["Counterevidence"],
                    evidence_ids=r["Evidence anchors"],
                    kill_criterion=r["Kill criterion"],
                    reopen_only_if=r["Tier change"],
                )
            )
    for drug, reason, eids, reopen in [
        (
            "Ataluren",
            "EU authorization not renewed; no exact BUB1B functional evidence",
            "E28;E18",
            "New status plus exact-allele functional/exposure evidence",
        ),
        (
            "Generic PARP-inhibitor selection",
            "CIN alone does not establish HR deficiency or tumor dependency",
            "E24",
            "Validated tumor-specific repair dependency and normal-cell window",
        ),
        (
            "Systemic NMD inhibition",
            "Broad surveillance effects; harmful products and no case functional rescue",
            "E01;E03;E18",
            "Allele substrate limitation plus demonstrated useful output and selective safety",
        ),
        (
            "Antiapoptotic systemic rescue",
            "Fly brain-size rescue does not justify preserving potentially malignant cells",
            "E44",
            "Tumor-safe context-specific mechanism; not a systemic clinical proposal",
        ),
        (
            "BUB1 kinase inhibitors as BUB1B replacement",
            "BUB1 and BUB1B are distinct; checkpoint inhibition is not scaffold restoration",
            "E04;E07",
            "Direct genotype-compatible suppressor evidence",
        ),
    ]:
        rejections.append(
            dict(
                candidate=drug,
                objective="Context-dependent; not admitted",
                reason=reason,
                evidence_ids=eids,
                kill_criterion=reason,
                reopen_only_if=reopen,
            )
        )
    write_csv("rejections.csv", rejections)
    sc = load("search_curation.json")
    searches = []
    for i, (direction, tool, query, outcome) in enumerate(sc["searches"], 1):
        searches.append(
            dict(
                search_id=f"S{i:03d}",
                date=c["cutoff"],
                direction=direction,
                source_or_tool=tool,
                query_or_url=query,
                outcome=outcome,
                coverage_note="Selected material discovery query; not a complete export of discovery-service internal subqueries",
            )
        )
    for e in evidence:
        searches.append(
            dict(
                search_id="R_" + e["evidence_id"],
                date=c["cutoff"],
                direction="source_retrieval",
                source_or_tool=e["source_type"],
                query_or_url=e["url"],
                outcome=e["access"] + "; " + e["locator"],
                coverage_note=e["read_scope"],
            )
        )
    write_csv("search_log.csv", searches)
    write_csv(
        "source_coverage.csv",
        mapped(["source", "status", "anchor", "limitations"], sc["coverage"]),
    )
    write_csv(
        "datasets.csv",
        mapped(
            [
                "accession",
                "url",
                "access_status",
                "model",
                "genotype",
                "treatment_time_concentration",
                "replicates",
                "analysis_status",
                "proposed_analysis",
                "confounders",
            ],
            sc["datasets"],
        ),
    )
    hpo = json.loads((REFERENCE / "hpo_selected.json").read_text(encoding="utf-8"))
    write_csv("phenotype_terms.csv", hpo["terms"], ["id", "name", "def"])
    cv = json.loads((REFERENCE / "clinvar_summary.json").read_text(encoding="utf-8"))
    populations = []
    for ident in ["533901", "4600147"]:
        rec = cv["result"][ident]
        for variant in rec["variation_set"]:
            populations.append(
                dict(
                    clinvar_id=ident,
                    title=rec["title"],
                    annotations=variant.get("allele_freq_set", []),
                    interpretation="ClinVar-carried population annotations only; release/ancestry/count/coverage not independently verified. Empty is unknown, not absent.",
                )
            )
    (DATA / "population_annotations.json").write_text(
        json.dumps(populations, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Readable directed causal graph",
        "",
        "Generated from the exported graph. Bracket labels are edge IDs, evidence state and claim IDs. ",
        "Measured means reported in the stated model, not demonstrated in this child. Inferred, hypothesized and contradicted links remain visible.",
        "Conservative unproven-link count includes every non-measured edge; model-to-case transfers are counted separately.",
        "",
    ]
    for p in paths:
        lines.extend(
            [
                "## " + p["path_id"] + ": " + p["hypothesis"],
                "",
                f"Objective {p['objective']}; {p['category']}; {p['hop_count']} links; {p['unproven_links']} unproven links; {p['context_transfer_count']} specified transfer events.",
                "",
                "```text",
            ]
        )
        for eid in p["edge_ids"].split(";"):
            e = edge_map[eid]
            lines.append(
                f"{node_map[e['source']]['label']} --{e['relationship']} [{eid}; {e['evidence_state']}; {e['evidence_ids']}]--> {node_map[e['target']]['label']}"
            )
        lines.extend(
            [
                "```",
                "",
                "**Weakest link:** " + p["weakest_link"],
                "**Decisive test:** " + p["decisive_test"],
                "**Competing effects:** " + p["feedback_timing_competing_effects"],
                "**Transfers:** " + p["transfer_events"],
                "",
            ]
        )
    lines.extend(
        [
            "## Complete edge adjacency",
            "",
            "| Edge | From | Relationship / sign | To | State | Evidence |",
            "|---|---|---|---|---|---|",
        ]
    )
    for e in edges:
        lines.append(
            "| "
            + " | ".join(
                [
                    e["edge_id"],
                    e["source"],
                    e["relationship"] + " / " + e["sign"],
                    e["target"],
                    e["evidence_state"],
                    e["evidence_ids"],
                ]
            )
            + " |"
        )
    (DATA / "graph.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = dict(
        evidence_claims=len(evidence),
        nodes=len(nodes),
        edges=len(edges),
        paths=len(paths),
        candidates=len(candidates),
        tiers=dict(Counter(r["Tier"] for r in candidates)),
        evidence_states=dict(Counter(e["evidence_state"] for e in edges)),
    )
    (DATA / "package_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
