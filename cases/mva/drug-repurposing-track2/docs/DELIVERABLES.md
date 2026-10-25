# Deliverable contract: public-evidence research assessment

Search cutoff: 2026-10-22. Historical and current literature included; priority is causal and genotype relevance, not recency. Difficulty: high. This is a bounded evidence synthesis, not an exhaustive systematic review. No patient data analysis, clinical dosing, or treatment recommendation.

## Checklist fixed before data extraction and graph construction

- `docs/report.md`: concise executive assessment, objectives A to E, variant mechanism before drug shortlist, conditional conclusion.
- `docs/appendices.md`: phenotype audit; axes A to I; all 13 what-if branches; paths and missing-link status; experiments and failure rules; up to three distinct submission concepts; access and saturation limitations.
- `data/nodes.csv`: UTF-8 CSV, unique `node_id`, `node_type`, `label`, `context`, `description`. Separate abundance, checkpoint function, normal cells and tumor cells. Include variant, allele, transcript, protein state/function, complex/interaction, process, cell/tissue phenotype, drug/metabolite/target, safety outcome, model and evidence source.
- `data/edges.csv`: UTF-8 CSV, `edge_id,source,target,relationship,sign,evidence_ids,evidence_state,species,tissue,genotype,cell_state,exposure,limitations,transfer_assumptions`. Referential integrity to nodes and evidence. States: measured/inferred/hypothesized/contradicted. A source establishes only the edge it actually supports. Hypothesized edges explicitly cite the bounding evidence and name the missing experiment.
- `data/evidence.csv`: claim-level ledger with identifiers/URLs, access level and read scope, locator, model and genotype, intervention/comparator, exposure/duration, endpoint, direction/magnitude, replication/statistical limitations, specificity/orthogonal controls, conflicts, case relevance, transfer assumptions, funding/conflicts, six separate evidence grades and study-family identifier. Unknown values have explicit reasons, never fabricated values.
- `data/candidates.csv`: required user columns: Candidate, Approval status/date, Objective, Axis, Mechanistic path, Variant condition, Evidence anchors, Counterevidence, Weakest link, Exposure, Key risks, Tier, Decisive test, Kill criterion. Additional tier-change and pediatric/context fields permitted. No Tier 1 without full gate satisfaction; tools/supplements listed separately.
- `data/paths.csv`: directed edge sequences, short/long/conditional/rejected, unsupported-edge count, species/tissue/genotype transfer count, weakest link, timing/feedback/conflicts, decisive test. Counts computed from final curated path definitions, not therapeutic scores.
- `data/rejections.csv`, `data/search_log.csv`, `data/source_coverage.csv`, `data/datasets.csv`: auditable reasons, queries, coverage/access status, dataset genotype/exposure/replicates/analysis status. Do not imply all listed databases were searched.
- `data/graph.md`: readable adjacency/path representation using evidence-linked directed text. A rendered illustration is unavailable; practical text graph is the deliverable and is generated from the CSVs.
- `scripts/build_package.py`, `scripts/validate_package.py`: standard-library reproducible generation and validation; input extraction in `curation/`; no network or interpreter-state dependency at rerun.
- `README.md`, `environment.txt`: reproducibility/GitHub layout, ranking rules, privacy and third-party redistribution boundaries. No public push requested.
- `data/validation.json`: final mechanical checks for exact files, nonempty schemas, IDs/references, no placeholder tokens, allowed states/types/tiers, path continuity and counts, coverage of all axes/conditions, and file hashes. Run from a fresh process on final files.

## Interpretation constraints

All clinical observations remain user-reported, not medically verified. Public database alleles with matching protein consequences are candidates for report reconciliation, not proven patient genotypes. Source-reported quantities retain their original units and precision. No weighted probability or efficacy score. No patient dose. No synthetic results. Negative search findings are coverage-limited.

Validation covers this static public-evidence package and graph consistency. It cannot validate patient diagnosis, actual phase, tissue expression, treatment interactions, or drug efficacy; no patient experiments or raw omics reanalysis were performed.
