# BUB1B-associated MVA: variant-aware drug-repurposing package

**Authors:** İlsu Duman, Özge Beğli, Hatice Mutlu and Kutsal Baran Ozkurt\*. Nexgene AI, https://nexgene.ai. Correspondence: kutsal.ozkurt@nexgene.ai

**MVA Hackathon 2026, Track 2.** Start with the manuscript in [`docs/report.md`](docs/report.md) (typeset PDF: [`pdf/NexgeneAI-MVA-Track2-Report.pdf`](pdf/NexgeneAI-MVA-Track2-Report.pdf)).

The conclusion is **no approved medicine with a complete demonstrated causal chain for constitutional correction of the reported genotype; Tier 1 is empty.** Every candidate is a repurposing *hypothesis for follow-up*, conditional on allele, tissue and tumour context, not evidence that a medicine works, and not medical advice.

## Mechanism-based strategies

The portfolio is organised into four independent strategies arising from the two reported *BUB1B* variants (`p.Leu737*` nonsense and `p.Asn1002Lys` missense):

1. **Nonsense suppression and BUBR1 restoration** (read-through; targets `p.Leu737*`): azithromycin, amikacin, gentamicin.
2. **Proteostasis and metabolic rescue** (raise or stabilise functional BUBR1): arimoclomol, metformin, trametinib, prescription niacin ER (complementary), hydroxychloroquine.
3. **Aneuploidy-associated cellular stress modulation (PP2A reactivation)**: fingolimod.
4. **RMS (tumour-directed)**: larotrectinib (conditional), vincristine, olaparib.

Considered and deprioritised (rejection ledger): potassium citrate ER, sirolimus, dasatinib, temsirolimus, bortezomib, somatropin.

## Repository layout

```text
docs/          manuscript, appendices, claim audit, deliverable contract
data/          generated CSV ledgers, causal graph, validation hashes
curation/      manual extraction JSON (edit these, not the CSVs)
reference/     public RefSeq, Ensembl, UniProt, ClinVar and HPO snapshots
scripts/       standard-library build, validate and optional PDF render
pdf/           typeset PDF copies of the Markdown deliverables
```

| Path | Role |
|---|---|
| [`docs/report.md`](docs/report.md) | Academic manuscript |
| [`docs/appendices.md`](docs/appendices.md) | Phenotype audit, axes A to I, 13-state decision matrix |
| [`docs/claim_audit.md`](docs/claim_audit.md) | Load-bearing claim-scope audit |
| [`data/evidence.csv`](data/evidence.csv) | Claim-level ledger with six-axis grades |
| [`data/candidates.csv`](data/candidates.csv), [`data/rejections.csv`](data/rejections.csv) | Shortlist and deprioritised ledger |
| [`data/nodes.csv`](data/nodes.csv), [`data/edges.csv`](data/edges.csv), [`data/paths.csv`](data/paths.csv), [`data/graph.md`](data/graph.md) | Directed causal graph |
| `curation/` | Manual extraction inputs that generate the ledgers |
| `reference/` | Public reference snapshots (not patient files) |
| `scripts/` | `normalize_references.py`, `build_package.py`, `validate_package.py`, `render_pdf.py` |

## Reproduce

Python 3.10+ and the standard library:

```text
python scripts/normalize_references.py
python scripts/build_package.py
python scripts/validate_package.py
```

No network, credentials or patient data are required. CSV ledgers and `data/graph.md` are generated from `curation/`, so corrections belong in the JSON inputs, not in the CSVs. Optional PDFs: `python scripts/render_pdf.py` (needs `markdown` and `xhtml2pdf`). Typesetting with Turkish characters is best if Times New Roman is placed in `assets/fonts/` locally; that folder is gitignored because the fonts are proprietary.

## Evidence conventions

Study families prevent duplicate papers, corrections and registry records from counting as independent replication. Source access is a read-scope label, not a quality score. Graph states are measured, inferred, hypothesised or contradicted; a conservative complete-chain check counts every non-measured edge as unproven and counts context transfers separately. No weighted efficacy probability is computed.

## Privacy and redistribution

A rare-feature phenotype combination can itself be identifying, so this report merits a privacy review before any public release. Keep controlled genomic data, clinical reports, pedigrees and identifiers outside any public repository. Full-text PDFs, copyrighted article text and restricted datasets are deliberately not redistributed. The package is licensed under [CC BY 4.0](LICENSE); third-party database records and labels remain under their original terms.

## AI-use disclosure

This analysis was produced with Nexgene AI's AI Co-Scientist, **Arf Machine**, and its underlying medical-reasoning engine, **ASA-W1**, operating as an autonomous, researcher-directed discovery pipeline; all AI-generated outputs were reviewed by the authors.
