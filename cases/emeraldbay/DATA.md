# Emerald Bay - Data Briefing

This document summarizes the [Tahoe Emerald Bay](https://huggingface.co/datasets/tahoebio/EmeraldBay) pharmacogenomic screen and how this case uses the official splits, feature artifacts, and proposed solution published with [Rhaister](https://huggingface.co/tahoebio/Rhaister).

## What Emerald Bay Is

Emerald Bay is a single-cell perturbation dataset from Tahoe Bio's MOSAIC high-throughput platform. It pairs two readouts for the same experiments:

1. **Transcriptional profiles** — 1,831,648 single cells with sparse raw counts
2. **Drug sensitivity** — scalar **growth rate** per (cell line, treatment) from five-day endpoint cell-count proportions

| Dimension | Count | Notes |
| --- | ---: | --- |
| Single cells | 1,831,648 | 116 parquet shards (~58 GB) |
| Cell lines | 52 | Cellosaurus IDs (e.g. `CVCL_1055`) |
| Treatments | 91 | 27 single drugs × doses, plus combinations |
| Summary-stat rows | 4,992 | Raw `(cell_line, condition)` table |
| Evaluation rows | 2,340 | After the Rhaister task filters (see below) |
| Genes (vocabulary) | 63,284 | Tahoe-100M tokens + 574 Emerald Bay–specific |

**License:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
**Source:** [tahoebio/EmeraldBay on Hugging Face](https://huggingface.co/datasets/tahoebio/EmeraldBay)

## Dataset Tables (Hugging Face)

The HF dataset exposes five configs. Only metadata is needed for the sensitivity target; expression shards are needed for perturbation features.

| Config | Rows | Role in this case |
| --- | ---: | --- |
| `expression_data` | 1.83M | Build mean-expression and DE features |
| `summary_statistics` | 4,992 | **Growth-rate target** (sensitivity readout) |
| `gene_metadata` | 63,284 | Map gene token IDs → symbols |
| `cell_line_metadata` | 487 | Driver mutations (1–51 rows per line) |
| `drug_metadata` | 27 | MoA, targets, SMILES (single-drug only) |

### `summary_statistics` (the prediction target)

| Column | Description |
| --- | --- |
| `cell_line` | Cellosaurus ID |
| `condition` | Compound × concentration, e.g. `[('Encorafenib', 0.1, 'uM')]` |
| `growth_rate` | Scalar drug response (higher = more growth under treatment) |

### `expression_data` (per-cell transcriptomics)

| Column | Description |
| --- | --- |
| `genes` | Integer token IDs (non-zero genes only) |
| `expressions` | Raw counts, aligned with `genes` |
| `drug` | Treatment name (`DMSO_TF` = vehicle) |
| `drugname_drugconc` | Condition key matching `summary_statistics.condition` |
| `cell_line` | Cellosaurus ID |
| `sample` | Replicate identifier |
| `BARCODE_SUB_LIB_ID` | Unique per cell |

## Prediction Task (This Case)

**Goal:** Predict held-out growth rates for unseen (cell line, drug) combinations.

**Evaluation:** Mean test **R²** across five official Emerald Bay splits (`EmeraldBay/split_0` ... `split_4`) used by [Rhaister](https://huggingface.co/tahoebio/Rhaister).

**Filtering** (from Rhaister's `dataset.toml`, same as the paper):

- Drop rows with missing `growth_rate`
- Drop multi-drug combination conditions (regex: `),(` in condition string)
- Drop `DMSO_T0` time-zero controls
- Mean-aggregate replicate `(cell_line, condition)` rows

Result: **2,340** `(cell_line, condition)` pairs × **52** cell lines.

**Split design:** RIFIVDU holdouts — each split holds out a disjoint set of cell lines and drugs so test pairs are genuinely unseen combinations. Splits live in `Rhaister/splits/EmeraldBay/split_{0..4}/split.toml`.

## Feature Blocks Used by the Discovered Predictor

| Block | Source | Dim | Description |
| --- | --- | ---: | --- |
| `mean_expr_2k` | Expression shards | 2,000 | Mean raw counts per (cell, condition) on Rhaister's static 2K gene list |
| `cell_eval` | Rhaister data prep | 2,000 | Log2 fold-change vs DMSO |
| `pdex` | Rhaister data prep | 2,000 | Perturbation DE expression |
| `pdex_pv` | Rhaister data prep | 2,000 | −log₁₀ p-values |
| `pdex_fdr` | Rhaister data prep | 2,000 | FDR-adjusted p-values |

Rhaister is the authors' proposed solution for this task and reports the following baselines on the same splits:

| Method | Mean R² (5 splits) | Source |
| --- | ---: | --- |
| Rhaister proposed solution | 0.26 | [Rhaister page](https://huggingface.co/tahoebio/Rhaister) |
| Rhaister proposed solution + features | 0.31 | Same |
| Arf Machine (this case) | **0.44** | `benchmark_interface.py` |

## Local Data Layout

This repository does not redistribute Emerald Bay data or Rhaister feature artifacts. The benchmark loader expects the following local files, either under the default ignored `data/` layout or through explicit environment-variable overrides in `env.example`:

```text
{repo}/
  data/
    emeraldbay/                          # EMERALDBAY_DATA_ROOT
      metadata/
        summary_statistics.parquet
        gene_metadata.parquet
      expression_data/train-*.parquet    # optional raw shards, ~58 GB
    rhaister_data/                       # RHAISTER_DATA_ROOT
      EmeraldBay/
        growth_rate_long.parquet         # filtered target (2340 rows)
        expression_means_2k.parquet      # mean expression features
        cell_eval/all_delta.parquet
        pdex/all_pdex.parquet
    Rhaister/                            # RHAISTER_REPO
      splits/EmeraldBay/
      rhaister/
```

The evaluation path requires these prepared artifacts:

- `metadata/summary_statistics.parquet`
- `metadata/gene_metadata.parquet`
- `rhaister_data/EmeraldBay/growth_rate_long.parquet`
- `rhaister_data/EmeraldBay/expression_means_2k.parquet`
- `rhaister_data/EmeraldBay/cell_eval/all_delta.parquet`
- `rhaister_data/EmeraldBay/pdex/all_pdex.parquet`
- `Rhaister/splits/EmeraldBay/dataset.toml`
- `Rhaister/splits/EmeraldBay/split_{0..4}/split.toml`

On Windows PowerShell, set `$env:PYTHONUTF8 = "1"` before running Rhaister-backed evaluation code.

## Citation

Emerald Bay dataset and Rhaister proposed solution:

```bibtex
@article{svensson2026back,
  title={Back to basics: Observed statistics are sufficient to predict drug responses},
  author={Svensson, Valentine and Khan, Umair and Heydari, Hamed and others},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.06.09.731197}
}
```

## Further Reading

- [Emerald Bay dataset card](https://huggingface.co/datasets/tahoebio/EmeraldBay) — schema, tutorials, Colab notebook
- [Rhaister page](https://huggingface.co/tahoebio/Rhaister) — proposed solution, baselines, splits, and reproduction
- [Loading tutorial notebook](https://huggingface.co/datasets/tahoebio/EmeraldBay/blob/main/tutorials/loading_data.ipynb)
