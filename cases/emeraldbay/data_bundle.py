"""Official EmeraldBay split loading plus read-only feature blocks."""

from __future__ import annotations

import importlib.util
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

DATASET = "EmeraldBay"
OFFICIAL_SPLITS = [f"EmeraldBay/split_{i}" for i in range(5)]
RHAISTER_BASELINE_R2 = 0.2568
RHAISTER_FEATURE_R2 = 0.3049
MEAN_EXPR_BLOCK = "mean_expr_2k"
RHAISTER_FEATURE_SOURCE_NAMES = ("cell_eval", "pdex", "pdex_pv", "pdex_fdr")


def _repo_root() -> Path:
    env = os.environ.get("ARF_MACHINE_EXPERIMENTS_ROOT")
    if env:
        path = Path(env)
        if path.is_dir():
            return path
    case_dir = Path(__file__).resolve().parent
    for base in [case_dir, *case_dir.parents]:
        if (base / "cases").is_dir() and (base / "README.md").is_file():
            return base
    raise FileNotFoundError(
        "Cannot find arf-machine-experiments repo root. "
        "Set ARF_MACHINE_EXPERIMENTS_ROOT or run from the repo."
    )


def _default_data_root() -> Path:
    return _repo_root() / "data" / "emeraldbay"


def _default_rhaister_data_root() -> Path:
    return _repo_root() / "data" / "rhaister_data"


def _default_expression_features() -> Path:
    return _default_rhaister_data_root() / "EmeraldBay" / "expression_means_2k.parquet"


def _require_abs_path(env_name: str, default: Path | None = None) -> Path:
    raw = os.environ.get(env_name)
    if not raw and default is not None:
        os.environ.setdefault(env_name, str(default))
        raw = os.environ.get(env_name)
    if not raw:
        raise EnvironmentError(
            f"Set {env_name} to an absolute path before running the benchmark."
        )
    path = Path(raw)
    if not path.is_absolute():
        raise EnvironmentError(f"{env_name} must be an absolute path, got {raw!r}")
    return path


def _rhaister_repo() -> Path:
    env = os.environ.get("RHAISTER_REPO")
    if env:
        path = Path(env)
        if path.is_dir():
            return path
        raise FileNotFoundError(f"RHAISTER_REPO is not a directory: {env!r}")
    candidate = _repo_root() / "data" / "Rhaister"
    if candidate.is_dir():
        return candidate
    raise FileNotFoundError(
        "Rhaister repo not found. Put tahoebio/Rhaister under data/Rhaister "
        "or set RHAISTER_REPO. See cases/emeraldbay/env.example."
    )


def _emeraldbay_data_root() -> Path:
    return _require_abs_path("EMERALDBAY_DATA_ROOT", _default_data_root())


def _rhaister_data_root() -> Path:
    return _require_abs_path("RHAISTER_DATA_ROOT", _default_rhaister_data_root())


def _expression_features_path() -> Path:
    return _require_abs_path(
        "EMERALDBAY_EXPRESSION_FEATURES", _default_expression_features()
    )


def _summary_stats_path() -> Path:
    return _emeraldbay_data_root() / "metadata" / "summary_statistics.parquet"


def _gene_metadata_path() -> Path:
    return _emeraldbay_data_root() / "metadata" / "gene_metadata.parquet"


def verify_required_data() -> None:
    root = _rhaister_data_root() / DATASET
    required = [
        _summary_stats_path(),
        _gene_metadata_path(),
        _expression_features_path(),
        root / "cell_eval" / "all_delta.parquet",
        root / "pdex" / "all_pdex.parquet",
        _rhaister_repo() / "splits" / DATASET / "dataset.toml",
    ]
    required.extend(
        _rhaister_repo() / "splits" / DATASET / f"split_{i}" / "split.toml"
        for i in range(5)
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing EmeraldBay/Rhaister artifacts:\n" + "\n".join(missing)
        )


def _load_expression_feature_map() -> (
    Tuple[Dict[Tuple[str, str], np.ndarray], List[str]]
):
    path = _expression_features_path()
    print(f"Loading {MEAN_EXPR_BLOCK} from {path}...")
    df = pd.read_parquet(path)
    meta = {"cell_line", "condition"}
    gene_cols = [c for c in df.columns if c not in meta]
    cells = df["cell_line"].astype(str).to_numpy()
    conditions = df["condition"].astype(str).to_numpy()
    vectors = df[gene_cols].to_numpy(dtype=np.float32)
    feature_map = {
        (str(cells[i]), str(conditions[i])): vectors[i] for i in range(len(cells))
    }
    print(f"  {len(feature_map)} (cell_line, condition) vectors, dim={len(gene_cols)}")
    return feature_map, gene_cols


def _align_features(
    cells: np.ndarray,
    treatments: np.ndarray,
    feature_map: Dict[Tuple[str, str], np.ndarray],
    n_features: int,
) -> Tuple[np.ndarray, int]:
    X = np.zeros((len(cells), n_features), dtype=np.float32)
    missing = 0
    for i, (cell, treatment) in enumerate(zip(cells, treatments)):
        vec = feature_map.get((str(cell), str(treatment)))
        if vec is None:
            missing += 1
        else:
            X[i] = vec
    return X, missing


def _build_mean_expression_block(
    data: dict,
    feature_map: Dict[Tuple[str, str], np.ndarray],
    gene_cols: List[str],
) -> dict:
    n_features = len(gene_cols)
    X_train, n_miss_tr = _align_features(
        data["train_cells"], data["train_treatments"], feature_map, n_features
    )
    X_test, n_miss_te = _align_features(
        data["test_cells"], data["test_treatments"], feature_map, n_features
    )
    train_mask = np.array(
        [
            (str(c), str(t)) in feature_map
            for c, t in zip(data["train_cells"], data["train_treatments"])
        ]
    )
    if int(train_mask.sum()) == 0:
        raise RuntimeError(
            f"No training rows have {MEAN_EXPR_BLOCK!r} features; cannot z-score"
        )
    mu = X_train[train_mask].mean(axis=0)
    sigma = X_train[train_mask].std(axis=0)
    sigma = np.where(sigma > 1e-6, sigma, 1.0)
    X_train = ((X_train - mu) / sigma).astype(np.float32)
    X_test = ((X_test - mu) / sigma).astype(np.float32)

    cell_to_idx = data["cell_to_idx"]
    treat_to_idx = data["treat_to_idx"]
    feat_mat = np.zeros(
        (data["n_cells"], data["n_treatments"], n_features), dtype=np.float32
    )
    n_filled = 0
    for (cell, treatment), vec in feature_map.items():
        ci = cell_to_idx.get(cell)
        ti = treat_to_idx.get(treatment)
        if ci is None or ti is None:
            continue
        feat_mat[ci, ti] = ((vec - mu) / sigma).astype(np.float32)
        n_filled += 1
    print(
        f"  block {MEAN_EXPR_BLOCK!r}: dim={n_features}, missing "
        f"{n_miss_tr}/{len(X_train)} train, {n_miss_te}/{len(X_test)} test rows; "
        f"feat_mat populated {n_filled}/{data['n_cells'] * data['n_treatments']}"
    )
    return {
        "X_train": X_train,
        "X_test": X_test,
        "feat_mat": feat_mat,
        "gene_cols": [f"{MEAN_EXPR_BLOCK}__{g}" for g in gene_cols],
    }


def _ensure_rhaister_import() -> Any:
    repo = _rhaister_repo()
    rhaister_pkg = repo / "rhaister"
    module_path = rhaister_pkg / "prepare_sensitivity.py"
    if not module_path.is_file():
        raise FileNotFoundError(f"Missing {module_path}")
    for entry in (str(repo), str(rhaister_pkg)):
        if entry not in sys.path:
            sys.path.insert(0, entry)
    spec = importlib.util.spec_from_file_location(
        "rhaister_prepare_sensitivity", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Rhaister module from {module_path}")
    ps = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ps)

    ps._cache_is_valid = lambda _split_name: False
    summary_path = _summary_stats_path()

    _aggregated_cache: Dict[str, Any] = {}

    def load_data(dataset: str):
        if dataset != DATASET:
            raise ValueError(f"Unexpected dataset {dataset!r}")
        if dataset in _aggregated_cache:
            return _aggregated_cache[dataset]
        if not summary_path.is_file():
            raise FileNotFoundError(f"Missing {summary_path}")
        df = pd.read_parquet(summary_path)
        cfg = ps._load_dataset_config(dataset)
        aggregated = ps._apply_filters_and_aggregate(df, ps._filters(cfg))
        _aggregated_cache[dataset] = aggregated
        return aggregated

    def prepare_all(split_name: str, with_features=None):
        dataset, _ = ps.parse_split_name(split_name)
        print(f"Loading {dataset} parquet (filters from dataset.toml)...")
        df = load_data(dataset)
        print(f"  {len(df)} (cell_line, condition) rows after filter+aggregate")
        print(
            f"  {df['cell_line'].nunique()} cell lines, {df['condition'].nunique()} conditions"
        )

        print("Parsing split...")
        split_info = ps.load_split(split_name)
        train_df, test_df = ps.make_splits(df, split_info)
        print(f"  Train: {len(train_df)} rows  |  Test: {len(test_df)} rows")

        train_cells, train_treatments, y_train = ps.to_arrays(train_df)
        test_cells, test_treatments, y_test = ps.to_arrays(test_df)

        all_cells = sorted(df["cell_line"].unique())
        all_treatments = sorted(df["condition"].unique())
        cell_to_idx = {c: i for i, c in enumerate(all_cells)}
        treat_to_idx = {t: i for i, t in enumerate(all_treatments)}

        data = {
            "train_cells": train_cells,
            "train_treatments": train_treatments,
            "y_train": y_train,
            "test_cells": test_cells,
            "test_treatments": test_treatments,
            "y_test": y_test,
            "cell_to_idx": cell_to_idx,
            "treat_to_idx": treat_to_idx,
            "n_cells": len(all_cells),
            "n_treatments": len(all_treatments),
            "evaluate_test": ps.make_evaluator(y_test),
        }
        if with_features:
            data = ps.attach_features(data, dataset, source=with_features)
        return data

    ps.load_data = load_data
    ps.prepare_all = prepare_all
    return ps


@dataclass
class SplitBundle:
    """Data passed to main.predict_sensitivity. Test labels stay inside evaluate_test."""

    split_name: str
    y_train: np.ndarray
    c_train: np.ndarray
    t_train: np.ndarray
    c_test: np.ndarray
    t_test: np.ndarray
    n_cells: int
    n_treatments: int
    cell_to_idx: Dict[str, int]
    treat_to_idx: Dict[str, int]
    idx_to_cell: List[str]
    idx_to_treat: List[str]
    feature_blocks: Dict[str, np.ndarray]
    X_train_feature_blocks: Dict[str, np.ndarray]
    X_test_feature_blocks: Dict[str, np.ndarray]
    feature_gene_cols: Dict[str, List[str]]
    feat_mat: Optional[np.ndarray]
    X_train_features: Optional[np.ndarray]
    X_test_features: Optional[np.ndarray]
    evaluate_test: Callable[[np.ndarray], Dict[str, float]]

    @classmethod
    def from_prepared(cls, data: dict, split_name: str) -> "SplitBundle":
        cell_to_idx = data["cell_to_idx"]
        treat_to_idx = data["treat_to_idx"]
        idx_to_cell = [
            cell for cell, _ in sorted(cell_to_idx.items(), key=lambda kv: kv[1])
        ]
        idx_to_treat = [
            treat for treat, _ in sorted(treat_to_idx.items(), key=lambda kv: kv[1])
        ]
        c_train = np.array(
            [cell_to_idx[c] for c in data["train_cells"]], dtype=np.int64
        )
        t_train = np.array(
            [treat_to_idx[t] for t in data["train_treatments"]], dtype=np.int64
        )
        c_test = np.array([cell_to_idx[c] for c in data["test_cells"]], dtype=np.int64)
        t_test = np.array(
            [treat_to_idx[t] for t in data["test_treatments"]], dtype=np.int64
        )
        blocks = data["feature_blocks"]
        row_train_blocks = data["X_train_feature_blocks"]
        row_test_blocks = data["X_test_feature_blocks"]
        concat_train = np.concatenate(
            [row_train_blocks[name] for name in blocks], axis=1
        ).astype(np.float32)
        concat_test = np.concatenate(
            [row_test_blocks[name] for name in blocks], axis=1
        ).astype(np.float32)
        concat_mat = np.concatenate([blocks[name] for name in blocks], axis=-1).astype(
            np.float32
        )
        return cls(
            split_name=split_name,
            y_train=np.asarray(data["y_train"], dtype=np.float64),
            c_train=c_train,
            t_train=t_train,
            c_test=c_test,
            t_test=t_test,
            n_cells=int(data["n_cells"]),
            n_treatments=int(data["n_treatments"]),
            cell_to_idx=cell_to_idx,
            treat_to_idx=treat_to_idx,
            idx_to_cell=idx_to_cell,
            idx_to_treat=idx_to_treat,
            feature_blocks=blocks,
            X_train_feature_blocks=row_train_blocks,
            X_test_feature_blocks=row_test_blocks,
            feature_gene_cols=data["feature_gene_cols"],
            feat_mat=concat_mat,
            X_train_features=concat_train,
            X_test_features=concat_test,
            evaluate_test=data["evaluate_test"],
        )


class SplitBundleCache:
    """In-process cache only; the evaluator does not write split artifacts."""

    def __init__(self) -> None:
        self._ps: Any = None
        self._mean_expr: Optional[
            Tuple[Dict[Tuple[str, str], np.ndarray], List[str]]
        ] = None
        self._rhaister_feature_maps: Optional[Dict[str, Tuple[dict, List[str]]]] = None
        self._bundles: Dict[str, SplitBundle] = {}

    def _prepare_sensitivity(self) -> Any:
        if self._ps is None:
            os.environ.setdefault("PYTHONUTF8", "1")
            verify_required_data()
            self._ps = _ensure_rhaister_import()
        return self._ps

    def _mean_expression(self) -> Tuple[Dict[Tuple[str, str], np.ndarray], List[str]]:
        if self._mean_expr is None:
            self._mean_expr = _load_expression_feature_map()
        return self._mean_expr

    def _rhaister_features(self, ps: Any) -> Dict[str, Tuple[dict, List[str]]]:
        if self._rhaister_feature_maps is None:
            self._rhaister_feature_maps = {
                source: ps._resolve_source(source, DATASET)
                for source in RHAISTER_FEATURE_SOURCE_NAMES
            }
        return self._rhaister_feature_maps

    def get_bundle(self, split_name: str) -> SplitBundle:
        if split_name not in self._bundles:
            ps = self._prepare_sensitivity()
            data = ps.prepare_all(split_name, with_features=None)

            feature_map, gene_cols = self._mean_expression()
            blocks = {
                MEAN_EXPR_BLOCK: _build_mean_expression_block(
                    data, feature_map, gene_cols
                )
            }

            for source, (source_map, source_gene_cols) in self._rhaister_features(
                ps
            ).items():
                blocks[source] = ps._build_block(
                    data, source_map, source_gene_cols, source
                )

            data = dict(data)
            data["feature_blocks"] = {
                name: block["feat_mat"] for name, block in blocks.items()
            }
            data["X_train_feature_blocks"] = {
                name: block["X_train"] for name, block in blocks.items()
            }
            data["X_test_feature_blocks"] = {
                name: block["X_test"] for name, block in blocks.items()
            }
            data["feature_gene_cols"] = {
                name: block["gene_cols"] for name, block in blocks.items()
            }
            self._bundles[split_name] = SplitBundle.from_prepared(data, split_name)
        return self._bundles[split_name]


_CACHE = SplitBundleCache()


def get_split_bundle(split_name: str) -> SplitBundle:
    return _CACHE.get_bundle(split_name)
