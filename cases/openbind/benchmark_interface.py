"""
OpenBind EV-A71 2A affinity benchmark harness.

This script loads the public OpenBind affinity reference table, calls a released
predictor module, aggregates structure-level predictions to compound level, and
reports benchmark metrics.
"""

from __future__ import annotations

import argparse
import importlib
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.stats


CASE_DIR = Path(__file__).resolve().parent
COVERAGE_FLOOR = 0.98
RMSE_NORMALIZER = 3.0
RMSE_WEIGHT = 0.15
MW_DELTA_WEIGHT = 0.50
COVERAGE_WEIGHT = 2.0
OFFICIAL_MW_SPEARMAN = 0.483
OFFICIAL_CLOGP_SPEARMAN = 0.174
OFFICIAL_GNINA_SPEARMAN = 0.453
OFFICIAL_BOLTZ2_RMSE = 1.091
DEFAULT_PREDICTOR_MODULE = "main_optimized"


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OpenBind EV-A71 2A affinity benchmark harness.",
    )
    parser.add_argument(
        "--benchmark-path",
        type=Path,
        required=True,
        help=(
            "Path to the cloned EV-A71_2A_benchmark repo "
            "(https://github.com/OpenBind-Consortium/EV-A71_2A_benchmark)."
        ),
    )
    parser.add_argument(
        "--predictor",
        default=DEFAULT_PREDICTOR_MODULE,
        help="Predictor module name (default: %(default)s).",
    )
    return parser.parse_args(argv)


def _load_predictor(module_name: str):
    if str(CASE_DIR) not in sys.path:
        sys.path.insert(0, str(CASE_DIR))
    module = importlib.import_module(module_name)
    return module.predict_affinity


def _read_reference(benchmark_dir: Path) -> pd.DataFrame:
    reference_path = (
        benchmark_dir / "affinity" / "reference" / "fragalysis_compound_reference.csv"
    )
    if not reference_path.exists():
        raise FileNotFoundError(f"OpenBind reference table not found: {reference_path}")

    reference = pd.read_csv(reference_path)
    required_cols = {"fragalysis_code", "smiles", "experimental_pKD"}
    missing = required_cols - set(reference.columns)
    if missing:
        raise ValueError(f"{reference_path} is missing columns: {sorted(missing)}")

    reference = reference.dropna(
        subset=["fragalysis_code", "smiles", "experimental_pKD"]
    ).copy()
    reference["fragalysis_code"] = reference["fragalysis_code"].astype(str)
    reference["smiles"] = reference["smiles"].astype(str)
    reference["experimental_pKD"] = pd.to_numeric(
        reference["experimental_pKD"], errors="coerce"
    )
    reference = reference.dropna(subset=["experimental_pKD"]).copy()
    return reference


def _public_compounds(reference: pd.DataFrame) -> pd.DataFrame:
    return reference[["fragalysis_code", "smiles"]].copy()


def _prediction_to_frame(predictions: Any, compounds: pd.DataFrame) -> pd.DataFrame:
    if isinstance(predictions, pd.DataFrame):
        frame = predictions.copy()
    elif isinstance(predictions, pd.Series):
        frame = pd.DataFrame(
            {
                "fragalysis_code": compounds["fragalysis_code"].to_numpy(),
                "predicted_affinity": predictions.to_numpy(),
            }
        )
    elif isinstance(predictions, dict):
        frame = pd.DataFrame(
            {
                "fragalysis_code": list(predictions.keys()),
                "predicted_affinity": list(predictions.values()),
            }
        )
    else:
        values = np.asarray(predictions, dtype=float)
        if len(values) != len(compounds):
            raise ValueError(
                "Array-like predictions must have the same length as the input compounds"
            )
        frame = pd.DataFrame(
            {
                "fragalysis_code": compounds["fragalysis_code"].to_numpy(),
                "predicted_affinity": values,
            }
        )

    required_cols = {"fragalysis_code", "predicted_affinity"}
    missing = required_cols - set(frame.columns)
    if missing:
        raise ValueError(f"Predictions are missing columns: {sorted(missing)}")

    frame = frame[["fragalysis_code", "predicted_affinity"]].copy()
    frame["fragalysis_code"] = frame["fragalysis_code"].astype(str)
    frame["predicted_affinity"] = pd.to_numeric(
        frame["predicted_affinity"], errors="coerce"
    )
    return frame


def _join_codes(values: pd.Series) -> str:
    return "; ".join(sorted(set(values.dropna().astype(str))))


def _compound_level(reference: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    merged = reference.merge(
        predictions[["fragalysis_code", "predicted_affinity"]],
        on="fragalysis_code",
        how="left",
    )
    compound = (
        merged.groupby("smiles", as_index=False)
        .agg(
            experimental_pKD=("experimental_pKD", "first"),
            predicted_affinity=("predicted_affinity", "mean"),
            n_structures_with_reference=("fragalysis_code", "nunique"),
            n_structures_with_prediction=("predicted_affinity", "count"),
            fragalysis_codes=("fragalysis_code", _join_codes),
        )
        .sort_values("smiles")
    )
    return compound


def _spearman(y_true: pd.Series, y_pred: pd.Series) -> float:
    if len(y_true) < 3 or y_pred.nunique(dropna=True) <= 1:
        return 0.0
    rho, _ = scipy.stats.spearmanr(y_true, y_pred, nan_policy="omit")
    if rho is None or not np.isfinite(rho):
        return 0.0
    return float(rho)


def _pearson(y_true: pd.Series, y_pred: pd.Series) -> float:
    if len(y_true) < 3 or y_pred.nunique(dropna=True) <= 1:
        return 0.0
    rho, _ = scipy.stats.pearsonr(y_true, y_pred)
    if rho is None or not np.isfinite(rho):
        return 0.0
    return float(rho)


def _score_compound_table(compound: pd.DataFrame) -> dict[str, float]:
    total = int(len(compound))
    scored = compound.dropna(subset=["experimental_pKD", "predicted_affinity"]).copy()
    n_scored = int(len(scored))
    coverage = float(n_scored / total) if total else 0.0

    if n_scored == 0:
        return {
            "spearman_rho": 0.0,
            "pearson_r": 0.0,
            "rmse": RMSE_NORMALIZER,
            "mae": RMSE_NORMALIZER,
            "coverage": 0.0,
            "n_compounds_total": float(total),
            "n_compounds_scored": 0.0,
            "n_missing_predictions": float(total),
        }

    y_true = scored["experimental_pKD"].astype(float)
    y_pred = scored["predicted_affinity"].astype(float)
    residual = y_pred - y_true
    rmse = float(np.sqrt(np.mean(np.square(residual))))
    mae = float(np.mean(np.abs(residual)))
    return {
        "spearman_rho": _spearman(y_true, y_pred),
        "pearson_r": _pearson(y_true, y_pred),
        "rmse": rmse if np.isfinite(rmse) else RMSE_NORMALIZER,
        "mae": mae if np.isfinite(mae) else RMSE_NORMALIZER,
        "coverage": coverage,
        "n_compounds_total": float(total),
        "n_compounds_scored": float(n_scored),
        "n_missing_predictions": float(total - n_scored),
    }


def _load_official_metrics(benchmark_dir: Path) -> dict[str, dict[str, float]]:
    metrics_path = benchmark_dir / "plotting" / "tables" / "affinity_metrics.csv"
    if not metrics_path.exists():
        return {}
    table = pd.read_csv(metrics_path)
    output: dict[str, dict[str, float]] = {}
    for _, row in table.iterrows():
        method = str(row["method"])
        output[method] = {
            "spearman": float(row["Spearman rho"]),
            "rmse": float(row["RMSE"]) if pd.notna(row.get("RMSE")) else 0.0,
        }
    return output


def _combined_score(metrics: dict[str, float]) -> float:
    rmse = min(max(metrics["rmse"], 0.0), RMSE_NORMALIZER)
    normalized_rmse = rmse / RMSE_NORMALIZER
    coverage_shortfall = max(0.0, COVERAGE_FLOOR - metrics["coverage"])
    return float(
        metrics["spearman_rho"]
        - RMSE_WEIGHT * normalized_rmse
        + MW_DELTA_WEIGHT * metrics["delta_vs_mw_spearman"]
        - COVERAGE_WEIGHT * coverage_shortfall
    )


def run_benchmark(
    benchmark_path: Path,
    predictor_module: str = DEFAULT_PREDICTOR_MODULE,
):
    try:
        predict_affinity = _load_predictor(predictor_module)
        reference = _read_reference(benchmark_path)
        compounds = _public_compounds(reference)

        raw_predictions = predict_affinity(compounds)
        predictions = _prediction_to_frame(raw_predictions, compounds)
        compound = _compound_level(reference, predictions)
        metrics = _score_compound_table(compound)

        official = _load_official_metrics(benchmark_path)
        mw_spearman = official.get("molecular_weight", {}).get(
            "spearman", OFFICIAL_MW_SPEARMAN
        )
        clogp_spearman = official.get("clogp", {}).get(
            "spearman", OFFICIAL_CLOGP_SPEARMAN
        )
        gnina_spearman = official.get("gnina", {}).get(
            "spearman", OFFICIAL_GNINA_SPEARMAN
        )
        boltz2_rmse = official.get("boltz-2", {}).get("rmse", OFFICIAL_BOLTZ2_RMSE)

        metrics["predictor_module"] = predictor_module
        metrics["baseline_mw_spearman"] = float(mw_spearman)
        metrics["baseline_clogp_spearman"] = float(clogp_spearman)
        metrics["baseline_gnina_spearman"] = float(gnina_spearman)
        metrics["baseline_boltz2_rmse"] = float(boltz2_rmse)
        metrics["delta_vs_mw_spearman"] = float(metrics["spearman_rho"] - mw_spearman)
        metrics["delta_vs_clogp_spearman"] = float(
            metrics["spearman_rho"] - clogp_spearman
        )
        metrics["delta_vs_gnina_spearman"] = float(
            metrics["spearman_rho"] - gnina_spearman
        )
        metrics["delta_vs_mw_rmse"] = 0.0
        metrics["delta_vs_boltz2_rmse"] = float(boltz2_rmse - metrics["rmse"])
        metrics["combined_score"] = _combined_score(metrics)

        return True, metrics
    except Exception:
        return False, traceback.format_exc()


if __name__ == "__main__":
    args = _parse_args()
    success, result = run_benchmark(
        benchmark_path=args.benchmark_path,
        predictor_module=args.predictor,
    )
    if success:
        for key, value in result.items():
            print(f"{key}: {value}")
    else:
        raise SystemExit(result)
