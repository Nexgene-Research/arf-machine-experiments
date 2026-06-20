"""Benchmark harness for the EmeraldBay drug sensitivity case."""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback

import numpy as np

from data_bundle import (
    OFFICIAL_SPLITS,
    RHAISTER_BASELINE_R2,
    RHAISTER_FEATURE_R2,
    get_split_bundle,
)
from main import predict_sensitivity


def evaluate() -> tuple[bool, dict[str, float] | str]:
    try:
        os.environ.setdefault("PYTHONUTF8", "1")

        split_metrics = []
        feature_shapes = {}
        for split_name in OFFICIAL_SPLITS:
            bundle = get_split_bundle(split_name)
            if not feature_shapes:
                feature_shapes = {
                    name: tuple(int(v) for v in block.shape)
                    for name, block in bundle.feature_blocks.items()
                }

            y_pred = predict_sensitivity(bundle)
            y_pred = np.asarray(y_pred, dtype=np.float64)
            if y_pred.ndim != 1:
                raise ValueError(
                    f"{split_name}: predictions must be 1D, got shape {y_pred.shape}"
                )
            if y_pred.shape[0] != bundle.c_test.shape[0]:
                raise ValueError(
                    f"{split_name}: expected {bundle.c_test.shape[0]} predictions, got {y_pred.shape[0]}"
                )
            if not np.all(np.isfinite(y_pred)):
                raise ValueError(f"{split_name}: predictions contain NaN or inf values")
            split_metrics.append(bundle.evaluate_test(y_pred))

        r2_values = [m["sensitivity/r2"] for m in split_metrics]
        pearson_values = [m["sensitivity/pearson"] for m in split_metrics]
        mse_values = [m["sensitivity/mse"] for m in split_metrics]
        mae_values = [m["sensitivity/mae"] for m in split_metrics]

        mean_r2 = float(np.mean(r2_values))
        metrics: dict[str, float] = {
            "combined_score": mean_r2,
            "mean_r2": mean_r2,
            "worst_r2": float(np.min(r2_values)),
            "mean_pearson": float(np.mean(pearson_values)),
            "mean_mse": float(np.mean(mse_values)),
            "mean_mae": float(np.mean(mae_values)),
            "delta_vs_rhaister_base": mean_r2 - RHAISTER_BASELINE_R2,
            "delta_vs_rhaister_features": mean_r2 - RHAISTER_FEATURE_R2,
            "feature_block_count": float(len(feature_shapes)),
        }
        for i, r2 in enumerate(r2_values):
            metrics[f"split_{i}_r2"] = float(r2)
            metrics[f"split_{i}_pearson"] = float(pearson_values[i])
            metrics[f"split_{i}_mse"] = float(mse_values[i])
            metrics[f"split_{i}_mae"] = float(mae_values[i])
        for name, shape in feature_shapes.items():
            metrics[f"feature_{name}_dim"] = float(shape[-1])

        return True, metrics
    except Exception:
        return False, traceback.format_exc()


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate EmeraldBay predictor.")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print metrics as JSON to stdout.",
    )
    args = parser.parse_args()

    ok, result = evaluate()
    if not ok:
        print(result, file=sys.stderr)
        return 1

    metrics = result
    if args.json:
        print(json.dumps(metrics, indent=2, sort_keys=True))
    else:
        print(f"mean_r2:      {metrics['mean_r2']:.4f}")
        print(f"mean_pearson: {metrics['mean_pearson']:.4f}")
        print(f"mean_mse:     {metrics['mean_mse']:.4f}")
        print(f"mean_mae:     {metrics['mean_mae']:.4f}")
        print(f"worst_r2:     {metrics['worst_r2']:.4f}")
        for i in range(5):
            print(
                f"split_{i}: r2={metrics[f'split_{i}_r2']:.4f} "
                f"pearson={metrics[f'split_{i}_pearson']:.4f} "
                f"mse={metrics[f'split_{i}_mse']:.4f} "
                f"mae={metrics[f'split_{i}_mae']:.4f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
