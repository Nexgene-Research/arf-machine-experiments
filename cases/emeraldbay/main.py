"""
EmeraldBay drug growth-rate predictor discovered by Arf Machine.

Single-fit LightGBM on concatenated feature blocks plus categorical cell/drug IDs,
with a compact ridge residual correction per official Rhaister split.
"""

from __future__ import annotations

import numpy as np

from data_bundle import SplitBundle

try:
    import lightgbm as lgb
except Exception as e:
    lgb = None
    print(
        f"[predict] warning: lightgbm unavailable, using legacy ensemble only: {e}",
        flush=True,
    )

HP_ALS = 30
HP_DREG_LAM = 0.1
HP_RIDGE_ALPHA = 30.0
HP_BIAS_LAM = 3.0
HP_D2D_ALPHA = 10.0
HP_D2D_GAMMA = 0.2
HP_LGB_ROUNDS = 900
HP_LGB_LEAVES = 192
HP_LGB_LR = 0.03
HP_LGB_MIN_DATA = 64
HP_LGB_FF = 0.75
HP_LGB_BF = 0.85
HP_LGB_L2 = 4.0
HP_RESID_ALPHA = 2.5


def _als_additive(y, c_idx, t_idx, n_cell, n_treat, n_iter):
    """Fit y ~= mu + cell_effect + treatment_effect on observed training rows."""
    mu = float(np.mean(y))
    cell_eff = np.zeros(n_cell, dtype=np.float64)
    treat_eff = np.zeros(n_treat, dtype=np.float64)

    cell_count = np.bincount(c_idx, minlength=n_cell).astype(np.float64)
    treat_count = np.bincount(t_idx, minlength=n_treat).astype(np.float64)
    safe_cell = np.maximum(cell_count, 1.0)
    safe_treat = np.maximum(treat_count, 1.0)

    for _ in range(n_iter):
        resid = y - mu - cell_eff[c_idx]
        treat_sum = np.bincount(t_idx, weights=resid, minlength=n_treat)
        treat_eff = np.where(treat_count > 0, treat_sum / safe_treat, 0.0)

        resid = y - mu - treat_eff[t_idx]
        cell_sum = np.bincount(c_idx, weights=resid, minlength=n_cell)
        cell_eff = np.where(cell_count > 0, cell_sum / safe_cell, 0.0)

    return mu, cell_eff, treat_eff


def _drug_transfer(y_imp, obs_mask, test_pairs, lam, holdout_set):
    n_cell, _ = y_imp.shape
    non_holdout = np.array(
        [c for c in range(n_cell) if obs_mask[c].any() and c not in holdout_set],
        dtype=np.int64,
    )
    if non_holdout.size == 0:
        return np.zeros(len(test_pairs), dtype=np.float64), np.zeros(
            len(test_pairs), dtype=bool
        )

    test_by_cell: dict[int, list[tuple[int, int]]] = {}
    for i, (c, t) in enumerate(test_pairs):
        test_by_cell.setdefault(int(c), []).append((i, int(t)))

    y_pred = np.zeros(len(test_pairs), dtype=np.float64)
    covered = np.zeros(len(test_pairs), dtype=bool)

    for heldout_cell, items in test_by_cell.items():
        observed_drugs = np.where(obs_mask[heldout_cell])[0]
        target_drugs = np.array(sorted({t for _, t in items}), dtype=np.int64)
        if observed_drugs.size == 0 or target_drugs.size == 0:
            continue

        X = y_imp[np.ix_(non_holdout, observed_drugs)]
        Y = y_imp[np.ix_(non_holdout, target_drugs)]
        x_mean = X.mean(axis=0)
        y_mean = Y.mean(axis=0)
        Xc = X - x_mean
        Yc = Y - y_mean
        gram = Xc @ Xc.T
        coef = np.linalg.solve(gram + lam * np.eye(non_holdout.size), Yc)
        pred = y_mean + coef.T @ (Xc @ (y_imp[heldout_cell, observed_drugs] - x_mean))

        target_pos = {int(t): j for j, t in enumerate(target_drugs)}
        for test_i, t in items:
            y_pred[test_i] = pred[target_pos[t]]
            covered[test_i] = True

    return y_pred, covered


def _build_cell_features(bundle: SplitBundle):
    feat_names = ["cell_eval", "pdex", "pdex_pv", "pdex_fdr", "mean_expr_2k"]
    Xtr_list = []
    Xte_list = []
    for k in feat_names:
        Xtr = np.asarray(bundle.X_train_feature_blocks[k], dtype=np.float64)
        Xte = np.asarray(bundle.X_test_feature_blocks[k], dtype=np.float64)
        if k in ("pdex_pv", "pdex_fdr"):
            Xtr = -np.log10(np.clip(np.abs(Xtr), 1e-12, None))
            Xte = -np.log10(np.clip(np.abs(Xte), 1e-12, None))
        Xtr_list.append(Xtr)
        Xte_list.append(Xte)
    return np.concatenate(Xtr_list, axis=1), np.concatenate(Xte_list, axis=1)


def _ridge_fit_predict(
    Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, alpha: float
) -> np.ndarray:
    Xtr = np.asarray(Xtr, dtype=np.float64)
    ytr = np.asarray(ytr, dtype=np.float64)
    Xte = np.asarray(Xte, dtype=np.float64)

    if Xtr.shape[0] == 0:
        base = float(ytr.mean()) if ytr.size else 0.0
        return np.full(Xte.shape[0], base, dtype=np.float64)

    xm = Xtr.mean(axis=0)
    xs = Xtr.std(axis=0)
    xs = np.where(xs > 1e-8, xs, 1.0)
    Xtrz = (Xtr - xm) / xs
    Xtez = (Xte - xm) / xs

    ym = float(ytr.mean())
    yc = ytr - ym
    p = Xtrz.shape[1]
    A = Xtrz.T @ Xtrz + float(alpha) * np.eye(p, dtype=np.float64)
    b = Xtrz.T @ yc
    try:
        w = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        w = np.linalg.lstsq(A, b, rcond=1e-10)[0]
    return (ym + Xtez @ w).astype(np.float64)


def _ridge_fit_weights(
    Xtr: np.ndarray, ytr: np.ndarray, alpha: float
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    Xtr = np.asarray(Xtr, dtype=np.float64)
    ytr = np.asarray(ytr, dtype=np.float64)
    if Xtr.shape[0] == 0:
        p = Xtr.shape[1]
        return (
            np.zeros(p, dtype=np.float64),
            np.ones(p, dtype=np.float64),
            float(ytr.mean()) if ytr.size else 0.0,
            np.zeros(p, dtype=np.float64),
        )
    xm = Xtr.mean(axis=0)
    xs = Xtr.std(axis=0)
    xs = np.where(xs > 1e-8, xs, 1.0)
    Xtrz = (Xtr - xm) / xs
    ym = float(ytr.mean())
    yc = ytr - ym
    p = Xtrz.shape[1]
    A = Xtrz.T @ Xtrz + float(alpha) * np.eye(p, dtype=np.float64)
    b = Xtrz.T @ yc
    try:
        w = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        w = np.linalg.lstsq(A, b, rcond=1e-10)[0]
    return xm, xs, ym, w


def _eb_drug_bias(resid: np.ndarray, t_idx: np.ndarray, n_treat: int, lam: float):
    t_cnt = np.bincount(t_idx, minlength=n_treat).astype(np.float64)
    t_sum = np.bincount(t_idx, weights=resid, minlength=n_treat).astype(np.float64)
    t_mean = np.where(t_cnt > 0, t_sum / np.maximum(t_cnt, 1.0), 0.0)
    shrink = t_cnt / (t_cnt + float(lam))
    return shrink * t_mean


def _blend_weights(train_preds: np.ndarray, y: np.ndarray) -> np.ndarray:
    C = np.asarray(train_preds, dtype=np.float64)
    A = C.T @ C + 1e-4 * np.eye(C.shape[1], dtype=np.float64)
    b = C.T @ y
    try:
        w = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        w = np.linalg.lstsq(A, b, rcond=1e-10)[0]
    w = np.maximum(w, 0.0)
    s = float(w.sum())
    if s > 1e-12:
        w /= s
    else:
        w = np.array([0.24, 0.16, 0.24, 0.2, 0.16], dtype=np.float64)
    return w


def predict_sensitivity(bundle: SplitBundle) -> np.ndarray:
    """Return predicted growth_rate for test rows of one official split."""
    print(
        f"[predict] {bundle.split_name} start n_train={bundle.y_train.shape[0]}",
        flush=True,
    )

    y = np.asarray(bundle.y_train, dtype=np.float64)
    c_tr = bundle.c_train
    t_tr = bundle.t_train
    c_te = bundle.c_test
    t_te = bundle.t_test

    if lgb is not None:
        try:
            feat_names = ["cell_eval", "pdex", "pdex_pv", "pdex_fdr", "mean_expr_2k"]
            Xtr_list = []
            Xte_list = []
            for k in feat_names:
                Xtr_k = np.asarray(bundle.X_train_feature_blocks[k], dtype=np.float32)
                Xte_k = np.asarray(bundle.X_test_feature_blocks[k], dtype=np.float32)
                if k in ("pdex_pv", "pdex_fdr"):
                    Xtr_k = -np.log10(np.clip(np.abs(Xtr_k), 1e-12, None)).astype(
                        np.float32, copy=False
                    )
                    Xte_k = -np.log10(np.clip(np.abs(Xte_k), 1e-12, None)).astype(
                        np.float32, copy=False
                    )
                Xtr_list.append(Xtr_k)
                Xte_list.append(Xte_k)

            Xtr = np.concatenate(Xtr_list, axis=1).astype(np.float32, copy=False)
            Xte = np.concatenate(Xte_list, axis=1).astype(np.float32, copy=False)

            Xtr_all = np.hstack(
                [
                    Xtr,
                    c_tr.astype(np.float32).reshape(-1, 1),
                    t_tr.astype(np.float32).reshape(-1, 1),
                ]
            ).astype(np.float32, copy=False)
            Xte_all = np.hstack(
                [
                    Xte,
                    c_te.astype(np.float32).reshape(-1, 1),
                    t_te.astype(np.float32).reshape(-1, 1),
                ]
            ).astype(np.float32, copy=False)
            cat_idx = [Xtr_all.shape[1] - 2, Xtr_all.shape[1] - 1]

            print(f"[predict] {bundle.split_name} train lgb", flush=True)
            dtrain = lgb.Dataset(
                Xtr_all,
                label=y,
                categorical_feature=cat_idx,
                free_raw_data=True,
            )
            params = {
                "objective": "regression",
                "metric": "l2",
                "learning_rate": HP_LGB_LR,
                "num_leaves": HP_LGB_LEAVES,
                "min_data_in_leaf": HP_LGB_MIN_DATA,
                "feature_fraction": HP_LGB_FF,
                "bagging_fraction": HP_LGB_BF,
                "bagging_freq": 1,
                "lambda_l1": 0.0,
                "lambda_l2": HP_LGB_L2,
                "max_bin": 255,
                "verbosity": -1,
                "num_threads": 10,
                "seed": 2025,
            }
            booster = lgb.train(params, dtrain, num_boost_round=HP_LGB_ROUNDS)
            pred_tr_lgb = booster.predict(Xtr_all).astype(np.float64)
            pred_te_lgb = booster.predict(Xte_all).astype(np.float64)

            print(f"[predict] {bundle.split_name} train residual ridge", flush=True)
            resid = y - pred_tr_lgb

            mu_l = float(np.mean(y))
            c_cnt = np.bincount(c_tr, minlength=bundle.n_cells).astype(np.float64)
            t_cnt = np.bincount(t_tr, minlength=bundle.n_treatments).astype(np.float64)
            c_sum = np.bincount(
                c_tr, weights=y - mu_l, minlength=bundle.n_cells
            ).astype(np.float64)
            t_sum = np.bincount(
                t_tr, weights=y - mu_l, minlength=bundle.n_treatments
            ).astype(np.float64)
            c_eff = np.where(c_cnt > 0, c_sum / np.maximum(c_cnt, 1.0), 0.0)
            t_eff = np.where(t_cnt > 0, t_sum / np.maximum(t_cnt, 1.0), 0.0)

            Xtr_resid = np.column_stack(
                [
                    c_tr.astype(np.float64),
                    t_tr.astype(np.float64),
                    c_eff[c_tr],
                    t_eff[t_tr],
                    pred_tr_lgb,
                ]
            )
            Xte_resid = np.column_stack(
                [
                    c_te.astype(np.float64),
                    t_te.astype(np.float64),
                    c_eff[c_te],
                    t_eff[t_te],
                    pred_te_lgb,
                ]
            )
            pred_te_resid = _ridge_fit_predict(
                Xtr_resid,
                resid,
                Xte_resid,
                alpha=HP_RESID_ALPHA,
            )
            y_pred = pred_te_lgb + pred_te_resid
            y_pred = np.clip(y_pred, -1.0, 1.0).astype(np.float64, copy=False)
            print(f"[predict] {bundle.split_name} done lgb+resid", flush=True)
            return y_pred
        except Exception as e:
            print(
                f"[predict] {bundle.split_name} lgb path failed, fallback legacy: {type(e).__name__}: {e}",
                flush=True,
            )

    mu, cell_eff, treat_eff = _als_additive(
        y, c_tr, t_tr, bundle.n_cells, bundle.n_treatments, HP_ALS
    )
    additive_te = mu + cell_eff[c_te] + treat_eff[t_te]
    additive_tr = mu + cell_eff[c_tr] + treat_eff[t_tr]

    obs_mask = np.zeros((bundle.n_cells, bundle.n_treatments), dtype=bool)
    obs_mask[c_tr, t_tr] = True
    y_imp = mu + cell_eff[:, None] + treat_eff[None, :]
    y_imp[c_tr, t_tr] = y

    holdout_set = {int(c) for c in c_te}
    train_pairs = list(zip(c_tr.tolist(), t_tr.tolist()))
    transfer_tr, covered_tr = _drug_transfer(
        y_imp, obs_mask, train_pairs, lam=HP_DREG_LAM, holdout_set=holdout_set
    )
    test_pairs = list(zip(c_te.tolist(), t_te.tolist()))
    transfer_te, covered_te = _drug_transfer(
        y_imp, obs_mask, test_pairs, lam=HP_DREG_LAM, holdout_set=holdout_set
    )
    transfer_tr = np.where(covered_tr, transfer_tr, additive_tr)
    transfer_te = np.where(covered_te, transfer_te, additive_te)

    Xtr_feat, Xte_feat = _build_cell_features(bundle)
    feat_tr = _ridge_fit_predict(Xtr_feat, y, Xtr_feat, alpha=HP_RIDGE_ALPHA)
    feat_te = _ridge_fit_predict(Xtr_feat, y, Xte_feat, alpha=HP_RIDGE_ALPHA)

    bias_vec = _eb_drug_bias(y - feat_tr, t_tr, bundle.n_treatments, HP_BIAS_LAM)
    bias_tr = feat_tr + bias_vec[t_tr]
    bias_te = feat_te + bias_vec[t_te]

    tr_cells = np.unique(c_tr)
    cell_to_row = {int(c): i for i, c in enumerate(tr_cells.tolist())}
    R_all = np.full((tr_cells.size, bundle.n_treatments), np.nan, dtype=np.float64)
    for i in range(y.shape[0]):
        rr = cell_to_row[int(c_tr[i])]
        R_all[rr, int(t_tr[i])] = y[i] - bias_tr[i]

    d2d_tr = bias_tr.copy()
    d2d_te_pred = bias_te.copy()
    te_drugs = np.unique(t_te)
    for d_i, d in enumerate(te_drugs, start=1):
        if d_i == 1 or d_i == len(te_drugs) or (d_i % 10 == 0):
            print(
                f"[predict] {bundle.split_name} d2d drug {d_i}/{len(te_drugs)}",
                flush=True,
            )
        d = int(d)
        rows_d = np.where(~np.isnan(R_all[:, d]))[0]
        if rows_d.size < 10:
            continue
        other = np.array(
            [j for j in range(bundle.n_treatments) if j != d], dtype=np.int64
        )
        Xd = R_all[np.ix_(rows_d, other)]
        yd = R_all[rows_d, d]
        valid_cols = np.where(np.all(np.isfinite(Xd), axis=0))[0]
        if valid_cols.size == 0:
            continue
        Xd2 = Xd[:, valid_cols]
        other2 = other[valid_cols]
        xm_d, xs_d, ym_d, w_d = _ridge_fit_weights(Xd2, yd, HP_D2D_ALPHA)

        tr_idx_d = np.where(t_tr == d)[0]
        for j in tr_idx_d:
            cc = int(c_tr[j])
            rr = cell_to_row.get(cc)
            if rr is None:
                continue
            rv = R_all[rr, other2]
            if not np.all(np.isfinite(rv)):
                continue
            z = (rv - xm_d) / xs_d
            d2d_tr[j] = bias_tr[j] + HP_D2D_GAMMA * (ym_d + float(z @ w_d))

        for j in np.where(t_te == d)[0]:
            cc = int(c_te[j])
            rr = cell_to_row.get(cc)
            if rr is None:
                continue
            rv = R_all[rr, other2]
            if not np.all(np.isfinite(rv)):
                continue
            z = (rv - xm_d) / xs_d
            d2d_te_pred[j] = bias_te[j] + HP_D2D_GAMMA * (ym_d + float(z @ w_d))

    P_tr = np.column_stack([additive_tr, transfer_tr, feat_tr, bias_tr, d2d_tr])
    w = _blend_weights(P_tr, y)

    y_pred = (
        w[0] * additive_te
        + w[1] * transfer_te
        + w[2] * feat_te
        + w[3] * bias_te
        + w[4] * d2d_te_pred
    )
    print(f"[predict] {bundle.split_name} done blend_w={np.round(w, 4)}", flush=True)
    return np.asarray(y_pred, dtype=np.float64)
