"""
Evolvable OpenBind EV-A71 2A affinity predictor.
"""

from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        result = float(value)
    except Exception:
        return default
    if not np.isfinite(result):
        return default
    return result


def _descriptor_row(smiles: str) -> dict[str, float]:
    smiles_text = str(smiles).strip() if isinstance(smiles, str) else ""
    canonical_smiles = smiles_text.split(" ")[0]
    if not canonical_smiles:
        return {
            "mw": 330.0,
            "clogp": 2.5,
            "tpsa": 70.0,
            "hbd": 1.0,
            "hba": 5.0,
            "rot": 4.0,
            "rings": 2.0,
            "aromatic_rings": 1.0,
            "heavy_atoms": 24.0,
            "hetero_atoms": 6.0,
            "fraction_csp3": 0.3,
        }

    atom_tokens = re.findall(r"Br|Cl|[B-IK-Za-ik-z]", canonical_smiles)
    atomic_weights = {
        "H": 1.008,
        "B": 10.81,
        "C": 12.011,
        "c": 12.011,
        "N": 14.007,
        "n": 14.007,
        "O": 15.999,
        "o": 15.999,
        "F": 18.998,
        "P": 30.974,
        "S": 32.06,
        "s": 32.06,
        "Cl": 35.45,
        "Br": 79.904,
        "I": 126.90,
    }
    atom_counts = {atom: atom_tokens.count(atom) for atom in set(atom_tokens)}
    heavy_atoms = max(1, len(atom_tokens))
    hetero_atoms = sum(
        count for atom, count in atom_counts.items() if atom not in {"C", "c", "H"}
    )
    aromatic_atoms = sum(count for atom, count in atom_counts.items() if atom.islower())
    carbon_atoms = atom_counts.get("C", 0) + atom_counts.get("c", 0)
    halogens = (
        atom_counts.get("F", 0)
        + atom_counts.get("Cl", 0)
        + atom_counts.get("Br", 0)
        + atom_counts.get("I", 0)
    )
    mw = sum(
        atomic_weights.get(atom, 12.0) * count for atom, count in atom_counts.items()
    )
    ring_digits = len(set(re.findall(r"\d", canonical_smiles)))
    branch_count = canonical_smiles.count("(")
    hbd = (
        canonical_smiles.count("[nH]")
        + canonical_smiles.count("N")
        + canonical_smiles.count("O")
    )
    hba = (
        canonical_smiles.count("N")
        + canonical_smiles.count("n")
        + canonical_smiles.count("O")
        + canonical_smiles.count("o")
        + canonical_smiles.count("S")
        + canonical_smiles.count("s")
    )
    rot = max(
        0.0,
        canonical_smiles.count("C")
        + canonical_smiles.count("N")
        - ring_digits
        - branch_count,
    )
    tpsa = 12.0 * hba + 8.0 * hbd
    clogp = (
        0.045 * carbon_atoms
        + 0.18 * halogens
        - 0.12 * hetero_atoms
        + 0.05 * aromatic_atoms
    )

    return {
        "mw": _safe_float(mw, 330.0),
        "clogp": _safe_float(clogp, 2.5),
        "tpsa": _safe_float(tpsa, 70.0),
        "hbd": _safe_float(hbd, 1.0),
        "hba": _safe_float(hba, 5.0),
        "rot": _safe_float(rot, 4.0),
        "rings": _safe_float(ring_digits, 2.0),
        "aromatic_rings": _safe_float(aromatic_atoms / 6.0, 1.0),
        "heavy_atoms": _safe_float(heavy_atoms, 24.0),
        "hetero_atoms": _safe_float(hetero_atoms, 6.0),
        "fraction_csp3": _safe_float(
            (carbon_atoms - aromatic_atoms) / max(carbon_atoms, 1), 0.3
        ),
    }


def _centered(values: Iterable[float], midpoint: float, scale: float) -> np.ndarray:
    array = np.asarray(list(values), dtype=float)
    return (array - midpoint) / max(scale, 1e-9)


def predict_affinity(compounds: pd.DataFrame) -> pd.DataFrame:
    """
    Return structure-level affinity predictions in pKD-like units.

    Required input columns:
    - fragalysis_code
    - smiles

    Required output columns:
    - fragalysis_code
    - predicted_affinity
    """
    if "fragalysis_code" not in compounds.columns or "smiles" not in compounds.columns:
        raise ValueError("compounds must contain fragalysis_code and smiles columns")

    features = pd.DataFrame([_descriptor_row(str(s)) for s in compounds["smiles"]])

    ### >>> ARF-CHECKPOINT-START: subsite-gated local SAR with buried-polar frustration correction
    smiles_series = compounds["smiles"].astype(str)

    mw = features["mw"].to_numpy(dtype=float)
    ha = np.maximum(features["heavy_atoms"].to_numpy(dtype=float), 1.0)
    logp = features["clogp"].to_numpy(dtype=float)
    tpsa = features["tpsa"].to_numpy(dtype=float)
    rot = np.maximum(features["rot"].to_numpy(dtype=float), 0.0)
    arom = features["aromatic_rings"].to_numpy(dtype=float)
    ring = features["rings"].to_numpy(dtype=float)
    hba = features["hba"].to_numpy(dtype=float)
    hbd = features["hbd"].to_numpy(dtype=float)
    hetero = features["hetero_atoms"].to_numpy(dtype=float)
    csp3 = np.clip(features["fraction_csp3"].to_numpy(dtype=float), 0.0, 1.0)

    mw_c = _centered(mw, midpoint=341.0, scale=121.0)
    heavy_c = _centered(ha, midpoint=24.0, scale=8.4)
    rot_c = _centered(rot, midpoint=4.0, scale=3.9)
    arom_c = _centered(arom, midpoint=1.28, scale=1.10)
    ring_c = _centered(ring, midpoint=2.2, scale=1.68)
    hba_c = _centered(hba, midpoint=5.7, scale=3.0)
    hbd_c = _centered(hbd, midpoint=1.0, scale=1.40)
    hetero_c = _centered(hetero, midpoint=6.0, scale=3.7)
    csp3_c = _centered(csp3, midpoint=0.34, scale=0.22)

    imhb_pattern_count = (
        smiles_series.str.count(
            r"NCCO|OCCN|NCCN|NCO|OCN|NC=O|O=CN|NCOC|COCN|NCC=O|O=CCN|NCCS|SCCN"
        ).to_numpy(dtype=float)
        + smiles_series.str.count(
            r"n1cccc(?:[OoNn])c1|[OoNn]1ccccn1|ncc[oOnN]|[OoNn]ccn"
        ).to_numpy(dtype=float)
        + smiles_series.str.count(r"N1CCO|O1CCN|N1CCN|O1CCO|N1COC|O1CNC").to_numpy(
            dtype=float
        )
    )
    imhb_capacity = np.maximum(np.minimum(hbd, hba), 0.0)
    imhb_pairs = np.minimum(1.18 * imhb_pattern_count, imhb_capacity)
    imhb_ratio = np.clip(imhb_pairs / np.maximum(hba + hbd, 1.0), 0.0, 0.90)

    acid_cooh = smiles_series.str.contains(r"C\(=O\)O|\[O-\]", regex=True).to_numpy(
        dtype=float
    )
    acid_sulfonamide = smiles_series.str.contains(
        r"S\(=O\)\(=O\)N", regex=True
    ).to_numpy(dtype=float)
    acid_imide_like = smiles_series.str.contains(
        r"N[Cc]?\(=O\)|\[nH\]", regex=True
    ).to_numpy(dtype=float)
    base_tertiary_amine = smiles_series.str.contains(
        r"N\([^)]*\)\([^)]*\)", regex=True
    ).to_numpy(dtype=float)
    base_pyridine = smiles_series.str.contains(
        r"n1ccccc1|c1ncccc1|c1ccncc1", regex=True
    ).to_numpy(dtype=float)
    base_imidazole = smiles_series.str.contains(
        r"n1c[nH]cc1|[nH]1cncc1|n1cncc1", regex=True
    ).to_numpy(dtype=float)
    base_aniline = smiles_series.str.contains(
        r"Nc1|c1ccc\(N\)cc1|c1cc\(N\)ccc1", regex=True
    ).to_numpy(dtype=float)
    base_guanidine = smiles_series.str.contains(
        r"N=C\(N\)N|NC\(=N\)N", regex=True
    ).to_numpy(dtype=float)
    quaternary_n = smiles_series.str.contains(r"\[N\+\]", regex=True).to_numpy(
        dtype=float
    )

    pka_acid = np.maximum.reduce(
        [4.5 * acid_cooh, 6.5 * acid_sulfonamide, 8.2 * acid_imide_like]
    )
    pka_base = np.maximum.reduce(
        [
            9.8 * base_tertiary_amine,
            6.8 * base_imidazole,
            5.2 * base_pyridine,
            5.0 * base_aniline,
            12.5 * base_guanidine,
            12.0 * quaternary_n,
        ]
    )
    acidic_factor = np.where(pka_acid > 0.0, np.power(10.0, pka_acid - 7.4), 0.0)
    basic_factor = np.where(pka_base > 0.0, np.power(10.0, 7.4 - pka_base), 0.0)
    frac_neutral = np.clip(1.0 / (1.0 + acidic_factor + basic_factor), 3e-4, 1.0)

    logd = logp + np.log10(frac_neutral)
    logd_c = _centered(logd, midpoint=2.20, scale=1.55)
    ion_delta = np.maximum(logp - logd, 0.0)
    ion_penalty = np.clip(0.12 * ion_delta, 0.0, 0.34)

    ring5_hetero = smiles_series.str.count(r"[n,o,s][1-9][^ ]{2,18}[1-9]").to_numpy(
        dtype=float
    )
    alpha_carbonyl_hetero = smiles_series.str.count(
        r"[NO]C\(=O\)|C\(=O\)[NO]"
    ).to_numpy(dtype=float)
    topo_da_prox = (
        smiles_series.str.count(r"NCCO|OCCN|NCCN|NCO|OCN|NCC\(=O\)|OCC\(=O\)").to_numpy(
            dtype=float
        )
        + smiles_series.str.count(r"[nH]1cccc[n,o,s]1|[n,o,s]1cccc[nH]1").to_numpy(
            dtype=float
        )
        + smiles_series.str.count(r"N1CCO|O1CCN|N1CCN|O1CCO").to_numpy(dtype=float)
    )
    maskable_hetero = np.minimum(
        np.maximum(
            0.52 * ring5_hetero + 0.24 * alpha_carbonyl_hetero + 0.24 * topo_da_prox,
            0.0,
        ),
        np.maximum(hba + hbd, 0.0),
    )

    eff_tpsa = tpsa * (1.0 - 0.35 * imhb_ratio)
    topo_mask = 10.5 * np.minimum(maskable_hetero, 3.1)
    imhb_mask = 18.0 * np.minimum(imhb_pairs, 2.3)
    exposed_psa = np.maximum(0.0, eff_tpsa - np.maximum(topo_mask, imhb_mask))
    exp_tpsa_c = _centered(exposed_psa, midpoint=73.0, scale=34.0)

    mw_sat = np.minimum(mw, 432.0)
    mw_backbone = 5.46 + 0.75 * _centered(mw_sat, midpoint=341.0, scale=118.0)

    mw_win = np.exp(-np.square((mw - 350.0) / 166.0))
    logd_win = np.exp(-np.square((logd - 2.80) / 1.50))
    etpsa_win = np.exp(-np.square((exposed_psa - 75.0) / 34.0))
    arom_frac = np.clip(arom / ha, 0.0, 1.0)
    arom_win = np.exp(-np.square((arom_frac - 0.19) / 0.11))
    fsp3_win = np.exp(-np.square((csp3 - 0.34) / 0.20))
    mpo = mw_win * logd_win * etpsa_win * (0.80 + 0.14 * arom_win + 0.06 * fsp3_win)
    sar_gate = 0.69 + 0.31 * (0.60 * mw_win + 0.40 * mpo)

    typed_hbond = np.maximum(0.94 * hbd + 0.75 * hba - 1.86, 0.0)
    typed_pi = np.maximum(arom - 0.30, 0.0)
    typed_hyd = np.maximum(arom + 0.30 * ring + 0.28 * np.maximum(logd_c, 0.0), 0.0)
    typed_cat = smiles_series.str.contains(r"\[N\+\]|\[nH\]|n", regex=True).to_numpy(
        dtype=float
    )

    charge_proxy = (
        0.23 * np.minimum(hba / ha, 0.56)
        + 0.14 * np.minimum(hbd / ha, 0.30)
        + 0.10 * np.minimum(hetero / ha, 0.52)
        + 0.03 * np.maximum(np.abs(logd_c) - 0.50, 0.0)
        + 0.06 * np.minimum(ion_delta, 2.0)
    )
    charge_balance = np.clip(1.0 - 0.13 * np.abs(hba_c - hbd_c), 0.80, 1.06)

    typed_score = (
        1.76 * typed_hbond + 1.18 * typed_pi + 1.03 * typed_hyd + 1.05 * typed_cat
    )
    typed_eff = (typed_score * (0.92 + 1.10 * charge_proxy) * charge_balance) / (
        ha * (1.0 + 0.11 * rot)
    )
    ce_term = np.clip(
        0.248
        * typed_eff
        * (0.58 + 0.42 * mpo)
        * (0.89 + 0.11 * (1.0 - 0.30 * imhb_ratio)),
        -0.015,
        0.136,
    )

    nitrile = smiles_series.str.contains(r"C#N", regex=True).to_numpy(dtype=float)
    aldehyde = smiles_series.str.contains(r"[CHch]\(=O\)|CHO", regex=True).to_numpy(
        dtype=float
    )
    ketone = smiles_series.str.contains(r"C\(=O\)C", regex=True).to_numpy(dtype=float)
    amide = smiles_series.str.contains(r"C\(=O\)N|NC\(=O\)", regex=True).to_numpy(
        dtype=float
    )
    acid = acid_cooh
    sulfonamide = acid_sulfonamide
    tertiary_amine = base_tertiary_amine
    halogen_count = smiles_series.str.count(r"F|Cl|Br|I").to_numpy(dtype=float)

    hetero_arom = np.minimum(
        arom,
        (
            smiles_series.str.count(r"n").to_numpy(dtype=float)
            + 0.55 * smiles_series.str.count(r"o|s").to_numpy(dtype=float)
        )
        / 1.9,
    )
    phenyl_like = np.maximum(arom - hetero_arom, 0.0)
    fused_proxy = np.maximum(ring - arom, 0.0)
    hrb = np.clip(
        0.007 * hetero_arom
        + 0.004 * np.minimum(fused_proxy, 1.8)
        - 0.005 * np.maximum(phenyl_like - 1.4, 0.0),
        -0.02,
        0.04,
    )

    tri = np.maximum(0.0, 1.25 * rot / np.sqrt(ha) - 0.14 * ring - 0.10 * imhb_pairs)
    tri_penalty = 0.013 * np.maximum(tri - 2.0, 0.0)

    covalent_warhead = np.clip(nitrile + aldehyde + 0.6 * ketone, 0.0, 1.0)
    chemotype_term = np.clip(
        0.033 * nitrile
        + 0.011 * amide
        + 0.006 * tertiary_amine
        + 0.009 * aldehyde
        + 0.010 * np.minimum(halogen_count, 2.0)
        - 0.012 * np.maximum(halogen_count - 2.0, 0.0)
        - 0.024 * acid
        - 0.011 * sulfonamide
        + 0.012 * nitrile * logd_win
        + 0.006 * amide * etpsa_win
        + 0.014 * covalent_warhead * (0.75 + 0.25 * mpo)
        + 0.35 * hrb,
        -0.18,
        0.19,
    )

    internal_satisfaction = 1.46 * imhb_pairs
    unsat_polar = np.maximum(
        (hba + 0.72 * hbd)
        - (
            typed_hbond
            + 0.58 * typed_cat
            + 0.22 * charge_proxy * ha
            + internal_satisfaction
        ),
        0.0,
    )
    unsat_penalty = (
        0.0049
        * unsat_polar
        * (0.85 + 0.15 * np.maximum(exp_tpsa_c, 0.0))
        * (1.0 + 0.06 * np.maximum(tri - 2.0, 0.0))
    )
    exposed_hbd = np.maximum(hbd - 0.88 * hba - 0.16 * imhb_pairs, 0.0)
    exposed_hbd_penalty = (
        0.0108 * exposed_hbd * (0.66 + 0.34 * np.maximum(exp_tpsa_c, 0.0))
    )

    cc_hyd = np.maximum(
        arom + 0.55 * np.maximum(logd, 0.0) + 0.25 * np.minimum(halogen_count, 2.0), 0.0
    )
    cc_pi = np.maximum(arom + 0.30 * hetero_arom, 0.0)
    cc_ion = np.maximum(
        tertiary_amine + quaternary_n + 0.8 * acid + 0.5 * base_imidazole, 0.0
    )
    cat_proxy = np.minimum(nitrile + 0.7 * aldehyde + 0.5 * ketone, 1.8)
    # DEBUG: define missing cc_hb term used by cc_score; this prevents NameError and keeps contact-efficiency logic active
    cc_hb = np.maximum(typed_hbond + 0.45 * np.maximum(imhb_pairs, 0.0), 0.0)

    cc_score = (
        1.18 * cc_hb * (0.68 + 0.32 * np.maximum(logd_win, 0.0))
        + 0.52 * cc_hyd * (0.72 + 0.28 * mw_win)
        + 0.70 * cc_pi * (0.70 + 0.30 * arom_win)
        + 0.98 * cc_ion * (0.62 + 0.38 * np.maximum(exp_tpsa_c + 0.4, 0.0))
        + 0.42 * cat_proxy
    )
    cc_c = _centered(cc_score, midpoint=2.6, scale=1.7)

    # DEBUG: moved buried-polar block after cc_c is defined to fix UnboundLocalError while preserving intended frustration correction
    ### >>> ARF-CHECKPOINT-START: strengthen buried-polar frustration with exposed-polar and ionization coupling
    buried_polar_sat = np.clip(
        (
            typed_hbond
            + 0.56 * typed_cat
            + 0.78 * imhb_pairs
            + 0.30 * np.maximum(cc_c + 0.22, 0.0)
            + 0.20 * np.maximum((hba * hbd) / np.maximum(ha, 1.0) - 0.10, 0.0)
        )
        / np.maximum(0.92 + 0.76 * hba + 0.60 * hbd, 1.0),
        0.0,
        1.25,
    )
    buried_unsat = np.clip(1.0 - buried_polar_sat, 0.0, 1.0)
    desolv_term = 0.026 * np.clip(buried_polar_sat - 0.60, -0.42, 0.40)
    buried_unsat_penalty = (
        0.024
        * buried_unsat
        * (0.74 + 0.26 * np.maximum(exp_tpsa_c + 0.18, 0.0))
        * (0.86 + 0.14 * np.maximum(ion_delta - 0.8, 0.0))
    )
    ### <<< ARF-CHECKPOINT-END

    residual_term = sar_gate * (
        0.095 * heavy_c
        + 0.072 * logd_c
        - 0.0122 * np.square(logd_c)
        - 0.0235 * np.abs(exp_tpsa_c)
        - 0.0041 * np.square(exp_tpsa_c)
        - 0.039 * np.maximum(rot_c, 0.0)
        + 0.061 * arom_c
        + 0.036 * ring_c
        + 0.016 * (arom_win - 0.57)
        + 0.041 * (0.72 * hba_c + 0.28 * hbd_c)
        - 0.0080 * np.abs(hba_c - hbd_c)
        - 0.0055 * np.abs(hetero_c)
        - 0.0028 * np.abs(csp3_c)
        + 0.026 * (logd_c * (0.72 - np.abs(exp_tpsa_c)))
        - 0.0047 * (mw_c * np.maximum(rot_c, 0.0))
        + 0.0086 * arom_c * logd_win
        + 0.0060 * fsp3_win * (0.62 - np.abs(arom_c))
        + 0.021 * etpsa_win
        + 0.005 * np.minimum(maskable_hetero, 2.0) * (0.6 + 0.4 * mpo)
        + 0.014 * cc_c
        - 0.007 * ion_penalty
        - tri_penalty
    )

    lipe_proxy = (
        5.35 + 0.52 * (0.095 * heavy_c) + ce_term + chemotype_term + 0.010 * cc_c
    ) - logd
    lipe_win = np.exp(-np.square((lipe_proxy - 5.0) / 1.50))
    lipe_term = 0.021 * lipe_win * (0.72 + 0.28 * mpo)

    s1_proxy = np.clip(
        0.84 * np.maximum(typed_hbond, 0.0)
        + 0.34 * amide
        + 0.31 * nitrile
        + 0.24 * np.minimum(hba, 7.0) / 7.0
        - 0.17 * acid
        - 0.11 * np.maximum(exposed_psa - 98.0, 0.0) / 44.0,
        0.0,
        3.5,
    )
    s2_proxy = np.clip(
        0.64 * np.maximum(typed_hyd, 0.0)
        + 0.58 * np.maximum(typed_pi, 0.0)
        + 0.22 * np.minimum(halogen_count, 2.0)
        + 0.15 * np.minimum(np.maximum(ring - arom, 0.0), 2.0)
        - 0.16 * np.maximum(rot - 8.0, 0.0) / 4.0,
        0.0,
        3.4,
    )
    s1p_proxy = np.clip(
        0.56 * np.maximum(typed_cat, 0.0)
        + 0.26 * base_imidazole
        + 0.18 * tertiary_amine
        + 0.14 * np.minimum(hbd, 3.0) / 3.0
        - 0.12 * np.maximum(np.abs(logd_c) - 1.8, 0.0),
        0.0,
        2.9,
    )
    ### >>> ARF-CHECKPOINT-START: subsite-gated local SAR with stronger engagement dynamic range and ion-lipophilicity coupling
    hyd_enclosure = np.clip(
        0.47 * np.maximum(typed_hyd - 0.94, 0.0)
        + 0.31 * np.maximum(typed_pi - 0.38, 0.0)
        + 0.20 * np.minimum(np.maximum(ring - arom, 0.0), 2.1)
        + 0.16 * np.minimum(halogen_count, 2.0),
        0.0,
        2.55,
    )
    polarity_fit = np.exp(-np.square((exposed_psa - 77.0) / 32.0))
    ion_lip_fit = np.exp(-np.square((logd - 2.72) / 1.18))
    ion_gate = np.clip(0.70 + 0.30 * frac_neutral, 0.68, 1.02)
    engage_shape = np.exp(-np.square((rot - 4.8) / 3.2))
    contact_richness = np.clip(
        1.03 * s1_proxy
        + 0.76 * s2_proxy
        + 0.57 * s1p_proxy
        + 0.61 * np.maximum(cc_c + 0.14, 0.0)
        + 0.45 * np.maximum(ce_term * 10.0, 0.0)
        + 0.12 * np.maximum(mpo - 0.20, 0.0)
        + 0.10 * covalent_warhead * (0.72 + 0.28 * np.maximum(cc_c, 0.0))
        + 0.10 * hyd_enclosure * polarity_fit
        + 0.06 * hyd_enclosure * ion_lip_fit * ion_gate
        + 0.05 * engage_shape * (0.70 + 0.30 * np.maximum(mpo, 0.0)),
        0.0,
        6.95,
    )
    contact_gate = 1.0 / (1.0 + np.exp(-(1.42 * (contact_richness - 1.95))))
    ### <<< ARF-CHECKPOINT-END

    ketoamide_like = smiles_series.str.contains(
        r"N[Cc]?\(=O\)C\(=O\)|C\(=O\)C\(=O\)N|NC\(=O\)C\(=O\)", regex=True
    ).to_numpy(dtype=float)
    heteroaryl_fused = np.clip(
        np.minimum(np.maximum(ring - arom, 0.0), 2.0), 0.0, 2.0
    ) * np.clip(
        (
            smiles_series.str.count(r"n").to_numpy(dtype=float)
            + 0.5 * smiles_series.str.count(r"o|s").to_numpy(dtype=float)
        )
        / 2.0,
        0.0,
        1.2,
    )
    aromatic_overload = np.clip((arom_frac - 0.34) / 0.22, 0.0, 1.0)
    flex_over = np.clip((rot - 6.9) / 3.1, 0.0, 1.0)
    neutral_balance = np.exp(-np.square((logd - 2.75) / 1.55)) * np.exp(
        -np.square((exposed_psa - 76.0) / 35.0)
    )

    ### >>> ARF-CHECKPOINT-START: local SAR fragment table with balanced aromatic-flexibility and polarity satisfaction gating
    polar_balance = np.exp(-np.square((hba - 5.25) / 2.55)) * np.exp(
        -np.square((hbd - 1.12) / 1.36)
    )
    unsat_gate = np.exp(-np.square((buried_unsat - 0.33) / 0.29))
    ring_flex_balance = np.clip((ring - 0.47 * rot) / 2.02, 0.0, 1.2)
    logd_frag_fit = np.exp(-np.square((logd - 2.72) / 1.22))
    amide_nitrile_synergy = (
        amide * nitrile * np.exp(-np.square((exposed_psa - 78.0) / 29.0))
    )
    halogen_eff = np.clip(
        np.minimum(halogen_count, 2.0) / np.maximum(np.sqrt(ha), 1.0), 0.0, 0.36
    )
    aromatic_drag = np.clip((arom_frac - 0.31) / 0.21, 0.0, 1.0)
    frag_raw = (
        0.248 * nitrile
        + 0.194 * ketoamide_like
        + 0.104 * amide * np.exp(-np.square((exposed_psa - 78.4) / 29.4))
        + 0.088 * np.minimum(hetero_arom, 2.2) / 2.2
        + 0.070 * np.minimum(heteroaryl_fused, 1.8)
        + 0.053 * ring_flex_balance
        + 0.050 * polar_balance
        + 0.033 * logd_frag_fit
        + 0.034 * amide_nitrile_synergy
        + 0.020 * halogen_eff
        - 0.084 * aromatic_overload
        - 0.074 * flex_over
        - 0.060 * acid
        - 0.045 * np.clip((exposed_psa - 101.0) / 38.0, 0.0, 1.0)
        - 0.015 * aromatic_drag
    )
    frag_bonus = np.clip(
        0.076
        * contact_gate
        * neutral_balance
        * unsat_gate
        * ion_gate
        * frag_raw
        * (0.91 + 0.09 * mpo),
        -0.020,
        0.046,
    )
    ### <<< ARF-CHECKPOINT-END

    ### >>> ARF-CHECKPOINT-START: RSAE final calibration with engagement uplift and softened anti-size penalty to protect rank signal
    rsae = (
        mw_backbone
        + residual_term
        + ce_term
        + chemotype_term
        + lipe_term
        + desolv_term
        - unsat_penalty
        - buried_unsat_penalty
        - exposed_hbd_penalty
        - ion_penalty
        - 0.0229 * np.maximum(logd - 2.93, 0.0)
    )
    engage_uplift = np.clip(contact_richness - 1.88, 0.0, 3.1) / 3.1
    lipophilicity_fit = np.exp(-np.square((logd - 2.72) / 1.24))
    size_anti_bias = np.clip((mw - 438.0) / 165.0, 0.0, 1.0)
    predicted = (
        rsae
        + 0.100 * cc_c
        + 0.0136 * hrb
        - 0.0056 * tri_penalty
        + 0.0024 * engage_uplift * ion_gate
        + frag_bonus
        + 0.0011
        * contact_gate
        * ion_gate
        * np.clip(lipophilicity_fit - 0.45, -0.18, 0.48)
        + 0.0060 * contact_gate * np.clip(polar_balance - 0.40, -0.21, 0.55)
        + 0.0036
        * contact_gate
        * ion_gate
        * np.clip(neutral_balance - 0.30, -0.15, 0.56)
        + 0.0018 * np.clip(hyd_enclosure - 0.82, 0.0, 1.52)
        + 0.0010 * contact_gate * np.clip(ring_flex_balance - 0.30, 0.0, 0.94)
        - 0.0030 * np.maximum(ion_delta - 1.05, 0.0)
        - 0.0024 * np.maximum(np.abs(exp_tpsa_c) - 1.70, 0.0)
        - 0.0022 * np.maximum(aromatic_overload - 0.51, 0.0)
        - 0.0013 * np.maximum(flex_over - 0.41, 0.0)
        - 0.0014 * np.maximum(buried_unsat - 0.66, 0.0)
        - 0.0007 * size_anti_bias
    )
    ### <<< ARF-CHECKPOINT-END
    predicted = np.clip(predicted, 3.0, 10.0)
    ### <<< ARF-CHECKPOINT-END

    return pd.DataFrame(
        {
            "fragalysis_code": compounds["fragalysis_code"].astype(str).to_numpy(),
            "predicted_affinity": np.asarray(predicted, dtype=float),
        }
    )
