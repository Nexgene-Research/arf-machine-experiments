"""OpenBind affinity predictor with Arf Machine refined parameters baked in.

Generated from the Arf Machine hyperparameter-refined predictor.
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
    ha = np.maximum(features["heavy_atoms"].to_numpy(dtype=float), 1.213485472977653)
    logp = features["clogp"].to_numpy(dtype=float)
    tpsa = features["tpsa"].to_numpy(dtype=float)
    rot = np.maximum(features["rot"].to_numpy(dtype=float), 0.023458560925294643)
    arom = features["aromatic_rings"].to_numpy(dtype=float)
    ring = features["rings"].to_numpy(dtype=float)
    hba = features["hba"].to_numpy(dtype=float)
    hbd = features["hbd"].to_numpy(dtype=float)
    hetero = features["hetero_atoms"].to_numpy(dtype=float)
    csp3 = np.clip(features["fraction_csp3"].to_numpy(dtype=float), -0.017266172360911533, 1.0450992152586298)

    mw_c = _centered(mw, midpoint=332.78452494556893, scale=105.210528433645)
    heavy_c = _centered(ha, midpoint=19.06691682521492, scale=8.583664964447252)
    rot_c = _centered(rot, midpoint=3.3178930173757104, scale=4.607883312239936)
    arom_c = _centered(arom, midpoint=1.427916903142194, scale=1.343949309267032)
    ring_c = _centered(ring, midpoint=2.454263829988568, scale=1.7658391610051711)
    hba_c = _centered(hba, midpoint=4.703092323886426, scale=3.596861826080691)
    hbd_c = _centered(hbd, midpoint=1.1538518523603327, scale=1.149770587644469)
    hetero_c = _centered(hetero, midpoint=6.035030699764558, scale=2.9959294549867947)
    csp3_c = _centered(csp3, midpoint=0.26795284843294037, scale=0.20408414369713382)

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
    imhb_capacity = np.maximum(np.minimum(hbd, hba), 0.03029532796781592)
    imhb_pairs = np.minimum(0.8986597855443015 * imhb_pattern_count, imhb_capacity)
    imhb_ratio = np.clip(imhb_pairs / np.maximum(hba + hbd, 1.2169666681524567), 0.03583994458091428, 0.7188563310434178)

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
        [4.6341658829276575 * acid_cooh, 8.12465708556217 * acid_sulfonamide, 9.993203180228482 * acid_imide_like]
    )
    pka_base = np.maximum.reduce(
        [
            9.566322809027575 * base_tertiary_amine,
            6.573555828974753 * base_imidazole,
            5.712397688976973 * base_pyridine,
            4.940568157132973 * base_aniline,
            9.606652016194644 * base_guanidine,
            9.412519682268162 * quaternary_n,
        ]
    )
    acidic_factor = np.where(pka_acid > -0.039293411074423516, np.power(8.366397055803425, pka_acid - 6.987141221609433), 0.017599015851836302)
    basic_factor = np.where(pka_base > 0.0040862326955286654, np.power(12.168742194543, 5.616528686199118 - pka_base), -0.04633819835936784)
    frac_neutral = np.clip(1.2230615979552493 / (1.122314927450877 + acidic_factor + basic_factor), 0.002077986634893862, 0.7894247903379735)

    logd = logp + np.log10(frac_neutral)
    logd_c = _centered(logd, midpoint=1.717882583640754, scale=1.5413802625629638)
    ion_delta = np.maximum(logp - logd, 0.008180780013542442)
    ion_penalty = np.clip(0.09958768905429563 * ion_delta, -0.04029392146311155, 0.27937164568384965)

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
            0.6192902370330099 * ring5_hetero + 0.21495790793412323 * alpha_carbonyl_hetero + 0.22406685771444115 * topo_da_prox,
            -0.014557613090313741,
        ),
        np.maximum(hba + hbd, -0.02339106902233775),
    )

    eff_tpsa = tpsa * (1.1464736188679026 - 0.4039960073116805 * imhb_ratio)
    topo_mask = 10.423182888140165 * np.minimum(maskable_hetero, 2.56131509127057)
    imhb_mask = 19.740394441391043 * np.minimum(imhb_pairs, 2.3096007429362926)
    exposed_psa = np.maximum(0.01243827716703979, eff_tpsa - np.maximum(topo_mask, imhb_mask))
    exp_tpsa_c = _centered(exposed_psa, midpoint=87.00816226939146, scale=33.083507074086775)

    mw_sat = np.minimum(mw, 346.2669664758953)
    mw_backbone = 5.263932075502313 + 0.5842765738340986 * _centered(mw_sat, midpoint=417.40511902580687, scale=138.51910342034583)

    mw_win = np.exp(-np.square((mw - 308.3270262061832) / 182.1551124653912))
    logd_win = np.exp(-np.square((logd - 3.305993135059309) / 1.1590443304107225))
    etpsa_win = np.exp(-np.square((exposed_psa - 78.81125662323397) / 40.5675831374602))
    arom_frac = np.clip(arom / ha, -0.04994179885478944, 1.2493090883265328)
    arom_win = np.exp(-np.square((arom_frac - 0.16767403780744264) / 0.10523324965182554))
    fsp3_win = np.exp(-np.square((csp3 - 0.26588036946119786) / 0.19199798491831505))
    mpo = mw_win * logd_win * etpsa_win * (0.6525457442164925 + 0.10695035272389673 * arom_win + 0.04826424482260582 * fsp3_win)
    sar_gate = 0.86241861425091 + 0.2590848017906366 * (0.5366340261057061 * mw_win + 0.43490064815147084 * mpo)

    typed_hbond = np.maximum(0.7128502847212432 * hbd + 0.8450081272171871 * hba - 2.0323109922683686, 0.009444071913323164)
    typed_pi = np.maximum(arom - 0.24657363336602284, -0.04049075094229701)
    typed_hyd = np.maximum(arom + 0.3551893022975167 * ring + 0.21173098767006082 * np.maximum(logd_c, 0.04818914463187869), 0.03714450859823539)
    typed_cat = smiles_series.str.contains(r"\[N\+\]|\[nH\]|n", regex=True).to_numpy(
        dtype=float
    )

    charge_proxy = (
        0.2170864003991932 * np.minimum(hba / ha, 0.4626611735739197)
        + 0.15679893756314053 * np.minimum(hbd / ha, 0.3658609087762253)
        + 0.09814652114750573 * np.minimum(hetero / ha, 0.6249555325248266)
        + 0.03223144191642066 * np.maximum(np.abs(logd_c) - 0.6076209104498775, -0.04084843851545422)
        + 0.06350588666910784 * np.minimum(ion_delta, 2.477556876019334)
    )
    charge_balance = np.clip(1.1659380273243152 - 0.11788955092177038 * np.abs(hba_c - hbd_c), 0.7480459131458389, 0.9924805986539535)

    typed_score = (
        2.012814569224099 * typed_hbond + 1.0397707487945582 * typed_pi + 0.8547524938205306 * typed_hyd + 1.14059971254951 * typed_cat
    )
    typed_eff = (typed_score * (0.8460122430477759 + 1.2259870607050503 * charge_proxy) * charge_balance) / (
        ha * (1.1667684756611683 + 0.13288754243711184 * rot)
    )
    ce_term = np.clip(
        0.2505131902645668
        * typed_eff
        * (0.5265304124762811 + 0.4391528481510463 * mpo)
        * (0.7456941724264656 + 0.08927185401019455 * (0.9694973081395826 - 0.3586405612524302 * imhb_ratio)),
        -0.011518589319220184,
        0.15039491758224335,
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
            + 0.6713793203078017 * smiles_series.str.count(r"o|s").to_numpy(dtype=float)
        )
        / 1.5468890164033469,
    )
    phenyl_like = np.maximum(arom - hetero_arom, -0.014068433091699355)
    fused_proxy = np.maximum(ring - arom, -0.007901911064907214)
    hrb = np.clip(
        0.006275883098021008 * hetero_arom
        + 0.00456286018199392 * np.minimum(fused_proxy, 1.4626169671956977)
        - 0.004106660905070166 * np.maximum(phenyl_like - 1.4197182659323597, 0.03531731041594069),
        -0.016968430620232202,
        0.03509729349967969,
    )

    tri = np.maximum(0.0359014469148906, 1.4701155091378562 * rot / np.sqrt(ha) - 0.13090812313287967 * ring - 0.09017300468545693 * imhb_pairs)
    tri_penalty = 0.012344277646542019 * np.maximum(tri - 1.6726692196743966, 0.004551554626992106)

    covalent_warhead = np.clip(nitrile + aldehyde + 0.5421016149073443 * ketone, 0.013526087032311333, 0.81417559688501)
    chemotype_term = np.clip(
        0.02606709584917596 * nitrile
        + 0.013431020866043185 * amide
        + 0.0058832697598735955 * tertiary_amine
        + 0.01035468305309982 * aldehyde
        + 0.007867759438174547 * np.minimum(halogen_count, 1.5331566255725935)
        - 0.01410508500095011 * np.maximum(halogen_count - 1.5591350876429277, -0.027810201640191767)
        - 0.029087846926020886 * acid
        - 0.013531604524912887 * sulfonamide
        + 0.009123275123149971 * nitrile * logd_win
        + 0.0071101979113678436 * amide * etpsa_win
        + 0.013586217912077349 * covalent_warhead * (0.8515548539066903 + 0.28476155309854645 * mpo)
        + 0.4270479408553301 * hrb,
        -0.14858705080989898,
        0.18364731970050674,
    )

    internal_satisfaction = 1.6593169253223004 * imhb_pairs
    unsat_polar = np.maximum(
        (hba + 0.651745549860765 * hbd)
        - (
            typed_hbond
            + 0.6960485260229327 * typed_cat
            + 0.26265831310446974 * charge_proxy * ha
            + internal_satisfaction
        ),
        0.03541945963263798,
    )
    unsat_penalty = (
        0.006373354103466073
        * unsat_polar
        * (0.8778698425322965 + 0.17331765502313506 * np.maximum(exp_tpsa_c, 0.023671153691879852))
        * (0.9732856232869661 + 0.057562979979275904 * np.maximum(tri - 2.2542848691207036, -0.0021848345081642894))
    )
    exposed_hbd = np.maximum(hbd - 0.6915302376487864 * hba - 0.165389713854596 * imhb_pairs, -0.04461564495618206)
    exposed_hbd_penalty = (
        0.011620573506062057 * exposed_hbd * (0.6619424057385573 + 0.37883204758781197 * np.maximum(exp_tpsa_c, 0.013743344941907877))
    )

    cc_hyd = np.maximum(
        arom + 0.42216445302489 * np.maximum(logd, -0.03818571142549132) + 0.19749481752818224 * np.minimum(halogen_count, 1.8093028803673574), 0.034549953559820275
    )
    cc_pi = np.maximum(arom + 0.31817548319879196 * hetero_arom, -0.04418877202246418)
    cc_ion = np.maximum(
        tertiary_amine + quaternary_n + 0.9747740351545217 * acid + 0.5557215145829879 * base_imidazole, 0.008930097113265222
    )
    cat_proxy = np.minimum(nitrile + 0.5297923610120187 * aldehyde + 0.39284309638611764 * ketone, 1.441382952972194)
    # DEBUG: define missing cc_hb term used by cc_score; this prevents NameError and keeps contact-efficiency logic active
    cc_hb = np.maximum(typed_hbond + 0.5442618572466351 * np.maximum(imhb_pairs, -0.000550898820494021), -0.045371648281047484)

    cc_score = (
        1.0617910771216377 * cc_hb * (0.8246856952868352 + 0.30474519143772066 * np.maximum(logd_win, -0.0006854838215988477))
        + 0.5821237703165969 * cc_hyd * (0.6246730217362527 + 0.3094962707109382 * mw_win)
        + 0.7358199645701134 * cc_pi * (0.7300357218918896 + 0.34536227994970153 * arom_win)
        + 0.8967116831872728 * cc_ion * (0.682332416960062 + 0.39421874013089486 * np.maximum(exp_tpsa_c + 0.3345624109389179, -0.04003313751430028))
        + 0.37376264347502547 * cat_proxy
    )
    cc_c = _centered(cc_score, midpoint=2.4287467498847723, scale=1.3348507295981509)

    # DEBUG: moved buried-polar block after cc_c is defined to fix UnboundLocalError while preserving intended frustration correction
    ### >>> ARF-CHECKPOINT-START: strengthen buried-polar frustration with exposed-polar and ionization coupling
    buried_polar_sat = np.clip(
        (
            typed_hbond
            + 0.5049851108123091 * typed_cat
            + 0.6959217790439601 * imhb_pairs
            + 0.30614052155821025 * np.maximum(cc_c + 0.17173490251448945, 0.04305895978232147)
            + 0.18479784297908955 * np.maximum((hba * hbd) / np.maximum(ha, 0.9084573266765367) - 0.09385868490663256, 0.045579875338174314)
        )
        / np.maximum(0.9338767841870539 + 0.6996842951965718 * hba + 0.6357871040119074 * hbd, 0.7680835824165518),
        0.025663694510065612,
        1.2122468233004027,
    )
    buried_unsat = np.clip(0.7781266657500666 - buried_polar_sat, 0.010641866983914826, 1.1364171299041133)
    desolv_term = 0.02402688511759296 * np.clip(buried_polar_sat - 0.45747065412577725, -0.44050637219753563, 0.44637589010639733)
    buried_unsat_penalty = (
        0.02258336677322769
        * buried_unsat
        * (0.6486923903269524 + 0.3167972120576375 * np.maximum(exp_tpsa_c + 0.1888402416766213, 0.008408862946931233))
        * (0.7626911030840314 + 0.1677283094082334 * np.maximum(ion_delta - 0.9863690509870462, 0.024081042674941867))
    )
    ### <<< ARF-CHECKPOINT-END

    residual_term = sar_gate * (
        0.1084116848600651 * heavy_c
        + 0.08332495648772786 * logd_c
        - 0.009806920682780919 * np.square(logd_c)
        - 0.018582403274103158 * np.abs(exp_tpsa_c)
        - 0.004706593156536247 * np.square(exp_tpsa_c)
        - 0.035826910698426896 * np.maximum(rot_c, 0.031894474226386)
        + 0.060209548905528286 * arom_c
        + 0.03030669573391039 * ring_c
        + 0.015590302256265069 * (arom_win - 0.4757216848195521)
        + 0.04973888999672312 * (0.8120596621962641 * hba_c + 0.2584867351406084 * hbd_c)
        - 0.0069976160727194345 * np.abs(hba_c - hbd_c)
        - 0.005398492055149903 * np.abs(hetero_c)
        - 0.0035145245277590654 * np.abs(csp3_c)
        + 0.020595197246117313 * (logd_c * (0.8925918551438613 - np.abs(exp_tpsa_c)))
        - 0.0055022916239465145 * (mw_c * np.maximum(rot_c, 0.04658620926721791))
        + 0.010415167511748213 * arom_c * logd_win
        + 0.004293869961837246 * fsp3_win * (0.4893814625254735 - np.abs(arom_c))
        + 0.022735304255101275 * etpsa_win
        + 0.004722765369516324 * np.minimum(maskable_hetero, 2.3683731596618456) * (0.6587486298235227 + 0.32372754878579707 * mpo)
        + 0.016627050015685076 * cc_c
        - 0.005858096297738261 * ion_penalty
        - tri_penalty
    )

    lipe_proxy = (
        4.723514504024911 + 0.5329568178429491 * (0.11434878787446688 * heavy_c) + ce_term + chemotype_term + 0.010431213506261722 * cc_c
    ) - logd
    lipe_win = np.exp(-np.square((lipe_proxy - 5.18532217499403) / 1.4729827867387133))
    lipe_term = 0.01648271987221793 * lipe_win * (0.7779849150406243 + 0.263823907256643 * mpo)

    s1_proxy = np.clip(
        0.6553671337169753 * np.maximum(typed_hbond, -0.02717688006518746)
        + 0.37950219572726074 * amide
        + 0.3250308591350252 * nitrile
        + 0.21043725122364595 * np.minimum(hba, 5.876094896002277) / 6.53395431707812
        - 0.17110134780524644 * acid
        - 0.09682656815130175 * np.maximum(exposed_psa - 111.47256408196, 0.027727142047139487) / 53.778017254645256,
        0.016216448845756584,
        3.5341606170018007,
    )
    s2_proxy = np.clip(
        0.7817502084117812 * np.maximum(typed_hyd, 0.013197589203640845)
        + 0.6795373100568969 * np.maximum(typed_pi, -0.03377982049087586)
        + 0.2353435548901897 * np.minimum(halogen_count, 1.550193778644782)
        + 0.17782061329974408 * np.minimum(np.maximum(ring - arom, 0.03894885351218009), 2.2484159111538973)
        - 0.15401253122423375 * np.maximum(rot - 8.138549356380153, 0.025133384140234288) / 3.1670217456235052,
        -0.043478605641654615,
        3.515372710181965,
    )
    s1p_proxy = np.clip(
        0.6349425592095317 * np.maximum(typed_cat, 0.01846342929949616)
        + 0.23013345183612216 * base_imidazole
        + 0.1811876340828188 * tertiary_amine
        + 0.16839419997509242 * np.minimum(hbd, 3.6885617536443713) / 2.7141464767043497
        - 0.1252869675786976 * np.maximum(np.abs(logd_c) - 1.8296324328185505, -0.007352433558570407),
        -0.006827545564286164,
        3.428801361861398,
    )
    ### >>> ARF-CHECKPOINT-START: subsite-gated local SAR with stronger engagement dynamic range and ion-lipophilicity coupling
    hyd_enclosure = np.clip(
        0.5598086932099993 * np.maximum(typed_hyd - 1.0572642410800503, 0.006280458907441622)
        + 0.24401687341834902 * np.maximum(typed_pi - 0.4595405732710814, 0.009212676049059688)
        + 0.1822007226365244 * np.minimum(np.maximum(ring - arom, 0.01230310801040705), 2.547800237655179)
        + 0.12594427094606947 * np.minimum(halogen_count, 2.200326554098997),
        0.041932289350169696,
        2.2589711667109564,
    )
    polarity_fit = np.exp(-np.square((exposed_psa - 94.72831872859658) / 29.045839936947562))
    ion_lip_fit = np.exp(-np.square((logd - 2.7795677458035195) / 1.4223638671769805))
    ion_gate = np.clip(0.5499545018339375 + 0.2878348681001922 * frac_neutral, 0.5235037476203575, 0.8204367366904187)
    engage_shape = np.exp(-np.square((rot - 3.87800059301733) / 3.8469561325321684))
    contact_richness = np.clip(
        1.0555871448972296 * s1_proxy
        + 0.7094980751590025 * s2_proxy
        + 0.4383409533177699 * s1p_proxy
        + 0.661906629442577 * np.maximum(cc_c + 0.16939775332528476, 0.03953759802133901)
        + 0.47923794708194417 * np.maximum(ce_term * 10.53210060880618, -0.033183530175840456)
        + 0.14891211569221974 * np.maximum(mpo - 0.23882399268773505, -0.014895028837676136)
        + 0.07642160532520285 * covalent_warhead * (0.7879799701692449 + 0.2678721354283585 * np.maximum(cc_c, 0.04742550418351377))
        + 0.08967211666012817 * hyd_enclosure * polarity_fit
        + 0.05309974105864441 * hyd_enclosure * ion_lip_fit * ion_gate
        + 0.061973967486129086 * engage_shape * (0.5676843486236731 + 0.3500239470549452 * np.maximum(mpo, -0.03661464684893224)),
        -0.016848151437579204,
        6.196670656273982,
    )
    contact_gate = 0.8325720189293653 / (1.200951545312691 + np.exp(-(1.7148219302127592 * (contact_richness - 2.1361185212912432))))
    ### <<< ARF-CHECKPOINT-END

    ketoamide_like = smiles_series.str.contains(
        r"N[Cc]?\(=O\)C\(=O\)|C\(=O\)C\(=O\)N|NC\(=O\)C\(=O\)", regex=True
    ).to_numpy(dtype=float)
    heteroaryl_fused = np.clip(
        np.minimum(np.maximum(ring - arom, 0.007349937614080968), 2.11270713628219), -0.007610714986905134, 2.3012753257865803
    ) * np.clip(
        (
            smiles_series.str.count(r"n").to_numpy(dtype=float)
            + 0.5373708121772427 * smiles_series.str.count(r"o|s").to_numpy(dtype=float)
        )
        / 2.034878613091185,
        -0.0499324014935805,
        1.3078833637476515,
    )
    aromatic_overload = np.clip((arom_frac - 0.3524160257968846) / 0.23299408474782768, 0.04708813241666712, 0.9147771652898965)
    flex_over = np.clip((rot - 6.204865328605439) / 2.4147868434972355, -0.022071223525611475, 1.2190048840318328)
    neutral_balance = np.exp(-np.square((logd - 2.825282563842893) / 1.4945787689101413)) * np.exp(
        -np.square((exposed_psa - 64.02419353049399) / 28.786208697141834)
    )

    ### >>> ARF-CHECKPOINT-START: local SAR fragment table with balanced aromatic-flexibility and polarity satisfaction gating
    polar_balance = np.exp(-np.square((hba - 5.999709951765213) / 1.9316628665589306)) * np.exp(
        -np.square((hbd - 1.3907789536008626) / 1.2488360645312162)
    )
    unsat_gate = np.exp(-np.square((buried_unsat - 0.33106163994886734) / 0.22005898019720257))
    ring_flex_balance = np.clip((ring - 0.5134527324562725 * rot) / 2.2109431539770967, -0.0057802445651884225, 0.9653124047124166)
    logd_frag_fit = np.exp(-np.square((logd - 3.052206956177743) / 1.3544902252226074))
    amide_nitrile_synergy = (
        amide * nitrile * np.exp(-np.square((exposed_psa - 68.09674181354724) / 28.976516588345827))
    )
    halogen_eff = np.clip(
        np.minimum(halogen_count, 1.5574140568496295) / np.maximum(np.sqrt(ha), 1.2106568519072323), 0.02627768375207614, 0.29372152893476705
    )
    aromatic_drag = np.clip((arom_frac - 0.38479245182311983) / 0.18749401252179373, 0.03676781443334313, 0.8126078227503206)
    frag_raw = (
        0.23054911576459106 * nitrile
        + 0.20372014781260425 * ketoamide_like
        + 0.12903055761521473 * amide * np.exp(-np.square((exposed_psa - 74.37362488944537) / 27.614846917081305))
        + 0.09025982183072616 * np.minimum(hetero_arom, 1.733517880691858) / 1.7527799838553806
        + 0.056375303881041575 * np.minimum(heteroaryl_fused, 1.8638373660302028)
        + 0.053459196567664624 * ring_flex_balance
        + 0.047678911689016126 * polar_balance
        + 0.031591401490096864 * logd_frag_fit
        + 0.03476459879696593 * amide_nitrile_synergy
        + 0.015883948508179126 * halogen_eff
        - 0.07728310626665032 * aromatic_overload
        - 0.07659934424003004 * flex_over
        - 0.05153657386523234 * acid
        - 0.04768467472498367 * np.clip((exposed_psa - 98.23989078900154) / 33.67420351290785, 0.0077780395611784575, 1.0759746991530978)
        - 0.014731352800261499 * aromatic_drag
    )
    frag_bonus = np.clip(
        0.06781445431642577
        * contact_gate
        * neutral_balance
        * unsat_gate
        * ion_gate
        * frag_raw
        * (0.9303592444950904 + 0.10066778036105872 * mpo),
        -0.021932399531611253,
        0.04901517037047339,
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
        - 0.018036509707494756 * np.maximum(logd - 2.705628346242082, 0.02366255023830288)
    )
    engage_uplift = np.clip(contact_richness - 1.8289544725081444, -0.04014479115011102, 2.828563215698701) / 3.1473015002965905
    lipophilicity_fit = np.exp(-np.square((logd - 2.4096880971039347) / 1.535279905782403))
    size_anti_bias = np.clip((mw - 382.80297441148474) / 132.41702046985603, -0.012035324274906499, 1.1058671799799145)
    predicted = (
        rsae
        + 0.12328223752332142 * cc_c
        + 0.012788874638191904 * hrb
        - 0.0037978921354413163 * tri_penalty
        + 0.0007072126255257129 * engage_uplift * ion_gate
        + frag_bonus
        + 0.0012077874454256576
        * contact_gate
        * ion_gate
        * np.clip(lipophilicity_fit - 0.4017460122534535, -0.21696192027516908, 0.5339565448254127)
        + 0.006483541169570466 * contact_gate * np.clip(polar_balance - 0.4847796299946965, -0.2530818983626887, 0.44892505236155694)
        + 0.005227656609205004
        * contact_gate
        * ion_gate
        * np.clip(neutral_balance - 0.3634364173550874, -0.17508936997176303, 0.4314003291706082)
        + 0.0010096966261512795 * np.clip(hyd_enclosure - 1.0088328948785168, 0.0005493134364064468, 1.2109346843427082)
        + 0.0001747269062799567 * contact_gate * np.clip(ring_flex_balance - 0.36061013372510586, -0.047879329001883124, 0.7895167262148661)
        - 0.0022661591898328363 * np.maximum(ion_delta - 1.1140830311731678, 0.029553853231209563)
        - 0.004324701291635535 * np.maximum(np.abs(exp_tpsa_c) - 1.9229168917845134, -0.013610998805350284)
        - 0.002602456607785064 * np.maximum(aromatic_overload - 0.6157452556939599, -0.00938240331754895)
        - 0.0006458819547469106 * np.maximum(flex_over - 0.36738766080790386, -0.03701071145673132)
        - 0.00019945443629609198 * np.maximum(buried_unsat - 0.6838139519303478, 0.044873326012558154)
        - 0.0009349640114683377 * size_anti_bias
    )
    ### <<< ARF-CHECKPOINT-END
    predicted = np.clip(predicted, 2.5618320893486173, 10.254929396605194)
    ### <<< ARF-CHECKPOINT-END

    return pd.DataFrame(
        {
            "fragalysis_code": compounds["fragalysis_code"].astype(str).to_numpy(),
            "predicted_affinity": np.asarray(predicted, dtype=float),
        }
    )
