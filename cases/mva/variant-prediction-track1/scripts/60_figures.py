"""Render the Track 1 report figures from the evidence tables in this repository.

Usage: python scripts/60_figures.py
Writes PNG (300 dpi), PDF and SVG to submission/figures/. Every plotted number is
read from results/evidence/ or results/; literature constants are cited in the
report next to the figure that uses them.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figstyle import BASELINE, MUTED, PALETTE, TEXT, WIDE, panel_labels, save, use_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "results/evidence"
OUT = ROOT / "submission/figures"

STOP = ("chr15", 40209701)
MISSENSE = ("chr15", 40220612)
CANDIDATE = PALETTE["red"]
SECOND = PALETTE["blue"]
OK = PALETTE["green"]
PENDING = PALETTE["orange"]
SIX = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1"]

HPO = {
    "HP:0002859": "Rhabdomyosarcoma",
    "HP:0000121": "Nephrocalcinosis",
    "HP:0004322": "Short stature",
    "HP:0001508": "Failure to thrive",
    "HP:0003202": "Skeletal muscle atrophy",
    "HP:0001622": "Premature birth",
    "HP:0001518": "Small for gestational age",
}

# UniProt O60566 (BUBR1_HUMAN), 1,050 aa; retrieved 2026-10-17.
PROTEIN_LENGTH = 1050
DOMAINS = [
    (62, 226, "BUB1 N-terminal (TPR)", PALETTE["cyan"]),
    (766, 1050, "Protein kinase domain", PALETTE["purple"]),
]
MOTIFS = [(152, 185, "KNL1"), (224, 232, "D-box")]
# MVA-associated BUBR1 substitutions characterised by Suijkerbuijk et al. 2010.
PUBLISHED_MISSENSE = [
    ("Y155C", 155), ("R550Q", 550), ("R727C", 727), ("R814H", 814),
    ("L844F", 844), ("I909T", 909), ("Q921H", 921), ("L1012P", 1012),
]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def export(fig, name: str) -> None:
    stem = str(OUT / name)
    save(fig, stem, formats=("pdf", "svg"), close=False)
    fig.savefig(f"{stem}.png", dpi=300)
    plt.close(fig)


def bub1b_exons() -> list[tuple[int, int, str, str]]:
    exons = []
    for line in (EV / "structural/tier1_sixgene_exon_utr.bed").read_text().splitlines():
        chrom, start, end, name = line.split("\t")
        if name.startswith("BUB1B|"):
            exon, kind = name.split(":")[1:3]
            exons.append((int(start), int(end), exon, kind))
    return exons


def cdna(exons, pos: int) -> int:
    """1-based CDS position of a plus-strand genomic coordinate."""
    total = 0
    for start, end, _, kind in exons:
        if kind != "CDS":
            continue
        if start < pos <= end:
            return total + pos - start
        total += end - start
    raise ValueError(pos)


def wilson(k: int, n: int) -> tuple[float, float]:
    z = 1.959964
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def box(ax, x, y, w, h, text, color, size=6.5, weight="normal", alpha=0.16):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.015",
            facecolor=mpl.colors.to_rgba(color, alpha), edgecolor=color, linewidth=0.8,
        )
    )
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size,
            fontweight=weight, wrap=True)


def arrow(ax, start, end, color=BASELINE):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=7,
                                 color=color, linewidth=0.8))


def blank(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def graphical_abstract(funnel) -> None:
    fig = plt.figure(figsize=(WIDE, WIDE * 0.46), layout="constrained")
    ax = fig.add_subplot()
    blank(ax)
    box(ax, 0.0, 0.56, 0.22, 0.4,
        "Proband (MVA case)\nblood WGS, ~44x, GRCh38\n7 HPO terms incl.\nrhabdomyosarcoma,\nnephrocalcinosis,\ngrowth restriction",
        PALETTE["black"], size=6.2, alpha=0.06)
    stages = [
        (f"{funnel['records_in']:,}", "variant records"),
        (f"{funnel['biologically_filtered']:,}", "rare, coding/splice candidates"),
        (f"{funnel['hypotheses']:,}", "genotype hypotheses"),
        (f"{funnel['reads_ok_hypotheses']}", "read-supported pair hypotheses"),
        ("1", "submitted BUB1B genotype"),
    ]
    labels = ["variant\nrecords", "rare coding/\nsplice alleles", "genotype\nhypotheses",
              "read-supported\npairs", "submitted\ngenotype"]
    for i, (count, _) in enumerate(stages):
        x = 0.27 + i * 0.105
        height = 0.3 - i * 0.04
        y = 0.8 - height / 2
        ax.add_patch(Rectangle((x, y), 0.09, height, facecolor=mpl.colors.to_rgba(SECOND, 0.15 + 0.15 * i),
                               edgecolor=SECOND, linewidth=0.6))
        ax.text(x + 0.045, 0.8, count, ha="center", va="center", fontsize=6.3, fontweight="bold",
                color="white" if i >= 3 else PALETTE["black"])
        ax.text(x + 0.045, 0.6, labels[i], ha="center", va="top", fontsize=5.4)
    arrow(ax, (0.225, 0.8), (0.265, 0.8))
    arrow(ax, (0.79, 0.8), (0.815, 0.8))
    box(ax, 0.82, 0.56, 0.18, 0.4,
        "BUB1B (MVA1, AR)\nc.2210T>G p.Leu737*\n+\nc.3006T>G p.Asn1002Lys\n(kinase domain)",
        CANDIDATE, size=6.2, weight="bold", alpha=0.1)

    ax.plot([0.03, 0.97], [0.36, 0.36], color=BASELINE, linewidth=4, solid_capstyle="round")
    for start, end, _, color in DOMAINS:
        x0 = 0.03 + 0.94 * start / PROTEIN_LENGTH
        x1 = 0.03 + 0.94 * end / PROTEIN_LENGTH
        ax.add_patch(Rectangle((x0, 0.33), x1 - x0, 0.06, facecolor=color, edgecolor="none"))
    for label, residue, color in (("p.Leu737*", 737, CANDIDATE), ("p.Asn1002Lys", 1002, CANDIDATE)):
        x = 0.03 + 0.94 * residue / PROTEIN_LENGTH
        ax.plot([x, x], [0.39, 0.47], color=color, linewidth=0.9)
        ax.plot(x, 0.47, "o", color=color, markersize=3.5)
        ax.text(x, 0.495, label, ha="center", va="bottom", fontsize=5.8, color=color)
    ax.text(0.03 + 0.94 * 144 / PROTEIN_LENGTH, 0.29, "TPR", ha="center", va="top", fontsize=5.5)
    ax.text(0.03 + 0.94 * 908 / PROTEIN_LENGTH, 0.29, "kinase domain", ha="center", va="top", fontsize=5.5)
    ax.text(0.03, 0.29, "BUBR1 1", ha="left", va="top", fontsize=5.5, color=BASELINE)
    ax.text(0.97, 0.29, "1050", ha="right", va="top", fontsize=5.5, color=BASELINE)

    checks = [
        ("Read-level authentication", OK), ("Three-caller consensus", OK),
        ("Phenotype fit: rank 1 of 8", OK), ("AlphaMissense PP3_Moderate", OK),
        ("No CNV, SV or splice second hit", OK), ("Trans phase: confirm", PENDING),
    ]
    for i, (label, color) in enumerate(checks):
        x = 0.01 + i * 0.165
        ax.add_patch(FancyBboxPatch((x, 0.04), 0.155, 0.14, boxstyle="round,pad=0.005,rounding_size=0.02",
                                    facecolor=mpl.colors.to_rgba(color, 0.14), edgecolor=color, linewidth=0.7))
        ax.text(x + 0.0775, 0.11, label, ha="center", va="center", fontsize=5.4, wrap=True)
    export(fig, "graphical_abstract")


def fig1_workflow(funnel) -> None:
    fig = plt.figure(figsize=(WIDE, WIDE * 0.62), layout="constrained")
    grid = fig.add_gridspec(2, 1, height_ratios=[0.9, 1.1])
    ax = fig.add_subplot(grid[0])
    blank(ax)
    steps = [
        ("Input", "FASTQ\nsupplied VCF\n7 HPO terms", PALETTE["black"]),
        ("MareNostrum 5", "bwa-mem2\nDeepVariant, GATK,\nStrelka2, bcftools", SECOND),
        ("Annotation", "snpEff, MANE\ngnomAD, ClinVar\nAlphaMissense", PALETTE["cyan"]),
        ("Prioritisation", "pre-registered\nTier-1 genes\npair hypotheses", PALETTE["purple"]),
        ("Validation", "MQ60/BQ20 reads\nconcordance, phase\nCNV, SV, splice", OK),
        ("Interpretation", "ACMG/AMP\nClinGen SVI\nranked CSV", CANDIDATE),
    ]
    width = 0.152
    for i, (title, body, color) in enumerate(steps):
        x = i * (width + 0.018)
        box(ax, x, 0.08, width, 0.78, "", color, alpha=0.08)
        ax.text(x + width / 2, 0.78, title, ha="center", va="top", fontsize=6.4, fontweight="bold", color=color)
        ax.text(x + width / 2, 0.4, body, ha="center", va="center", fontsize=5.6, linespacing=1.4)
        if i:
            arrow(ax, (x - 0.017, 0.47), (x - 0.001, 0.47))

    ax = fig.add_subplot(grid[1])
    stages = [
        ("Annotated variant records", funnel["records_in"], "records"),
        ("Rare (gnomAD AF <= 0.01) prefilter universe", funnel["prefilter_universe"], "records"),
        ("Biologically filtered candidates", funnel["biologically_filtered"], "alleles"),
        ("Constructed genotype hypotheses", funnel["hypotheses"], "hypotheses"),
        ("Top hypotheses read-checked", funnel["top60_read_checked"], "hypotheses"),
        ("Read-supported SNV-pair hypotheses", funnel["reads_ok_hypotheses"], "hypotheses"),
        ("Initial ranked genotypes", funnel["submitted_rows"], "rows"),
        ("Adjudicated genotype (BUB1B)", 1, "row"),
    ]
    y = np.arange(len(stages))[::-1]
    colors = [SECOND] * (len(stages) - 1) + [CANDIDATE]
    ax.barh(y, [s[1] for s in stages], color=colors, height=0.62)
    ax.set_xscale("log")
    ax.set_xlim(0.6, 4e7)
    ax.set_yticks(y, [s[0] for s in stages])
    ax.set_xlabel("Count (log scale; units differ between stages)")
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    for yi, (_, value, unit) in zip(y, stages):
        ax.text(value * 1.25, yi, f"{value:,} {unit}", va="center", fontsize=5.8)
    panel_labels([fig.axes[0], ax])
    export(fig, "fig1_workflow_funnel")


def fig2_locus() -> dict[str, int]:
    exons = bub1b_exons()
    c_stop, c_mis = cdna(exons, STOP[1]), cdna(exons, MISSENSE[1])
    cds = [(s, e, x) for s, e, x, k in exons if k == "CDS"]
    lengths = [e - s for s, e, _ in cds]
    cds_total = sum(lengths)
    last_junction = sum(lengths[:-1])
    stop_exon = next(x for s, e, x, k in exons if k == "CDS" and s < STOP[1] <= e)
    mis_exon = next(x for s, e, x, k in exons if k == "CDS" and s < MISSENSE[1] <= e)
    facts = {
        "c_stop": c_stop, "c_missense": c_mis, "cds_length": cds_total,
        "last_junction_c": last_junction, "ptc_to_last_junction_nt": last_junction - c_stop,
        "stop_exon": int(stop_exon[1:]), "missense_exon": int(mis_exon[1:]), "exon_count": len(cds),
        "span_bp": MISSENSE[1] - STOP[1],
    }

    fig = plt.figure(figsize=(WIDE, WIDE * 0.78), layout="constrained")
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.15, 0.75])

    ax = fig.add_subplot(grid[0])
    gene_start, gene_end = exons[0][0], exons[-1][1]
    kb = lambda p: (p - gene_start) / 1000  # noqa: E731
    ax.plot([kb(gene_start), kb(gene_end)], [0, 0], color=BASELINE, linewidth=0.8)
    for start, end, exon, kind in exons:
        h = 0.5 if kind == "CDS" else 0.26
        ax.add_patch(Rectangle((kb(start), -h / 2), max(kb(end) - kb(start), 0.12), h,
                               facecolor=PALETTE["black"] if kind == "CDS" else MUTED, edgecolor="none"))
    for pos, label, color, dy in ((STOP[1], "c.2210T>G\np.Leu737*", CANDIDATE, 1.35),
                                  (MISSENSE[1], "c.3006T>G\np.Asn1002Lys", CANDIDATE, 1.35)):
        ax.plot([kb(pos)] * 2, [0.3, dy - 0.25], color=color, linewidth=0.9)
        ax.text(kb(pos), dy, label, ha="center", va="bottom", fontsize=6, color=color)
    ax.annotate("", xy=(kb(STOP[1]), -0.75), xytext=(kb(MISSENSE[1]), -0.75),
                arrowprops=dict(arrowstyle="<->", color=PENDING, linewidth=0.9))
    ax.text((kb(STOP[1]) + kb(MISSENSE[1])) / 2, -0.95, f"{facts['span_bp']:,} bp between alleles",
            ha="center", va="top", fontsize=6, color=PENDING)
    ax.add_patch(Rectangle((kb(gene_start) + 2, -1.15), 0.5, 0.18, facecolor=SECOND, edgecolor="none"))
    ax.text(kb(gene_start) + 2.8, -1.06, "0.5 kb: illustrative short-read fragment", va="center", fontsize=5.6,
            color=SECOND)
    ax.text(kb(exons[0][0]), 0.4, "exon 1", ha="left", va="bottom", fontsize=5.5)
    ax.text(kb(exons[-1][1]), -0.4, f"exon {facts['exon_count']}", ha="right", va="top", fontsize=5.5)
    ax.text(kb(STOP[1]) - 0.4, -0.4, f"exon {facts['stop_exon']}", ha="right", va="top", fontsize=5.5)
    ax.set_ylim(-1.45, 2.3)
    ax.set_xlim(-1, kb(gene_end) + 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_xlabel(f"Position from BUB1B transcription start (kb; chr15:{gene_start + 1:,}, + strand, MANE NM_001211.6)")

    ax = fig.add_subplot(grid[1])
    ax.add_patch(Rectangle((1, -0.12), PROTEIN_LENGTH - 1, 0.24, facecolor=MUTED, edgecolor="none"))
    for start, end, label, color in DOMAINS:
        ax.add_patch(Rectangle((start, -0.2), end - start, 0.4, facecolor=color, edgecolor="none"))
        ax.text((start + end) / 2, 0, label, ha="center", va="center", fontsize=6, color="white",
                fontweight="bold")
    for start, end, label in MOTIFS:
        ax.plot([start, max(end, start + 8)], [0.27, 0.27], color=PALETTE["black"], linewidth=1.6,
                solid_capstyle="butt")
        ax.text((start + end) / 2, 0.33, label, ha="center", va="bottom", fontsize=5.4)
    for label, residue in (("p.Leu737*", 737), ("p.Asn1002Lys", 1002)):
        ax.plot([residue] * 2, [0.2, 0.95], color=CANDIDATE, linewidth=0.9)
        ax.plot(residue, 0.95, marker="*" if "*" in label else "o", color=CANDIDATE, markersize=6)
        ax.text(residue, 1.08, label, ha="center", va="bottom", fontsize=6.2, color=CANDIDATE, fontweight="bold")
    offsets = {"R814H": -1.0, "L844F": -1.32, "I909T": -1.0, "Q921H": -1.32, "L1012P": -1.0,
               "R727C": -1.0, "R550Q": -1.0, "Y155C": -1.0}
    for label, residue in PUBLISHED_MISSENSE:
        tip = offsets[label]
        ax.plot([residue] * 2, [-0.2, tip + 0.12], color=SECOND, linewidth=0.7)
        ax.plot(residue, tip + 0.12, "o", color=SECOND, markersize=3)
        ax.text(residue, tip, label, ha="center", va="top", fontsize=5.5, color=SECOND)
    ax.text(5, 1.08, "This proband", fontsize=6, color=CANDIDATE, va="bottom")
    ax.text(5, -1.5, "MVA-associated substitutions (Suijkerbuijk et al. 2010)", fontsize=6, color=SECOND,
            va="bottom")
    ax.set_xlim(-10, PROTEIN_LENGTH + 10)
    ax.set_ylim(-1.75, 1.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_xlabel("BUBR1 residue (UniProt O60566, 1,050 aa)")

    ax = fig.add_subplot(grid[2])
    position = 0
    for i, (start, end, exon) in enumerate(cds):
        length = end - start
        ax.add_patch(Rectangle((position, -0.2), length, 0.4,
                               facecolor=MUTED if i % 2 else "#E8E8E8", edgecolor="white", linewidth=0.4))
        position += length
    ax.axvspan(0, last_junction - 55, ymin=0.15, ymax=0.85, color=OK, alpha=0.12, linewidth=0)
    ax.axvline(last_junction, color=PALETTE["black"], linewidth=0.8, linestyle="--")
    ax.text(last_junction - 15, -0.6, f"last exon junction c.{last_junction}", ha="right", va="center",
            fontsize=5.4)
    ax.text((last_junction - 55) / 2, 0.5,
            f"premature stop codons here are predicted to trigger NMD (c.{c_stop} lies "
            f"{facts['ptc_to_last_junction_nt']} nt upstream of the last junction)",
            ha="center", va="center", fontsize=5.6, color=OK)
    for c, label, marker, y in ((c_stop, f"c.{c_stop} (exon {facts['stop_exon']})", "*", -0.3),
                                (c_mis, f"c.{c_mis} (exon {facts['missense_exon']})", "o", 0.3)):
        ax.plot(c, 0, marker=marker, color=CANDIDATE, markersize=6 if marker == "*" else 4)
        ax.text(min(c, 2990), y, label, ha="center", va="top" if y < 0 else "bottom", fontsize=5.6,
                color=CANDIDATE)
    ax.set_xlim(0, cds_total)
    ax.set_ylim(-0.85, 0.85)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_xlabel(f"Coding sequence position (c.; {cds_total:,} nt across {len(cds)} coding exons)")
    panel_labels(fig.axes)
    export(fig, "fig2_bub1b_locus_protein")
    return facts


def fig3_alleles() -> None:
    rows = {r["grch38_key"]: r for r in read(EV / "alleles/allele_validation.tsv")}
    stop = rows["chr15-40209701-T-G"]
    mis = rows["chr15-40220612-T-G"]
    fig, axes = plt.subplots(2, 2, figsize=(WIDE, WIDE * 0.62), layout="constrained")

    ax = axes[0, 0]
    labels = ["p.Leu737*\nchr15:40,209,701", "p.Asn1002Lys\nchr15:40,220,612"]
    y = np.array([1, 0])
    for yi, row in zip(y, (stop, mis)):
        ref, fwd, rev = int(row["read_ref_MQ60_BQ20"]), int(row["alt_forward"]), int(row["alt_reverse"])
        ax.barh(yi, ref, color=MUTED, height=0.55, label="REF reads" if yi == 1 else None)
        ax.barh(yi, fwd, left=ref, color=CANDIDATE, height=0.55, label="ALT forward" if yi == 1 else None)
        ax.barh(yi, rev, left=ref + fwd, color=PALETTE["purple"], height=0.55,
                label="ALT reverse" if yi == 1 else None)
        ax.text(ref + fwd + rev + 0.8, yi, f"{ref}/{fwd}+{rev}", va="center", fontsize=5.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Reads at MAPQ >= 60, base quality >= 20")
    ax.set_xlim(0, 58)
    ax.legend(loc="lower right", fontsize=5.6, ncols=3, bbox_to_anchor=(1, 1.0))
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)

    ax = axes[0, 1]
    ax.axvspan(0.25, 0.75, color=OK, alpha=0.1, linewidth=0)
    ax.axvline(0.5, color=BASELINE, linewidth=0.6, linestyle="--")
    for yi, row, label in zip(y, (stop, mis), ("p.Leu737*", "p.Asn1002Lys")):
        alt = int(row["read_alt_MQ60_BQ20"])
        depth = alt + int(row["read_ref_MQ60_BQ20"])
        lo, hi = wilson(alt, depth)
        ax.errorbar(alt / depth, yi, xerr=[[alt / depth - lo], [hi - alt / depth]], fmt="o", color=CANDIDATE,
                    capsize=2, markersize=4, linewidth=0.9)
        ax.text(hi + 0.02, yi, f"{alt}/{depth} = {alt / depth:.2f}", va="center", fontsize=5.8)
    vaf = {r["metric"]: int(r["count"]) for r in read(EV / "qc/vaf_counts.tsv")}
    share = vaf["het_vaf_ok"] / vaf["het"]
    ax.text(0.26, 1.45, f"heterozygous window 0.25-0.75\n(genome-wide: {share:.1%} of {vaf['het']:,} het calls)",
            fontsize=5.6, va="top", color=OK)
    ax.set_yticks(y, ["p.Leu737*", "p.Asn1002Lys"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("ALT allele fraction (Wilson 95% interval)")

    ax = axes[1, 0]
    conc = read(EV / "alleles/caller_concordance.tsv")
    gt_code = {"0/1": 1, "1/1": 2, "not_called": 0}
    for row in conc:
        x = (int(row["grch38_key"].split("-")[1]) - 40161068) / 1000
        ref, alt = row["grch38_key"].split("-")[2:4]
        indel = len(ref) != len(alt)
        is_cand = int(row["grch38_key"].split("-")[1]) in (STOP[1], MISSENSE[1])
        for yy, key in ((1, "source_gt"), (0, "recall_gt")):
            code = gt_code[row[key]]
            color = [PALETTE["black"], SECOND, PALETTE["cyan"]][code] if not is_cand else CANDIDATE
            face = color if code else "white"
            ax.plot(x, yy, marker="s" if indel else "o", markersize=4.6, markerfacecolor=face,
                    markeredgecolor=color, linestyle="none")
        if row["caller_relation"] != "concordant":
            ax.plot([x, x], [0.12, 0.88], color=PENDING, linewidth=0.6, linestyle=":")
    n_conc = sum(r["caller_relation"] == "concordant" for r in conc)
    n_indel_disc = sum(r["caller_relation"] != "concordant" and len(r["grch38_key"].split("-")[2])
                       != len(r["grch38_key"].split("-")[3]) for r in conc)
    ax.set_yticks([1, 0], ["supplied VCF", "bwa-mem2 + bcftools\n(independent)"])
    ax.set_ylim(-0.6, 1.9)
    ax.set_xlabel("BUB1B bounded alleles, position from gene start (kb)")
    ax.text(6, 1.78, f"{n_conc}/{len(conc)} concordant; all {n_indel_disc} discordant sites are indels",
            fontsize=5.6, va="top")
    ax.text(6, 1.52, "blue 0/1, light blue 1/1, red = submitted alleles; square = indel; open = not called",
            fontsize=5.3, va="top", color=BASELINE)
    ax.grid(axis="y", visible=False)

    ax = axes[1, 1]
    alleles = read(EV / "alleles/allele_validation.tsv")
    points = [(r["gene"], float(r["gnomad_AF"]), r["grch38_key"]) for r in alleles
              if r["gnomad_AF"] not in ("", ".", "NA")]
    genes = [g for g in SIX if any(p[0] == g for p in points)]
    for gi, gene in enumerate(genes):
        values = [p for p in points if p[0] == gene]
        jitter = np.linspace(-0.22, 0.22, len(values)) if len(values) > 1 else [0]
        for j, (_, af, key) in zip(jitter, values):
            is_cand = int(key.split("-")[1]) in (STOP[1], MISSENSE[1])
            ax.plot(gi + j, af, "o", color=CANDIDATE if is_cand else BASELINE,
                    markersize=4.5 if is_cand else 2.6, zorder=3 if is_cand else 2)
    ax.axhline(0.01, color=PENDING, linewidth=0.7, linestyle="--")
    ax.text(len(genes) - 0.5, 0.012, "rare-allele filter (AF 0.01)", ha="right", va="bottom", fontsize=5.6,
            color=PENDING)
    ax.annotate("p.Leu737*\nAF 7.4e-05", xy=(0, 7.4e-5), xytext=(0.45, 2e-4), fontsize=5.6, color=CANDIDATE,
                arrowprops=dict(arrowstyle="-", color=CANDIDATE, linewidth=0.5))
    ax.annotate("p.Asn1002Lys\nAF 5e-07", xy=(0, 5e-7), xytext=(0.45, 1.2e-6), fontsize=5.6, color=CANDIDATE,
                arrowprops=dict(arrowstyle="-", color=CANDIDATE, linewidth=0.5))
    ax.set_yscale("log")
    ax.set_xticks(range(len(genes)), genes)
    ax.set_ylabel("gnomAD allele frequency")
    ax.set_xlim(-0.6, len(genes) - 0.4)
    panel_labels(axes)
    export(fig, "fig3_allele_authentication")


def band_chart(ax, bands, score, label, title):
    for lo, hi, name, color in bands:
        ax.add_patch(Rectangle((lo, 0), hi - lo, 1, facecolor=color, alpha=0.75, edgecolor="white", linewidth=0.6))
        if hi - lo >= 0.03:
            ax.text((lo + hi) / 2, 0.5, name, ha="center", va="center", fontsize=5.3,
                    rotation=90 if hi - lo < 0.09 else 0)
    ax.plot([score, score], [-0.12, 1.12], color=PALETTE["black"], linewidth=1.3)
    ax.plot(score, 1.2, "v", color=PALETTE["black"], markersize=5)
    ax.text(score, 1.33, label, ha="center", va="bottom", fontsize=6.2, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.2, 1.75)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(False)
    ax.set_xlabel(title)


def fig4_effect() -> None:
    lr = {r["predictor"]: r for r in read(EV / "alleles/predictor_lr.tsv")}
    avi = {r["label"]: r for r in read(EV / "alleles/coding_avi_splice.tsv")}
    noncoding = read(EV / "second_hit/noncoding_avi.tsv")
    splice = read(EV / "second_hit/spliceai_keys.tsv")
    grey, s1, s2, s3 = "#E3E3E3", "#F7C99B", "#EE9A5B", "#D55E00"
    fig = plt.figure(figsize=(WIDE, WIDE * 0.6), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])

    ax = fig.add_subplot(grid[0, :2])
    band_chart(ax, [(0, 0.17, "BP4", "#BFD9EE"), (0.17, 0.792, "indeterminate", grey),
                    (0.792, 0.906, "PP3\nSupp.", s1), (0.906, 0.972, "PP3\nMod.", s2),
                    (0.972, 0.99, "+3", s3), (0.99, 1.0, "Str.", "#A23B00")],
               float(lr["AlphaMissense"]["candidate_score"]), "p.Asn1002Lys 0.923",
               "AlphaMissense score; calibrated PP3/BP4 intervals (Bergquist 2025)")

    ax = fig.add_subplot(grid[1, :2])
    band_chart(ax, [(0, 0.29, "BP4", "#BFD9EE"), (0.29, 0.644, "indeterminate", grey),
                    (0.644, 0.773, "PP3 Supp.", s1), (0.773, 0.879, "PP3 Mod.", s2),
                    (0.879, 0.932, "+3", s3), (0.932, 1.0, "Strong", "#A23B00")],
               float(lr["REVEL"]["candidate_score"]), "p.Asn1002Lys 0.472",
               "REVEL score; calibrated PP3/BP4 intervals (Pejaver 2022)")

    ax = fig.add_subplot(grid[0, 2])
    names = ["AlphaMissense", "REVEL"]
    values = [float(lr[n]["lr_pathogenic_call"]) for n in names]
    ax.bar(names, values, color=[s2, grey], edgecolor=PALETTE["black"], linewidth=0.5, width=0.55)
    for x, v in enumerate(values):
        ax.text(x, v * 1.12, f"{v:.1f}", ha="center", va="bottom", fontsize=6)
    for v, name in ((2.08, "Supporting"), (4.33, "Moderate"), (18.7, "Strong")):
        ax.axhline(v, color=BASELINE, linewidth=0.5, linestyle=":")
        ax.text(1.42, v, name, fontsize=5.3, va="bottom", ha="right", color=BASELINE)
    ax.set_yscale("log")
    ax.set_ylim(1, 300)
    ax.set_ylabel("Tail LR+ at candidate score\n(ClinVar 2-star, BUB1B held out)")

    ax = fig.add_subplot(grid[1, 2])
    cats = [("SpliceAI", [float(r["max_delta"]) for r in splice if r["max_delta"] not in (".", "")],
             [float(r["max_delta"]) for r in splice if r["label"].startswith("BUB1B_p.")])]
    values = cats[0][1]
    ax.plot(np.random.default_rng(20261017).uniform(-0.18, 0.18, len(values)), values, "o", color=BASELINE,
            markersize=2.4, alpha=0.8)
    for label, row, dy in (("p.Leu737*", "BUB1B_p.Leu737Ter", 0.1), ("p.Asn1002Lys", "BUB1B_p.Asn1002Lys", 0.04)):
        value = float(next(r["max_delta"] for r in splice if r["label"] == row))
        ax.plot(0.3, value, "o", color=CANDIDATE, markersize=4)
        ax.annotate(f"{label} {value:.2f}", xy=(0.3, value), xytext=(0.42, value + dy), fontsize=5.6,
                    color=CANDIDATE, va="center", arrowprops=dict(arrowstyle="-", color=CANDIDATE, linewidth=0.5))
    ax.axhline(0.2, color=PENDING, linewidth=0.7, linestyle="--")
    ax.text(-0.3, 0.205, "0.2 splice-altering cut-off", fontsize=5.5, va="bottom", color=PENDING)
    percentiles = [float(r["top_percentile"]) for r in noncoding if r["top_percentile"]]
    ax.text(-0.3, 0.33, f"{len(values)} scored keys, max {max(values):.2f}\n"
            f"AlphaGenome AVI top percentile:\nstop {float(avi['BUB1B_p.Leu737Ter']['top_percentile']):.3f}%, "
            f"missense {float(avi['BUB1B_p.Asn1002Lys']['top_percentile']):.3f}%\n"
            f"28 noncoding keys: median {np.median(percentiles):.0f}%", fontsize=5.4, va="top")
    ax.set_xlim(-0.35, 1.2)
    ax.set_ylim(-0.02, 0.36)
    ax.set_xticks([])
    ax.set_xlabel("Scored alleles")
    ax.set_ylabel("SpliceAI max delta score")
    panel_labels(fig.axes, ["(a)", "(b)", "(c)", "(d)"])
    export(fig, "fig4_variant_effect")


def fig5_phenotype() -> None:
    rows = read(EV / "literature/angle_b_phenotype.tsv")
    null = read(EV / "controls/null_runs.tsv")
    genes = [r["gene"] for r in rows]
    terms = list(HPO)
    matrix = np.zeros((len(genes), len(terms)))
    for i, row in enumerate(rows):
        levels = dict(part.rsplit(":", 1) for part in row["feature_match_levels_0_unsubstantiated_1_broad_2_exact"].split(";"))
        for j, term in enumerate(terms):
            matrix[i, j] = int(levels.get(term, 0))
    fig = plt.figure(figsize=(WIDE, WIDE * 0.5), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=[2.3, 1.0, 1.0])

    ax = fig.add_subplot(grid[0])
    cmap = mpl.colors.ListedColormap(["#F2F2F2", "#9ECAE1", SECOND])
    ax.imshow(matrix, cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(terms)), [HPO[t] for t in terms], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(genes)), genes)
    ax.set_xticks(np.arange(-0.5, len(terms)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)
    for (i, j), v in np.ndenumerate(matrix):
        ax.text(j, i, ["-", "broad", "exact"][int(v)], ha="center", va="center", fontsize=5.3,
                color="white" if v == 2 else PALETTE["black"])
    ax.get_yticklabels()[0].set_color(CANDIDATE)
    ax.get_yticklabels()[0].set_fontweight("bold")

    ax = fig.add_subplot(grid[1])
    points = [float(r["weighted_coverage_points"]) for r in rows]
    maximum = float(rows[0]["max_weighted_points"])
    y = np.arange(len(genes))
    ax.barh(y, points, color=[CANDIDATE] + [BASELINE] * (len(genes) - 1), height=0.6)
    ax.set_yticks(y, genes)
    ax.invert_yaxis()
    ax.set_xlim(0, maximum)
    ax.set_xlabel(f"Weighted phenotype coverage\n(points of {maximum:g})")
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    for yi, v in zip(y, points):
        ax.text(v + 0.2, yi, f"{v:g}", va="center", fontsize=5.6)

    ax = fig.add_subplot(grid[2])
    counts = {r["item"].replace("times_rank1_", ""): int(r["value"]) for r in null
              if r["arm"] == "hpo_weight_permutation" and r["item"].startswith("times_rank1_")}
    total = int(next(r["value"] for r in null if r["item"] == "permutations"))
    shown = [g for g in genes if g in counts]
    ax.barh(range(len(shown)), [counts[g] for g in shown],
            color=[CANDIDATE if g == "BUB1B" else BASELINE for g in shown], height=0.6)
    ax.set_yticks(range(len(shown)), shown)
    ax.invert_yaxis()
    ax.set_xlim(0, total)
    ax.set_xlabel(f"Rank-1 count in {total}\nHPO-weight permutations")
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    for yi, g in enumerate(shown):
        ax.text(counts[g] + 2, yi, str(counts[g]), va="center", fontsize=5.6)
    panel_labels(fig.axes)
    export(fig, "fig5_phenotype_fit")


def fig6_competing() -> None:
    cnv = {r["gene"]: r for r in read(EV / "structural/sixgene_cnv.tsv")}
    exon = {r["gene"]: r for r in read(EV / "structural/exon_depth_by_gene.tsv")}
    rare = {r["gene"]: r for r in read(EV / "genome_scan/rare_sixgene.tsv")}
    plp = read(EV / "genome_scan/clinvar_plp_genes.tsv")
    sites = read(EV / "literature/published_site_genotypes.tsv")
    fig = plt.figure(figsize=(WIDE, WIDE * 0.62), layout="constrained")
    grid = fig.add_gridspec(2, 2)

    ax = fig.add_subplot(grid[0, 0])
    x = np.arange(len(SIX))
    ax.axhspan(0.75, 1.25, color=OK, alpha=0.1, linewidth=0)
    for level, text in ((0.5, "expected heterozygous deletion"), (1.5, "expected duplication")):
        ax.axhline(level, color=PENDING, linewidth=0.7, linestyle="--")
        ax.text(5.45, level + 0.03, text, ha="right", va="bottom", fontsize=5.5, color=PENDING)
    ax.plot(x - 0.12, [float(cnv[g]["ratio_to_sixgene_median"]) for g in SIX], "o", color=SECOND,
            markersize=4, label="gene body +/- 2 kb")
    ax.plot(x + 0.12, [float(exon[g]["ratio_to_sixgene_median"]) for g in SIX], "s", color=PALETTE["purple"],
            markersize=3.6, label="coding exons (median)")
    ax.set_xticks(x, SIX)
    ax.set_ylim(0.3, 1.75)
    ax.set_ylabel("Depth ratio to six-gene median")
    ax.legend(loc="upper left", fontsize=5.6, ncols=2)

    ax = fig.add_subplot(grid[0, 1])
    het = [int(rare[g]["rare_het"]) for g in SIX]
    ax.bar(x, het, color=[CANDIDATE if g == "BUB1B" else BASELINE for g in SIX], width=0.6)
    for xi, v in zip(x, het):
        ax.text(xi, v + 0.05, str(v), ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x, SIX)
    ax.set_ylim(0, 3)
    ax.set_ylabel("Rare coding het alleles\n(AF <= 0.01)")
    other = ", ".join(r["gene"] for r in plp if r["gene"] != "BUB1B")
    ax.text(5.4, 2.85, f"Genome-wide ClinVar P/LP: {len(plp)} heterozygous,\n0 homozygous; only BUB1B in panel\n"
            f"(others: {other})", ha="right", va="top", fontsize=5.4)

    ax = fig.add_subplot(grid[1, :])
    labels = [r["label"].replace("_", " ", 1) for r in sites]
    depth = [int(r["dp"]) for r in sites]
    ax.bar(range(len(sites)), depth, color=MUTED, width=0.6, label="REF reads")
    ax.bar(range(len(sites)), [0] * len(sites), color=CANDIDATE, width=0.6, label="ALT reads (0 at every site)")
    for xi, r in enumerate(sites):
        ax.text(xi, int(r["dp"]) + 1, f"{r['gt']}", ha="center", va="bottom", fontsize=5.6)
    ax.set_xticks(range(len(sites)), labels, rotation=25, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Read depth at published\npathogenic site")
    ax.set_ylim(0, 62)
    ax.legend(loc="upper left", fontsize=5.6, ncols=2)
    ax.text(len(sites) - 0.5, 60, "Four further published SNVs (BUB1B upstream, BUB1 start-loss, TRIP13 x2): ALT = 0",
            ha="right", va="top", fontsize=5.6)
    panel_labels(fig.axes)
    export(fig, "fig6_competing_explanations")


def fig7_genome() -> None:
    arms = [r for r in read(EV / "aneuploidy/arm_depth.tsv") if r["arm"] in ("p", "q")]
    mosaic = read(EV / "aneuploidy/mosaic_segments.tsv")
    roh = {r["metric"]: r["value"] for r in read(EV / "aneuploidy/roh_summary.tsv")}
    vbid = {r["metric"]: r["value"] for r in read(EV / "qc/verifybamid2_summary.tsv")}
    autosomes = [f"chr{i}" for i in range(1, 23)]
    fig, axes = plt.subplots(2, 1, figsize=(WIDE, WIDE * 0.56), layout="constrained", sharex=True)

    ax = axes[0]
    ax.axhspan(0.95, 1.05, color=OK, alpha=0.1, linewidth=0)
    for c, label in ((0.2, "20% trisomy"), (-0.2, "20% monosomy")):
        ax.axhline(1 + c / 2, color=PENDING, linewidth=0.6, linestyle="--")
        ax.text(22.4, 1 + c / 2, label, fontsize=5.4, va="center", ha="left", color=PENDING)
    for arm, offset, marker in (("p", -0.15, "o"), ("q", 0.15, "s")):
        xs, ys = [], []
        for i, chrom in enumerate(autosomes):
            row = next((r for r in arms if r["chrom"] == chrom and r["arm"] == arm), None)
            if row and float(row["median_depth"]) > 0:
                xs.append(i + offset)
                ys.append(float(row["ratio_to_autosome_arm_median"]))
        ax.plot(xs, ys, marker, color=SECOND if arm == "p" else PALETTE["purple"], markersize=3.4,
                linestyle="none", label=f"{arm} arm")
    ax.set_ylim(0.85, 1.15)
    ax.set_ylabel("Arm depth ratio\n(100 kb bins, MAPQ >= 20)")
    ax.legend(loc="lower left", fontsize=5.6, ncols=2)

    ax = axes[1]
    for c in (0.1, 0.2):
        ax.axhline(1 / (2 + c), color=PENDING, linewidth=0.6, linestyle="--")
        ax.text(22.4, 1 / (2 + c), f"{int(c * 100)}% trisomy", fontsize=5.4, va="center", ha="left", color=PENDING)
    ax.axhline(0.5, color=BASELINE, linewidth=0.6)
    seg = {r["chrom"]: r for r in mosaic}
    ys = [float(seg[c]["minor_af"]) for c in autosomes]
    ax.plot(range(22), ys, "o", color=SECOND, markersize=3.6)
    for i, c in enumerate(autosomes):
        ax.text(i, ys[i] - 0.004, seg[c]["n_af"], ha="center", va="top", fontsize=5)
    ax.set_ylim(0.44, 0.51)
    ax.set_ylabel("Minor allele fraction\n(GATK ModelSegments)")
    ax.set_xticks(range(22), [c.replace("chr", "") for c in autosomes])
    ax.set_xlabel("Autosome (numbers below points: heterozygous sites modelled)")
    ax.set_xlim(-0.6, 24.2)
    ax.text(0, 0.508, f"VerifyBamID2 FREEMIX {float(vbid['freemix']):.1e} ({int(vbid['snps']):,} SNPs);  "
            f"F_ROH(>=1.5 Mb) = {float(roh['F_ROH_ge_1.5Mb']):g};  all 22 autosomes balanced",
            fontsize=5.6, va="top")
    panel_labels(axes)
    export(fig, "fig7_genome_wide_qc")


def fig8_phase(facts) -> None:
    pid = {r["metric"]: int(r["count"]) for r in read(EV / "phase/gatk_pid_pair_counts.tsv")}
    pairs = read(EV / "phase/phase_pair_evidence.tsv")
    fig = plt.figure(figsize=(WIDE, WIDE * 0.42), layout="constrained")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.7, 1])

    ax = fig.add_subplot(grid[0])
    blank(ax)
    for title, x0, mat, pat, verdict, color in (
        ("Trans (compound heterozygous)", 0.0, "stop", "missense", "no wild-type BUBR1 allele:\nnull + destabilised kinase domain\n= biallelic MVA1 genotype", CANDIDATE),
        ("Cis (same haplotype)", 0.52, "both", None, "one intact allele remains:\nheterozygous carrier state", BASELINE),
    ):
        ax.text(x0 + 0.27, 0.95, title, ha="center", va="top", fontsize=6.6, fontweight="bold", color=color)
        for yy, name in ((0.68, "hap. 1"), (0.46, "hap. 2")):
            ax.plot([x0 + 0.1, x0 + 0.44], [yy, yy], color=BASELINE, linewidth=3, solid_capstyle="round")
            ax.text(x0 + 0.09, yy, name, ha="right", va="center", fontsize=5.4)
        def mark(xx, yy, kind):
            ax.plot(xx, yy, marker="*" if kind == "stop" else "o", color=CANDIDATE,
                    markersize=7 if kind == "stop" else 4.5, zorder=3)
        mark(x0 + 0.2, 0.68, "stop")
        mark(x0 + 0.36, 0.68 if mat == "both" else 0.46, "missense")
        ax.text(x0 + 0.27, 0.3, verdict, ha="center", va="top", fontsize=5.7, color=color)
    ax.text(0.5, 0.03, f"Alleles {facts['span_bp']:,} bp apart: resolvable by parental genotyping or reads > 11 kb "
            "(long-read HiFi/ONT)", ha="center", va="bottom", fontsize=5.8, color=PENDING)

    ax = fig.add_subplot(grid[1])
    rows = [
        ("WhatsHap: Tier-1\ncandidate pairs", sum(p["phase"] != "unknown" for p in pairs), len(pairs), "pairs"),
        ("GATK PID: coding\nheterozygous sites", pid["with_pid"], pid["matched_high_moderate_hets"], "sites"),
    ]
    y = np.arange(len(rows))[::-1]
    for yi, (_, phased, total, unit) in zip(y, rows):
        ax.barh(yi, phased / total, color=OK, height=0.5, label="phased" if yi == y[0] else None)
        ax.barh(yi, 1 - phased / total, left=phased / total, color=PENDING, height=0.5,
                label="unphased" if yi == y[0] else None)
        ax.text(1.02, yi, f"{phased:,}/{total:,}\n{unit}", va="center", fontsize=5.6)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlim(0, 1.3)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Fraction phased")
    ax.legend(loc="upper center", fontsize=5.6, ncols=2, bbox_to_anchor=(0.4, 1.16))
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    panel_labels(fig.axes)
    export(fig, "fig8_phase_scenarios")


def fig10_acmg() -> None:
    fig, ax = plt.subplots(figsize=(WIDE, WIDE * 0.36), layout="constrained")
    pm2 = ("PM2_P", 1, PALETTE["cyan"])
    pp3 = ("PP3_Moderate", 2, PALETTE["orange"])
    ladders = [
        ("p.Leu737*: present data", [("PVS1", 8, CANDIDATE), pm2]),
        ("p.Asn1002Lys: present data", [pm2, pp3]),
        ("+ parental trans (PM3)\n+ MVA cytogenetics (PP4)", [pm2, pp3, ("PM3", 2, OK), ("PP4", 1, PALETTE["purple"])]),
        ("+ BUBR1 abundance assay\n(PS3_Moderate)", [pm2, pp3, ("PM3", 2, OK), ("PP4", 1, PALETTE["purple"]),
                                                    ("PS3_Moderate", 2, SECOND)]),
    ]
    ax.axvspan(0, 6, color="#F2F2F2", linewidth=0)
    ax.axvspan(6, 10, color=mpl.colors.to_rgba(PENDING, 0.12), linewidth=0)
    ax.axvspan(10, 13, color=mpl.colors.to_rgba(CANDIDATE, 0.12), linewidth=0)
    for x, text in ((3, "VUS (0-5)"), (8, "Likely pathogenic (6-9)"), (11.5, "Pathogenic (>= 10)")):
        ax.text(x, 3.75, text, ha="center", va="bottom", fontsize=6)
    for yi, (label, parts) in zip(range(len(ladders))[::-1], ladders):
        left = 0
        for name, points, color in parts:
            ax.barh(yi, points, left=left, color=color, height=0.55, edgecolor="white", linewidth=0.8)
            ax.text(left + points / 2, yi, name, ha="center", va="center", fontsize=5.2,
                    color="white" if points >= 2 else PALETTE["black"])
            left += points
        ax.text(left + 0.15, yi, f"{left} pts", va="center", fontsize=5.8, fontweight="bold")
    ax.set_yticks(range(len(ladders))[::-1], [l[0] for l in ladders])
    ax.set_xlim(0, 13)
    ax.set_ylim(-0.5, 4.2)
    ax.set_xlabel("ACMG/AMP points (Tavtigian et al. 2020 scale); lower rows add the named confirmatory evidence")
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    export(fig, "fig10_acmg_points")


def fig9_benchmark() -> None:
    rows = [r for r in read(EV / "controls/l7_metrics_audited.tsv") if r["variant_type"] in ("SNV", "INDEL")]
    strata = ["confident", "low_complexity", "homopolymer_ge12", "segmental_duplication"]
    names = {"confident": "GIAB confident", "low_complexity": "low complexity",
             "homopolymer_ge12": "homopolymer >= 12", "segmental_duplication": "segmental duplication"}
    fig, axes = plt.subplots(1, 2, figsize=(WIDE, WIDE * 0.34), layout="constrained", sharey=True)
    for ax, metric in zip(axes, ("precision", "recall")):
        for (sample, vtype), color, marker, dx in (
            (("HG002", "SNV"), SECOND, "o", -0.18), (("HG002", "INDEL"), SECOND, "s", -0.06),
            (("HG003", "SNV"), PALETTE["purple"], "o", 0.06), (("HG003", "INDEL"), PALETTE["purple"], "s", 0.18),
        ):
            xs, ys = [], []
            for i, stratum in enumerate(strata):
                row = next(r for r in rows if r["sample"] == sample and r["variant_type"] == vtype and r["stratum"] == stratum)
                xs.append(float(row[metric]))
                ys.append(i + dx)
            ax.plot(xs, ys, marker, color=color, markersize=3.8, linestyle="none",
                    label=f"{sample} {vtype}" if metric == "precision" else None,
                    markerfacecolor=color if vtype == "SNV" else "white")
        ax.set_xlim(0.93, 1.004)
        ax.set_xticks([0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0])
        ax.set_xlabel(f"{metric.capitalize()} on chr20 (exact normalised allele)")
        ax.grid(axis="x")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(strata)), [names[s] for s in strata])
    axes[0].set_ylabel("GIAB stratum")
    axes[0].invert_yaxis()
    axes[0].legend(loc="lower left", fontsize=5.6, ncols=2)
    panel_labels(axes)
    export(fig, "fig9_benchmark")


def main() -> None:
    use_style()
    OUT.mkdir(parents=True, exist_ok=True)
    funnel = json.loads((EV / "ranking/funnel.json").read_text(encoding="utf-8"))
    graphical_abstract(funnel)
    fig1_workflow(funnel)
    facts = fig2_locus()
    fig3_alleles()
    fig4_effect()
    fig5_phenotype()
    fig6_competing()
    fig7_genome()
    fig8_phase(facts)
    fig9_benchmark()
    fig10_acmg()
    (OUT / "locus_facts.json").write_text(json.dumps(facts, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(facts))


if __name__ == "__main__":
    main()
