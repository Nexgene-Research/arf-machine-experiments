"""Aggregate plot of measured Phase C pair statuses across Tier-1 genes."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from figstyle import PALETTE, TEXT, figure, save, use_style

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/track1/phase_pair_evidence.tsv"
BLOCKS = ROOT / "results/track1/phase_blocks_summary.tsv"
STEM = ROOT / "results/track1/phase_aggregate"


def main() -> None:
    if any(STEM.with_suffix(ext).exists() for ext in (".pdf", ".svg", ".png")):
        raise FileExistsError("phase figure already saved")
    with SOURCE.open(newline="", encoding="utf-8") as fh:
        pairs = list(csv.DictReader(fh, delimiter="\t"))
    with BLOCKS.open(newline="", encoding="utf-8") as fh:
        block_rows = list(csv.DictReader(fh, delimiter="\t"))
    labels = ["BUB1", "BUB1B"]
    counts = Counter((r["gene"], r["phase"]) for r in pairs)
    if any(r["phase"] not in {"unknown", "trans", "cis"} for r in pairs):
        raise ValueError("unrecognized status")
    use_style()
    fig, ax = figure(width=TEXT, ratio=0.52)
    left = [0, 0]
    for status, color in (
        ("trans", PALETTE["blue"]),
        ("cis", PALETTE["red"]),
        ("unknown", PALETTE["orange"]),
    ):
        value = [counts[(gene, status)] for gene in labels]
        ax.barh(
            labels,
            value,
            left=left,
            color=color,
            edgecolor="#333333",
            linewidth=0.5,
            label=status,
        )
        left = [x + y for x, y in zip(left, value)]
    ax.set_xlim(0, max(left) + 1)
    ax.set_xlabel("Candidate pairs assessed (count)")
    ax.set_ylabel("Tier-1 gene with tested pairs")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    ax.text(
        0.98,
        0.98,
        f"{len(block_rows)} regional phase blocks; no pair linked",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
    )
    for y, total in enumerate(left):
        ax.text(total + 0.08, y, str(total), va="center", fontsize=7)
    save(fig, str(STEM), close=False)
    fig.savefig(STEM.with_suffix(".png"), dpi=300)
    import matplotlib.pyplot as plt

    plt.close(fig)
    print(
        f"plotted candidate pairs={len(pairs)}; regional blocks={len(block_rows)}; unknown={sum(r['phase'] == 'unknown' for r in pairs)}"
    )


if __name__ == "__main__":
    main()
