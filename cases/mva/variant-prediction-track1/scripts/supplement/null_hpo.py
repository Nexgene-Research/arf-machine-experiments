"""HPO-weight permutation null for the frozen angle-B matrix.

Does not reread the phenotype document. Permutes the seven term weights onto
the existing per-gene match levels. Seed 20261017, 100 permutations.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

DEST = Path(__file__).resolve().parents[1]
TABLE = DEST / "compiled" / "literature" / "angle_b_phenotype.tsv"
WEIGHTS = {
    "HP:0002859": 5.0,
    "HP:0000121": 3.0,
    "HP:0004322": 2.0,
    "HP:0001508": 1.0,
    "HP:0003202": 1.0,
    "HP:0001622": 0.5,
    "HP:0001518": 1.0,
}
TERMS = list(WEIGHTS)
SEED = 20261017
N = 100


def parse_levels(text: str) -> dict[str, int]:
    out = {}
    for item in text.split(";"):
        term, value = item.rsplit(":", 1)
        out[term] = int(value)
    return out


def score(levels: dict[str, int], weights: dict[str, float]) -> float:
    return sum(weights[term] * levels.get(term, 0) / 2 for term in TERMS)


def main() -> None:
    with TABLE.open(encoding="utf-8", newline="") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    observed = []
    for row in genes:
        levels = parse_levels(row["feature_match_levels_0_unsubstantiated_1_broad_2_exact"])
        observed.append((score(levels, WEIGHTS), row["gene"], levels))
    observed.sort(reverse=True)
    top = observed[0][1]
    rng = random.Random(SEED)
    top_counts = {row["gene"]: 0 for row in genes}
    beat = 0
    for _ in range(N):
        shuffled = TERMS[:]
        rng.shuffle(shuffled)
        remap = dict(zip(TERMS, shuffled))
        weights = {term: WEIGHTS[remap[term]] for term in TERMS}
        ranked = sorted(
            ((score(levels, weights), gene) for _s, gene, levels in observed),
            reverse=True,
        )
        top_counts[ranked[0][1]] += 1
        if ranked[0][1] != top:
            beat += 1
    dest = DEST / "computed" / "null" / "null_runs.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "arm\titem\tvalue\tnote\n",
        f"hpo_weight_permutation\tobserved_top\t{top}\tangle_B weighted coverage\n",
        f"hpo_weight_permutation\tpermutations\t{N}\tseed={SEED}\n",
        f"hpo_weight_permutation\ttimes_top_changes\t{beat}\tother gene outranks {top}\n",
    ]
    for gene, count in sorted(top_counts.items(), key=lambda item: -item[1]):
        lines.append(
            f"hpo_weight_permutation\ttimes_rank1_{gene}\t{count}\t"
            "weight permutation of seven HPO terms; not a genome-wide null\n"
        )
    dest.write_text("".join(lines), encoding="utf-8")
    print(f"null hpo top={top} changed={beat}/{N}")


if __name__ == "__main__":
    main()
