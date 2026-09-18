#!/usr/bin/env python3
"""Test target-direction coherence against an artifact-shaped null.

Reversal directions are permuted within each regulator, preserving every
regulator's lost/gained counts. The production run uses 20,000 permutations
and belongs on a Slurm cluster; its compact output is committed in the repo.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--census",
        type=Path,
        default=ROOT / "data/census/census_master_edges.csv.gz",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "data/census/rna_content_confound.json",
    )
    parser.add_argument("--permutations", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20_260_911)
    return parser.parse_args()


def group_coherence(
    reversals: pd.DataFrame, key: str
) -> tuple[int, int, float]:
    grouped = reversals.groupby(key)["direction"].agg(
        n_dir="nunique", n_rev="size"
    )
    grouped = grouped[grouped["n_rev"] >= 2]
    single = int(grouped["n_dir"].eq(1).sum())
    return int(len(grouped)), single, float(single / len(grouped))


def main() -> int:
    args = parse_args()
    census = pd.read_csv(args.census)
    trans = census[~census["is_self"]].copy()
    reversals = trans[trans["kind"].eq("rev")].copy()
    if len(trans) != 83_489 or len(reversals) != 4_379:
        raise ValueError("input differs from the frozen trans-regulatory census")

    n_targets, target_single, target_fraction = group_coherence(
        reversals, "target"
    )
    n_regulators, regulator_single, regulator_fraction = group_coherence(
        reversals, "regulator"
    )

    rng = np.random.default_rng(args.seed)
    work = reversals[["regulator", "target", "direction"]].copy()
    directions = work["direction"].to_numpy()
    regulator_indices = work.groupby("regulator").indices
    targets = work["target"].to_numpy()
    draws = np.empty(args.permutations)
    n_as_or_more_extreme = 0
    for index in range(args.permutations):
        permuted = directions.copy()
        for positions in regulator_indices.values():
            permuted[positions] = rng.permutation(directions[positions])
        grouped = pd.DataFrame(
            {"target": targets, "direction": permuted}
        ).groupby("target")["direction"].agg(n_dir="nunique", n_rev="size")
        grouped = grouped[grouped["n_rev"] >= 2]
        fraction = float(grouped["n_dir"].eq(1).mean())
        draws[index] = fraction
        n_as_or_more_extreme += fraction >= target_fraction

    result = {
        "_script": "scripts/analysis/rna_content_direction_coherence.py",
        "_source": "data/census/census_master_edges.csv.gz",
        "_motivation": (
            "test whether per-regulator global RNA-content shifts can explain "
            "the target-level organisation of reversal directions"
        ),
        "_claim_boundary": (
            "rules out a per-regulator global RNA-content artifact as the "
            "source of reversal organisation; it does not measure the "
            "mu-parallel response itself"
        ),
        "census": {
            "trans_edges": int(len(trans)),
            "reversals": int(len(reversals)),
            "regulators": int(trans["regulator"].nunique()),
            "targets": int(trans["target"].nunique()),
        },
        "direction_coherence": {
            "_what": (
                "reversal direction is coherent within targets more than "
                "expected after preserving each regulator's direction mix"
            ),
            "min_reversals_per_group": 2,
            "targets": {
                "n_groups": n_targets,
                "single_orientation": target_single,
                "fraction": target_fraction,
            },
            "regulators": {
                "n_groups": n_regulators,
                "single_orientation": regulator_single,
                "fraction": regulator_fraction,
            },
            "within_regulator_permutation_null": {
                "_design": (
                    "direction shuffled within each regulator while holding "
                    "its lost/gained counts fixed"
                ),
                "n_draws": args.permutations,
                "seed": args.seed,
                "null_mean": float(draws.mean()),
                "null_sd": float(draws.std()),
                "observed": target_fraction,
                "plus_one_p": float(
                    (n_as_or_more_extreme + 1) / (args.permutations + 1)
                ),
            },
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
