#!/usr/bin/env python3
"""Test whether state-dependent total-RNA shifts explain census reversals.

The input is the full-atlas per-cell total-UMI table with columns
``is_ntc``, ``condition``, ``n_umi`` and ``perturbed_gene``. The table is
derived from the source atlas and is not redistributed because of its size.
The script writes the compact summary and per-regulator table used by
Supplementary Figure S5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
STATE_A = "Rest"
STATE_B = "Stim48hr"


def partial_spearman(
    x: pd.Series, y: pd.Series, covariates: list[pd.Series]
) -> tuple[float, float, int]:
    """Rank-residual partial Spearman correlation."""
    n = len(x)
    k = len(covariates)
    ranked_x = stats.rankdata(x)
    ranked_y = stats.rankdata(y)
    design = np.column_stack(
        [np.ones(n)] + [stats.rankdata(covariate) for covariate in covariates]
    )
    residual_x = ranked_x - design @ np.linalg.lstsq(
        design, ranked_x, rcond=None
    )[0]
    residual_y = ranked_y - design @ np.linalg.lstsq(
        design, ranked_y, rcond=None
    )[0]
    rho = float(np.corrcoef(residual_x, residual_y)[0, 1])
    degrees_freedom = n - 2 - k
    test_statistic = rho * np.sqrt(
        degrees_freedom / max(1e-12, 1 - rho**2)
    )
    p_value = float(2 * stats.t.sf(abs(test_statistic), degrees_freedom))
    return rho, p_value, int(n)


def redacted_path(path: Path) -> str:
    text = str(path.resolve())
    home = str(Path.home())
    return text.replace(home, "$HOME", 1) if text.startswith(home) else text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--percell",
        type=Path,
        required=True,
        help="full-atlas per-cell total-UMI CSV or CSV.GZ",
    )
    parser.add_argument(
        "--census",
        type=Path,
        default=ROOT / "data/census/census_master_edges.csv.gz",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=ROOT / "data/census/rna_content_umi_shift_full_atlas.json",
    )
    parser.add_argument(
        "--regulator-out",
        type=Path,
        default=ROOT / "data/census/rna_content_per_regulator.tsv",
    )
    parser.add_argument("--min-cells", type=int, default=25)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cells = pd.read_csv(args.percell)
    census = pd.read_csv(args.census)

    required_cell_columns = {"is_ntc", "condition", "n_umi", "perturbed_gene"}
    missing = required_cell_columns.difference(cells.columns)
    if missing:
        raise ValueError(f"per-cell table lacks columns: {sorted(missing)}")

    if len(census) != 84_044:
        raise ValueError(f"census has {len(census)} edges; expected 84,044")
    trans = census[~census["is_self"]].copy()
    if len(trans) != 83_489:
        raise ValueError(f"census has {len(trans)} trans edges; expected 83,489")
    if int(trans["kind"].eq("rev").sum()) != 4_379:
        raise ValueError("reversal count differs from the frozen value 4,379")

    ntc = (
        cells[cells["is_ntc"]]
        .groupby("condition")["n_umi"]
        .median()
        .rename("ntc_umi")
    )
    perturbed = (
        cells[~cells["is_ntc"]]
        .groupby(["perturbed_gene", "condition"])["n_umi"]
        .agg(["median", "size"])
        .rename(columns={"median": "pert_umi", "size": "n_cells"})
        .reset_index()
    )
    perturbed = perturbed.join(ntc, on="condition")
    perturbed["log2_shift"] = np.log2(
        perturbed["pert_umi"] / perturbed["ntc_umi"]
    )

    shifts = perturbed.pivot(
        index="perturbed_gene", columns="condition", values="log2_shift"
    )
    counts = perturbed.pivot(
        index="perturbed_gene", columns="condition", values="n_cells"
    )
    for state in (STATE_A, STATE_B):
        if state not in shifts.columns:
            raise ValueError(f"{state} is missing from the per-cell table")
    shifts = shifts[[STATE_A, STATE_B]].dropna()
    counts = counts.reindex(shifts.index)[[STATE_A, STATE_B]]
    n_tested = len(shifts)
    keep = (counts >= args.min_cells).all(axis=1)
    shifts = shifts.loc[keep].copy()
    counts = counts.loc[keep]
    shifts["opposite_sign"] = np.sign(shifts[STATE_A]) != np.sign(
        shifts[STATE_B]
    )
    shifts["shift_gap"] = (shifts[STATE_B] - shifts[STATE_A]).abs()
    shifts["n_cells"] = counts[STATE_A] + counts[STATE_B]

    burden = (
        trans.assign(rev=trans["kind"].eq("rev"))
        .groupby("regulator")
        .agg(n_edges=("rev", "size"), n_rev=("rev", "sum"))
    )
    burden["rev_rate"] = burden["n_rev"] / burden["n_edges"]
    symbols = (
        trans[["regulator", "regulator_ensg"]]
        .drop_duplicates()
        .set_index("regulator_ensg")["regulator"]
    )
    shifts["regulator"] = shifts.index.map(symbols)
    tested = (
        shifts.dropna(subset=["regulator"])
        .join(burden, on="regulator")
        .dropna(subset=["rev_rate"])
    )

    if len(tested) < 20:
        raise ValueError(
            f"only {len(tested)} census regulators passed the filters; "
            "the full-atlas input appears incomplete"
        )

    opposite = tested[tested["opposite_sign"]]
    stable = tested[~tested["opposite_sign"]]
    mann_whitney = stats.mannwhitneyu(
        opposite["rev_rate"], stable["rev_rate"], alternative="two-sided"
    )
    cliffs_delta = (
        2 * mann_whitney.statistic / (len(opposite) * len(stable)) - 1
    )

    carrier = tested["n_rev"] > 0
    carrier_table = np.array(
        [
            [
                int((carrier & tested["opposite_sign"]).sum()),
                int((~carrier & tested["opposite_sign"]).sum()),
            ],
            [
                int((carrier & ~tested["opposite_sign"]).sum()),
                int((~carrier & ~tested["opposite_sign"]).sum()),
            ],
        ]
    )
    fisher_odds, fisher_p = stats.fisher_exact(carrier_table)

    shift_gap_spearman = stats.spearmanr(
        tested["shift_gap"], tested["rev_rate"]
    )
    partial_results: dict[str, dict[str, float | int]] = {}
    for name, covariates in (
        ("n_edges", [tested["n_edges"]]),
        ("n_cells", [tested["n_cells"]]),
        ("n_edges+n_cells", [tested["n_edges"], tested["n_cells"]]),
    ):
        rho, p_value, n = partial_spearman(
            tested["shift_gap"], tested["rev_rate"], covariates
        )
        partial_results[name] = {"rho": rho, "p": p_value, "n": n}

    covariate_results: dict[str, dict[str, float]] = {}
    for name in ("n_edges", "n_cells"):
        correlation = stats.spearmanr(tested["shift_gap"], tested[name])
        covariate_results[name] = {
            "rho": float(correlation.statistic),
            "p": float(correlation.pvalue),
        }

    quartile = pd.qcut(
        tested["n_edges"], 4, labels=False, duplicates="drop"
    )
    numerator = 0.0
    denominator = 0.0
    strata: list[dict[str, int]] = []
    for level in sorted(pd.unique(quartile)):
        in_stratum = quartile == level
        t11 = int((carrier & tested["opposite_sign"] & in_stratum).sum())
        t12 = int((~carrier & tested["opposite_sign"] & in_stratum).sum())
        t21 = int((carrier & ~tested["opposite_sign"] & in_stratum).sum())
        t22 = int((~carrier & ~tested["opposite_sign"] & in_stratum).sum())
        total = t11 + t12 + t21 + t22
        if total == 0:
            continue
        numerator += t11 * t22 / total
        denominator += t12 * t21 / total
        strata.append(
            {
                "edge_quartile": int(level),
                "opp_carrier": t11,
                "opp_noncarrier": t12,
                "same_carrier": t21,
                "same_noncarrier": t22,
            }
        )

    reversal_edges = trans[trans["kind"].eq("rev")].copy()
    delta_shift = (tested[STATE_B] - tested[STATE_A]).rename("delta_shift")
    regulator_delta = pd.Series(
        delta_shift.to_numpy(), index=tested["regulator"].to_numpy()
    )
    regulator_delta = regulator_delta[~regulator_delta.index.duplicated()]
    reversal_edges["delta_shift"] = reversal_edges["regulator"].map(
        regulator_delta
    )
    scored = reversal_edges.dropna(subset=["delta_shift"])
    per_regulator = (
        scored.assign(
            gained=scored["direction"].eq("gained_on_activation")
        )
        .groupby("regulator")
        .agg(frac_gained=("gained", "mean"), n_rev=("gained", "size"))
    )
    per_regulator["delta_shift"] = per_regulator.index.map(regulator_delta)
    direction_tested = per_regulator[per_regulator["n_rev"] >= 3]
    direction_spearman = stats.spearmanr(
        direction_tested["delta_shift"], direction_tested["frac_gained"]
    )

    gained = scored.loc[
        scored["direction"].eq("gained_on_activation"), "delta_shift"
    ]
    lost = scored.loc[
        scored["direction"].eq("lost_on_activation"), "delta_shift"
    ]
    edge_test = stats.mannwhitneyu(gained, lost, alternative="two-sided")
    edge_cliffs = 2 * edge_test.statistic / (len(gained) * len(lost)) - 1

    reversal_by_target = reversal_edges.groupby("target")["direction"].agg(
        ["nunique", "size"]
    )
    multi_target = reversal_by_target[reversal_by_target["size"] >= 2]

    summary = {
        "percell_table": redacted_path(args.percell),
        "label": "full_atlas",
        "min_cells_per_state": args.min_cells,
        "n_cells": int(len(cells)),
        "ntc_median_umi": {key: float(value) for key, value in ntc.items()},
        "ntc_umi_fold_stim_over_rest": float(ntc[STATE_B] / ntc[STATE_A]),
        "n_perturbations_tested": int(n_tested),
        "n_perturbations_passing_cell_floor": int(len(shifts)),
        "n_regulators_matched_to_census": int(len(tested)),
        "analysis_status": "RAN",
        "opposite_sign": {
            "n": int(len(opposite)),
            "median_reversal_rate": float(opposite["rev_rate"].median()),
        },
        "same_sign": {
            "n": int(len(stable)),
            "median_reversal_rate": float(stable["rev_rate"].median()),
        },
        "mannwhitney": {
            "U": float(mann_whitney.statistic),
            "p": float(mann_whitney.pvalue),
            "cliffs_delta": float(cliffs_delta),
        },
        "carrier_fisher": {
            "table_opp_carrier_noncarrier_same_carrier_noncarrier": [
                int(value) for value in carrier_table.ravel()
            ],
            "odds_ratio": float(fisher_odds),
            "p": float(fisher_p),
        },
        "spearman_shiftgap_vs_revrate": {
            "rho": float(shift_gap_spearman.statistic),
            "p": float(shift_gap_spearman.pvalue),
            "n": int(len(tested)),
        },
        "partial_spearman": partial_results,
        "shiftgap_vs_covariate": covariate_results,
        "carrier_cmh_by_edge_quartile": {
            "odds_ratio": float(numerator / denominator),
            "strata": strata,
        },
        "direction_vs_rna_shift": {
            "n_reversal_edges_scored": int(len(scored)),
            "n_regulators_scored": int(scored["regulator"].nunique()),
            "regulator_level_spearman": {
                "rho": float(direction_spearman.statistic),
                "p": float(direction_spearman.pvalue),
                "n": int(len(direction_tested)),
                "min_reversals_per_regulator": 3,
            },
            "edge_level_mannwhitney": {
                "median_delta_gained": float(gained.median()),
                "median_delta_lost": float(lost.median()),
                "p": float(edge_test.pvalue),
                "cliffs_delta": float(edge_cliffs),
                "caveat": "edges are clustered within regulators; not independent",
            },
        },
        "target_direction_coherence": {
            "n_targets_with_ge2_reversals": int(len(multi_target)),
            "single_orientation": int(
                multi_target["nunique"].eq(1).sum()
            ),
            "fraction": float(multi_target["nunique"].eq(1).mean()),
        },
    }

    per_regulator["shift_Rest"] = per_regulator.index.map(
        tested.set_index("regulator")[STATE_A]
    )
    per_regulator["shift_Stim48hr"] = per_regulator.index.map(
        tested.set_index("regulator")[STATE_B]
    )
    per_regulator["n_edges"] = per_regulator.index.map(
        tested.set_index("regulator")["n_edges"]
    )

    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.regulator_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    per_regulator.to_csv(args.regulator_out, sep="\t")
    print(
        f"wrote {args.summary_out} and {args.regulator_out} "
        f"for {len(tested)} census regulators"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
