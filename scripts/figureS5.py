#!/usr/bin/env python3
"""Render Supplementary Figure S5 from the committed compact outputs."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import figstyle as style  # noqa: E402


def main() -> int:
    coherence = json.loads(
        (ROOT / "data/census/rna_content_confound.json").read_text()
    )
    full_atlas = json.loads(
        (
            ROOT
            / "data/census/rna_content_umi_shift_full_atlas.json"
        ).read_text()
    )
    regulators = pd.read_csv(
        ROOT / "data/census/rna_content_per_regulator.tsv", sep="\t"
    )
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)

    direction = coherence["direction_coherence"]
    null = direction["within_regulator_permutation_null"]
    direct_test = full_atlas["direction_vs_rna_shift"][
        "regulator_level_spearman"
    ]

    assert direction["targets"]["n_groups"] == 913
    assert direction["targets"]["single_orientation"] == 633
    assert direction["regulators"]["n_groups"] == 227
    assert null["plus_one_p"] < 1e-3
    assert full_atlas["label"] == "full_atlas"
    assert full_atlas["n_cells"] == 14_853_023
    assert full_atlas["n_regulators_matched_to_census"] == 616
    assert full_atlas["ntc_median_umi"]["Rest"] == 7_674.0
    assert full_atlas["ntc_median_umi"]["Stim48hr"] == 11_105.0
    assert direct_test["n"] == 187
    assert len(regulators) == full_atlas["direction_vs_rna_shift"][
        "n_regulators_scored"
    ]

    style.setup()
    figure = plt.figure(figsize=(7.2, 5.2), facecolor="white")
    grid = figure.add_gridspec(
        2,
        2,
        wspace=0.42,
        hspace=0.62,
        left=0.11,
        right=0.985,
        top=0.89,
        bottom=0.10,
    )

    axis = figure.add_subplot(grid[0, 0])
    umi = full_atlas["ntc_median_umi"]
    states = ["Rest", "Stim48hr"]
    values = [umi[state] for state in states]
    axis.bar(
        range(2),
        values,
        color=[style.PAL["REST"], style.PAL["STIM"]],
        width=0.62,
    )
    for index, value in enumerate(values):
        axis.text(
            index,
            value * 1.02,
            f"{value:,.0f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    axis.set_xticks(range(2))
    axis.set_xticklabels(["Rest", "Stim 48 h"])
    axis.set_ylabel("median total UMI\n(non-targeting cells)")
    axis.set_ylim(0, max(values) * 1.22)
    axis.set_title(
        "Global RNA content rises\n"
        f"{full_atlas['ntc_umi_fold_stim_over_rest']:.2f}-fold on activation",
        fontsize=8,
    )
    style.panel_tag(axis, "a", dx=-0.26, dy=1.16)

    axis = figure.add_subplot(grid[0, 1])
    fractions = [
        direction["targets"]["fraction"],
        direction["regulators"]["fraction"],
    ]
    axis.bar(
        range(2),
        [fraction * 100 for fraction in fractions],
        color=[style.PAL["CORE"], style.PAL["GREY"]],
        width=0.62,
    )
    for index, (fraction, key) in enumerate(
        zip(fractions, ["targets", "regulators"])
    ):
        values_for_key = direction[key]
        axis.text(
            index,
            fraction * 100 + 1.8,
            f"{values_for_key['single_orientation']}/"
            f"{values_for_key['n_groups']}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    axis.set_xticks(range(2))
    axis.set_xticklabels(["by target", "by regulator"])
    axis.set_ylabel("single-orientation (%)")
    axis.set_ylim(0, 100)
    axis.set_title(
        "Direction is coherent within\ntargets, not within regulators",
        fontsize=8,
    )
    style.panel_tag(axis, "b", dx=-0.26, dy=1.16)

    axis = figure.add_subplot(grid[1, 0])
    null_mean = null["null_mean"]
    null_sd = null["null_sd"]
    x_values = np.linspace(
        null_mean - 5 * null_sd,
        max(null_mean + 5 * null_sd, null["observed"] + null_sd),
        400,
    )
    axis.fill_between(
        x_values,
        np.exp(-0.5 * ((x_values - null_mean) / null_sd) ** 2),
        color=style.PAL["CONC"],
        linewidth=0,
        label=f"within-regulator null\n({null['n_draws']:,} draws)",
    )
    axis.axvline(null["observed"], color=style.PAL["LOST"], linewidth=1.6)
    axis.text(
        null["observed"] - 0.006,
        0.55,
        f"observed\n{null['observed']:.1%}",
        ha="right",
        va="center",
        fontsize=7,
        color=style.PAL["LOST"],
    )
    axis.set_xlabel("targets with a single reversal orientation")
    axis.set_ylabel("null density")
    axis.set_yticks([])
    axis.set_xlim(x_values.min(), max(x_values.max(), null["observed"] * 1.02))
    axis.set_ylim(0, 1.18)
    axis.set_title(
        f"$p$ = {null['plus_one_p']:.0e} against a\nper-regulator artifact",
        fontsize=8,
    )
    axis.legend(
        loc="upper left",
        fontsize=6,
        frameon=False,
        handlelength=1.1,
        borderaxespad=0.2,
    )
    style.panel_tag(axis, "c", dx=-0.20, dy=1.16)

    axis = figure.add_subplot(grid[1, 1])
    subset = regulators[
        regulators["n_rev"] >= direct_test["min_reversals_per_regulator"]
    ]
    assert len(subset) == direct_test["n"]
    axis.axhline(50, color=style.PAL["GREY"], linewidth=0.7, linestyle=":")
    axis.axvline(0, color=style.PAL["GREY"], linewidth=0.7, linestyle=":")
    axis.scatter(
        subset["delta_shift"],
        subset["frac_gained"] * 100,
        s=np.clip(subset["n_rev"] * 1.6, 6, 60),
        facecolor=style.PAL["CORE"],
        edgecolor="white",
        linewidth=0.35,
        alpha=0.75,
        zorder=3,
    )
    axis.set_xlabel(
        "change in the regulator's global RNA output\n"
        r"$\log_2$(perturbed/control), Stim 48 h minus Rest"
    )
    axis.set_ylabel("reversals running\n'gained on activation' (%)")
    axis.set_ylim(-6, 106)
    axis.set_title(
        "Reversal direction does not follow\nthe RNA-content shift "
        + rf"($\rho$ = {direct_test['rho']:+.3f}, $p$ = {direct_test['p']:.2f})",
        fontsize=8,
    )
    axis.text(
        0.02,
        0.03,
        f"{direct_test['n']} regulators with "
        + r"$\geq$"
        + "3 reversals\n"
        + f"full atlas, {full_atlas['n_cells'] / 1e6:.1f}M cells",
        transform=axis.transAxes,
        fontsize=6,
        va="bottom",
        ha="left",
        color=style.PAL["GREY"],
    )
    style.panel_tag(axis, "d", dx=-0.22, dy=1.16)

    for extension in ("pdf", "png"):
        figure.savefig(
            output / f"FigS5.{extension}",
            dpi=300,
            facecolor="white",
            bbox_inches="tight",
        )
    plt.close(figure)
    print(f"wrote {output / 'FigS5.pdf'} and {output / 'FigS5.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
