#!/usr/bin/env python3
"""Supplementary Figure S3 - single self-contained builder.

S3 carries two analyses the manuscript keeps in one supplementary figure:

  a-f  two-donor Arce protein-channel composition pilot, built here
  g-i  reconstruction geometry and estimator robustness, demoted from Figure 4

Both halves are built here and merged as vectors, so the script is standalone: a-f
into figures/FigS3_base.pdf, g-i into figures/Fig4_supp.pdf. scripts/figure4.py writes
a byte-identical Fig4_supp.pdf from the same frozen bundle, but neither script needs
the other to have run.

No statistic is recomputed. The g-i half reads the sealed v5.1 production tables in
data/response_shape/production_v5_1/ and checks each one against the SHA-256 recorded
in their manifest; the full custody/validation chain stays in scripts/figure4.py, which
remains the auditor of that bundle. This is a renderer.

Inputs previously read from a directory outside the repository are now committed
under data/state_decomposition/.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pypdf import PdfReader, PdfWriter, Transformation, PageObject

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'data/state_decomposition'
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
GAP = 12.0  # pt of whitespace between the two source figures

# ----------------------------- shared style (was scripts/figstyle.py) ---------
PAL = dict(
    REST="#2c6fbb", STIM="#e07b39",
    LOST="#c1272d",     # up in Rest -> down in Stim48  (edge lost on stimulation)
    GAINED="#2a9d8f",   # down in Rest -> up in Stim48  (edge gained on stimulation)
    GREY="#9aa0a6", CONC="#c7ccd1", CORE="#6a51a3",
    INK="#22262b", HAIR="#dfe2e6", HL="#f2c14e",
)


def setup(scale=1.0):
    def s(v):
        return round(v * scale, 2)
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "DejaVu Sans", "font.size": s(9),
        "axes.titlesize": s(10), "axes.titleweight": "bold", "axes.labelsize": s(9),
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#565b61", "axes.linewidth": 0.9,
        "xtick.labelsize": s(8), "ytick.labelsize": s(8), "xtick.color": "#565b61",
        "ytick.color": "#565b61", "legend.fontsize": s(7.5), "legend.frameon": False,
        "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
    })


def panel_tag(ax, s, dx=-0.04, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="right", color="#22262b")


# ----------------------------- a-f: state decomposition -----------------------
#!/usr/bin/env python
"""Corrected state-decomposition figure. Writes Fig5.*
(the v1 figure is left untouched).

WHAT CHANGED AND WHY -- every panel here replaces one that sat on its own null:

v1 panel a/b (reconstructed vs published, rho=0.956/0.969)  -> DROPPED to Methods.
    total_log2fc is MATHEMATICALLY INVARIANT to the state model (f@mu is the pooled
    mean by the law of total expectation), so this check returns rho=0.956 even under
    randomly shuffled state labels. It validates plumbing, not the decomposition.
v1 panel e ("Most reversals survive composition adjustment", 115/116 and 54/54) -> DELETED.
    Under a control-vs-control null where true composition is exactly zero, expected
    retention is 115.0/116 = 99.1%. Observed 99.1%. Only 10/116 flips are even eligible
    for removal; for the bolded 54/54 bar, 0/54 are eligible. It is arithmetic.
v1 panel f (MED12 relays) -> DELETED. MED12->KLF2 has the same sign in both states, so
    the -sign(med12 x intermediate) factor cancels exactly and relay_predicts_flip
    reduces to sign(inter_rest) != sign(inter_stim). It is a reversal-overlap test
    mislabelled as a relay test.

What replaces them is the test that actually discriminates: composition adjustment is
only informative where it has DYNAMIC RANGE, and where it does (multidim, 30.2% of flips
reachable) it destroys reversals and non-reversal sign flips at the SAME rate --
OR=0.93, p=0.87. Simpson's paradox requires specificity. This is its refutation.

All numbers are read from decomposition_null_results.json (decomp_null_summary.py);
nothing is hardcoded. log_fc is log2 (DESeq2). Usage: fig8_state_decomposition_v2.py
"""




setup()
C_REV, C_STABLE = PAL["CORE"], PAL["CONC"]
C_REST, C_STIM = PAL["REST"], PAL["STIM"]
C_INK, C_GREY, C_HAIR = PAL["INK"], PAL["GREY"], PAL["HAIR"]
C_NULL = PAL["LOST"]

N = json.load(open(AUDIT / "decomposition_null_results.json"))
primary = json.load(open(AUDIT / "decomposition_results.json"))
occ = pd.read_csv(AUDIT / "occupancy_effects.csv")

# Substrate panels (a-c): UMAP of the ADT protein channel these states are read from.
# Built by scripts/arce_adt_umap.py, which drops the four HTO hashtags (the state label is their
# argmax -- leaving them in makes the embedding circular) and the isotype controls.
UM = pd.read_csv(ROOT / "data/arce/arce_adt_umap_coords.csv.gz")
UST = json.load(open(ROOT / "data/arce/arce_adt_umap_stats.json"))
_rng = np.random.default_rng(0)


def _umap(ax, key, colours, title, order, tag):
    """One UMAP panel, coloured by `key`. Points drawn in random order so no group is buried."""
    idx = _rng.permutation(len(UM))
    g = UM[key].to_numpy()[idx]
    ax.scatter(UM.u1.to_numpy()[idx], UM.u2.to_numpy()[idx], s=0.6,
               c=[colours[v] for v in g], linewidths=0, alpha=0.35, rasterized=True)
    for k in order:
        ax.scatter([], [], s=18, c=colours[k], label=k, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_box_aspect(1)
    ax.set_title(title, pad=8)
    ax.legend(loc="lower right", handletextpad=0.4, borderpad=0.3)
    panel_tag(ax, tag, dx=-0.16, dy=1.08)


def umap_state(ax):
    _umap(ax, "act", {"Rest": C_REST, "Stim": C_STIM},
          f"The protein activation axis is real\n(AUROC {UST['auroc_activation_all']:.3f}, "
          f"{UST['n_features_used']} surface proteins)", ["Rest", "Stim"], "a")


def umap_donor(ax):
    d = UST["auroc_activation_by_donor"]
    _umap(ax, "donor", {"A": C_REV, "B": PAL["GAINED"]},
          "and both donors span it\n(" + ", ".join(f"{k} {v:.3f}" for k, v in sorted(d.items())) + ")",
          sorted(d), "b")


def umap_subset(ax):
    _umap(ax, "subset", {"Teff": C_INK, "Treg": PAL["HL"]},
          "Teff and Treg separate within state\n(d-f use the "
          f"{UST['n_by_subset']['Teff']:,} Teff cells)", ["Teff", "Treg"], "c")


def occupancy(ax):
    """Composition is REAL: perturbations move ADT state occupancy, reproducibly."""
    for cond, colour, marker in (("Resting", C_REST, "o"), ("Stimulated", C_STIM, "s")):
        p = occ[occ.condition.eq(cond)].pivot(index="regulator", columns="donor",
                                              values="score_shift").dropna()
        d = N["occupancy_donor_reproducibility"][cond]
        ax.scatter(p.iloc[:, 0], p.iloc[:, 1], s=28, marker=marker, color=colour,
                   alpha=0.78, edgecolor="#565b61", lw=0.45, zorder=3,
                   label=rf"{cond.replace('ing','')}: $\rho$={d['spearman_rho']:.2f} "
                         rf"[{d['boot95'][0]:.2f}, {d['boot95'][1]:.2f}]")
        if "MED12" in p.index:
            x, y = p.loc["MED12"].to_numpy()[:2]
            ax.annotate("MED12", (x, y), xytext=(5, 5 if cond == "Resting" else -10),
                        textcoords="offset points", fontsize=6.6, fontweight="bold",
                        color=colour, zorder=5)
    lo, hi = -0.09, 0.13
    ax.plot([lo, hi], [lo, hi], color=C_INK, lw=0.9, ls="--", zorder=1)
    ax.axhline(0, color=C_HAIR, lw=0.8, zorder=0)
    ax.axvline(0, color=C_HAIR, lw=0.8, zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_box_aspect(1)  # square box, not set_aspect -- keeps the title level with b/c
    ax.set_xlabel("donor A ADT activation-score shift")
    ax.set_ylabel("donor B ADT activation-score shift")
    ax.set_title("State redistribution is real\nand reproduces across donors", pad=8)
    ax.legend(loc="upper left", handletextpad=0.5, labelspacing=0.5)
    n = N["occupancy_donor_reproducibility"]["Resting"]["n_regulators"]
    ax.text(0.97, 0.04, f"n = {n} perturbations", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=6.4, color="#737980")
    panel_tag(ax, "d", dx=-0.16, dy=1.08)


def matched_contrast(ax):
    """The contrast -- read against its PERMUTATION NULL CENTRE, which is not zero."""
    m = primary["matched_primary"]
    per = pd.Series(m["per_regulator"]).sort_values(ascending=False)
    est = N["matched_estimator"]
    centre = est["median_stat_null_centre"]
    # summary diamond at the BOTTOM (forest-plot convention); regulators above it
    y = np.arange(len(per)) + 1.6
    ax.scatter(per.to_numpy(), y, s=38, color=[C_REV if v > 0 else C_GREY for v in per],
               edgecolor="#565b61", lw=0.5, zorder=3)
    # the null centre is the honest reference line; zero is NOT the null
    ax.axvspan(centre - est["median_stat_null_sd"], centre + est["median_stat_null_sd"],
               color=C_NULL, alpha=0.10, zorder=0)
    ax.axvline(centre, color=C_NULL, lw=1.4, ls="-", zorder=2)
    ax.axvline(0, color=C_HAIR, lw=1.0, ls=":", zorder=1)
    lo, hi = m["cluster_bootstrap_ci"]
    obs = m["mean_regulator_median_difference"]
    ax.axhline(0.9, color=C_HAIR, lw=0.9)
    ax.errorbar(obs, 0.2, xerr=[[obs - lo], [hi - obs]], fmt="D", ms=6, color=C_INK,
                ecolor=C_INK, elinewidth=2.0, capsize=4, capthick=1.8, zorder=4)
    ax.set_yticks([0.2] + list(y))
    ax.set_yticklabels(["regulator\nsummary"] + list(per.index))
    ax.get_yticklabels()[0].set_fontweight("bold")
    ax.set_ylim(-1.5, len(per) + 2.2)
    ax.set_xlim(-0.058, 0.062)
    ax.set_box_aspect(1)
    ax.set_xlabel("reversal minus matched stable composition fraction")
    ax.set_title("Reversal edges show no detectable\nexcess composition", pad=8)
    # label the reference lines in place -- a legend box collides with the summary row
    ax.annotate(f"permutation null\ncentre ({centre:+.3f})", (centre, len(per) + 1.5),
                xytext=(-4, 0), textcoords="offset points", fontsize=6.1, color=C_NULL,
                ha="right", va="center", fontweight="bold")
    ax.annotate("zero", (0, len(per) + 1.5), xytext=(5, 0), textcoords="offset points",
                fontsize=6.1, color="#737980", ha="left", va="center")
    ax.annotate(f"centred effect {est['median_stat_centred_effect']:+.3f},  "
                f"permutation $p$ = {est['median_stat_p_one_sided']:.2f}",
                (0.002, -1.0), fontsize=6.5, color="#737980", ha="center", va="center")
    panel_tag(ax, "e", dx=-0.16, dy=1.08)


def flip_loss(ax):
    """THE RESULT: where the adjustment has range, it is indifferent to reversal."""
    models = [("k3", "1-D activation\nscore (k=3)"), ("multidim", "10-PC ADT\n(K=3)")]
    x = np.arange(len(models))
    width = 0.34
    for off, kk, colour, lab in ((-width / 2, "rev", C_REV, "reversal"),
                                 (+width / 2, "stable", C_STABLE, "stable (comparator)")):
        rates = [N["flip_loss_specificity"][m][f"{kk}_rate"] for m, _ in models]
        bars = ax.bar(x + off, rates, width, color=colour, edgecolor="#565b61", lw=0.65,
                      label=lab, zorder=2)
        for i, (bar, (m, _)) in enumerate(zip(bars, models)):
            d = N["flip_loss_specificity"][m]
            n = d[f"{kk}_retained"] + d[f"{kk}_lost"]
            ax.text(bar.get_x() + bar.get_width() / 2, d[f"{kk}_rate"] + 0.004,
                    f"{d[f'{kk}_retained']}/{n}", ha="center", va="bottom", fontsize=6.8,
                    color=C_INK)
            # control-vs-control null band: what retention would be if composition were ZERO
            nb = N["flip_retention_vs_null"][m][kk]
            ax.add_patch(plt.Rectangle((bar.get_x(), nb["null_rate_95"][0]), width,
                                       nb["null_rate_95"][1] - nb["null_rate_95"][0],
                                       facecolor=C_NULL, alpha=0.16, edgecolor="none",
                                       zorder=3))
            ax.plot([bar.get_x(), bar.get_x() + width], [nb["null_rate"]] * 2,
                    color=C_NULL, lw=1.3, zorder=4)
    for i, (m, _) in enumerate(models):
        d = N["flip_loss_specificity"][m]
        # the eligibility line is the point: at k=3 the test has no dynamic range
        ax.text(i, 1.032, f"OR {d['fisher_or']:.2f}, $p$ = {d['fisher_p']:.3f}",
                ha="center", va="bottom", fontsize=7.0, color=C_INK,
                fontweight="bold" if m == "multidim" else "normal")
        ax.text(i, 1.014, f"only {d['eligible_frac']:.0%} of flips reachable"
                if m == "k3" else f"{d['eligible_frac']:.0%} of flips reachable",
                ha="center", va="bottom", fontsize=6.3,
                color=C_NULL if m == "k3" else "#737980")
    ax.set_ylim(0.855, 1.055)
    ax.set_yticks(np.arange(0.86, 1.001, 0.04))
    ax.set_yticklabels([f"{v:.0%}" for v in np.arange(0.86, 1.001, 0.04)])
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in models])
    ax.set_box_aspect(1)
    ax.set_ylabel("sign flips retained after\nremoving all measured composition")
    ax.set_title("Composition erodes reversals and\nstable flips at the same rate", pad=26)
    h, _ = ax.get_legend_handles_labels()
    h.append(Patch(facecolor=C_NULL, alpha=0.16, edgecolor=C_NULL,
                   label="null: composition $\\equiv$ 0 (95%)"))
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              handletextpad=0.5, columnspacing=1.1, fontsize=6.8)
    panel_tag(ax, "f", dx=-0.16, dy=1.08)


# 3 square panels across 11.0in with wspace=0.62 -> each ~2.6in wide, so each row is ~2.6in
# tall; height is set to that plus room for the two-line titles, panel f's legend and the footnote.
fig = plt.figure(figsize=(11.0, 7.1))
grid = fig.add_gridspec(2, 3, wspace=0.62, hspace=0.30)
# row 1 (a-c): the substrate -- the ADT protein space the states are read from
umap_state(fig.add_subplot(grid[0, 0]))
umap_donor(fig.add_subplot(grid[0, 1]))
umap_subset(fig.add_subplot(grid[0, 2]))
# row 2 (d-f): the composition pilot computed on that space
occupancy(fig.add_subplot(grid[1, 0]))
matched_contrast(fig.add_subplot(grid[1, 1]))
flip_loss(fig.add_subplot(grid[1, 2]))
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"FigS3_base.{ext}", dpi=300, bbox_inches="tight",
                facecolor="white")
print(f"wrote {OUT}/FigS3_base.pdf")


# ----------------------------- g-i: response-shape geometry -------------------
# The three panels demoted from Figure 4, kept in full as Supplementary Figure S3g-i.
# Rest-versus-48h reconstruction repeats the geometry of Figure 1b in different units and on
# a 43% subset; physical-guide and leave-one-donor agreement repeat the guide and donor
# robustness Figure 1g,h already carries. They are quality control on the per-cell estimator
# rather than independent evidence, so they belong here, but nothing is discarded.
#
# Geometry, limits and labels below are the v5.1 presentation layer verbatim (originally
# scripts/response_shape/fig4_presentation_v5_1.py). The panels are rendered at the width of
# the Arce pilot above so the two halves composite without either being rescaled.
BUNDLE = ROOT / "data/response_shape/production_v5_1"
MANIFEST = BUNDLE / "figure4_production_manifest_v5_1.json"
DONORS = ("D1", "D2", "D3", "D4")
CONDITIONS = ("Rest", "Stim48hr")
TEXT_PT = 7.0          # the frozen builder's minimum final text size
DISPLAY_QUANTILE = 0.995   # cut the scatter axes here, not at the max, so a sub-1% tail
DISPLAY_PAD = 1.10         # cannot flatten the cloud onto the origin


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(key: str) -> pd.DataFrame:
    """One sealed production table, refused unless it matches its manifest SHA-256."""
    entry = json.loads(MANIFEST.read_text())["artifacts"][key]
    path = BUNDLE / entry["filename"]
    digest = _sha256(path)
    if digest != entry["sha256"]:
        raise ValueError(f"{path.name}: sha256 {digest} does not match the sealed manifest")
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")


unit_effects = _artifact("unit_effects")
edge_summary = _artifact("edge_continuous_summary")
cal_envelopes = _artifact("cal_envelopes")
census = _artifact("reversal_census")


def _reversal_summary() -> pd.DataFrame:
    return edge_summary.loc[edge_summary["analysis_role"].eq("REVERSAL")].copy()


def _envelope(metric: str, condition: str) -> float | None:
    scope = "cross_condition" if condition == "delta" else "condition"
    row = cal_envelopes.set_index(["metric", "scope", "condition"]).loc[(metric, scope, condition)]
    return float(row["absolute_q95"]) if bool(row["available"]) else None


def _canonical_donors(value: object) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and np.isnan(value)) or str(value).strip() == "":
        return ()
    donors = tuple(str(value).split(","))
    if len(donors) != len(set(donors)) or not set(donors).issubset(DONORS) \
            or donors != tuple(donor for donor in DONORS if donor in donors):
        raise ValueError(f"common donor IDs are malformed or out of frozen order: {value!r}")
    return donors


def _stability() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Two-arm guide agreement, and the leave-one-donor recomputation of the same edges."""
    reversal = _reversal_summary()
    supported = reversal.loc[reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")].set_index("analysis_edge_id")
    guides = supported[["oriented_delta__guide_A_guide_total_effect",
                        "oriented_delta__guide_B_guide_total_effect"]].rename(
        columns={"oriented_delta__guide_A_guide_total_effect": "A",
                 "oriented_delta__guide_B_guide_total_effect": "B"}).dropna()
    units = unit_effects.loc[unit_effects["analysis_role"].eq("REVERSAL")
                             & unit_effects["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
    rows = []
    for donor in DONORS:
        estimates: list[tuple[str, float]] = []
        for edge_id, block in units.groupby("analysis_edge_id"):
            if edge_id not in supported.index:
                continue
            frozen_common = _canonical_donors(supported.loc[edge_id, "common_donor_ids"])
            # only four-donor edges are evaluable: dropping one must still leave a balanced three
            if donor not in frozen_common or len(frozen_common) < 4:
                continue
            retained_common = tuple(item for item in frozen_common if item != donor)
            retained = block.loc[block["donor"].isin(retained_common)]
            if any(set(retained.loc[retained["condition"].eq(condition), "donor"]) != set(retained_common)
                   for condition in CONDITIONS):
                continue
            means = retained.groupby("condition")["total_effect"].mean()
            estimates.append((edge_id, (means["Stim48hr"] - means["Rest"])
                              * int(supported.loc[edge_id, "orientation"])))
        if estimates:
            ids = [edge for edge, _ in estimates]
            values = np.asarray([value for _, value in estimates])
            reference = supported.loc[ids, "oriented_delta__total_effect"].to_numpy(float)
            agreement = float(np.mean(np.sign(values) == np.sign(reference)))
            median_abs = float(np.median(np.abs(values - reference)))
        else:
            agreement = np.nan
            median_abs = np.nan
        rows.append({"donor": donor, "n_evaluable": len(estimates),
                     "sign_agreement": agreement, "median_absolute_difference": median_abs})
    return guides, pd.DataFrame(rows)


def display_limit(pairs, quantile=DISPLAY_QUANTILE, pad=DISPLAY_PAD) -> tuple[float, int]:
    """Symmetric limit at a robust quantile, with the count of points outside it.

    `pairs` are the (x, y) series actually drawn. A point counts as off-scale when either
    coordinate leaves the box, which is what the reader sees.
    """
    series: list[tuple[np.ndarray, np.ndarray]] = []
    for x_values, y_values in pairs:
        x = np.asarray(list(x_values), dtype=float)
        y = np.asarray(list(y_values), dtype=float)
        if x.size != y.size:
            raise ValueError("display_limit requires paired x and y series")
        finite = np.isfinite(x) & np.isfinite(y)
        if finite.any():
            series.append((x[finite], y[finite]))
    if not series:
        return 1e-6, 0
    pooled = np.concatenate([np.concatenate((np.abs(x), np.abs(y))) for x, y in series])
    limit = max(1e-6, float(np.quantile(pooled, quantile)) * pad)
    off_scale = sum(int(((np.abs(x) > limit) | (np.abs(y) > limit)).sum()) for x, y in series)
    return limit, off_scale


def note_off_scale(ax, off_scale: int, limit: float, total: int, corner="lower right") -> None:
    """State the clipped tail in a free corner, in DATA coordinates.

    These panels are set_aspect("equal", adjustable="box"), which shrinks the drawn box inside
    the axes rectangle while transAxes still spans the full rectangle -- so an axes-fraction
    annotation at 0.98 lands in the margin outside the plot. Data coordinates stay on the box.
    """
    if off_scale <= 0:
        return
    top = corner.startswith("upper")
    left = corner.endswith("left")
    ax.text(-limit * 0.96 if left else limit * 0.96,
            limit * 0.96 if top else -limit * 0.96,
            f"{off_scale:,} of {total:,} beyond axis",
            ha="left" if left else "right", va="top" if top else "bottom",
            fontsize=TEXT_PT, color="#666B70",
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none"))


def panel_rest_stim(ax):
    reversal = _reversal_summary()
    table = reversal.loc[reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
    census_direction = census.set_index("reversal_edge_id")["direction"]
    gained = table["analysis_edge_id"].map(census_direction).eq("gained_on_activation")
    for mask, colour, marker, label in ((gained, "#2A9D8F", "o", "gained"),
                                        (~gained, "#C1272D", "^", "lost")):
        block = table.loc[mask]
        ax.scatter(block["Rest__total_effect"], block["Stim48hr__total_effect"], s=7, alpha=0.30,
                   color=colour, marker=marker, edgecolors="none", rasterized=True,
                   label=f"{label} (n={len(block):,})")
    limit, off_scale = display_limit([(table["Rest__total_effect"], table["Stim48hr__total_effect"])])
    ax.plot([-limit, limit], [-limit, limit], color="#C5C9CD", lw=0.8)
    ax.axhline(0, color="#9AA0A6", lw=0.7)
    ax.axvline(0, color="#9AA0A6", lw=0.7)
    for condition, orientation in (("Rest", "vertical"), ("Stim48hr", "horizontal")):
        threshold = _envelope("total_effect", condition)
        if threshold is not None:
            for sign in (-1, 1):
                (ax.axvline if orientation == "vertical" else ax.axhline)(
                    sign * threshold, color="#767D84", lw=0.7, ls="--")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    clipped = f"\n{off_scale} of {len(table):,} beyond axis" if off_scale else ""
    ax.set_xlabel(f"Rest donor-mean total effect{clipped}")
    ax.set_ylabel("48 h donor-mean total effect")
    ax.set_title("Rest vs 48 h total effect", loc="left", pad=7)
    # Direct labels rather than a legend box: gained clusters upper-left and lost lower-right, so
    # both diagonal corners carry data and any framed legend would sit on points.
    counts = {"gained": int(gained.sum()), "lost": int((~gained).sum())}
    ax.text(-limit * 0.92, limit * 0.86, f"gained (n={counts['gained']:,})", color="#2A9D8F",
            fontsize=TEXT_PT, fontweight="bold", ha="left", va="center")
    ax.text(limit * 0.92, -limit * 0.86, f"lost (n={counts['lost']:,})", color="#C1272D",
            fontsize=TEXT_PT, fontweight="bold", ha="right", va="center")


def panel_stability(fig, spec):
    grid = spec.subgridspec(1, 2, wspace=0.42)
    ax_guide = fig.add_subplot(grid[0, 0])
    ax_lodo = fig.add_subplot(grid[0, 1])
    guides, lodo = _stability()

    ax_guide.scatter(guides["A"], guides["B"], s=7, alpha=0.28, color="#457B9D",
                     edgecolors="none", rasterized=True)
    limit, off_scale = display_limit([(guides["A"], guides["B"])])
    ax_guide.plot([-limit, limit], [-limit, limit], color="#B8BDC2", lw=0.8)
    ax_guide.axhline(0, color="#9AA0A6", lw=0.6)
    ax_guide.axvline(0, color="#9AA0A6", lw=0.6)
    ax_guide.set(xlim=(-limit, limit), ylim=(-limit, limit),
                 xlabel="guide A oriented change", ylabel="guide B oriented change")
    ax_guide.set_aspect("equal", adjustable="box")
    agreement = float(np.mean(np.sign(guides["A"]) == np.sign(guides["B"]))) if len(guides) else np.nan
    ax_guide.set_title("Physical-guide agreement", loc="left", pad=7)
    ax_guide.text(-limit * 0.96, limit * 0.96, f"n={len(guides):,}; sign={agreement:.1%}",
                  ha="left", va="top", fontsize=TEXT_PT, color="#565B61")
    note_off_scale(ax_guide, off_scale, limit, len(guides), corner="lower right")

    x = np.arange(len(lodo))
    ax_lodo.bar(x, lodo["sign_agreement"], color="#6A8EAE", width=0.62)
    # Donor identity alone on the ticks, with the operation named once as the axis label.
    ax_lodo.set_xticks(x, list(lodo["donor"]))
    # The baseline stays at zero: these are proportions, and truncating the axis would blow up
    # the single point separating 99.8% from 100.0% into a visible step.
    ax_lodo.set_ylim(0, 1.10)
    ax_lodo.set_ylabel("direction agreement")
    evaluable = {int(row.n_evaluable) for row in lodo.itertuples(index=False)}
    shared_evaluable = evaluable.pop() if len(evaluable) == 1 else None
    for xx, row in enumerate(lodo.itertuples(index=False)):
        if row.n_evaluable == 0:
            label = "n=0"
        elif shared_evaluable is not None:
            label = f"{row.sign_agreement:.1%}"
        else:
            label = f"{row.sign_agreement:.1%}\nn={row.n_evaluable:,}"
        ax_lodo.text(xx, 0.02 if row.n_evaluable == 0 else row.sign_agreement + 0.02, label,
                     ha="center", va="bottom", fontsize=TEXT_PT)
    # The shared n rides on the axis label rather than floating inside a bar, where dark fill
    # left it unreadable.
    ax_lodo.set_xlabel(
        f"donor left out (n = {shared_evaluable:,} four-donor edges)"
        if shared_evaluable is not None else "donor left out"
    )
    ax_lodo.set_title("Leave-one-donor agreement", loc="left", pad=7)
    return [ax_guide, ax_lodo]


setup(1.0)
# One row of three, sized to the width of the Arce pilot render above, so the two halves of
# Supplementary Figure S3 share a page width without rescaling either.
sfig = plt.figure(figsize=(10.24, 3.55), facecolor="white")
sgrid = sfig.add_gridspec(1, 3, left=0.070, right=0.985, top=0.790, bottom=0.175, wspace=0.40)
ax_g = sfig.add_subplot(sgrid[0, 0])
panel_rest_stim(ax_g)
stability_axes = panel_stability(sfig, sgrid[0, 1:])
panel_tag(ax_g, "g", dx=-0.055, dy=1.09)
panel_tag(stability_axes[0], "h", dx=-0.16, dy=1.09)
panel_tag(stability_axes[1], "i", dx=-0.16, dy=1.09)
sfig.suptitle("Response-shape reconstruction geometry and estimator robustness",
              x=0.070, y=0.975, ha="left", fontsize=11.0, fontweight="bold", color="#22262B")
sfig.savefig(OUT / "Fig4_supp.pdf", dpi=600, facecolor="white")
sfig.savefig(OUT / "Fig4_supp.png", dpi=300, facecolor="white")
plt.close(sfig)
print(f"wrote {OUT}/Fig4_supp.pdf")


# ----------------------------- compose a-f with g-i ---------------------------
def compose() -> Path:
    top = PdfReader(OUT / 'FigS3_base.pdf').pages[0]
    bottom = PdfReader(OUT / 'Fig4_supp.pdf').pages[0]
    tw, th = float(top.mediabox.width), float(top.mediabox.height)
    bw, bh = float(bottom.mediabox.width), float(bottom.mediabox.height)
    width = max(tw, bw)
    height = th + bh + GAP
    page = PageObject.create_blank_page(width=width, height=height)
    # sources are centred rather than scaled: scaling one to the other's width would
    # change its font size relative to the other half of the same figure.
    page.merge_transformed_page(bottom, Transformation().translate((width - bw) / 2, 0))
    page.merge_transformed_page(top, Transformation().translate((width - tw) / 2, bh + GAP))
    writer = PdfWriter()
    writer.add_page(page)
    out = OUT / 'FigS3.pdf'
    with out.open('wb') as handle:
        writer.write(handle)
    subprocess.run(['pdftoppm', '-png', '-r', '300', '-singlefile', str(out),
                    str(OUT / 'FigS3')], check=True)
    (OUT / 'FigS3_base.pdf').unlink(missing_ok=True)
    (OUT / 'FigS3_base.png').unlink(missing_ok=True)
    (OUT / 'Fig4_supp.pdf').unlink(missing_ok=True)
    (OUT / 'Fig4_supp.png').unlink(missing_ok=True)
    return out


if __name__ == '__main__':
    path = compose()
    reader = PdfReader(path)
    box = reader.pages[0].mediabox
    print(f'wrote {path}  ({float(box.width):.1f} x {float(box.height):.1f} pt)')
