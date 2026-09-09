#!/usr/bin/env python3
"""Supplementary Figure S2: separate-dataset recurrence and the composition falsification."""

import os
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# ----------------- STYLE CONSTANTS -----------------
PAL = dict(
    REST="#2c6fbb", STIM="#e07b39",
    LOST="#c1272d",     # up in Rest -> down in Stim48  (edge lost on stimulation)
    GAINED="#2a9d8f",   # down in Rest -> up in Stim48  (edge gained on stimulation)
    GREY="#9aa0a6", CONC="#c7ccd1", CORE="#6a51a3",
    INK="#22262b", HAIR="#dfe2e6", HL="#f2c14e",
)
A4_PORTRAIT = (8.27, 11.69)
SOURCE_WIDTH_IN = A4_PORTRAIT[0]
SOURCE_HEIGHT_IN = A4_PORTRAIT[1] * 0.60
FONT_FLOOR_PT = 6.70
PRINT_WIDTHS_MM = (165.1, 180.0)

def setup(scale=1.0):
    def s(v): return round(v * scale, 2)
    matplotlib.rcParams.update({
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

# ----------------- ARCE STATS FUNCTIONS -----------------

def boot_rate(df, col_kind, n_boot=10000, seed=3):
    rng = np.random.default_rng(seed)
    regs = df.regulator.unique()
    by = {r: g[g.our_kind == col_kind].arce_flipped.values for r, g in df.groupby("regulator")}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(regs, size=len(regs), replace=True)
        v = np.concatenate([by[r] for r in pick if len(by[r])])
        if len(v): out.append(v.mean())
    return np.percentile(out, [2.5, 97.5])

def strata_tables(df, reg_col="regulator", kind_col="our_kind", flip_col="arce_flipped"):
    out = []
    for reg, d in df.groupby(reg_col):
        r, s = d[d[kind_col] == "rev"], d[d[kind_col] == "stable"]
        a = int(r[flip_col].sum()); b = len(r) - a
        c = int(s[flip_col].sum()); dd = len(s) - c
        out.append((reg, a, b, c, dd))
    return out

def mh_odds_ratio(tabs):
    num = den = 0.0
    for _, a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0: continue
        num += a * d / n
        den += b * c / n
    if den == 0: return np.inf if num > 0 else np.nan
    return num / den

def cluster_bootstrap_or(df, n_boot=10000, seed=0, reg_col="regulator", flip_col="arce_flipped"):
    rng = np.random.default_rng(seed)
    regs = df[reg_col].unique()
    by = {r: d for r, d in df.groupby(reg_col)}
    ors = []
    for _ in range(n_boot):
        pick = rng.choice(regs, size=len(regs), replace=True)
        tabs = []
        for i, r in enumerate(pick):
            d = by[r]
            rv, st = d[d.our_kind == "rev"], d[d.our_kind == "stable"]
            a = int(rv[flip_col].sum()); b = len(rv) - a
            c = int(st[flip_col].sum()); dd = len(st) - c
            tabs.append((f"{r}_{i}", a, b, c, dd))
        o = mh_odds_ratio(tabs)
        if np.isfinite(o): ors.append(o)
    return np.array(ors)

def per_regulator_rates(df, reg_col="regulator", flip_col="arce_flipped"):
    rows = []
    for reg, d in df.groupby(reg_col):
        r, s = d[d.our_kind == "rev"], d[d.our_kind == "stable"]
        a = int(r[flip_col].sum()); b = len(r) - a
        c = int(s[flip_col].sum()); dd = len(s) - c
        or_k = ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({
            "regulator": reg, "n_rev": len(r), "n_stable": len(s),
            "rev_flip": a, "rev_flip_rate": a / len(r) if len(r) else np.nan,
            "stable_flip": c, "stable_flip_rate": c / len(s) if len(s) else np.nan,
            "or_haldane": or_k,
        })
    return rows

def leave_one_out(df, reg_col="regulator", flip_col="arce_flipped"):
    rows = []
    for r in sorted(df[reg_col].unique()):
        sub = df[df[reg_col] != r]
        tabs = strata_tables(sub, reg_col=reg_col, flip_col=flip_col)
        o = mh_odds_ratio(tabs)
        rows.append({"dropped": r, "or_mh": o})
    return rows

# ----------------- DATA LOADING -----------------

# Atlas Composition Data
AC_d = pd.read_csv(ROOT / "data/atlas/atlas_edge_composition.csv.gz")
AC_rev = AC_d[AC_d.kind.eq("rev")].dropna(subset=["comp_real", "comp_zi", "minmag"]).copy()
AC_n_rev = len(AC_rev)
AC_n_flip = int((AC_rev.comp_real > AC_rev.minmag).sum())
AC_med = dict(floor=AC_rev.comp_zi.median(), real=AC_rev.comp_real.median(), end=AC_rev.minmag.median())
C_REAL, C_FLOOR, C_END = PAL["CORE"], PAL["GREY"], PAL["LOST"]
C_INK, C_HAIR = PAL["INK"], PAL["HAIR"]

# Arce Data
AR_E = pd.read_csv(ROOT / "data/arce/arce_validation_edges.csv.gz")
AR_PS = AR_E[(AR_E.arce_dataset == "perturbseq")].copy()
AR_BK = AR_E[(AR_E.arce_dataset == "bulk")].copy()
with open(ROOT / "data/arce/arce_validation_results.json") as f:
    AR_R = json.load(f)
AR_rev, AR_stab = AR_PS[AR_PS.our_kind == "rev"], AR_PS[AR_PS.our_kind == "stable"]
AR_rev_rate, AR_stab_rate = AR_rev.arce_flipped.mean(), AR_stab.arce_flipped.mean()
AR_tabs = strata_tables(AR_PS)
AR_OR_MH = mh_odds_ratio(AR_tabs)
AR_boot = cluster_bootstrap_or(AR_PS, n_boot=10000, seed=1)
AR_OR_LO, AR_OR_HI = np.percentile(AR_boot, [2.5, 97.5])
AR_pr = pd.DataFrame(per_regulator_rates(AR_PS)).sort_values("n_rev")
AR_loo = pd.DataFrame(leave_one_out(AR_PS))
AR_LOO_LO, AR_LOO_HI = AR_loo.or_mh.min(), AR_loo.or_mh.max()

# Target Axis Data
with open(ROOT / "data/census/regulator_target_structure.json") as f:
    TA_STRUCT = json.load(f)
with open(ROOT / "data/census/target_theory_falsification.json") as f:
    TA_FALS = json.load(f)
C_OBS = "#33475b"

# ----------------- PANEL FUNCTIONS -----------------
def panel_a(ax, tag="a"):
    xs = {"stable": 0.0, "rev": 1.0}
    C_REV, C_STAB = PAL["CORE"], PAL["CONC"]
    rev_ci, stab_ci = boot_rate(AR_PS, "rev"), boot_rate(AR_PS, "stable")

    # per-regulator paired lines, point size ~ number of reversal edges
    for _, r in AR_pr.iterrows():
        sz = 12 + 46 * np.sqrt(r.n_rev / AR_pr.n_rev.max())
        up = r.rev_flip_rate > r.stable_flip_rate
        ax.plot([xs["stable"], xs["rev"]], [r.stable_flip_rate, r.rev_flip_rate],
                "-", color=(PAL["GREY"] if up else PAL["LOST"]),
                lw=0.9, alpha=0.85 if up else 1.0, zorder=2 if up else 3)
        ax.scatter([xs["stable"]], [r.stable_flip_rate], s=sz, color=C_STAB,
                   edgecolor="#565b61", lw=0.6, zorder=4)
        ax.scatter([xs["rev"]], [r.rev_flip_rate], s=sz,
                   color=(C_REV if up else PAL["LOST"]),
                   edgecolor="#565b61", lw=0.6, zorder=4)

    # pooled estimates with CLUSTER-bootstrap CIs
    for xv, rate, ci, col in [(xs["stable"], AR_stab_rate, stab_ci, "#3f4448"),
                              (xs["rev"], AR_rev_rate, rev_ci, C_REV)]:
        ax.errorbar([xv], [rate], yerr=[[rate - ci[0]], [ci[1] - rate]], fmt="_",
                    ms=30, mew=2.6, color=col, ecolor=col, elinewidth=2.0,
                    capsize=5, capthick=2.0, zorder=6)

    ax.annotate(f"{AR_rev_rate:.1%}", (xs["rev"], AR_rev_rate), xytext=(22, -12),
                textcoords="offset points", fontsize=9.5, fontweight="bold",
                color=C_REV, va="center")
    ax.annotate(f"{AR_stab_rate:.1%}", (xs["stable"], AR_stab_rate), xytext=(-24, 0),
                textcoords="offset points", fontsize=9.5, fontweight="bold",
                color="#3f4448", va="center", ha="right")
    ax.set_xlim(-0.62, 2.28); ax.set_ylim(-0.04, 1.08)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"stable edges\n(n = {len(AR_stab):,})", f"reversing edges\n(n = {len(AR_rev)})"])
    ax.set_yticks(np.arange(0, 1.01, 0.25))
    ax.set_yticklabels([f"{int(v*100)}%" for v in np.arange(0, 1.01, 0.25)])
    ax.set_ylabel("edges that re-reverse in Arce et al.")
    ax.set_title("Reversals recur; the same regulators'\nstable edges do not", pad=10)
    ax.legend(handles=[
        Line2D([], [], color=PAL["GREY"], lw=1.0, label="regulator stratum (11)"),
        Line2D([], [], color=PAL["LOST"], lw=1.0, label="against the effect (2)"),
        Line2D([], [], color=PAL["CORE"], marker="_", ls="none", ms=13, mew=2.4, label="pooled, 95% cluster boot."),
    ], loc="lower right", bbox_to_anchor=(1.03, 0.07), handletextpad=0.7, borderaxespad=0.0)
    panel_tag(ax, tag, dx=-0.13, dy=1.12)

def panel_concordance(ax, tag="b"):
    b2 = TA_STRUCT["B_target_determines_direction"]["B2_sign_concordance_all_edges"]
    pooled_obs, pooled_null = b2["observed"], b2["null_mean"]
    def excesses(key):
        return [v["excess"] for k, v in TA_FALS[key].items() if isinstance(v, dict) and not v.get("insufficient") and not k.startswith("_")]
    bm = excesses("A4_B2_within_baseMean_decile")
    idg = excesses("A5_B2_within_indegree_decile")
    pooled_excess = pooled_obs - pooled_null
    ax.axhline(0, color="#565b61", lw=0.9)
    ax.axhline(pooled_excess, ls="--", lw=1.0, color=C_OBS)
    ax.text(2.45, pooled_excess + 0.004, f"pooled +{pooled_excess:.3f}", fontsize=6.2, color=C_OBS, va="bottom", ha="center")
    for x, vals, lab, c in [(1, bm, "baseline\nexpression", "#2c6fbb"), (2, idg, "target\nin-degree", "#3f8f7d")]:
        jitter = (np.arange(len(vals)) - (len(vals) - 1) / 2) * 0.035
        ax.scatter(np.full(len(vals), x) + jitter, vals, s=22, color=c, zorder=3, edgecolor="white", lw=0.4)
        ax.text(x, -0.006, f"{len(vals)}/{len(vals)} deciles\nexcess > 0", ha="center", fontsize=6.0, color=c, va="top")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["baseline expression", "target in-degree"], fontsize=7)
    ax.set_xlim(0.5, 2.9); ax.set_ylim(-0.02, max(max(bm), max(idg)) * 1.25)
    ax.set_ylabel("excess same-sign edge pairs\nover the within-regulator null")
    ax.set_title(f"Edges on one target share sign:\n{pooled_obs:.1%} against {pooled_null:.1%};"
                 " excess remains positive\nacross both target stratifications", fontsize=8)
    panel_tag(ax, tag)

def panel_ladder(ax, tag="c"):
    data = [AC_rev.comp_zi.clip(1e-4), AC_rev.comp_real.clip(1e-4), AC_rev.minmag.clip(1e-4)]
    labels = ["random-state\nfloor", "real\ncomposition", "reversal\nendpoint"]
    colors = [C_FLOOR, C_REAL, C_END]
    parts = ax.violinplot([np.log10(x) for x in data], positions=[0, 1, 2], showmedians=True, widths=0.8)
    for pc, c in zip(parts["bodies"], colors):
        pc.set_facecolor(c); pc.set_alpha(0.55); pc.set_edgecolor("#565b61"); pc.set_linewidth(0.6)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color(C_INK); parts[key].set_linewidth(1.0)
    for i, (m, c) in enumerate(zip([AC_med["floor"], AC_med["real"], AC_med["end"]], colors)):
        ax.text(i, np.log10(m) + 0.12, f"{m:.3f}", ha="center", fontsize=7, fontweight="bold", color=C_INK)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(labels)
    yt = np.array([0.001, 0.01, 0.1, 1.0])
    ax.set_yticks(np.log10(yt)); ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.set_ylabel("|log2 fold change|  (log scale)")
    ax.set_title("Composition is real but sits far below\nthe endpoint it would need to flip", pad=8)
    ax.annotate("", xy=(1, np.log10(AC_med["end"])), xytext=(1, np.log10(AC_med["real"])),
                arrowprops=dict(arrowstyle="<->", color="#737980", lw=1.0))
    # Place the ratio at the geometric midpoint on the log scale.  Summing the
    # two log coordinates puts the label far below the interval it annotates.
    midpoint = np.log10(np.sqrt(AC_med["real"] * AC_med["end"]))
    ax.text(1.06, midpoint, f"{AC_med['end']/AC_med['real']:.0f}x",
            fontsize=7.5, color="#737980", va="center", fontweight="bold")
    panel_tag(ax, tag, dx=-0.17, dy=1.08)

def panel_scatter(ax, tag="d"):
    x = AC_rev.minmag.clip(1e-3).to_numpy()
    y = AC_rev.comp_real.clip(1e-4).to_numpy()
    ax.scatter(x, y, s=7, color=C_REAL, alpha=0.35, edgecolors="none", rasterized=True, zorder=3)
    lim_lo, lim_hi = 3e-3, 3.0
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color=C_END, lw=1.4, ls="--", zorder=4)
    ax.text(0.5, 0.9, "flip threshold (y = x):\ncomposition would have to\nreach the endpoint",
            transform=ax.transAxes, fontsize=6.6, color=C_END, ha="left", va="top")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim_lo, lim_hi); ax.set_ylim(1e-4, lim_hi)
    ax.set_xlabel("reversal endpoint  |log2FC| (min of the two)")
    ax.set_ylabel("|composition| at that endpoint")
    ax.set_title(f"None of {AC_n_rev:,} reversals could be\nflipped by composition ({AC_n_flip} above the line)", pad=8)
    panel_tag(ax, tag, dx=-0.17, dy=1.08)

def panel_c(ax, tag="e"):
    groups = []
    for ds, d, lab in [("perturbseq", AR_PS, "Perturb-CITE-seq\nCRISPRi, 2 donors"),
                       ("bulk", AR_BK, "bulk RNA-seq\nCRISPR KO, 3 donors")]:
        m12 = d[d.regulator == "MED12"]
        groups.append((lab,
                       m12[m12.our_kind == "stable"].arce_flipped.mean(),
                       m12[m12.our_kind == "rev"].arce_flipped.mean(),
                       int((m12.our_kind == "rev").sum()),
                       int((m12.our_kind == "stable").sum()),
                       AR_R[ds]["med12"]["fisher_or"]))
    x = np.arange(len(groups))
    w = 0.33
    ax.bar(x - w / 2, [g[1] for g in groups], w, color=PAL["CONC"], edgecolor="#565b61", lw=0.7, label="MED12 stable edges")
    ax.bar(x + w / 2, [g[2] for g in groups], w, color=PAL["CORE"], edgecolor="#565b61", lw=0.7, label="MED12 reversing edges")
    for xi, g in zip(x, groups):
        ax.text(xi - w / 2, g[1] + 0.016, f"{g[1]:.1%}", ha="center", fontsize=7.4, color="#3f4448")
        ax.text(xi + w / 2, g[2] + 0.016, f"{g[2]:.1%}", ha="center", fontsize=7.4, fontweight="bold", color=PAL["CORE"])
        ax.text(xi, 0.523, f"OR {g[5]:.2f}", ha="center", fontsize=7.6, fontweight="bold", color=PAL["INK"])
        ax.text(xi, 0.492, f"{g[3]} rev / {g[4]:,} stable", ha="center", fontsize=6.4, color="#8a9096")
    ax.set_xticks(x); ax.set_xticklabels([g[0] for g in groups], fontsize=7.6)
    ax.set_ylim(0, 0.80); ax.set_yticks(np.arange(0, 0.61, 0.2))
    ax.set_yticklabels([f"{int(v*100)}%" for v in np.arange(0, 0.61, 0.2)])
    ax.set_ylabel("edges that re-reverse in Arce et al.")
    ax.set_title("MED12: a baseline reversal rate (Figure 2b)\nand real reversals are both true", pad=10)
    ax.legend(loc="upper center", ncol=2, handletextpad=0.6, columnspacing=1.4, bbox_to_anchor=(0.5, 1.005))
    panel_tag(ax, tag, dx=-0.13, dy=1.12)

def panel_d(ax, tag="f"):
    cb = AR_R["cbfb_il2ra"]
    ps_cb = cb["perturbseq"]
    sc = cb["screen_ko_surface_protein"]
    flip = sc["sign_relative_to_knockdown_effect"]
    series = [
        ("our atlas\nCRISPRi, transcript", ps_cb["our_rest"], ps_cb["our_stim"], PAL["CORE"], "o"),
        ("Arce Perturb-seq\nCRISPRi, transcript", ps_cb["arce_rest"], ps_cb["arce_stim"], PAL["REST"], "s"),
        ("Arce screen\nCRISPR KO, surface protein", flip * sc["lfc_rest_lowhigh"], flip * sc["lfc_stim_lowhigh"], PAL["STIM"], "^"),
    ]
    for (lab, a, b, col, mk), dx in zip(series, (-0.035, 0.0, 0.035)):
        ax.plot([dx, 1 + dx], [a, b], "-", color=col, lw=1.6, marker=mk, ms=7, mec="#565b61", mew=0.6, label=lab, zorder=3)
    ax.axhline(0, ls="-", lw=1.0, color="#565b61", zorder=1)
    ax.set_xlim(-0.55, 2.05); ax.set_ylim(-1.35, 3.55)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Rest", "Stim 48 h"])
    ax.set_ylabel("effect of losing CBFB on IL2RA\n(log2 fold change)")
    ax.set_title("CBFB $\\rightarrow$ IL2RA: the one donor-reproduced\nreversal Arce covers", pad=10)
    ax.legend(loc="upper right", bbox_to_anchor=(1.03, 1.02), handletextpad=0.7, labelspacing=0.85, borderaxespad=0.0)
    ax.text(-0.52, 1.55, "raises\nIL2RA", fontsize=6.6, color="#565b61", style="italic", va="center")
    ax.text(-0.52, -0.75, "lowers\nIL2RA", fontsize=6.6, color="#565b61", style="italic", va="center")
    panel_tag(ax, tag, dx=-0.16, dy=1.12)

# ----------------- RETUNE & BUILD -----------------
def _all_axis_text(ax):
    out = [ax.title, ax.xaxis.label, ax.yaxis.label]
    out += list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.texts)
    legend = ax.get_legend()
    if legend is not None:
        out += list(legend.get_texts())
        if legend.get_title() is not None: out.append(legend.get_title())
    return [t for t in out if t is not None and t.get_text()]

def _enforce_font_floor(axes, floor=FONT_FLOOR_PT):
    for ax in axes:
        for artist in _all_axis_text(ax):
            if artist.get_fontsize() < floor:
                artist.set_fontsize(floor)

def _retune(ax_a, ax_b, ax_c, ax_d, ax_e, ax_f):
    ax_a.set_title("Reversals recur;\nstable edges do not")
    ax_c.set_title("Composition is real\nbut far too small")
    ax_d.set_title("No reversal flipped\nby composition")
    ax_e.set_title("MED12: baseline rate\nand real reversals")
    ax_f.set_title("CBFB$\\rightarrow$IL2RA\nin three assays")
    ax_c.set_xticklabels(["random\nstate", "real", "reversal\nendpoint"])
    ax_d.set_xlabel("reversal endpoint |log2FC|\n(min of the two)")
    for ax, ncol in ((ax_e, 2), (ax_f, 1)):
        legend = ax.get_legend()
        if legend is None: continue
        handles = list(legend.legend_handles)
        labels = [txt.get_text() for txt in legend.get_texts()]
        legend.remove()
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.32),
                  ncol=ncol, frameon=False, borderaxespad=0.0, fontsize=5.8,
                  labelspacing=0.16, handletextpad=0.50, columnspacing=0.9)
    ax_e.set_xticklabels(["Perturb-CITE-seq\nCRISPRi\n2 donors", "bulk RNA-seq\nCRISPR KO\n3 donors"])
    for artist in ax_e.texts:
        text = artist.get_text()
        if "rev /" in text:
            x, _ = artist.get_position()
            artist.set_text(text.replace(" / ", "\n"))
            artist.set_position((0.26 if x < 0.5 else 0.76, 0.905))
            artist.set_transform(ax_e.transAxes)
            artist.set_ha("center"); artist.set_va("top"); artist.set_fontsize(5.8)
        elif text.startswith("OR "):
            x, _ = artist.get_position()
            artist.set_position((0.26 if x < 0.5 else 0.76, 1.005))
            artist.set_transform(ax_e.transAxes)
            artist.set_ha("center"); artist.set_va("top")
    ax_a.set_xticklabels([f"stable\nn = {len(AR_stab):,}", f"reversing\nn = {len(AR_rev):,}"])
    for artist in ax_a.texts:
        if artist.get_text() == "22.0%":
            artist.set_position((10, 14)); artist.set_ha("left"); artist.set_va("bottom")
    legend = ax_a.get_legend()
    if legend is not None:
        handles = list(legend.legend_handles)
        legend.remove()
        ax_a.legend(handles, ["regulator strata (11)", "against effect (2)", "pooled, 95% cluster bootstrap"],
                    loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=1,
                    borderaxespad=0.0, labelspacing=0.16, handletextpad=0.50, fontsize=5.6)
    ax_b.set_title("Edges on one target share sign\n65.0% observed vs 55.0% null\npositive excess in both strata")
    ax_b.set_xticklabels(["baseline expr.", "in-degree"])
    ax_b.yaxis.labelpad = 6.0
    for artist in ax_b.texts:
        # The tag is drawn as "b" by panel_concordance; this branch used to test
        # for "c" (a leftover from an earlier panel order), so the reposition
        # never fired and the tag sat on the last line of the 3-line title.
        # This title is three lines (one taller than panels a and c), so the tag
        # is top-aligned with the first title line and set left of the title
        # block, where the y-axis label does not reach.
        if artist.get_text() == "b":
            artist.set_position((-0.30, 1.32))
        elif "deciles" in artist.get_text():
            x, _ = artist.get_position()
            artist.set_position((x, 0.002)); artist.set_va("bottom")
        elif artist.get_text().startswith("pooled +"):
            artist.set_position((2.84, 0.126))
            artist.set_ha("right")
            artist.set_bbox(dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.86))

def build():
    setup()
    fig = plt.figure(figsize=(SOURCE_WIDTH_IN, SOURCE_HEIGHT_IN), facecolor="white")
    grid = fig.add_gridspec(2, 3, wspace=0.62, hspace=0.70, left=0.118, right=0.980, top=0.910, bottom=0.205)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[1, 0])
    ax_e = fig.add_subplot(grid[1, 1])
    ax_f = fig.add_subplot(grid[1, 2])
    panel_a(ax_a, "a")
    panel_concordance(ax_b, "b")
    panel_ladder(ax_c, "c")
    panel_scatter(ax_d, "d")
    panel_c(ax_e, "e")
    panel_d(ax_f, "f")
    axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]
    for ax in axes:
        ax.set_box_aspect(1.0)
        ax.set_anchor("C")
    _retune(*axes)
    _enforce_font_floor(axes)
    fig.canvas.draw()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"FigS2.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote FigS2 (six panels; composition n={AC_n_rev:,})")

if __name__ == "__main__":
    build()
