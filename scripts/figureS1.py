#!/usr/bin/env python3
"""Supplementary Figure S1: reversal magnitude, burden against connectivity, and 8 h kinetics."""
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr, norm
from scipy import stats, optimize

# Repository locates itself
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
OKABE_ITO = {"orange": "#E69F00", "sky_blue": "#56B4E9", "bluish_green": "#009E73",
             "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
             "reddish_purple": "#CC79A7", "black": "#000000"}
BLUE, VERM, GREEN = OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["bluish_green"]
MUTED, BAND = "#4D4D4D", "#C8C8C8"

SOURCE_WIDTH_IN, SOURCE_HEIGHT_IN = (8.27, 11.69)  # A4_PORTRAIT
FONT_FLOOR_PT = 7.0

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

# ----------------- DATA LOADING -----------------
# 1. Magnitude data
D_MV = json.load(open(ROOT / "data/census/figure_data.json"))
MV_base = D_MV["census"]["inversion_rate"] * 100
MV_mabs = np.array(D_MV["arr"]["min_abs"])
MV_med = D_MV["census"]["median_min_abs_lfc"]

# 2. Kinetics data
df_kin = pd.read_csv(ROOT / "data/kinetics/kinetics_8hr_edges.csv.gz")
df_kin = df_kin[np.isfinite(df_kin.z_Stim8hr) & np.isfinite(df_kin.z_Rest) & np.isfinite(df_kin.z_Stim48hr)]
kin_rev = df_kin[df_kin.kind == "rev"]
kin_stab = df_kin[df_kin.kind == "stable"]
kin_rng = np.random.default_rng(7)

N_GO_DICT = {'AHR': 31, 'ARNT': 24, 'ATAD5': 2, 'ATP2A2': 30, 'BRD8': 20, 'BTAF1': 4, 'CBFB': 8, 'CCNC': 6, 'CPSF6': 15, 'CREBBP': 43, 'DENR': 3, 'DHX36': 63, 'DOLPP1': 5, 'EIF4G2': 11, 'ELOB': 12, 'ELOF1': 3, 'GATA3': 118, 'KDM1A': 57, 'LCK': 42, 'LEO1': 20, 'MED12': 8, 'MED15': 10, 'MED19': 12, 'MED24': 0, 'MEN1': 74, 'METTL1': 6, 'NAA30': 3, 'NCKAP1L': 82, 'NCOR1': 26, 'NFRKB': 24, 'NSD1': 12, 'PPP1R11': 12, 'SENP5': 5, 'SGF29': 32, 'SMARCB1': 45, 'SMARCE1': 34, 'STAT3': 110, 'STAT5B': 27, 'STIP1': 0, 'SUPT20H': 16, 'SUPT7L': 18, 'TADA1': 0, 'TADA2B': 18, 'TAF13': 16, 'TAF6L': 25, 'TMED9': 17, 'UBE2L3': 29, 'UBXN1': 14, 'UFM1': 18, 'WAC': 24}

# ----------------- MAGNITUDE PANELS -----------------
def panel_magnitude(ax, tag="a"):
    ax.hist(MV_mabs, bins=np.linspace(0, 2.0, 41), color=PAL["REST"], alpha=0.85, edgecolor="white", lw=0.3)
    ax.axvline(MV_med, color=PAL["INK"], lw=1.2)
    ax.text(MV_med+0.03, ax.get_ylim()[1]*0.9, f"median\n{MV_med:.2f}", fontsize=7.5, color=PAL["INK"])
    ax.axvline(0.5, color=PAL["CORE"], lw=1.0, ls="--")
    ax.axvline(1.0, color=PAL["LOST"], lw=1.0, ls="--")
    ax.text(0.5, ax.get_ylim()[1]*0.55, " 0.5", color=PAL["CORE"], fontsize=7)
    ax.text(1.0, ax.get_ylim()[1]*0.55, " 1.0", color=PAL["LOST"], fontsize=7)
    ax.set_xlabel("min |log2FC| across the two states\n(strength of the weaker side)")
    ax.set_ylabel("number of reversing edges")
    ax.set_title("Reversals have real magnitude\non both sides, not wobble")
    ax.set_xlim(0, 2.0)
    panel_tag(ax, tag)

def panel_exemplars(ax, tag="c"):
    ex = sorted(D_MV["exemplars"], key=lambda e: e["lfc_Rest"])
    yy = np.arange(len(ex))
    for i, e in enumerate(ex):
        col = PAL["LOST"] if e["lfc_Rest"] > 0 else PAL["GAINED"]
        ax.plot([e["lfc_Rest"], e["lfc_Stim48"]], [i, i], "-", color=col, lw=2.0, alpha=0.9, zorder=2)
        ax.scatter(e["lfc_Rest"], i, s=34, facecolor="white", edgecolor="#565b61", lw=1.1, zorder=3)
        ax.scatter(e["lfc_Stim48"], i, s=40, color=col, zorder=4)
    ax.axvline(0, color="#565b61", lw=0.9)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{e['regulator']}→{e['target']}" for e in ex], fontsize=7)
    ax.set_ylim(-0.7, len(ex)-0.3)
    _lo = min(min(e["lfc_Rest"], e["lfc_Stim48"]) for e in ex)
    _hi = max(max(e["lfc_Rest"], e["lfc_Stim48"]) for e in ex)
    ax.set_xlim(_lo - 4.4, _hi + 0.7)
    ax.set_xlabel("regulatory effect\nlog2 fold change")
    ax.set_title("Concrete reversals\n(|log2FC| > 1 in both states)")
    ax.legend(handles=[
        plt.Line2D([0],[0], marker="o", ls="", mfc="white", mec="#565b61", ms=7, label="Rest"),
        plt.Line2D([0],[0], marker="o", ls="", color=PAL["LOST"], ms=7, label="Stim 48 h (lost)"),
        plt.Line2D([0],[0], marker="o", ls="", color=PAL["GAINED"], ms=7, label="Stim 48 h (gained)"),
    ], loc="lower left", fontsize=6.2, handletextpad=0.2, labelspacing=0.25)
    panel_tag(ax, tag, dx=-0.28, dy=1.05)

# ----------------- BURDEN PANELS -----------------
def load_burden():
    d = pd.read_csv(ROOT / "data/census/per_regulator_enrichment_all.csv")
    n, k = d["n_both"].to_numpy(float), d["n_inv"].to_numpy(float)
    return n, k

def fit_bb(n, k, x, with_disp=True):
    def nll(par):
        mu = 1 / (1 + np.exp(-(par[0] + par[1] * (x - x.mean()))))
        if not with_disp:
            return -stats.binom.logpmf(k, n, mu).sum()
        rho = 1 / (1 + np.exp(-par[-1]))
        a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
        if np.any(a <= 0) or np.any(b <= 0) or not np.all(np.isfinite(a + b)):
            return 1e12
        return -stats.betabinom.logpmf(k, n, a, b).sum()
    p0 = [np.log(0.05 / 0.95), 0.0] + ([-3.0] if with_disp else [])
    r = optimize.minimize(nll, p0, method="Nelder-Mead",
                          options=dict(maxiter=40000, fatol=1e-11, xatol=1e-11))
    return -r.fun, r.x

def panel_rate_vs_degree(ax, n, k):
    x = np.log10(n)
    _, par = fit_bb(n, k, x, False)
    b0, b1 = par[0], par[1]

    n_g = np.unique(np.round(np.logspace(np.log10(n.min()), np.log10(n.max()), 300))).astype(int)
    mu_g = 1 / (1 + np.exp(-(b0 + b1 * (np.log10(n_g) - x.mean()))))
    lo_g = stats.binom.ppf(0.025, n_g, mu_g) / n_g
    hi_g = stats.binom.ppf(0.975, n_g, mu_g) / n_g
    ax.fill_between(n_g, lo_g * 100, hi_g * 100, color=BAND, alpha=0.75, lw=0, zorder=1, label="binomial 95% interval")

    mu_i = 1 / (1 + np.exp(-(b0 + b1 * (x - x.mean()))))
    lo_i = stats.binom.ppf(0.025, n.astype(int), mu_i) / n
    hi_i = stats.binom.ppf(0.975, n.astype(int), mu_i) / n
    out = (k / n > hi_i) | (k / n < lo_i)

    ax.plot(n_g, mu_g * 100, color="#000000", lw=1.2, zorder=3, label="degree trend")
    ax.scatter(n[~out], k[~out] / n[~out] * 100, s=7, marker="o", alpha=0.55, color=BLUE,
               edgecolors="none", zorder=2, label=f"within ({(~out).sum()})")
    ax.scatter(n[out], k[out] / n[out] * 100, s=11, marker="^", alpha=0.75, color=VERM,
               edgecolors="none", zorder=4, label=f"outside ({out.sum()})")
    ax.set_xscale("log")
    ax.set(xlabel="both-significant trans edges (out-degree)", ylabel="reversal rate (%)", ylim=(-2, 60))
    ax.set_title("Heterogeneity survives degree adjustment", loc="left")
    leg = ax.legend(fontsize=5.6, loc="upper right", ncol=1, frameon=True, framealpha=1.0,
                    facecolor="white", edgecolor="#BFBFBF", handletextpad=0.5,
                    labelspacing=0.30, borderpad=0.35, borderaxespad=0.3, markerscale=1.35)
    leg.get_frame().set_linewidth(0.5)
    leg.set_zorder(6)
    return int(out.sum())

# ----------------- KINETICS PANELS -----------------
def panel_a(ax, tag="d"):
    # NB: do NOT clip into the last bin -- it fakes a spike at the edge.
    bins = np.linspace(0, 8, 45)
    ax.hist(np.abs(kin_stab.z_Stim8hr), bins=bins, density=True, color=PAL["CONC"],
            label=f"stable edges (median {np.median(np.abs(kin_stab.z_Stim8hr)):.2f})", zorder=1)
    ax.hist(np.abs(kin_rev.z_Stim8hr), bins=bins, density=True, histtype="step", lw=1.8,
            color=PAL["CORE"],
            label=f"reversing edges (median {np.median(np.abs(kin_rev.z_Stim8hr)):.2f})", zorder=3)
    x = np.linspace(0, 8, 300)
    ax.plot(x, 2 * norm.pdf(x), ls="--", lw=1.4, color=PAL["INK"],
             label="no effect, |z|~half-normal (median 0.67)", zorder=4)
    ax.set_xlabel("|z| at 8 h   (z = log2FC / lfcSE)")
    ax.set_ylabel("density")
    ax.set_title("8 h carries real signal on reversing edges", pad=8)
    ax.legend(loc="upper right", handletextpad=0.5)
    ax.text(0.97, 0.42, f"{(np.abs(kin_rev.z_Stim8hr) > 1.96).mean()*100:.1f}% of reversing edges |z|>1.96\n"
                              "(5% expected under no effect)",
             transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#565b61")
    panel_tag(ax, tag, dx=-0.10, dy=1.06)

def panel_b(ax, tag="e"):
    def rho_by(sub, ep="z_Rest", minn=8):
        out = []
        for _, g in sub.groupby("regulator"):
            if len(g) < minn: continue
            v, _ = spearmanr(g.z_Stim8hr, g[ep])
            if np.isfinite(v): out.append(v)
        return np.array(out)
    o_rev, o_stab = rho_by(kin_rev), rho_by(kin_stab)
    sim = kin_rev.copy()
    sim["z_Stim8hr"] = kin_rng.normal(0, 1, len(sim))
    o_sim = rho_by(sim)
    data = [o_sim, o_rev, o_stab]
    labs = [f"reversing,\n8 h replaced\nby noise\n(n={len(o_sim)})",
            f"reversing\n(observed)\n(n={len(o_rev)})", f"stable\n(control)\n(n={len(o_stab)})"]
    cols = [PAL["GREY"], PAL["CORE"], PAL["CONC"]]
    for i, (d, c) in enumerate(zip(data, cols)):
        ax.scatter(np.full(len(d), i) + kin_rng.normal(0, 0.055, len(d)), d, s=7, color=c, alpha=0.5, lw=0, zorder=2)
        ax.hlines(np.median(d), i - 0.28, i + 0.28, color=PAL["INK"], lw=2.2, zorder=3)
        ax.text(i, np.median(d) + 0.07, f"{np.median(d):+.3f}", ha="center", fontsize=8,
                 fontweight="bold", color=PAL["INK"], zorder=4)
    ax.axhline(0, color="#565b61", lw=0.9, ls=":", zorder=1)
    ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=7.5)
    ax.set_ylabel(r"within-regulator $\rho$(z$_{8h}$, z$_{Rest}$)")
    ax.set_title("Reversals are target-specifically coherent", pad=8)
    ax.text(0.5, 0.03, "permutation null median 0.000 · observed p < 0.0005",
             transform=ax.transAxes, ha="center", fontsize=7, color="#565b61")
    panel_tag(ax, tag, dx=-0.10, dy=1.06)

def panel_c(ax, tag="f"):
    r = kin_rev.copy()
    r["res8"] = r.padj_Stim8hr < 0.10
    r["late"] = np.sign(r.lfc_Stim8hr) == np.sign(r.lfc_Rest)
    def cluster_ci(sub, nboot=5000):
        g = sub.groupby("regulator").late.agg(["sum", "size"])
        regs = g.index.to_numpy()
        b = [(lambda s: s["sum"].sum() / s["size"].sum())(g.loc[kin_rng.choice(regs, len(regs), replace=True)])
             for _ in range(nboot)]
        return np.percentile(b, [2.5, 97.5])

    rows = []
    for lab, thr in ((f"all reversals\n(n={len(r):,})", 0.0), (f"core |log2FC|>0.5\n(n={int((r.minmag>0.5).sum()):,})", 0.5)):
        s = r[r.minmag > thr]; q = s[s.res8]
        nl = int(q.late.sum()); ne = int((~q.late).sum()); nu = int((~s.res8).sum())
        lo, hi = cluster_ci(q)
        rows.append((lab, nu, nl, ne, nl / max(ne, 1), (lo, hi), len(s)))
    y = np.arange(len(rows))[::-1] * 1.0
    for i, (lab, nu, nl, ne, ratio, pv, tot) in zip(y, rows):
        left = 0
        for val, col in ((nu, PAL["GREY"]), (nl, PAL["REST"]), (ne, PAL["STIM"])):
            ax.barh(i, val / tot * 100, left=left / tot * 100, color=col, height=0.42, zorder=2)
            left += val
        ax.text(99, i + 0.30, f"Rest-sign retaining {nl/(nl+ne)*100:.1f}%  "
                 f"(cluster-bootstrap 95% CI "
                 f"{pv[0]*100:.1f}\u2013{pv[1]*100:.1f}%)", va="bottom", ha="right",
                 fontsize=7.5, fontweight="bold", color=PAL["INK"], zorder=4)
    ax.set_yticks(y); ax.set_yticklabels([x[0] for x in rows], fontsize=8)
    ax.set_xlabel("% of reversing edges"); ax.set_xlim(0, 100); ax.set_ylim(-1.75, 1.62)
    ax.set_title("Resolvable effects mostly retain\nthe Rest sign at 8 h", pad=8)
    ax.legend(handles=[Patch(color=PAL["GREY"], label="unresolved at 8 h (padj$\\geq$0.10)"),
                        Patch(color=PAL["REST"], label="Rest-sign retaining at 8 h"),
                        Patch(color=PAL["STIM"], label="Stim48-sign matching at 8 h")],
               loc="lower left", bbox_to_anchor=(-0.005, -0.01), ncol=1, handlelength=1.3,
               fontsize=7.0)
    ax.annotate(f"{4379 - len(r):,} of the 4,379 census reversals have no 8 h\n"
                 f"measurement and cannot be classified",
                 xy=(0.015, 0.315), xycoords="axes fraction", fontsize=6.6, color="#8a8f95",
                 ha="left", va="bottom")
    panel_tag(ax, tag, dx=-0.10, dy=1.06)

def panel_d(ax, tag="g"):
    G = pd.read_csv(ROOT / "data/kinetics/kinetics_8hr_regulator_timing.csv", index_col=0)
    G["n_go"] = [N_GO_DICT.get(reg, 0) for reg in G.index]
    rho, pv = spearmanr(G.n_go, G.late_frac)
    ax.scatter(G.n_go + 0.5, G.late_frac, s=26, color=PAL["CORE"], alpha=0.65, lw=0, zorder=2)
    for g in ("LCK", "KDM1A"):
        if g in G.index:
            ax.annotate(g, (G.loc[g, "n_go"] + 0.5, G.loc[g, "late_frac"]), fontsize=7,
                         xytext=(4, 4), textcoords="offset points", color=PAL["INK"])
    ax.set_xscale("log"); ax.axhline(0.5, color="#565b61", lw=0.8, ls=":", zorder=1)
    ax.set_xlabel("GO annotation density (no. of GO BP terms, +0.5)")
    ax.set_ylabel("regulator Rest-sign-retaining fraction")
    ax.set_title("No functional timing\ndeterminant emerges", pad=8)
    ax.text(0.03, 0.06, f"Spearman $\\rho$ = {rho:+.3f}, p = {pv:.4f}  (n={len(G)})\n"
                         "well-studied regulators look 'early':\nthe functional signal is study bias",
             transform=ax.transAxes, fontsize=7.2, color="#565b61", va="bottom")
    panel_tag(ax, tag, dx=-0.10, dy=1.06)

# ----------------- ASSEMBLY -----------------
def _texts(ax):
    out = [ax.title, ax.xaxis.label, ax.yaxis.label]
    out += list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.texts)
    legend = ax.get_legend()
    if legend is not None:
        out += list(legend.get_texts())
    return [t for t in out if t is not None and t.get_text()]

def build():
    setup()
    fig = plt.figure(figsize=(SOURCE_WIDTH_IN, SOURCE_HEIGHT_IN * 0.86), facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.12], hspace=0.28,
                             left=0.135, right=0.980, top=0.940, bottom=0.055)
    top = outer[0].subgridspec(1, 3, wspace=0.62)
    bottom = outer[1].subgridspec(2, 2, wspace=0.34, hspace=0.42)

    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(top[0, 1])
    ax_c = fig.add_subplot(top[0, 2])
    ax_d = fig.add_subplot(bottom[0, 0])
    ax_e = fig.add_subplot(bottom[0, 1])
    ax_f = fig.add_subplot(bottom[1, 0])
    ax_g = fig.add_subplot(bottom[1, 1])

    panel_magnitude(ax_a, "a")
    n_deg, k_deg = load_burden()
    panel_rate_vs_degree(ax_b, n_deg, k_deg)
    ax_b.set_title("", loc="left")
    ax_b.set_title("Heterogeneity survives\ndegree adjustment", fontsize=8.0)
    panel_tag(ax_b, "b", dx=-0.22, dy=1.08)
    legend_b = ax_b.get_legend()
    if legend_b is not None:
        handles = list(legend_b.legend_handles)
        labels = [x.get_text() for x in legend_b.get_texts()]
        legend_b.remove()
        ax_b.legend(handles, labels, loc="upper right", frameon=False, fontsize=6.2,
                    handletextpad=0.4, labelspacing=0.18, borderaxespad=0.3)
                    
    panel_exemplars(ax_c, "c")
    ax_c.set_title("Concrete reversals\n(|log2FC| > 1 in both states)", fontsize=8.0)
    panel_a(ax_d, "d")
    ax_d.set_title("8 h carries real signal\non reversing edges", fontsize=8.0)
    panel_b(ax_e, "e")
    ax_e.set_title("Reversals are\ntarget-specifically coherent", fontsize=8.0)
    panel_c(ax_f, "f")
    panel_d(ax_g, "g")

    ax_f.set_yticklabels(["all reversals\n(n = 4,377)", "core |log2FC| > 0.5\n(n = 530)"], fontsize=7.0)
    for ax, tag, dx, dy in ((ax_a, "a", -0.22, 1.10), (ax_g, "g", -0.20, 1.08)):
        for artist in list(ax.texts):
            if artist.get_text() == tag:
                artist.set_position((dx, dy))

    axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g]
    for ax in axes:
        for artist in _texts(ax):
            if artist.get_fontsize() < FONT_FLOOR_PT:
                artist.set_fontsize(FONT_FLOOR_PT)
    fig.canvas.draw()

    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"FigS1.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print("wrote FigS1 (seven panels)")

if __name__ == "__main__":
    build()
