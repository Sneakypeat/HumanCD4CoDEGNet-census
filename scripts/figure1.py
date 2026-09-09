import os, sys, math, hashlib, json, numpy as np, pandas as pd, matplotlib, matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import optimize
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# --- figstyle.py ---
"""Shared publication style + palette for the sign-reversal paper figures."""


PAL = dict(
    REST="#2c6fbb", STIM="#e07b39",
    LOST="#c1272d",     # up in Rest -> down in Stim48  (edge lost on stimulation)
    GAINED="#2a9d8f",   # down in Rest -> up in Stim48  (edge gained on stimulation)
    GREY="#9aa0a6", CONC="#c7ccd1", CORE="#6a51a3",
    INK="#22262b", HAIR="#dfe2e6", HL="#f2c14e",
)

_SCALE = [1.0]   # current text scale, so panel_tag matches setup()


def setup(scale=1.0):
    """scale < 1 shrinks every text element, for dense multi-panel figures that must
    fit a fixed page size (e.g. the A4 composite Figure 1)."""
    _SCALE[0] = scale

    def s(v):
        return round(v * scale, 2)
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


# A4 portrait canvas in inches. Figure 1 is saved at exactly this size, so it must NOT be
# written with bbox_inches="tight" (that crops to content and silently discards the page size).
A4_PORTRAIT = (8.27, 11.69)


os.makedirs(OUT, exist_ok=True)

def panel_tag(ax, s, dx=-0.04, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=12 * _SCALE[0], fontweight="bold",
            va="top", ha="right", color="#22262b")

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight",
                    facecolor="white")
    print("  saved", name + ".pdf/.png")


# --- fig1_phenomenon.py ---
#!/usr/bin/env python
"""Figure 1 - The sign-reversal phenomenon.

Refactored into panel functions (panel_scatter / panel_census / panel_direction) so it
renders standalone AND composes into the combined phenomenon+reproducibility+magnitude figure."""







D = json.load(open(ROOT / "data/census/figure_data.json"))
c = D["census"]
assert c["edges_sig_both"] == 83_489 and c["same_sign"] == 79_110
assert c["inverting"] == 4_379 and c["n_inverting_regulators"] == 324
x = np.array(D["arr"]["lfc_Rest"]); y = np.array(D["arr"]["lfc_Stim48"])
lost = (x > 0) & (y < 0)     # up in Rest -> down in Stim  (2298)
gain = (x < 0) & (y > 0)     # down in Rest -> up in Stim  (2081)


def panel_scatter(ax, tag="a"):
    lim = 3.0
    ax.add_patch(Rectangle((-lim, 0), lim, lim, color=PAL["GAINED"], alpha=0.06, zorder=0))
    ax.add_patch(Rectangle((0, -lim), lim, lim, color=PAL["LOST"], alpha=0.06, zorder=0))
    ax.axhline(0, color="#565b61", lw=0.9); ax.axvline(0, color="#565b61", lw=0.9)
    ax.plot([-lim, lim], [lim, -lim], ls="--", lw=0.9, color="#8a8f98", zorder=1)
    inside = (np.abs(x) <= lim) & (np.abs(y) <= lim)
    ax.scatter(x[lost & inside], y[lost & inside], s=7, c=PAL["LOST"], alpha=0.55, lw=0,
               label=f"lost on stimulation  (+→−, n={c['pos_to_neg']:,})")
    ax.scatter(x[gain & inside], y[gain & inside], s=7, c=PAL["GAINED"], alpha=0.55, lw=0,
               label=f"gained on stimulation  (−→+, n={c['neg_to_pos']:,})")
    # 15 of the 4,379 edges fall outside the +-3 box (x reaches 4.47, y reaches -6.61), and
    # Matplotlib drops out-of-range points silently rather than pinning them. They are drawn
    # on the frame as open markers with clip_on=False, and the end tick labels carry the
    # inequality. Widening the limits instead would make this equal-aspect square 13 units
    # across and flatten the y = -x structure the panel exists to show.
    for _mask, _colour in ((lost, PAL["LOST"]), (gain, PAL["GAINED"])):
        _m = _mask & ~inside
        if not _m.any():
            continue
        # pinned just inside the frame, not on it: drawn on the spine with clip_on=False the
        # top-edge marker ran into the two-line panel title.
        _e = lim - 0.10
        ax.scatter(np.clip(x[_m], -_e, _e), np.clip(y[_m], -_e, _e), s=11, marker="o",
                   facecolors="none", edgecolors=_colour, lw=0.7, zorder=4)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    _ticks = [-3, -2, -1, 0, 1, 2, 3]
    _labels = ["<\u22123", "\u22122", "\u22121", "0", "1", "2", ">+3"]
    ax.set_xticks(_ticks); ax.set_xticklabels(_labels)
    ax.set_yticks(_ticks); ax.set_yticklabels(_labels)
    ax.set_xlabel("regulatory effect in Rest\nlog2 fold change")
    ax.set_ylabel("regulatory effect in Stim 48 h\nlog2 fold change")
    ax.set_title("Edges that reverse sign between states", pad=8)
    ax.legend(loc="upper left", markerscale=1.6, handletextpad=0.3, borderpad=0.2)
    # upper right: the only quadrant with neither a point cloud nor the direction key. At
    # lower right it sat on the lost-edge cloud it was describing.
    ax.text(0.97, 0.93, "reversals track the\nanti-diagonal (y ≈ −x):\nstrong one way,\nstrong the other",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#565b61")
    panel_tag(ax, tag, dx=-0.16, dy=1.05)


def panel_census(ax, tag="b"):
    tot = c["edges_sig_both"]; conc = c["same_sign"]; inv = c["inverting"]
    ax.barh(0, conc, color=PAL["CONC"], edgecolor="white", height=0.55)
    ax.barh(0, inv, left=conc, color=PAL["INK"], edgecolor="white", height=0.55)
    ax.text(conc/2, 0, f"same sign  {conc:,}", ha="center", va="center", fontsize=7.5, color="#3d4348")
    # The leader points at the middle of the black sliver and stops short of the text. Anchored
    # centre at 0.80*tot it ran diagonally through "(5.2%)" on the callout's second line.
    # va="bottom" anchors the leader at the text's BOTTOM edge, so it runs up from the bar and
    # stops there. With va="top" the anchor sits above the text and the line either crosses it
    # (striking through "(5.2%)") or, once shrunk clear, degenerates into a stub.
    ax.annotate(f"reverse sign\n{inv:,}  ({c['inversion_rate']*100:.1f}%)",
                xy=(conc + inv * 0.5, 0.33), xytext=(tot * 1.00, 0.60),
                fontsize=7.6, fontweight="bold", ha="right", va="bottom", color=PAL["INK"],
                arrowprops=dict(arrowstyle="-", lw=0.8, color="#565b61",
                                shrinkA=0, shrinkB=2))
    ax.set_ylim(-0.5, 1.25); ax.set_xlim(0, tot*1.02)
    ax.set_yticks([]); ax.set_xticks([0, 40000, 80000])
    ax.set_xlabel("edges significant in both states", fontsize=8)
    # pad clears the two-line callout, which used to sit on the title
    ax.set_title("A rare minority of shared edges", pad=10)
    ax.spines["left"].set_visible(False); ax.tick_params(left=False)
    panel_tag(ax, tag, dx=-0.16, dy=1.30)


def panel_direction(ax, tag=None):
    inv = c["inverting"]
    ax.barh(0, c["pos_to_neg"], color=PAL["LOST"], edgecolor="white", height=0.55)
    ax.barh(0, c["neg_to_pos"], left=c["pos_to_neg"], color=PAL["GAINED"], edgecolor="white", height=0.55)
    ax.text(c["pos_to_neg"]/2, 0, f"lost\n{c['pos_to_neg']:,}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.text(c["pos_to_neg"]+c["neg_to_pos"]/2, 0, f"gained\n{c['neg_to_pos']:,}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.set_ylim(-0.6, 0.6); ax.set_xlim(0, inv)
    ax.set_yticks([]); ax.set_xticks([0, 2000, 4000])
    ax.set_xlabel("the 4,379 reversals, by direction", fontsize=8)
    ax.set_title("Balanced in direction", pad=6)
    ax.spines["left"].set_visible(False); ax.tick_params(left=False)
    if tag:
        panel_tag(ax, tag, dx=-0.10, dy=1.18)




# --- fig2_reproducibility.py ---
#!/usr/bin/env python
"""Figure 2 - Sign reversals reproduce in a disjoint donor-pair (and cross-guide) edge-level holdout.
Unique-edge reproduction with CLUSTER-BOOTSTRAP 95% CIs (resampling REGULATORS, not edges --
edges within a regulator share a perturbation and are not independent; an edge-level Wilson
interval understates the width by ~a third and, for the guide split, is spuriously informative),
all magnitude bins, against a within-regulator
target-label permutation null and the stable-edge ceiling. NOT leave-one-donor-out (by_donors is
pairwise): reversals are defined in one donor pair and tested in the complementary disjoint pair."""







D = json.load(open(ROOT / "data/holdout/holdout_donors_results.json"))
G = json.load(open(ROOT / "data/holdout/holdout_guide_results.json"))
REPRO, STABLE, NULL = "#2a9d8f", "#9aa0a6", PAL["LOST"]
TIERS = [">0", ">0.25", ">0.5", ">1.0"]

def bars(ax, res, title, tag, tag_dx=-0.17):
    rates, los, his, ns, nregs, nulls, repro = [], [], [], [], [], [], []
    for t in TIERS:
        b = res["by_threshold"][t]; r = b["unique_repro"] or 0.0
        lo, hi = b["unique_cluster_boot95"]
        rates.append(r*100); los.append((r-lo)*100); his.append((hi-r)*100)
        ns.append(b["n_unique_edges"]); nregs.append(b["n_regulators"])
        nulls.append(b["perm_null_mean"]*100)
        repro.append(round(r*b["n_unique_edges"]))
    x = np.arange(len(TIERS))
    ceiling = np.mean([res["by_threshold"][t]["stable_repro"] for t in TIERS])*100
    ax.axhline(ceiling, ls="--", lw=1.1, color=STABLE)
    # parked in headroom above every bar/value label so nothing collides
    ax.text(-0.45, 125, f"– – stable-edge ceiling (~{ceiling:.0f}%)", fontsize=6.3,
            color="#6b7178", ha="left", va="top")
    ax.text(-0.45, 117,
            "– – within-regulator permutation null, per bin, left to right: "
            + ", ".join(f"{v:.2f}" for v in nulls) + "%",
            fontsize=6.3, color=NULL, ha="left", va="top")
    ax.bar(x, rates, 0.6, color=REPRO, yerr=[los, his], capsize=3,
           error_kw=dict(lw=1.1, ecolor="#1f6f63"))
    # The null is NOT constant across bins (2.33, 2.26, 2.54, 4.56%): a single averaged line sat
    # nearest the >0.5 bin and misrepresented the headline one, so each bar carries its own.
    for i, nl in enumerate(nulls):
        ax.plot([i-0.33, i+0.33], [nl, nl], ls="--", lw=1.2, color=NULL, zorder=5,
                solid_capstyle="butt")
    for i, (r, n, g, k) in enumerate(zip(rates, ns, nregs, repro)):
        ax.text(i, min(r+his[i]+3, 104), f"{r:.1f}%", ha="center", fontsize=8,
                fontweight="bold", color="#1f6f63")
    # The fraction explains the bar height; the regulator count explains the whisker, because
    # the interval is bootstrapped over regulators. n = edges alone explains neither. Both sit
    # under the axis rather than inside the bar: drawn white-on-teal at the enforced 7 pt floor,
    # "127/166" is wider than its own bar on a narrow panel and the overhang went invisible
    # against the page.
    ax.set_xticks(x)
    ax.set_xticklabels([f"min|log2FC|\n{t}\n{k}/{n}, {g} reg."
                        for t, n, g, k in zip(TIERS, ns, nregs, repro)], fontsize=7.5)
    # Published on the axes so an alternative composer can relabel without recomputing:
    # the compact Figure 1 arrangement carries its own x-axis label and drops the prefix.
    ax.tier_counts = list(zip(TIERS, ns, nregs, repro))
    # The k/n label is centred on the bar and is wider than the bar itself, so on a narrow
    # panel (the compact Figure 1 arrangement) the leftmost one runs past the default
    # x-limit and is clipped. Pad the limits rather than shrinking the type.
    ax.set_xlim(x[0] - 0.66, x[-1] + 0.62)
    ax.set_ylim(0, 127); ax.set_ylabel("unique reversing edges reproduced (%)")
    ax.set_title(title)
    # tag_dx is exposed because in the A4 Figure 1 this panel sits hard against the left
    # page margin, where the default offset put the letter off-canvas and it vanished
    panel_tag(ax, tag, dx=tag_dx, dy=1.04)

def panel_guide(ax, tag="b"):
    """Cross-guide holdout (honestly underpowered)."""
    b0 = G["by_threshold"][">0"]
    r0 = b0["unique_repro"]*100; lo0, hi0 = [v*100 for v in b0["unique_cluster_boot95"]]
    ceil = b0["stable_repro"]*100; nl = b0["perm_null_mean"]*100
    ax.axhline(ceil, ls="--", lw=1.1, color=STABLE); ax.axhline(nl, ls="--", lw=1.1, color=NULL)
    ax.text(-0.65, 97, f"– – stable ceiling (~{ceil:.0f}%)", fontsize=6.5, color="#6b7178", ha="left", va="top")
    ax.text(-0.65, 91, f"– – permutation null ({nl:.2f}%)", fontsize=6.5, color=NULL, ha="left", va="top")
    ax.bar([0], [r0], 0.5, color=REPRO, yerr=[[r0-lo0],[hi0-r0]], capsize=4, error_kw=dict(lw=1.2, ecolor="#1f6f63"))
    ax.text(0.5, 0.52, "cluster-bootstrap 95% CI spans the\nentire range [0, 100%]: only 4 regulator\nclusters, so this panel carries no\nusable precision. Directional only.",
            transform=ax.transAxes, ha="center", va="center", fontsize=7.5, color="#8c2f2f",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fdf3f3", ec="#e0b4b4", lw=0.8))
    ax.text(0, min(hi0+3,104), f"{r0:.1f}%", ha="center", fontsize=9, fontweight="bold", color="#1f6f63")
    # Same reason as the tier panel: white-on-teal at the enforced 7 pt floor is wider than a
    # single narrow bar, so the counts sit under the axis where they cannot go invisible.
    _k = round(b0['unique_repro'] * b0['n_unique_edges'])
    ax.set_xticks([0])
    ax.set_xticklabels([f"all reversing\nedges (min|log2FC|>0)\n"
                        f"{_k}/{b0['n_unique_edges']}, {b0['n_regulators']} reg."], fontsize=7.5)
    ax.set_xlim(-0.7, 0.9); ax.set_ylim(0, 112); ax.set_ylabel("reproduced in held-out guide (%)")
    ax.set_title("Across guides: uninformative\n(8 edges, 4 regulator clusters)")
    ax.text(0.98, 0.20, "single-guide DE is much noisier;\ndirectionally consistent, not decisive",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.6, color="#565b61")
    panel_tag(ax, tag, dx=-0.2, dy=1.04)


DONOR_TITLE = "Reversals reproduce across disjoint donor pairs\n(cluster-bootstrap 95% CI over regulators)"



# --- fig1_transition.py ---
#!/usr/bin/env python
"""Fig 1 panel: sign reversal is a temporal transition (Rest -> 8 h -> 48 h trajectory).

This is the *visualization* that the flip develops over time: the two arms start on
opposite signs, sit near zero at the intermediate 8 h timepoint (the flip is not yet
complete), and reach their reversed signs by 48 h. The x-axis is categorical (the three
measured timepoints, evenly spaced) for legibility. The rigorous demonstration that the
8 h signal is real and not noise -- target-specific coherence and the kinetic partition --
is Supplementary Figure S2; this panel points there rather than resting on 8 h magnitude.
"""





_D = pd.read_csv(ROOT / "data/kinetics/kinetics_8hr_classified.csv.gz")


def panel_transition(ax, tag="c"):
    d = _D
    lost = d[(d["lfc_Rest"] > 0) & (d["lfc_Stim48hr"] < 0)]   # + at rest -> - at 48 h
    gain = d[(d["lfc_Rest"] < 0) & (d["lfc_Stim48hr"] > 0)]   # - at rest -> + at 48 h
    tp = ["lfc_Rest", "lfc_Stim8hr", "lfc_Stim48hr"]
    # x is REAL time in hours, not evenly spaced: 8 h sits near the left where it belongs
    # (only 8 h after rest), and the long 8 -> 48 h stretch carries most of the inversion.
    # An evenly spaced axis wrongly put 8 h at the midpoint and made the flip look linear.
    x = [0, 8, 48]

    def med(df):
        return ([df[c].median() for c in tp],
                [df[c].quantile(.25) for c in tp],
                [df[c].quantile(.75) for c in tp])
    lm, l25, l75 = med(lost)
    gm, g25, g75 = med(gain)

    ax.axhline(0, color="#565b61", lw=0.9)
    ax.axvline(8, color="#9aa0a6", ls=":", lw=1)
    ax.fill_between(x, l25, l75, color=PAL["LOST"], alpha=0.15, lw=0)
    ax.fill_between(x, g25, g75, color=PAL["GAINED"], alpha=0.15, lw=0)
    ax.plot(x, lm, "-o", color=PAL["LOST"], lw=2.3, ms=5,
            label="lost on stimulation  (+→−)")
    ax.plot(x, gm, "-o", color=PAL["GAINED"], lw=2.3, ms=5,
            label="gained on stimulation  (−→+)")

    # Scoped to the MEDIAN deliberately: both arms are still on their resting sign at 8 h
    # (lost +0.055, gained -0.013), but individual edges vary and 35.6% of resolvable ones
    # already match their 48 h sign, so a per-edge phrasing would overclaim.
    ax.set_title("Median effects at 8 h\nnot reversed", loc="left")
    ax.set_xlim(-2.0, 50.0)
    ax.set_ylim(-0.95, 0.95)
    ax.set_xticks(x)
    ax.set_xticklabels(["Rest\n0 h", "8 h", "48 h"])
    ax.set_ylabel("target log2 FC (median, IQR)")
    ax.set_xlabel("time after stimulation (h, to scale)")
    ax.annotate("still on the original sign at 8 h", xy=(8, 0.03), xytext=(20, 0.62),
                fontsize=6.8, color="#565b61", ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#9aa0a6", lw=0.9))
    ax.legend(loc="lower left", fontsize=6.6)
    # 2,296 not the 2,298 of panel a: two KAT7 reversals carry no 8 h measurement, so they
    # cannot appear on a three-timepoint trace. Stated here so the two counts do not read as
    # a contradiction. See data/census/census_master_manifest.json.
    if len(lost) != 2298:
        ax.annotate(f"{2298 - len(lost)} KAT7 reversals lack 8 h data",
                    xy=(0.98, 0.29), xycoords="axes fraction", fontsize=5.3,
                    color="#8a8f95", ha="right", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.82))
    panel_tag(ax, tag)


# --- fig_burden_heterogeneity.py ---
#!/usr/bin/env python3
"""Reversal-burden heterogeneity across regulators.

WHY THIS EXISTS. This module answers a different question from `degree_distribution_test.py`,
which fits the census degree distributions and supplies Figure 1f. That fit describes the shape
of n, how many targets a regulator has, and carries two scope limits: the census gate admitted a
regulator only with at least 50 significant targets in each state, and sampled subnets do not in
general inherit a parent network's degree form (Stumpf, Wiuf & May, PNAS 2005). Neither limit
touches the question here, because burden conditions on each regulator's own denominator.

The estimable question is how unevenly reversal burden is distributed and whether connectivity
explains it. Each regulator contributes k reversals out of n both-significant edges, so the
object is binomial overdispersion, which conditions on each regulator's own denominator: the gate
decides which n are observed, not the k-given-n relationship being tested.

  (a) Descending concentration/top-share curve. Assumption-free description of the observed
      catalogue; makes no generative claim, so the selection cannot invalidate it.
  (b) Rate against out-degree with pointwise central binomial outcome intervals around a
      coherent binomial degree trend. The positive degree slope is retained explicitly: this
      panel asks whether degree adjustment eliminates the remaining heterogeneity, not whether
      connectivity has no association with rate.
  (c) Residual overdispersion as the largest carriers are removed.

FIGURE CONVENTIONS (scientific-visualization skill).
  * Okabe-Ito on white. The audit puts blue/vermillion at grayscale delta-L* 8.2, below the 10.0
    screening threshold, so inside/outside band is encoded by marker SHAPE as well as colour and
    stays readable in grayscale and for colour-vision deficiency.
  * Bars in (c) are measured from zero, because bar length is read from a baseline.
  * Physical width is 180 mm (double column). Constrained layout, and no bbox_inches="tight",
    so the exported page keeps that width.
  * The band in (b) is a binomial 95% interval, named as such; the unit of replication is the
    regulator, n = 617. Out-degree is >= 1 throughout, so the log axis drops no observations.

TECHNICAL-POWER SENSITIVITY. The companion covariate closure found that measured edge
precision, on-target knockdown efficacy, nonlinear degree and available cell support did not
collapse rho. This excludes those measured explanations, not unmeasured technical variation,
and does not make the residual a regulator-intrinsic or biological mechanism.
"""






matplotlib.use("Agg")



# Repo-relative so the scripts run from any checkout, not only this machine.

# The scientific-visualization skill supplies a style sheet, the Okabe-Ito values and an export
# helper. It lives outside the repo under ~/.claude, so it is treated as OPTIONAL: without it the
# figure still builds from matplotlib defaults and the literal palette below. Nothing from the
# skill is vendored in.
SKILL = Path(os.path.expanduser("~/.claude/skills/scientific-visualization"))
HAVE_SKILL = (SKILL / "assets/publication.mplstyle").is_file()
if HAVE_SKILL:
    sys.path.insert(0, str(SKILL / "assets"))
    sys.path.insert(0, str(SKILL / "scripts"))



# Okabe-Ito is a published colour-vision-safe palette, not a skill asset; the literals keep the
# figure reproducible when the skill is absent.
OKABE_ITO = {"orange": "#E69F00", "sky_blue": "#56B4E9", "bluish_green": "#009E73",
             "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
             "reddish_purple": "#CC79A7", "black": "#000000"}


def _style():
    """Skill style sheet when present, otherwise a self-contained equivalent."""
    if HAVE_SKILL:
        return plt.style.context(str(SKILL / "assets/publication.mplstyle"))
    return plt.rc_context({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "axes.titleweight": "bold",
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "lines.linewidth": 1.5, "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": 300, "pdf.fonttype": 42, "figure.facecolor": "white",
    })


def _export(fig, stem, provenance):
    """Skill exporter when present, otherwise a plain savefig plus the same JSON manifest."""
    if HAVE_SKILL:
        from figure_export import export_figure
        return export_figure(fig, stem, formats=["pdf", "png"], dpi=600, bbox_inches=None,
                             provenance=provenance, write_manifest=True, overwrite=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=600, facecolor="white")
    Path(f"{stem}.export.json").write_text(json.dumps(
        {"provenance": provenance, "exporter": "fallback (skill absent)"}, indent=2) + "\n")
    return {"provenance": provenance}


SRC = ROOT / "data/census/per_regulator_enrichment_all.csv"
OUTSTEM = ROOT / "data/figure_provenance/Fig_burden_heterogeneity"
WIDTH_MM, HEIGHT_MM = 180.0, 72.0

BLUE, VERM, GREEN = OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["bluish_green"]
INK, MUTED, BAND = "#000000", "#4D4D4D", "#C8C8C8"


def load():
    d = pd.read_csv(SRC)
    assert len(d) == 617, f"expected the 617-regulator census family, got {len(d)}"
    assert d["gene"].notna().all() and d["gene"].is_unique, "regulator identities must be unique"
    n, k = d["n_both"].to_numpy(float), d["n_inv"].to_numpy(float)
    assert (k <= n).all() and (n >= 1).all(), "log axis assumes out-degree >= 1"
    assert int(n.sum()) == 83_489 and int(k.sum()) == 4_379
    assert int((k > 0).sum()) == 324
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
    assert r.success and np.isfinite(r.fun) and np.all(np.isfinite(r.x)), \
        f"burden mean-model optimisation failed: {r.message}"
    return -r.fun, r.x


def panel_concentration(ax, k):
    srt = np.sort(k)[::-1]
    frac = np.arange(1, len(k) + 1) / len(k)
    share = srt.cumsum() / k.sum()
    v = np.sort(k)
    g = (2 * (np.arange(1, len(v) + 1) * v).sum()) / (len(v) * v.sum()) - (len(v) + 1) / len(v)
    ax.plot([0, 1], [0, 1], ls=(0, (4, 3)), lw=1.0, color=MUTED)
    ax.plot(np.r_[0, frac], np.r_[0, share], lw=1.6, color=BLUE, solid_capstyle="round")
    ax.fill_between(np.r_[0, frac], np.r_[0, frac], np.r_[0, share], color=BLUE, alpha=0.12, lw=0)
    # Leaders live INSIDE the shaded region, which is the one placement that cannot cross
    # anything. That region is bounded above by a concave concentration curve and below by a straight
    # diagonal, so it is convex: a straight segment between two interior points stays interior.
    # Both earlier attempts put the labels in the wedge BELOW the diagonal, where every leader
    # must cross the diagonal to reach its point, and one also let the curve run through a label
    # at 1,708 sampled points. Label heights increase with point height, so leaders cannot cross
    # each other either.
    # Right-aligned, and each anchor sits left of its own height. A left-aligned label runs
    # rightward at fixed y until x passes y, which puts it straight across the diagonal; the
    # gate caught exactly that on the first attempt.
    anchors = {0.05: (0.52, 0.545), 0.10: (0.67, 0.695), 0.20: (0.82, 0.845)}
    for f, (tx, ty) in anchors.items():
        i = max(1, int(round(f * len(k)))) - 1
        ax.plot([frac[i]], [share[i]], marker="o", ms=3.6, mfc="white", mec=BLUE, mew=1.1, zorder=5)
        ax.annotate(f"top {int(f*100)}%: {share[i]*100:.0f}%",
                    xy=(frac[i], share[i]), xytext=(tx, ty), textcoords=ax.transAxes,
                    fontsize=6.3, color=INK, va="center", ha="right",
                    arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.55, alpha=0.85,
                                    shrinkA=2, shrinkB=3))
    ax.text(0.96, 0.06, f"Gini {g:.3f}", transform=ax.transAxes, ha="right",
            fontsize=7.5, fontweight="bold", color=BLUE)
    # Down-left along the diagonal, clear of the labels and their leaders; the opaque box keeps
    # the dashes off the glyphs.
    ax.text(0.275, 0.275, "equal burden", transform=ax.transAxes, rotation=41,
            fontsize=6.0, color=MUTED, ha="center", va="center", rotation_mode="anchor",
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white", edgecolor="none"))
    ax.set(xlim=(0, 1), ylim=(0, 1.03),
           xlabel="regulators ranked by reversal count\n(cumulative fraction)",
           ylabel="cumulative share of reversals")
    ax.set_title("Burden is concentrated", loc="left")
    return g


def panel_rate_vs_degree(ax, n, k):
    x = np.log10(n)
    # Fit the mean with the same binomial model used to construct the pointwise outcome
    # intervals.  Mixing a beta-binomial mean fit with binomial intervals is a hybrid visual
    # diagnostic and changed the outside count by one regulator.  The beta-binomial remains
    # the separate residual-dispersion model; this panel is the coherent binomial mean/null
    # display only.
    _, par = fit_bb(n, k, x, False)
    b0, b1 = par[0], par[1]

    # Envelope drawn on an integer grid: binom.ppf returns nan for non-integer n, which would
    # silently erase the band.
    n_g = np.unique(np.round(np.logspace(np.log10(n.min()), np.log10(n.max()), 300))).astype(int)
    mu_g = 1 / (1 + np.exp(-(b0 + b1 * (np.log10(n_g) - x.mean()))))
    lo_g = stats.binom.ppf(0.025, n_g, mu_g) / n_g
    hi_g = stats.binom.ppf(0.975, n_g, mu_g) / n_g
    assert np.all(np.isfinite(lo_g)) and np.all(np.isfinite(hi_g)), "binomial envelope is nan"
    ax.fill_between(n_g, lo_g * 100, hi_g * 100, color=BAND, alpha=0.75, lw=0, zorder=1,
                    label="binomial 95% interval")

    # Classify each regulator against a band built from its OWN n. Interpolating the plotting
    # grid mis-sorted regulators falling between grid points and inflated the count to 156.
    mu_i = 1 / (1 + np.exp(-(b0 + b1 * (x - x.mean()))))
    lo_i = stats.binom.ppf(0.025, n.astype(int), mu_i) / n
    hi_i = stats.binom.ppf(0.975, n.astype(int), mu_i) / n
    out = (k / n > hi_i) | (k / n < lo_i)

    ax.plot(n_g, mu_g * 100, color=INK, lw=1.2, zorder=3, label="degree trend")
    ax.scatter(n[~out], k[~out] / n[~out] * 100, s=7, marker="o", alpha=0.55, color=BLUE,
               edgecolors="none", zorder=2, label=f"within ({(~out).sum()})")
    ax.scatter(n[out], k[out] / n[out] * 100, s=11, marker="^", alpha=0.75, color=VERM,
               edgecolors="none", zorder=4, label=f"outside ({out.sum()})")
    ax.set_xscale("log")
    # ylim must clear BOTH the data maximum (57.1%, UBE2N at 4/7) and the widest point of the
    # binomial interval (50.0% at out-degree 2). An earlier cap of 42% silently pushed four
    # regulators off the panel, which is data suppression, not a cosmetic choice.
    # Headroom above the data maximum (57.1%) and the widest point of the interval (50.0%) so the
    # key has clearance under any style sheet. At ylim 60 the legend cleared the points only with
    # the skill's font metrics and covered two of them with matplotlib defaults.
    assert (k / n * 100).max() < 60 and (hi_g * 100).max() < 60, "axis would clip data or interval"
    ax.set(xlabel="both-significant trans edges (out-degree)",
           ylabel="reversal rate (%)", ylim=(-2, 60))
    # About the RESIDUAL, not about degree having no effect: the degree slope is positive and
    # significant, so "connectivity does not explain the rate" would be false.
    ax.set_title("Heterogeneity survives degree adjustment", loc="left")
    # Key goes in the upper-right, which is empty because every regulator above 42% has an
    # out-degree under 21. The opaque frame marks it as a key so its glyphs cannot be misread
    # as plotted regulators.
    # Single column on the right, with short labels. Four stacked rows make a tall block, so the
    # axis carries headroom above the data maximum (57.1%) for it to clear every point. Measured
    # Nothing is clipped at 60: the data maximum is 57.1% (UBE2N, 4/7) and the widest point of
    # the interval is 50.0% (out-degree 2). Fitting the key inside that ceiling instead needs the
    # short first label; a swept comparison showed the two-line "pointwise central 95% binomial
    # outcome interval" only clears the points if the legend drops to 5.0 pt, which is smaller
    # than every other annotation here. Verified under BOTH styled and default font metrics.
    leg = ax.legend(fontsize=5.6, loc="upper right", ncol=1, frameon=True, framealpha=1.0,
                    facecolor="white", edgecolor="#BFBFBF", handletextpad=0.5,
                    labelspacing=0.30, borderpad=0.35, borderaxespad=0.3, markerscale=1.35)
    leg.get_frame().set_linewidth(0.5)
    leg.set_zorder(6)
    return int(out.sum())


def panel_robustness(ax, n, k):
    rows = []
    for drop in (0, 5, 10, 20, 50):
        o = np.argsort(k)[::-1][drop:]
        xx = np.log10(n[o])
        a, _ = fit_bb(n[o], k[o], xx, False)
        b, pp = fit_bb(n[o], k[o], xx, True)
        rows.append((drop, float(1 / (1 + np.exp(-pp[-1]))), float(2 * (b - a))))
    xs = np.arange(len(rows))
    rho = [r[1] for r in rows]
    ax.bar(xs, rho, color=GREEN, width=0.62, zorder=2)          # baseline at zero, see docstring
    for i, (_d, r_, c_) in enumerate(rows):
        ax.text(i, r_ + 0.0022, f"{r_:.3f}", ha="center", fontsize=6.2, fontweight="bold")
        ax.text(i, 0.0045, f"{c_:,.0f}", ha="center", fontsize=5.6, color="white")
    ax.set_xticks(xs, [str(r[0]) for r in rows])
    ax.set(xlabel="largest carriers removed",
           ylabel=r"residual overdispersion $\rho$", ylim=(0, max(rho) * 1.30))
    ax.text(0.5, 0.965, r"white numerals: dispersion LRT $\chi^2$", transform=ax.transAxes,
            ha="center", va="top", fontsize=5.8, color=MUTED)
    ax.set_title("Not a few regulators", loc="left")
    return rows


def main():
    n, k = load()
    with _style():
        fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), layout="constrained")
        axes = fig.subplots(1, 3, width_ratios=[1.0, 1.15, 0.78])
        g = panel_concentration(axes[0], k)
        n_out = panel_rate_vs_degree(axes[1], n, k)
        rows = panel_robustness(axes[2], n, k)
        for lab, ax in zip("abc", axes):
            ax.text(-0.24, 1.06, lab, transform=ax.transAxes, fontsize=9,
                    fontweight="bold", va="bottom", ha="left", color=INK)

        # Measured on the real canvas before export. Eyeballing a downscaled PNG missed a
        # concentration curve running through a label at 1,708 sampled points.
        assert_no_overlaps(
            fig, axes=axes,
            intentional=lambda kind, detail: (
                # labels the diagonal it names; its opaque bbox masks the dashes
                (kind == "text-line" and detail.startswith("'equal burden'"))
                # panel letters are deliberately set outside the axes rectangle
                or (kind == "text-outside" and detail in "abc" and len(detail) == 1)))

        report = _export(fig, str(OUTSTEM), provenance={
                "raw_data": str(SRC.relative_to(ROOT)),
                "analysis_script": "scripts/burden_heterogeneity.py",
                "figure_script": "scripts/fig_burden_heterogeneity.py",
                "unit_of_replication": "regulator (n = 617)",
                "uncertainty": "binomial 95% interval about the fitted degree trend in panel b; "
                               "no resampling, no random seed used",
                "transformations": "none; k and n are census counts. log10 out-degree axis; "
                                   "out-degree >= 1 so no observations are dropped",
                "missing_data": "none; all 617 census regulators are plotted",
                "palette": "Okabe-Ito; inside/outside encoded by marker shape as well as colour",
                "alt_text": (
                    f"Three panels on reversal-burden heterogeneity across 617 regulators. "
                    f"(a) Descending concentration curve far above the equal-burden diagonal, "
                    f"Gini {g:.3f}; the most "
                    f"burdened 5% of regulators carry 62% of the 4,379 reversals. (b) Reversal "
                    f"rate against out-degree on a log axis; the fitted trend rises only from "
                    f"about 2% to 8% across three orders of magnitude while regulators at the "
                    f"same out-degree span 0% to 57%, and {n_out} of 617 fall outside their "
                    f"pointwise central 95% binomial outcome intervals. This interval count is "
                    f"descriptive, not a calibrated family-level test. (c) Residual overdispersion "
                    f"stays between 0.072 and 0.057 as the 5, 10, 20 and 50 largest carriers are "
                    f"removed."),
            })
        plt.close(fig)
    print(f"wrote {OUTSTEM.name}: Gini={g:.3f}  outside={n_out}/617  rho={rows[0][1]:.4f}")
    return report




# --- fig1f_degree_ccdf.py ---
#!/usr/bin/env python3
"""Figure 1f: regulator out-degree CCDF against a random network of the same size.

DESIGN. Follows Alam et al. 2017 (Nat Commun 8:16018) Figure 1c: one log-log complementary
cumulative plot carrying the observed network (filled), a random network of identical size
(open), the fitted tail models, a fitted exponent for each network, and a second right-hand
axis giving the cumulative percentage of nodes.

PROVENANCE. Nothing is fitted here. Observed-network tail parameters are read from
data/census/degree_distribution_test.regulator.json (produced by
degree_distribution_test.py, the same fits rendered as repository figure R6) and the random
comparator from data/census/degree_random_comparator.json. Survival curves are drawn
with degree_distribution_test.model_survival, so the lines are the fitted models themselves
and not a redrawing of them.

SCOPE. The fits describe the 617 regulators represented in this signal-gated census. The gate
admits a regulator only with at least 50 significant targets in each state, so this is the
observed census range and not the complete T-cell regulatory network; sampled subnets do not
in general inherit a parent network's degree form (Stumpf, Wiuf & May, PNAS 2005). Stated in
the legend and Methods, not implied away here.
"""













REG_FIT = ROOT / "data/census/degree_distribution_test.regulator.json"
RANDOM = ROOT / "data/census/degree_random_comparator.json"
POISSON = ROOT / "data/census/degree_poisson_fit.json"

INK, GREY, HAIR = PAL["INK"], PAL["GREY"], PAL["HAIR"]
PL, LN, EX, CO, PO = "#6a51a3", "#2c6fbb", "#d9a441", "#2a9d8f", "#c0392b"
CUM = "#2c6fbb"
ALPHA_BOUND = 15.0        # fit_alpha's upper bound; a fit landing here has no power-law tail


def _ccdf(v):
    v = np.sort(np.asarray(v))
    return v, 1.0 - np.arange(len(v)) / len(v)


def _exp_label(alpha):
    return r"$\gamma > %.0f$" % ALPHA_BOUND if alpha >= ALPHA_BOUND - 1e-3 else r"$\gamma = %.2f$" % alpha


def panel_degree_ccdf(ax, tag=None, fs=1.0, compact=False):
    fit = json.loads(REG_FIT.read_text())
    rnd = json.loads(RANDOM.read_text())
    pois = json.loads(POISSON.read_text())

    obs = np.asarray(rnd["observed"]["degrees"], float)
    ran = np.asarray(rnd["random"]["degrees"], float)
    xmin, alpha = int(fit["xmin"]), float(fit["alpha"])
    mc = fit["model_comparison"]

    # --- empirical curves -------------------------------------------------------------
    xo, yo = _ccdf(obs)
    xr, yr = _ccdf(ran)
    ax.plot(xr, yr, ls="none", marker="o", ms=2.6 * fs, mfc="none", mec=GREY, mew=0.7,
            label=f"random network ({_exp_label(rnd['random']['alpha'])[1:-1].split('=')[0].strip()})"
            if False else None, zorder=2)
    ax.plot(xo, yo, ls="none", marker="o", ms=2.6 * fs, mfc=INK, mec="none", zorder=3)

    # --- fitted tail models, rescaled onto the unconditional CCDF ----------------------
    # model_survival returns P(K >= k | K >= xmin); multiply by the empirical P(K >= xmin)
    # so the curves lie on the same axis as the data rather than floating an order up.
    p_tail = float((obs >= xmin).mean())
    grid = np.unique(np.round(np.logspace(np.log10(xmin), np.log10(obs.max()), 260)).astype(int))
    for name, colour, style, lw in (
        ("power_law", PL, "-", 1.5),
        ("lognormal", LN, (0, (5, 2)), 1.1),
        ("exponential", EX, (0, (1.4, 1.6)), 1.1),
        ("cutoff_power_law", CO, (0, (6, 2, 1.4, 2)), 1.1),
    ):
        key = {"power_law": None}.get(name, name + "_params")
        params = {"alpha": alpha} if name == "power_law" else mc[key]
        ax.plot(grid, p_tail * model_survival(grid, name, params, xmin),
                color=colour, ls=style, lw=lw * fs, zorder=4, solid_capstyle="round")
    # Poisson: the no-structure reference, fitted on the same tail (degree_poisson_fit.py).
    _po = p_tail * poisson_survival(grid, pois["poisson_lambda"], xmin)
    _floor = 1.0 / (len(obs) * 2.2)
    _po = np.where(_po >= _floor, _po, np.nan)
    ax.plot(grid, _po, color=PO, ls=(0, (2.6, 1.6)), lw=1.1 * fs, zorder=4)

    # --- right axis: cumulative percentage of regulators -------------------------------
    # Compact mode drops this axis. It is 1 - CCDF, i.e. the same series read the other way,
    # and at one third of a page row it consumes the only two label-free corners the panel
    # has. The full version keeps it for fidelity to the Alam-2017 design.
    k90 = float(np.percentile(obs, 90))
    if not compact:
        ax2 = ax.twinx()
        ax2.plot(xo, 100.0 * (1.0 - yo + 1.0 / len(xo)), lw=1.0 * fs, color=CUM, alpha=0.75,
                 solid_capstyle="round", zorder=1)
        ax2.set_ylim(0, 100)
        ax2.set_ylabel("cumulative % of regulators", color=CUM, fontsize=8 * fs, labelpad=1.0)
        ax2.tick_params(axis="y", colors=CUM, labelsize=8 * fs)
        ax2.spines["right"].set_visible(True)
        ax2.spines["right"].set_color(CUM)
        ax2.spines["top"].set_visible(False)
        ax2.axhline(90, color=CUM, lw=0.7 * fs, ls=(0, (3, 3)), zorder=0)
        ax2.axvline(k90, color=CUM, lw=0.7 * fs, ls=(0, (3, 3)), zorder=0)

    # --- axes -------------------------------------------------------------------------
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.8, obs.max() * 1.6)
    ax.set_ylim(1.0 / (len(obs) * 2.2), 1.4)
    ax.set_xlabel("regulator out-degree $k$", fontsize=9 * fs)
    ax.set_ylabel(r"$P(K \geq k)$", fontsize=9 * fs)
    ax.tick_params(labelsize=8 * fs)

    # Compact mode labels the two empirical series in place. A six-entry key is physically
    # wider than a third-of-a-row panel and covers the low-k plateau wherever it is anchored;
    # the four fitted curves are named in the legend text instead.
    if compact:
        ax.annotate("census", xy=(obs[len(obs) // 2], 0.5), xytext=(2.0, 0.30),
                    fontsize=6.6 * fs, color=INK,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=INK, shrinkA=1, shrinkB=2))
        ax.annotate("random\nnetwork", xy=(float(np.median(ran)), 0.5),
                    xytext=(float(np.median(ran)) * 1.25, 0.85), fontsize=6.6 * fs, color=GREY,
                    ha="left", va="top", linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY, shrinkA=1, shrinkB=2))

    note = (f"census  {_exp_label(alpha)}, GOF $p$={fit['gof_p']:.2f}\n"
            f"random  {_exp_label(rnd['random']['alpha'])}, GOF $p$ < 0.001")
    if not compact:
        obs_gamma = f"γ > {ALPHA_BOUND:.0f}" if alpha >= ALPHA_BOUND - 1e-3 else f"γ = {alpha:.2f}"
        random_alpha = float(rnd["random"]["alpha"])
        random_gamma = (f"γ > {ALPHA_BOUND:.0f}" if random_alpha >= ALPHA_BOUND - 1e-3
                        else f"γ = {random_alpha:.2f}")
        note = (f"census  {obs_gamma}, tail k ≥ {xmin}, GOF p={fit['gof_p']:.2f}\n"
                f"random  {random_gamma}, GOF p < 0.001\n"
                f"90% of regulators have k ≤ {k90:.0f}")
    ax.text(0.03, 0.03, note, transform=ax.transAxes, fontsize=(6.6 if compact else 7) * fs,
            color=INK, va="bottom", ha="left", linespacing=1.5)

    if tag:
        panel_tag(ax, tag, dx=-0.18, dy=1.08)
    return dict(alpha=alpha, xmin=xmin, gof_p=fit["gof_p"],
                random_alpha=rnd["random"]["alpha"], random_gof_p=rnd["random"]["gof_p"],
                k90=k90)


def legend_handles(fs=1.0):
    from matplotlib.lines import Line2D
    return [
        Line2D([], [], ls="none", marker="o", ms=3.2 * fs, mfc=INK, mec="none", label="census"),
        Line2D([], [], ls="none", marker="o", ms=3.2 * fs, mfc="none", mec=GREY, mew=0.7,
               label="random network"),
        Line2D([], [], color=PL, lw=1.5 * fs, label="power law"),
        Line2D([], [], color=LN, lw=1.1 * fs, ls=(0, (5, 2)), label="lognormal"),
        Line2D([], [], color=EX, lw=1.1 * fs, ls=(0, (1.4, 1.6)), label="exponential"),
        Line2D([], [], color=CO, lw=1.1 * fs, ls=(0, (6, 2, 1.4, 2)), label="power law + cutoff"),
        Line2D([], [], color=PO, lw=1.1 * fs, ls=(0, (2.6, 1.6)), label="Poisson"),
    ]




# --- fig_census_network_full.py ---
#!/usr/bin/env python3
"""The whole sign-reversal census as one picture: all 4,379 reversing edges.

Figure 2c draws the 532-edge high-magnitude core because that is the largest subset in
which individual genes can still be labelled. This panel is the opposite trade: nothing is
labelled except the largest hubs, and in exchange the reader sees the census at full size.
It is a scale-and-shape panel, not a lookup table, and it exists so a reader can see that
the core Figure 2c dissects is a slice of something much bigger.

Same visual grammar as Figure 2c on purpose -- bipartite two-ring, regulators on the inner
ring, targets on the outer, red for effects lost on stimulation and teal for gained, gold for
the T-cell effector genes -- so the two figures read as the same object at two zoom levels.

Layout note: a force-directed layout collapses this graph (it is bipartite and strongly
convergent, so springs pull it into a blob). The two-ring with barycentric ordering is what
makes the convergence visible: targets pulled by the same regulators end up adjacent, so
convergence shows as bundles of chords landing together rather than as uniform hairball.

Reads  data/census/sign_inverting_edges.csv   (canonical 4,379)
Writes figures/fig_census_network_full.{png,pdf}
"""








matplotlib.use("Agg")










EFFECTOR = {"GZMA", "GZMB", "GNLY", "PRF1", "IFNG", "TNF", "NKG7", "CCL5", "CCL4", "CCL3",
            "CCR5", "CCR1", "CCR7", "CCR2", "IL2", "IL2RA", "IL2RB", "CD69", "CTLA4",
            "IL2RG", "SELL", "TBX21", "KLRB1", "XCL1", "XCL2"}

BARYCENTRE_ITERS = 60
R_IN, R_OUT = 0.58, 1.00

# Frozen label set shown by the named complete-census figure.  The renderer first takes the
# 14 largest reversal out-degrees, then suppresses labels closer than 0.075 radians to an
# already retained label.  Keeping the resolved tuple explicit prevents a later layout or
# degree-order change from silently changing the labels in Figure 1a.
COMPLETE_CENSUS_REGULATOR_LABELS = (
    "CCNC", "TADA2B", "ARNT", "SGF29", "AHR",
    "NCKAP1L", "ATAD5", "ARHGAP30", "ELOF1", "SENP5",
)

EDGE_SOURCE = ROOT / "data/census/sign_inverting_edges.csv"
EXPECTED_EDGE_SOURCE_SHA256 = "48430579f0974c9ff36d5e7d089710cfba97c5b5313e6c2838469de035c2b484"
with open(EDGE_SOURCE, "rb") as _source_handle:
    _source_sha256 = hashlib.sha256(_source_handle.read()).hexdigest()
if _source_sha256 != EXPECTED_EDGE_SOURCE_SHA256:
    raise RuntimeError(
        "complete-census edge source changed: "
        f"expected {EXPECTED_EDGE_SOURCE_SHA256}, observed {_source_sha256}"
    )
d = pd.read_csv(EDGE_SOURCE)
d["lost"] = d["lfc_Rest"] > 0          # positive at rest -> effect lost on stimulation
d["minmag"] = np.minimum(d["lfc_Rest"].abs(), d["lfc_Stim48"].abs())
assert len(d) == 4_379
assert not d[["regulator", "target"]].duplicated().any()
assert int(d["lost"].sum()) == 2_298

regs_l = sorted(set(d["regulator"]))
tars_l = sorted(set(d["target"]))
outdeg = d["regulator"].value_counts().to_dict()
indeg = d["target"].value_counts().to_dict()

# --- bipartite two-ring with barycentric crossing reduction (as Figure 2c) ---
r_ang = {g: 2 * np.pi * i / len(regs_l) for i, g in enumerate(regs_l)}
t_ang = {t: 0.0 for t in tars_l}
reg_targets, tar_regs = {g: [] for g in regs_l}, {t: [] for t in tars_l}
for u, v in zip(d["regulator"], d["target"]):
    reg_targets[u].append(v)
    tar_regs[v].append(u)


def cmean(angs):
    return float(np.angle(np.mean([np.exp(1j * a) for a in angs])))


for _ in range(BARYCENTRE_ITERS):
    for t in tars_l:
        t_ang[t] = cmean([r_ang[u] for u in tar_regs[t]])
    for g in regs_l:
        r_ang[g] = cmean([t_ang[v] for v in reg_targets[g]])

reg_order = sorted(regs_l, key=lambda g: r_ang[g])
tar_order = sorted(tars_l, key=lambda t: t_ang[t])
r_ang = {g: 2 * np.pi * i / len(reg_order) for i, g in enumerate(reg_order)}
t_ang = {t: 2 * np.pi * i / len(tar_order) for i, t in enumerate(tar_order)}
# 99 genes are BOTH regulator and target. Keep a separate coordinate per ROLE, so a
# dual-role gene appears once on each ring rather than having its regulator position
# silently overwritten by its target position (which would launch its out-edges from the
# wrong ring).
rpos = {g: (R_IN * np.cos(r_ang[g]), R_IN * np.sin(r_ang[g])) for g in regs_l}
tpos = {t: (R_OUT * np.cos(t_ang[t]), R_OUT * np.sin(t_ang[t])) for t in tars_l}


def panel_census_network(ax, tag: str = "", alpha: float = 0.16, lw: float = 0.30,
                         n_reg_labels: int = 14, fs: float = 1.0,
                         legend: bool = True, labels: bool = True,
                         reg_uniform: bool = False, short_legend: bool = False,
                         lab_pitch: float = 0.028,
                         hi_mag: float | None = None, hi_alpha: float = 0.50,
                         hi_lw: float = 1.7, lo_scale: float = 0.90,
                         dual_role: bool = False) -> None:
    """Draw the census network.

    alpha/lw MUST come down when the panel is drawn smaller. The same 4,379 chords cross far
    more densely per unit area on a small axes, so transparency that reads as tone at full
    page stacks into opaque colour: at 4.2 in the default settings turn the whole graph into a
    solid red/teal ring and the convergence bundles disappear. Roughly halve alpha when
    halving the panel's linear size. fs scales the label type with the panel.

    hi_mag marks the high-magnitude subgraph inside the complete census instead of redrawing
    it as a separate panel. Chords whose weaker endpoint exceeds hi_mag log2 units are drawn at
    hi_alpha and hi_lw x lw; every other chord keeps its colour but drops to lo_scale x alpha.
    Emphasis therefore rides on opacity alone, leaving the red/teal channel to carry direction
    as it does everywhere else in the paper. Passing None restores the flat census view.

    dual_role colours BOTH copies of the genes that are both regulator and target. Such a gene
    already has a coordinate on each ring (see the rpos/tpos note above), so marking only the
    outer copy would leave its inner square looking unrelated. The outer circle takes a purple
    fill, the inner square a purple edge, in the palette colour reserved for node attributes so
    the red/teal direction channel is untouched. Effector targets keep their gold fill and take
    the purple edge instead, because effector status is the stronger label on that node.

    labels=False drops every gene label, regulator and target alike. reg_uniform=True draws
    all regulator markers at one size instead of scaling them by out-degree. Both are for the
    Figure 1 panel, which exists to show the shape of the census rather than to be read off:
    out-degree there would be a second encoding of the same quantity the chord density already
    shows, and no individual gene is being pointed at.
    """
    strong = (d["minmag"] > hi_mag).to_numpy() if hi_mag is not None else np.zeros(len(d), bool)
    ts = np.linspace(0, 1, 14)
    # Faint pass first, emphasised pass second: matplotlib draws in call order within a zorder,
    # and the 532 marked chords must not be buried under the 3,847 they are drawn among.
    for emphasise in (False, True) if hi_mag is not None else (False,):
        a = hi_alpha if emphasise else (alpha * lo_scale if hi_mag is not None else alpha)
        w = lw * hi_lw if emphasise else lw
        z = 1.6 if emphasise else 1
        for u, v, lost, hot in zip(d["regulator"], d["target"], d["lost"], strong):
            if hi_mag is not None and hot != emphasise:
                continue
            x1, y1 = rpos[u]
            x2, y2 = tpos[v]
            mx, my = (x1 + x2) * 0.30, (y1 + y2) * 0.30
            bx = (1 - ts) ** 2 * x1 + 2 * (1 - ts) * ts * mx + ts ** 2 * x2
            by = (1 - ts) ** 2 * y1 + 2 * (1 - ts) * ts * my + ts ** 2 * y2
            ax.plot(bx, by, "-", lw=w, alpha=a,
                    color=(PAL["LOST"] if lost else PAL["GAINED"]), zorder=z)

    # Targets are drawn in three passes so no class can be hidden by a neighbour landing on
    # top of it: plain grey first, then the gold effector landmarks, then the dual-role nodes
    # last of all. Order matters because adjacent ring positions overlap at this density.
    dual = set(regs_l) & set(tars_l) if dual_role else set()
    tsize = lambda g: (min(60, 12 + 3.2 * indeg[g]) if g in EFFECTOR
                       else min(26, 1.6 + 1.5 * indeg[g]))
    plain = [g for g in tars_l if g not in EFFECTOR and g not in dual]
    teff = [g for g in tars_l if g in EFFECTOR and g not in dual]
    ax.scatter([tpos[g][0] for g in plain], [tpos[g][1] for g in plain],
               s=[tsize(g) for g in plain], color="#9aa3ad", lw=0, zorder=2)
    ax.scatter([tpos[g][0] for g in teff], [tpos[g][1] for g in teff],
               s=[tsize(g) for g in teff], color=PAL["HL"],
               edgecolor="#7a5b12", lw=0.6, zorder=4)
    if dual:
        # A dual-role effector keeps its gold fill and takes the purple edge, because effector
        # status is the stronger label on that node; every other dual-role node is filled purple.
        dl = [g for g in tars_l if g in dual]
        ax.scatter([tpos[g][0] for g in dl], [tpos[g][1] for g in dl],
                   s=[tsize(g) for g in dl],
                   color=[PAL["HL"] if g in EFFECTOR else PAL["CORE"] for g in dl],
                   edgecolor=PAL["CORE"],
                   lw=[1.1 if g in EFFECTOR else 0.0 for g in dl], zorder=6)

    reg_s = ([14.0] * len(regs_l) if reg_uniform
             else [min(110, 6 + 2.4 * outdeg[g]) for g in regs_l])
    ax.scatter([rpos[g][0] for g in regs_l], [rpos[g][1] for g in regs_l],
               s=reg_s, marker="s", color="#33475b", edgecolor="white", lw=0.35, zorder=3)
    if dual:
        # Purple sits immediately outside the white stroke, flush against it, so the square keeps
        # the neighbour separation the white edge gives it and the two read as one border.
        # Marker `s` is an AREA in points squared, so the flush size has to be computed in linear
        # units and squared back: both strokes are centred on their path, so the purple ring's
        # side must exceed the square's by the white stroke plus the purple stroke.
        RW, PW = 0.35, 0.55                       # white and purple stroke widths, points
        dr = [g for g in regs_l if g in dual]
        ax.scatter([rpos[g][0] for g in dr], [rpos[g][1] for g in dr],
                   s=[(math.sqrt(min(110, 6 + 2.4 * outdeg[g])) + RW + PW) ** 2 for g in dr],
                   marker="s", facecolor="none", edgecolor=PAL["CORE"], lw=PW, zorder=3.5)

    if not labels:
        _finish(ax, alpha, lw, fs, legend, tag, short_legend, hi_mag, dual_role)
        return

    # label only what stays legible: the biggest regulators, and the effector genes
    shown_ang = []
    shown_regulators = []
    for g in sorted(regs_l, key=lambda x: -outdeg[x])[:n_reg_labels]:
        a = r_ang[g]
        if any(abs(a - b) < 0.075 or abs(abs(a - b) - 2 * np.pi) < 0.075 for b in shown_ang):
            continue          # SAGA subunits sit adjacent; label the largest, skip the pile
        shown_ang.append(a)
        shown_regulators.append(g)
        # Placed in the annulus between the rings and anchored to grow OUTWARD. Two earlier
        # placements failed: ha="center" on the ring put half of each name over the dark squares
        # (CCNC read as "CNC"), and growing inward from the ring filled the centre the legend
        # occupies. The annulus carries only faint chords, so bold type stays readable on it.
        ax.text((R_IN + 0.055) * np.cos(a), (R_IN + 0.055) * np.sin(a), g,
                fontsize=5.6 * fs, ha=("left" if np.cos(a) >= 0 else "right"),
                va="center", color="#33475b", fontweight="bold",
                rotation=np.degrees(a) + (180 if np.cos(a) < 0 else 0),
                rotation_mode="anchor", zorder=6)
    if n_reg_labels == 14:
        assert tuple(shown_regulators) == COMPLETE_CENSUS_REGULATOR_LABELS, (
            "complete-census regulator-label contract changed: "
            f"{tuple(shown_regulators)!r}"
        )
    # Crowded effector labels are separated in y, at their own radius. An earlier version
    # staggered them radially over three rings 0.16 apart, which pushed the axes limits out
    # to +-1.34 and shrank the ring itself: the labels were robbing the circle of space.
    # Horizontal type instead of radial is the other half of that: near the horizontal axis
    # the two orientations cost the same, but a radial label at the TOP of the ring spends
    # its whole length on the y extent (CTLA4 alone forced +-1.21), whereas horizontal type
    # costs one line height there. lab_pitch is that line height in DATA units, which only
    # the caller knows -- it depends on the axes size and the limits it is about to set.
    rr = R_OUT + 0.030
    pos = {t: (rr * np.cos(t_ang[t]), rr * np.sin(t_ang[t])) for t in teff}
    for side in (1, -1):
        col = sorted([t for t in teff if np.sign(pos[t][0]) == side], key=lambda t: pos[t][1])
        i = 0
        while i < len(col):
            j = i
            while j + 1 < len(col) and pos[col[j + 1]][1] - pos[col[j]][1] < lab_pitch:
                j += 1
            if j > i:          # spread the run about its own centroid, vertical order kept
                grp = col[i:j + 1]
                c = float(np.mean([pos[t][1] for t in grp]))
                for k, t in enumerate(grp):
                    pos[t] = (pos[t][0], c + (k - (len(grp) - 1) / 2) * lab_pitch)
            i = j + 1
    for t in teff:
        x, y = pos[t]
        ax.text(x, y, t, fontsize=5.4 * fs, ha=("left" if x >= 0 else "right"),
                va="center", color="#7a5b12", fontweight="bold", zorder=6)

    _finish(ax, alpha, lw, fs, legend, tag, short_legend, hi_mag, dual_role)


def _finish(ax, alpha, lw, fs, legend, tag, short_legend=False, hi_mag=None,
            dual_role=False) -> None:
    """Axes limits, legend and panel tag. Shared by the labelled and unlabelled paths."""
    ax.set_aspect("equal")
    ax.set_xlim(-1.30, 1.30)
    ax.set_ylim(-1.24, 1.24)
    ax.axis("off")
    n_lost = int(d["lost"].sum())
    if not legend:
        if tag:
            panel_tag(ax, tag, dx=-0.02, dy=1.00)
        return
    # The key sits inside the inner ring, so its box has to fit within R_IN. In the A4
    # composite the 7 pt physical floor widens every string, and the full labels then run out
    # over the regulator names. short_legend trims them to what still fits.
    if short_legend:
        texts = (f"lost on stimulation ({n_lost:,})", f"gained ({len(d)-n_lost:,})",
                 f"regulator ({len(regs_l)})", f"target ({len(tars_l):,})",
                 "fixed effector landmark")
    else:
        texts = (f"effect lost on stimulation ({n_lost:,})", f"effect gained ({len(d)-n_lost:,})",
                 f"regulator ({len(regs_l)})", f"target ({len(tars_l):,})",
                 "T-cell effector / activation gene")
    handles = [
        Line2D([0], [0], color=PAL["LOST"], lw=2.4, label=texts[0]),
        Line2D([0], [0], color=PAL["GAINED"], lw=2.4, label=texts[1]),
        Line2D([0], [0], marker="s", ls="", mfc="#33475b", mec="white", ms=7, label=texts[2]),
        Line2D([0], [0], marker="o", ls="", mfc="#9aa3ad", mec="none", ms=6, label=texts[3]),
        Line2D([0], [0], marker="o", ls="", mfc=PAL["HL"], mec="#7a5b12", ms=7, label=texts[4]),
    ]
    if hi_mag is not None:
        # Counted from the data, never a literal: the emphasised set must match what was drawn.
        n_hi = int((d["minmag"] > hi_mag).sum())
        handles.insert(2, Line2D([0], [0], color="#5b6570", lw=2.6,
                                 label=(f"bold: |log2FC| > {hi_mag:g} at both endpoints "
                                        f"({n_hi:,})" if not short_legend else
                                        f"bold: > {hi_mag:g} both states ({n_hi:,})")))
    if dual_role:
        n_dual = len(set(regs_l) & set(tars_l))          # counted, never a literal
        handles.append(Line2D([0], [0], marker="o", ls="", mfc=PAL["CORE"],
                              mec=PAL["CORE"], ms=6,
                              label=(f"both regulator and target ({n_dual})"
                                     if short_legend else
                                     f"gene that is both regulator and target ({n_dual})")))
    leg = ax.legend(handles=handles, loc="center", fontsize=6.4 * fs,
        frameon=hi_mag is not None,
        labelspacing=0.55 if short_legend else 0.75,
        handletextpad=0.5 if short_legend else 0.8,
        bbox_to_anchor=(0.5, 0.5))
    if hi_mag is not None:
        # Emphasised chords cross the inner ring, so the key needs its own ground to stay
        # legible; without hi_mag the faint chords never obscure it and the frame stays off.
        leg.get_frame().set(facecolor="white", edgecolor="none", alpha=0.58)
        leg.set_zorder(6)
    if tag:
        panel_tag(ax, tag, dx=-0.02, dy=1.00)




# --- Main Layout ---
#!/usr/bin/env python3
"""Versioned Figure 1 layout candidate with a larger census ring.

Relative to the promoted compact arrangement, panel a is enlarged while panels b-e
are narrowed.  Panels b-e use exactly 5 pt internal text and exactly 7 pt titles/tags.
All data, panel functions, statistics, denominators, axis limits and panels f-h are retained.
"""











ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
STEM = "Fig1"



def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Fail closed before rendering if any authoritative dependency, input, predecessor
# artifact or canonical predecessor has drifted.
setup(scale=0.82)

assert c["edges_sig_both"] == 83_489
# not by removing panel content or uniformly shrinking an A4 composition.

"""Is the regulator out-degree distribution heavy-tailed, and is 'power law' defensible?

WHICH NETWORK. Not the drawn ones. Repository figure R5 is the 532-edge high-magnitude core, whose
degrees are set by the |log2FC| > 0.5 threshold. Figure 1d is the reversing subgraph, where a
regulator's degree is HOW MANY OF ITS EDGES REVERSED, i.e. the outcome variable rather than
its connectivity. The object with actual regulator connectivity is the FULL CENSUS edge table
after the paper-wide Ensembl-ID self-edge mask: 83,489 both-significant trans edges over 617
regulators,
of which both drawn networks are subgraphs. The reversing subgraph is fitted too, but only as
a comparison.

WHAT THIS CAN AND CANNOT CLAIM. The 617 census regulators are not all perturbed genes. An
upstream gate retained 618 perturbations with significant on-target knockdown and at least 50
significant downstream targets in EACH state separately; CARD6 contributed no edge significant
in both endpoint states, leaving 617 regulators represented in the shared-edge census. Thus the
distribution is conditioned on signal-based selection, and a conclusion here is about the shape
over the observed census range, never about the T-cell regulatory network being scale-free.

METHOD (Clauset, Shalizi & Newman 2009), because a straight line on a log-log plot is not
evidence and most "scale-free" claims in biology do not survive this procedure:

  1. discrete power-law MLE for alpha at each candidate xmin, p(k) proportional to
     k^-alpha / zeta(alpha, xmin)
  2. xmin chosen to minimise the KS distance between the empirical and fitted CDFs
  3. GOODNESS OF FIT by semi-parametric bootstrap: synthesise datasets that are power law
     above xmin and resampled-empirical below it, refit each end to end, and take
     p = fraction of synthetic KS >= observed KS. p < 0.1 REJECTS the power law. A large p
     means "not ruled out", never "confirmed".
  4. RELATIVE FIT on the same selected tail against a discretised, tail-conditioned
     lognormal; a shifted geometric (the discrete exponential); and a power law with an
     exponential cutoff. Vuong tests are used only for the two non-nested comparisons.
     The cutoff model contains the pure power law at a boundary, so it is compared
     descriptively by AICc rather than passed through an invalid ordinary Vuong test.

Reads   data/census/census_master_edges.csv.gz  (validated master; 83,489 trans edges)
Writes  data/census/degree_distribution_test.json
        figures/FigS_degree_distribution.{pdf,png}
"""

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp, zeta
from scipy.stats import lognorm, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S

PAPER = Path(__file__).resolve().parents[1]
RD = PAPER / "data"
OUT = PAPER / "figures"

SEED = 20260813
N_BOOT = 5000          # synthetic datasets for the goodness-of-fit p
ALPHA_BOUNDS = (1.01, 15.0)
ANALYSIS_VERSION = 4       # v4 reports package-compatible raw R separately from Vuong z


# ----------------------------------------------------------- discrete power law

def _nll_powerlaw(alpha: float, x: np.ndarray, xmin: int) -> float:
    if not (ALPHA_BOUNDS[0] < alpha < ALPHA_BOUNDS[1]):
        return np.inf
    return len(x) * np.log(zeta(alpha, xmin)) + alpha * np.sum(np.log(x))


def fit_alpha(x: np.ndarray, xmin: int) -> float:
    r = minimize_scalar(_nll_powerlaw, bounds=ALPHA_BOUNDS, args=(x, xmin), method="bounded")
    return float(r.x)


def pl_cdf(k: np.ndarray, alpha: float, xmin: int) -> np.ndarray:
    """P(X <= k) for the discrete power law, via the Hurwitz zeta tail."""
    return 1.0 - zeta(alpha, k + 1) / zeta(alpha, xmin)


def ks_distance(x: np.ndarray, alpha: float, xmin: int) -> float:
    """Exact discrete KS distance, evaluated before and after every tied empirical jump."""
    vals, counts = np.unique(x, return_counts=True)
    empirical_after = np.cumsum(counts) / len(x)
    empirical_before = (np.cumsum(counts) - counts) / len(x)
    fitted_after = pl_cdf(vals, alpha, xmin)
    fitted_before = pl_cdf(vals - 1, alpha, xmin)
    return float(max(np.max(np.abs(empirical_after - fitted_after)),
                     np.max(np.abs(empirical_before - fitted_before))))


def choose_xmin(data: np.ndarray) -> tuple[int, float, float, int]:
    cand = np.unique(data)
    cand = cand[cand >= 1]
    cand = cand[:-1]
    best = (None, None, np.inf, 0)
    for xm in cand:
        tail = data[data >= xm]
        if len(tail) < 50:
            continue
        a = fit_alpha(tail, int(xm))
        d = ks_distance(tail, a, int(xm))
        if d < best[2]:
            best = (int(xm), a, d, len(tail))
    return best


def sample_powerlaw(n: int, alpha: float, xmin: int, rng) -> np.ndarray:
    """Exact unbounded discrete-power-law rejection sampler.

    The proposal rounds a continuous Pareto variate with lower bound xmin - 1/2. Convexity
    supplies a finite rejection bound, avoiding the silent finite-tail truncation produced by
    sampling from an arbitrarily capped probability table.
    """
    if n == 0:
        return np.empty(0, dtype=np.int64)
    lower = xmin - 0.5
    normalizer = zeta(alpha, xmin)
    bound = 1.0 / (normalizer * (alpha - 1.0) * lower ** (alpha - 1.0))
    chunks, have = [], 0
    while have < n:
        m = max(64, int((n - have) * 1.05) + 8)
        u = rng.random(m)
        y = lower * np.exp(-np.log1p(-u) / (alpha - 1.0))
        k = np.floor(y + 0.5).astype(np.int64)
        lo = k - 0.5
        q = (lower ** (alpha - 1.0) * lo ** (1.0 - alpha)
             * (-np.expm1((1.0 - alpha) * np.log((k + 0.5) / lo))))
        p = k.astype(float) ** (-alpha) / normalizer
        keep = rng.random(m) < np.minimum(1.0, p / (bound * q))
        accepted = k[keep]
        chunks.append(accepted)
        have += len(accepted)
    return np.concatenate(chunks)[:n]


def gof_pvalue(data: np.ndarray, alpha: float, xmin: int, d_obs: float, rng,
               n_boot: int = N_BOOT) -> tuple[float, int, int]:
    """Semi-parametric bootstrap. Each synthetic set is refitted END TO END, xmin included,
    which is what makes the p-value account for xmin having been chosen from the data."""
    below = data[data < xmin]
    n, n_tail = len(data), int((data >= xmin).sum())
    p_tail = n_tail / n
    worse, valid = 0, 0
    for _ in range(n_boot):
        n_t = rng.binomial(n, p_tail)
        tail = sample_powerlaw(n_t, alpha, xmin, rng)
        rest = rng.choice(below, size=n - n_t, replace=True) if len(below) else np.array([], int)
        syn = np.concatenate([tail, rest]).astype(int)
        xm_s, a_s, d_s, _ = choose_xmin(syn)
        if xm_s is None:
            continue
        valid += 1
        worse += int(d_s >= d_obs)
    # Finite-Monte-Carlo correction: a simulated p-value is never reported as exactly zero.
    return (worse + 1) / (valid + 1), worse, valid


# ------------------------------------------------------- alternative tail models

def loglik_powerlaw(x, alpha, xmin):
    return -alpha * np.log(x) - np.log(zeta(alpha, xmin))


def loglik_lognormal(x, mu, sig, xmin):
    """Integer-bin lognormal mass, conditioned on K >= xmin."""
    upper = lognorm.cdf(x + 0.5, s=sig, scale=np.exp(mu))
    lower = lognorm.cdf(np.maximum(x - 0.5, 1e-12), s=sig, scale=np.exp(mu))
    norm_tail = lognorm.sf(xmin - 0.5, s=sig, scale=np.exp(mu))
    return np.log(np.maximum((upper - lower) / max(norm_tail, 1e-300), 1e-300))


def loglik_exponential(x, q, xmin):
    """Discrete exponential, i.e. geometric mass on K-xmin = 0,1,..."""
    return np.log1p(-q) + (x - xmin) * np.log(q)


def _cutoff_log_normalizer(alpha: float, lam: float, xmin: int, max_observed: int) -> float:
    """Log sum of k^-alpha exp(-lambda*k), k=xmin..infinity, to numerical precision."""
    # exp(-lambda*k) makes a finite numerical grid safe. Thirty-five exponential e-folds
    # beyond xmin leaves < 1e-15 of an already upper-bounded geometric remainder. The
    # observed maximum is always retained, and the cap is merely a guard against pathological
    # optimiser proposals near lambda=0 (the pure PL handles that boundary separately).
    upper = max(max_observed * 2, xmin + 100, int(xmin + 35.0 / max(lam, 1e-9)))
    upper = min(upper, 5_000_000)
    support = np.arange(xmin, upper + 1, dtype=float)
    return float(logsumexp(-alpha * np.log(support) - lam * support))


def loglik_cutoff_powerlaw(x, alpha, lam, xmin):
    log_z = _cutoff_log_normalizer(alpha, lam, xmin, int(np.max(x)))
    return -alpha * np.log(x) - lam * x - log_z


def fit_cutoff_powerlaw(x: np.ndarray, xmin: int) -> tuple[float, float, np.ndarray]:
    """Fit exact discrete p(k) proportional to k^-alpha exp(-lambda*k)."""
    x = np.asarray(x, float)

    def nll(theta):
        alpha, lam = np.exp(theta)
        return -float(loglik_cutoff_powerlaw(x, alpha, lam, xmin).sum())

    starts = [(1.5, 1e-3), (0.05, 0.1), (2.0, 1e-2)]
    fits = [minimize(nll, np.log(start), method="Nelder-Mead",
                     options={"maxiter": 5000, "xatol": 1e-10, "fatol": 1e-9})
            for start in starts]
    opt = min(fits, key=lambda z: z.fun)
    if not opt.success:
        raise RuntimeError(f"cutoff power-law fit failed: {opt.message}")
    alpha, lam = np.exp(opt.x)
    ll = loglik_cutoff_powerlaw(x, float(alpha), float(lam), xmin)
    return float(alpha), float(lam), ll


def aicc(log_likelihood: float, n: int, n_parameters: int) -> float:
    return float(-2 * log_likelihood + 2 * n_parameters
                 + 2 * n_parameters * (n_parameters + 1) / (n - n_parameters - 1))


def fit_tail_models(x: np.ndarray, alpha: float, xmin: int) -> dict:
    """Fit all alternatives on the same integer support and tail conditioning as the PL."""
    x = np.asarray(x, float)
    lx = np.log(x)

    def nll_ln(theta):
        mu, log_sig = theta
        return -float(loglik_lognormal(x, mu, np.exp(log_sig), xmin).sum())

    opt = minimize(nll_ln, [float(lx.mean()), np.log(max(float(lx.std(ddof=1)), 0.1))],
                   method="Nelder-Mead", options={"maxiter": 5000, "xatol": 1e-10,
                                                   "fatol": 1e-10})
    if not opt.success:
        raise RuntimeError(f"discrete lognormal fit failed: {opt.message}")
    mu, sig = float(opt.x[0]), float(np.exp(opt.x[1]))

    mean_excess = float((x - xmin).mean())
    q = mean_excess / (1.0 + mean_excess)

    ll_pl = loglik_powerlaw(x, alpha, xmin)
    ll_ln = loglik_lognormal(x, mu, sig, xmin)
    ll_ex = loglik_exponential(x, q, xmin)
    cutoff_alpha, cutoff_lambda, ll_cutoff = fit_cutoff_powerlaw(x, xmin)
    log_likelihoods = {
        "power_law": float(ll_pl.sum()),
        "lognormal": float(ll_ln.sum()),
        "exponential": float(ll_ex.sum()),
        "cutoff_power_law": float(ll_cutoff.sum()),
    }
    aicc_values = {
        "power_law": aicc(log_likelihoods["power_law"], len(x), 1),
        "lognormal": aicc(log_likelihoods["lognormal"], len(x), 2),
        "exponential": aicc(log_likelihoods["exponential"], len(x), 1),
        "cutoff_power_law": aicc(log_likelihoods["cutoff_power_law"], len(x), 2),
    }
    best_aicc = min(aicc_values.values())
    out = {
        "lognormal_params": {"mu": mu, "sigma": sig},
        "exponential_params": {"q": q, "lambda": -float(np.log(q))},
        "cutoff_power_law_params": {"alpha": cutoff_alpha, "lambda": cutoff_lambda,
                                    "alpha_at_lower_bound": bool(cutoff_alpha < 1e-5)},
        "log_likelihood": log_likelihoods,
        "aicc": aicc_values,
        "delta_aicc": {k: float(v - best_aicc) for k, v in aicc_values.items()},
    }
    for name, ll in [("lognormal", ll_ln), ("exponential", ll_ex)]:
        raw_R, z, p = vuong(ll_pl, ll)
        out[f"power_law_vs_{name}"] = {
            # powerlaw.Fit.distribution_compare() returns this raw summed log-likelihood
            # ratio by default. The normalised statistic is reported separately because it
            # is the quantity used to calculate the non-nested two-sided p-value.
            "R": float(raw_R), "vuong_z": float(z), "p": float(p),
            "favours": ("power law" if raw_R > 0 else name)
            if p < 0.05 else "indistinguishable"}
    out["power_law_vs_cutoff_power_law"] = {
        "delta_log_likelihood": float(ll_pl.sum() - ll_cutoff.sum()),
        "test": ("not tested: cutoff power law contains the pure power law at lambda=0; "
                 "relative support is reported by AICc"),
    }
    return out


def add_global_vuong_fdr(results: dict) -> None:
    """BH-adjust the six non-nested PL-vs-lognormal/geometric comparisons together."""
    tests = []
    keys = ["regulator_outdegree_full_census", "target_indegree_full_census",
            "regulator_outdegree_reversing_subgraph"]
    for result_key in keys:
        for comparison in ("power_law_vs_lognormal", "power_law_vs_exponential"):
            item = results[result_key]["model_comparison"][comparison]
            tests.append((result_key, comparison, float(item["p"])))
    order = np.argsort([x[2] for x in tests])
    adjusted = np.empty(len(tests), float)
    running = 1.0
    for rank_index in range(len(tests) - 1, -1, -1):
        idx = int(order[rank_index])
        rank = rank_index + 1
        running = min(running, tests[idx][2] * len(tests) / rank)
        adjusted[idx] = running
    for (result_key, comparison, _), q in zip(tests, adjusted):
        results[result_key]["model_comparison"][comparison]["q_bh_across_six"] = float(q)


def upgrade_comparison_schema(results: dict) -> None:
    """Upgrade cached v3 comparisons that mislabeled the normalised z statistic as raw R."""
    keys = ["regulator_outdegree_full_census", "target_indegree_full_census",
            "regulator_outdegree_reversing_subgraph"]
    for result_key in keys:
        result = results[result_key]
        comparison = result["model_comparison"]
        n = int(result["n_tail"])
        ll_pl = float(comparison["log_likelihood"]["power_law"])
        for alternative in ("lognormal", "exponential"):
            item = comparison[f"power_law_vs_{alternative}"]
            if "vuong_z" in item:
                continue
            old_sample_sd_z = float(item["R"])
            raw_R = ll_pl - float(comparison["log_likelihood"][alternative])
            z = old_sample_sd_z * np.sqrt(n / (n - 1)) if n > 1 else 0.0
            p = float(2 * norm.sf(abs(z)))
            item.update({
                "R": float(raw_R),
                "vuong_z": float(z),
                "p": p,
                "favours": (("power law" if raw_R > 0 else alternative)
                             if p < 0.05 else "indistinguishable"),
            })
    results["analysis_version"] = ANALYSIS_VERSION


def model_survival(support: np.ndarray, name: str, params: dict, xmin: int) -> np.ndarray:
    """Conditional P(K >= k | K >= xmin) for each fitted discrete-tail model."""
    support = np.asarray(support, int)
    if name == "power_law":
        return zeta(params["alpha"], support) / zeta(params["alpha"], xmin)
    if name == "lognormal":
        return (lognorm.sf(support - 0.5, s=params["sigma"], scale=np.exp(params["mu"])) /
                lognorm.sf(xmin - 0.5, s=params["sigma"], scale=np.exp(params["mu"])))
    if name == "exponential":
        return params["q"] ** (support - xmin)
    if name == "cutoff_power_law":
        alpha, lam = params["alpha"], params["lambda"]
        upper = max(int(support.max()) * 2, xmin + 100,
                    int(xmin + 35.0 / max(lam, 1e-9)))
        upper = min(upper, 5_000_000)
        grid = np.arange(xmin, upper + 1, dtype=float)
        mass = np.exp(-alpha * np.log(grid) - lam * grid
                      - logsumexp(-alpha * np.log(grid) - lam * grid))
        survival = np.cumsum(mass[::-1])[::-1]
        return survival[support - xmin]
    raise ValueError(name)


def vuong(ll1: np.ndarray, ll2: np.ndarray) -> tuple[float, float, float]:
    """Raw log-likelihood ratio, normalised Vuong z, and two-sided p.

    Raw R > 0 favours model 1. This matches the default return value of the Python ``powerlaw``
    package. The p-value is calculated from z = R / sqrt(n * Var[d]); a large p means the two
    models are not distinguishable on this data, which is a real and common answer.
    """
    d = ll1 - ll2
    n = len(d)
    raw_R = float(d.sum())
    # Match powerlaw 1.5/Vuong's population variance convention for the pointwise
    # log-likelihood differences. This changes only the normalising statistic, not raw R.
    variance = float(np.mean((d - d.mean()) ** 2))
    if variance == 0:
        return raw_R, 0.0, 1.0
    z = raw_R / np.sqrt(n * variance)
    return raw_R, float(z), float(2 * norm.sf(abs(z)))


def compare_tail(x: np.ndarray, alpha: float, xmin: int) -> dict:
    return fit_tail_models(x, alpha, xmin)


def analyse(data: np.ndarray, label: str, rng, do_boot: bool = True) -> dict:
    xmin, alpha, d, n_tail = choose_xmin(data)
    at_bound = abs(alpha - ALPHA_BOUNDS[1]) < 1e-3 or abs(alpha - ALPHA_BOUNDS[0]) < 1e-3
    res = {"alpha_at_optimiser_bound": bool(at_bound), "label": label, "n": int(len(data)), "min": int(data.min()),
           "max": int(data.max()), "median": float(np.median(data)),
           "xmin": xmin, "alpha": float(alpha), "n_tail": n_tail,
           "ks_distance": float(d),
           "tail_fraction": round(n_tail / len(data), 4)}
    if do_boot:
        p, exceedances, valid_bootstraps = gof_pvalue(data, alpha, xmin, d, rng)
        bootstrap_mode = "serial; shared deterministic RNG stream"
        res["gof_p"] = round(p, 4)
        res["gof_exceedances"] = int(exceedances)
        res["gof_valid_bootstraps"] = int(valid_bootstraps)
        res["gof_execution"] = bootstrap_mode
        res["power_law_verdict"] = ("REJECTED (p < 0.1)" if p < 0.1
                                    else "not ruled out (this is NOT confirmation)")
    res["model_comparison"] = compare_tail(data[data >= xmin], alpha, xmin)

    # COMBINED VERDICT. Absolute goodness of fit and relative model comparison answer
    # different questions and can disagree, so both are reported.
    if do_boot:
        if res["gof_p"] < 0.1:
            res["verdict"] = "power law REJECTED by goodness of fit"
        else:
            res["verdict"] = "power-law tail compatible by GOF; this is not confirmation"
    return res


def ccdf(x: np.ndarray):
    """Tie-aware empirical survival P(X >= k) at each observed integer value."""
    xs, counts = np.unique(np.asarray(x, int), return_counts=True)
    return xs, np.cumsum(counts[::-1])[::-1] / counts.sum()


def main(figure_only: bool = False) -> None:
    """figure_only reuses data/census/degree_distribution_test.json so that layout
    edits do not re-run the 5,000-draw goodness-of-fit bootstrap."""
    rng = np.random.default_rng(SEED)
    # Everything comes from ONE file now: census_master_edges.csv.gz, built and validated by
    # build_census_master.py against all five published totals (84,044 edges, 617 regulators,
    # 4,379 reversals over 324 regulators, 79,665 same-sign, 2,298/2,081 lost/gained). An
    # earlier version of this figure read a stale kinetics_8hr_edges.csv.gz subset that omitted
    # MCAT and reported 322 reversing regulators. The corrected 8 h subset reports 323 because
    # KAT7 has no finite 8 h value; neither subset is admissible for the Rest--Stim48 topology.
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "direction", "is_self"])
    e = e[~e.is_self].copy()
    assert len(e) == 83_489 and int((e.kind == "rev").sum()) == 4_379
    out_all = e.groupby("regulator").size().to_numpy()
    in_all = e.groupby("target").size().to_numpy()
    rev = e[e.kind == "rev"]
    out_rev = rev.groupby("regulator").size().to_numpy()
    N_REV_REG = int(rev.regulator.nunique())            # 324
    no_outdeg = []                                      # nothing is missing any more

    cache = RD / "census/degree_distribution_test.json"
    if figure_only and cache.exists():
        candidate = json.loads(cache.read_text())
        if (candidate.get("analysis_version") == ANALYSIS_VERSION
                and candidate.get("n_bootstrap") == N_BOOT):
            res = candidate
            print("reusing cached fits (figure-only mode)")
        else:
            res = None
            print("cache predates the corrected discrete-tail analysis; recomputing")
    else:
        res = None
    if res is None:
        parts = {
            "regulator_outdegree_full_census": RD / "census/degree_distribution_test.regulator.json",
            "target_indegree_full_census": RD / "census/degree_distribution_test.target.json",
            "regulator_outdegree_reversing_subgraph": RD / "census/degree_distribution_test.reversal.json",
        }
        if figure_only and all(path.exists() for path in parts.values()):
            res = {"analysis_version": ANALYSIS_VERSION, "seed": SEED,
                   "n_bootstrap": N_BOOT,
                   "_source": ("data/census/census_master_edges.csv.gz, is_self=false "
                               "by endpoint Ensembl identity (83,489 trans census edges)"),
                   "_scope": ("618 perturbations passed the upstream on-target and >=50 "
                              "significant-target gate in each endpoint state; CARD6 contributed "
                              "no shared significant edge, leaving 617 represented census "
                              "regulators. This signal-selected census is not a sample of all "
                              "perturbed genes, so it cannot establish a whole-network "
                              "scale-free claim")}
            res.update({key: json.loads(path.read_text()) for key, path in parts.items()})
            print("assembled corrected single-family bootstrap results")
        else:
            res = {"analysis_version": ANALYSIS_VERSION, "seed": SEED,
                   "n_bootstrap": N_BOOT,
                   "_source": ("data/census/census_master_edges.csv.gz, is_self=false "
                               "by endpoint Ensembl identity (83,489 trans census edges)"),
                   "_scope": ("618 perturbations passed the upstream on-target and >=50 "
                              "significant-target gate in each endpoint state; CARD6 contributed "
                              "no shared significant edge, leaving 617 represented census "
                              "regulators. This signal-selected census is not a sample of all "
                              "perturbed genes, so it cannot establish a whole-network "
                              "scale-free claim"),
                   "regulator_outdegree_full_census": analyse(
                       out_all, "regulator out-degree", rng),
                   "target_indegree_full_census": analyse(
                       in_all, "target in-degree", rng),
                   "regulator_outdegree_reversing_subgraph": analyse(
                       out_rev, "reversals per regulator", rng)}
    upgrade_comparison_schema(res)
    add_global_vuong_fdr(res)
    # ------------------------------------------------------------------ figure
    S.setup()
    # Two plots above a full-width result strip remain legible when embedded at 6.5 inches.
    fig = plt.figure(figsize=(8.2, 7.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.00, 0.82],
                          left=0.085, right=0.985, bottom=0.125, top=0.94,
                          wspace=0.24, hspace=0.34)

    # ---- panel a: which regulators reverse at all, over the same degree axis ----
    # A regulator's out-degree is fixed by the census, so it cannot be split by direction the
    # way targets can. What CAN be shown is which regulators carry any reversal at all.
    ax = fig.add_subplot(gs[0, 0])
    r = res["regulator_outdegree_full_census"]
    rev_regs = set(rev.regulator.unique())     # canonical 324
    od = e.groupby("regulator").size()
    od_rev = od[od.index.isin(rev_regs)].to_numpy()
    od_non = od[~od.index.isin(rev_regs)].to_numpy()
    assert len(od_rev) == N_REV_REG, f"{len(od_rev)} != canonical {N_REV_REG}"
    for dat, col, lab in [(od_non, "#c7ccd1", f"no reversal ({len(od_non)} regulators)"),
                          (od_rev, "#33475b",
                           f"carries $\\geq$1 reversal ({len(od_rev)} regulators)")]:
        xs, cc = ccdf(dat)
        ax.loglog(xs, cc, "o", ms=2.6, color=col, alpha=0.75, label=lab)
    # The reversal-bearing split is descriptive: higher-degree regulators mechanically have
    # more opportunities to carry at least one reversal, so no rank-test p-value is reported.
    ax.set_xlabel("regulator out-degree $k$"); ax.set_ylabel("P(K $\\geq$ k)")
    ax.set_title("Regulator out-degree", fontsize=9.2)
    ax.legend(fontsize=7.0, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, -0.02))
    ax.text(0.02, 0.98,
            f"median $k$: {np.median(od_rev):.0f} (\u22651 reversal) vs "
            f"{np.median(od_non):.0f} (none)",
            transform=ax.transAxes, fontsize=7.4, va="top", color="#3d4348")
    S.panel_tag(ax, "a")

    # ---- panel b: targets split by the DIRECTION of the reversals they receive ----
    ax = fig.add_subplot(gs[0, 1])
    # labels straight from the master file, so the split is the SAME 2,298 / 2,081 shown in
    # Figure 1a rather than a re-derivation that could drift from it
    lost = rev[rev.direction == "lost_on_activation"].groupby("target").size()
    gain = rev[rev.direction == "gained_on_activation"].groupby("target").size()
    xs, cc = ccdf(in_all)
    ax.loglog(xs, cc, "o", ms=2.4, color="#c7ccd1", alpha=0.7,
              label=f"all census edges ({len(in_all):,} targets)")
    # the legend carries the edge counts, which are the SAME 2,298 / 2,081 as Figure 1a, so
    # the match is checkable at a glance without lengthening the title
    for dat, col, lab in [(lost.to_numpy(), S.PAL["LOST"],
                           f"effect LOST on stimulation ({int(lost.sum()):,} edges "
                           f"over {len(lost):,} targets)"),
                          (gain.to_numpy(), S.PAL["GAINED"],
                           f"effect GAINED on stimulation ({int(gain.sum()):,} edges "
                           f"over {len(gain):,} targets)")]:
        xs2, cc2 = ccdf(dat)
        ax.loglog(xs2, cc2, "o", ms=2.8, color=col, alpha=0.8, label=lab)
    both = len(set(lost.index) & set(gain.index))
    ax.set_xlabel("target in-degree $k$"); ax.set_ylabel("P(K $\\geq$ k)")
    ax.set_title("Target in-degree by reversal direction", fontsize=9.2)
    ax.legend(fontsize=6.9, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, -0.02))
    ax.text(0.02, 0.99,
            f"both directions: {both:,}/"
            f"{len(set(lost.index) | set(gain.index)):,} targets "
            f"({both / len(set(lost.index) | set(gain.index)):.1%})",
            transform=ax.transAxes, fontsize=7.4, va="top", color="#3d4348")
    S.panel_tag(ax, "b")

    res["direction_split"] = {
        "n_targets_lost": int(len(lost)), "n_targets_gained": int(len(gain)),
        "n_targets_both": int(both),
        "median_indegree_lost": float(np.median(lost)),
        "median_indegree_gained": float(np.median(gain)),
        "regulator_median_outdeg_reversing": float(np.median(od_rev)),
        "regulator_median_outdeg_nonreversing": float(np.median(od_non)),
        "_caveat": ("the regulator split is guaranteed by sampling: a regulator with more "
                    "edges is more likely to carry at least one reversal, so panel a is "
                    "descriptive only"),
    }

    # ---- panel c: empirical tails and all four fitted discrete candidate models ----
    # This is a graphical fit diagnostic, not a correlation plot. Ordered empirical and
    # theoretical quantiles can correlate near 1 even when a formal GOF test rejects a model.
    carrier = fig.add_subplot(gs[1, :]); carrier.axis("off")
    S.panel_tag(carrier, "c")
    inner = gs[1, :].subgridspec(1, 3, wspace=0.25)
    entries = [("regulator_outdegree_full_census", "regulator out-degree"),
               ("target_indegree_full_census", "target in-degree"),
               ("regulator_outdegree_reversing_subgraph", "reversals per regulator")]
    raw_data = [out_all, in_all, out_rev]
    styles = {
        "power_law": ("#7b3294", "-", "power law"),
        "lognormal": ("#2c7fb8", "--", "lognormal"),
        "exponential": ("#d99000", ":", "geometric"),
        "cutoff_power_law": ("#31966f", "-.", "power law + cutoff"),
    }
    shared_tail_ymin = max(0.4 / max(res[key]["n_tail"] for key, _ in entries), 1e-3)
    for j, (key, nice) in enumerate(entries):
        r = res[key]
        dat = np.asarray(raw_data[j])
        tail_values = dat[dat >= r["xmin"]]
        support = np.arange(r["xmin"], int(tail_values.max()) + 1)
        _, counts = np.unique(tail_values, return_counts=True)
        empirical = np.cumsum(counts[::-1])[::-1] / len(tail_values)
        empirical_x = np.unique(tail_values)
        ax = fig.add_subplot(inner[0, j])
        ax.step(empirical_x, empirical, where="post", color="#222222", lw=1.35,
                label="empirical")
        comp = r["model_comparison"]
        params = {
            "power_law": {"alpha": r["alpha"]},
            "lognormal": comp["lognormal_params"],
            "exponential": comp["exponential_params"],
            "cutoff_power_law": comp["cutoff_power_law_params"],
        }
        for model, (colour, line, label) in styles.items():
            fitted = model_survival(support, model, params[model], r["xmin"])
            ax.plot(support, fitted, color=colour, ls=line, lw=1.2, label=label)
        ax.set_xscale("log"); ax.set_yscale("log")
        tick_values = np.unique(np.rint(np.geomspace(support.min(), support.max(), 3)).astype(int))
        ax.xaxis.set_major_locator(FixedLocator(tick_values))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{int(value):,}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(shared_tail_ymin, 1.08)
        ax.set_xlabel("tail degree $k$" if j < 2 else "reversals $k$")
        if j == 0:
            ax.set_ylabel(r"$P(K\geq k\mid K\geq k_{min})$")
        else:
            ax.tick_params(labelleft=False)
        gof = r["gof_p"]
        verdict = "compatible" if gof >= 0.1 else "rejected"
        ax.set_title(nice, fontsize=8.0, pad=4)
        ax.text(0.03, 0.04,
                f"$k_{{min}}$={r['xmin']}; $n_{{tail}}$={r['n_tail']}\n"
                f"PL GOF $p$={gof:.3f} ({verdict})",
                transform=ax.transAxes, fontsize=6.7, va="bottom", color="#30363b")
    handles, labels = ax.get_legend_handles_labels()
    carrier.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.31),
                   ncol=5, fontsize=6.8, frameon=False, handlelength=2.4,
                   columnspacing=1.2)

    # Include the exact descriptive split behind panels a/b in the machine-readable result.
    cache.write_text(json.dumps(res, indent=2) + "\n")
    part_paths = {
        "regulator_outdegree_full_census": RD / "census/degree_distribution_test.regulator.json",
        "target_indegree_full_census": RD / "census/degree_distribution_test.target.json",
        "regulator_outdegree_reversing_subgraph": RD / "census/degree_distribution_test.reversal.json",
    }
    for result_key, path in part_paths.items():
        path.write_text(json.dumps(res[result_key], indent=2) + "\n")

    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"FigS_degree_distribution.{ext}", dpi=300, facecolor="white")

    for key in res:
        if not isinstance(res[key], dict) or "alpha" not in res[key]:
            continue
        r = res[key]
        print(f"{r['label']:28s} n={r['n']:5d}  kmin={r['xmin']:4d}  alpha={r['alpha']:.3f}  "
              f"KS={r['ks_distance']:.4f}  GOF p={r.get('gof_p', float('nan')):.3f}\n"
              f"      -> {r.get('verdict', '')}")
        for k, v in r["model_comparison"].items():
            if not k.startswith("power_law_vs_"):
                continue
            if "R" in v:
                print(f"      {k:32s} raw R={v['R']:+.3f} z={v['vuong_z']:+.3f} "
                      f"p={v['p']:.3g}  {v['favours']}")
            else:
                print(f"      {k:32s} deltaLL={v['delta_log_likelihood']:+.3f}  "
                      f"{v['test']}")


def run_single_analysis(which: str) -> None:
    """Run one fit/bootstrap family; useful for bounded-memory sequential orchestration."""
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "is_self"])
    e = e[~e.is_self]
    specs = {
        "regulator": (e.groupby("regulator").size().to_numpy(), "regulator out-degree", 0),
        "target": (e.groupby("target").size().to_numpy(), "target in-degree", 1),
        "reversal": (e[e.kind == "rev"].groupby("regulator").size().to_numpy(),
                     "reversals per regulator", 2),
    }
    if which not in specs:
        raise ValueError(f"unknown analysis {which!r}; choose one of {sorted(specs)}")
    data, label, seed_index = specs[which]
    # Advance through earlier families so a standalone family uses the same segment of the
    # deterministic RNG stream as a full run.
    rng = np.random.default_rng(SEED)
    ordered = [specs["regulator"], specs["target"], specs["reversal"]]
    result = None
    for j in range(seed_index + 1):
        dat_j, label_j, _ = ordered[j]
        result = analyse(dat_j, label_j, rng)
    assert result is not None and label_j == label
    out = RD / "census" / f"degree_distribution_test.{which}.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{label}: kmin={result['xmin']}, alpha={result['alpha']}, "
          f"KS={result['ks_distance']}, GOF p={result['gof_p']} -> {out}")




# --- Poisson Fit ---
#!/usr/bin/env python3
"""Tail-conditioned Poisson fit for the regulator out-degree panel.

WHY SEPARATE. Alam et al. 2017 Figure 1c compares the observed degree distribution against a
power law, a log-normal, an exponential AND a Poisson. `degree_distribution_test.py` fits the
first three plus an exponentially cut off power law, and its published output is quoted in the
manuscript as six non-nested Vuong comparisons and in repository figure R6. Adding a fifth model
inside that script would change that comparison set and silently invalidate those numbers, so the
Poisson is fitted here instead and reported alongside.

The Poisson is the reference for "no structure": it is the degree distribution of an
Erdos-Renyi graph, which is exactly the comparator drawn in Figure 1f. Conditioned on the same
k_min as every other tail model so the likelihoods are comparable.

Writes data/census/degree_poisson_fit.json.
"""












EDGES = ROOT / "data/census/census_master_edges.csv.gz"
FIT = ROOT / "data/census/degree_distribution_test.regulator.json"



def loglik_poisson(x: np.ndarray, lam: float, xmin: int) -> np.ndarray:
    """Per-observation log P(K = k | K >= xmin) for a Poisson truncated below xmin."""
    return stats.poisson.logpmf(x, lam) - np.log(stats.poisson.sf(xmin - 1, lam))


def fit_poisson(x: np.ndarray, xmin: int) -> float:
    def nll(lam):
        if lam <= 0:
            return np.inf
        v = loglik_poisson(x, float(lam), xmin)
        return np.inf if not np.all(np.isfinite(v)) else -v.sum()
    # The truncated MLE sits below the tail mean; bracket generously around it.
    r = optimize.minimize_scalar(nll, bounds=(1e-6, float(x.mean()) * 2.0), method="bounded",
                                 options=dict(xatol=1e-8))
    return float(r.x)


def poisson_survival(support: np.ndarray, lam: float, xmin: int) -> np.ndarray:
    """P(K >= k | K >= xmin)."""
    return stats.poisson.sf(np.asarray(support, int) - 1, lam) / stats.poisson.sf(xmin - 1, lam)


def main() -> None:
    fit = json.loads(FIT.read_text())
    xmin, alpha = int(fit["xmin"]), float(fit["alpha"])

    e = pd.read_csv(EDGES, usecols=["regulator", "is_self"])
    deg = e[~e.is_self].groupby("regulator").size().to_numpy().astype(np.int64)
    tail = deg[deg >= xmin]

    lam = fit_poisson(tail, xmin)
    ll_pois = loglik_poisson(tail, lam, xmin)
    ll_pl = loglik_powerlaw(tail, alpha, xmin)
    R, z, p = vuong(ll_pl, ll_pois)

    out = {
        "_what": "Tail-conditioned Poisson fit for regulator out-degree, Alam-2017 comparison set.",
        "_why": ("kept out of degree_distribution_test.py so the published six-comparison Vuong "
                 "set and repository figure R6 are unchanged"),
        "_script": "scripts/degree_poisson_fit.py",
        "edges_sha256": hashlib.sha256(EDGES.read_bytes()).hexdigest(),
        "xmin": xmin, "n_tail": int(tail.size),
        "poisson_lambda": lam,
        "log_likelihood": {"power_law": float(ll_pl.sum()), "poisson": float(ll_pois.sum())},
        "aicc": {"power_law": aicc(float(ll_pl.sum()), tail.size, 1),
                 "poisson": aicc(float(ll_pois.sum()), tail.size, 1)},
        "power_law_vs_poisson": {
            "R": float(R), "vuong_z": float(z), "p": float(p),
            "favours": ("power law" if z > 0 and p < 0.05 else
                        "poisson" if z < 0 and p < 0.05 else "indistinguishable"),
        },
    }
    out["delta_aicc_poisson_minus_power_law"] = out["aicc"]["poisson"] - out["aicc"]["power_law"]
    OUT.write_text(json.dumps(out, indent=2))

    print(f"xmin={xmin} n_tail={tail.size} lambda={lam:.1f}")
    print(f"  logL  power law {ll_pl.sum():10.2f}   poisson {ll_pois.sum():10.2f}")
    print(f"  AICc  power law {out['aicc']['power_law']:10.2f}   "
          f"poisson {out['aicc']['poisson']:10.2f}   "
          f"(delta {out['delta_aicc_poisson_minus_power_law']:+.1f})")
    print(f"  Vuong power law vs poisson: z={z:+.2f}  p={p:.3e}  -> "
          f"{out['power_law_vs_poisson']['favours']}")



def build():
    PAGE_SIZE = (8.27, 8.75)
    fig = plt.figure(figsize=PAGE_SIZE)
    outer = fig.add_gridspec(
        2, 1,
        height_ratios=[1.90, 1.0],
        hspace=0.30,
        left=0.060,
        # A little more page-side breathing room keeps the full c/h titles and
        # right-column axis descriptors inside the fixed 595.44-pt MediaBox.
        right=0.965,
        top=0.965,
        bottom=0.115,
    )
    
    # ---------------------------------------------------------------- upper block
    # Panel a spans the two upper tiers at left; b/c sit above d/e at right.
    upper = outer[0].subgridspec(1, 2, width_ratios=[1.28, 0.92], wspace=0.080)
    ax_a = fig.add_subplot(upper[0, 0])
    # hi_mag marks the 532-edge high-magnitude subgraph inside the complete census, the subgraph
    # Figure 3a redrew as its own panel in this same two-ring grammar. Emphasis rides on opacity
    # and linewidth only, leaving red/teal to carry direction as it does elsewhere.
    panel_census_network(
        ax_a, "a", alpha=0.085, lw=0.20, fs=1.30, n_reg_labels=14,
        labels=True, reg_uniform=False, legend=True, short_legend=True,
        lab_pitch=0.066, hi_mag=0.5, dual_role=True,
    )
    # The reference composition gives the census ring slightly more weight.  The
    # matched x/y spans enlarge it by about 2.5%; the positive x-centre moves the
    # ring itself left, clearing the neighbouring d axis without moving the title.
    ax_a.set_xlim(-1.05, 1.29)
    ax_a.set_ylim(-1.03, 1.03)
    # Separate the two nearly vertical inner-ring labels without changing their
    # identities, angles, or underlying node positions.
    for _text in ax_a.texts:
        if _text.get_text() == "ATAD5":
            _x, _y = _text.get_position()
            _text.set_position((_x + 0.025, _y))
        elif _text.get_text() == "CCNC":
            _x, _y = _text.get_position()
            _text.set_position((_x - 0.025, _y))
        elif _text.get_text() == "TADA2B":
            _x, _y = _text.get_position()
            _text.set_position((_x + 0.085, _y))
        elif _text.get_text() == "CCL4":
            _x, _y = _text.get_position()
            _text.set_position((_x, _y + 0.022))
        elif _text.get_text() == "CCR5":
            _x, _y = _text.get_position()
            _text.set_position((_x, _y - 0.006))
        elif _text.get_text() == "IFNG":
            _x, _y = _text.get_position()
            _text.set_position((_x, _y - 0.015))
    ax_a.set_title(
        f"Complete sign-reversal census of {c['inverting']:,} edges\n"
        "324 regulators · 2,577 targets",
        loc="left", pad=9,
    )
    
    # Right side follows the supplied reference exactly: b/c side-by-side, d/e below.
    right = upper[0, 1].subgridspec(
        2, 2,
        width_ratios=[1.08, 1.0],
        height_ratios=[1.0, 1.0],
        wspace=0.42,
        hspace=0.45,
    )
    
    ax_b = fig.add_subplot(right[0, 0])
    panel_scatter(ax_b, "b")
    ax_b.set_title("Edges that reverse sign\nbetween states", y=0.98, pad=0)
    ax_b.xaxis.labelpad = 1
    for _text in ax_b.texts:
        if _text.get_text() == "b":
            _text.set_position((-0.28, 1.12))
    _bh, _ = ax_b.get_legend_handles_labels()
    _leg = ax_b.legend(
        _bh, ["lost (+→−)", "gained (−→+)"], loc="lower left",
        frameon=True, framealpha=0.96, facecolor="white", edgecolor="#d0d3d6",
        handletextpad=0.30, borderpad=0.22, labelspacing=0.20, fontsize=5.0,
    )
    _leg.get_frame().set_linewidth(0.5)
    for _text in ax_b.texts:
        if _text.get_text().startswith("reversals track"):
            _text.set_text("reversals track the\nanti-diagonal (y ≈ −x)")
            _text.set_fontsize(5.0)
            _text.set_bbox(dict(boxstyle="square,pad=0.12", facecolor="white",
                                edgecolor="none", alpha=0.92))
    
    # Canonical census and direction panels, nested in c's cell.  No label or denominator is dropped.
    right_c = right[0, 1].subgridspec(2, 1, height_ratios=[1.00, 0.78], hspace=1.02)
    ax_c_top = fig.add_subplot(right_c[0, 0])
    panel_census(ax_c_top, "c")
    ax_c_bottom = fig.add_subplot(right_c[1, 0])
    panel_direction(ax_c_bottom)
    ax_c_top.xaxis.labelpad = 2
    ax_c_top.set_title("A rare minority\nof shared edges", pad=2)
    ax_c_bottom.set_title(ax_c_bottom.get_title(), pad=1)
    
    ax_d = fig.add_subplot(right[1, 0])
    panel_transition(ax_d, "d")
    # The reference leaves a clean gutter between the network labels and the two
    # left-hand right-column panels.  Shift b/d into the otherwise generous b-c
    # gutter; no panel or data artist is resized.
    for _ax in (ax_b, ax_d):
        _p = _ax.get_position()
        _ax.set_position([_p.x0 + 0.016, _p.y0, _p.width, _p.height])
    _dleg = ax_d.get_legend()
    if _dleg is not None:
        _dleg.set_bbox_to_anchor((0.0, 0.0))
        _dleg._loc = 3
    for _text in ax_d.texts:
        if _text.get_text().startswith("still on the original sign"):
            _text.set_text("still on the\noriginal sign at 8 h")
            _text.set_position((14, 0.61))
            _text.set_ha("left")
            _text.set_fontsize(6.0)
    _n, _k = load()
    ax_e = fig.add_subplot(right[1, 1])
    _gini = panel_concentration(ax_e, _k)
    ax_e.yaxis.labelpad = 0
    ax_e.set_xlabel("regulators ranked by\nreversal count\n(cumulative fraction)")
    panel_tag(ax_e, "e", dx=-0.18, dy=1.08)
    
    # ---------------------------------------------------------------- bottom row
    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.45, 1.60, 1.0], wspace=0.42)
    
    ax_f = fig.add_subplot(bottom[0, 0])
    _dd = panel_degree_ccdf(ax_f, compact=False, fs=0.85)
    ax_f_right = fig.axes[-1]
    ax_f_right.set_ylabel("cumulative %\nof regulators", labelpad=0)
    ax_f_right.yaxis.set_label_coords(1.12, 0.50)
    ax_f.set_title("Out-degree is compatible with\na power-law tail", loc="left", pad=3)
    panel_tag(ax_f, "f", dx=-0.18, dy=1.08)
    ax_f.legend(
        handles=legend_handles(fs=0.85),
        loc="upper center", bbox_to_anchor=(0.50, -0.22),
        ncol=3, fontsize=5.6, handlelength=1.7,
        columnspacing=0.75, handletextpad=0.35,
        labelspacing=0.30, borderpad=0.15, frameon=False,
    )
    
    ax_g = fig.add_subplot(bottom[0, 1])
    bars(ax_g, D, DONOR_TITLE, "g", tag_dx=-0.055)
    # Use the available f-g gutter rather than stacking the two vertical y labels.
    _gp = ax_g.get_position()
    ax_g.set_position([_gp.x0 + 0.015, _gp.y0, _gp.width, _gp.height])
    ax_g.set_xticklabels([f"{t}\n{k}/{n}\n{g} reg."
                          for t, n, g, k in ax_g.tier_counts])
    ax_g.set_xlabel("minimum |log2FC| threshold")
    ax_g.yaxis.set_label_coords(-0.14, 0.50)
    for _text in ax_g.texts:
        if _text.get_text().startswith("– – stable-edge ceiling"):
            _text.set_text("– – stable-edge ceiling (~100%)")
            _text.set_transform(ax_g.transAxes)
            _text.set_position((0.02, 0.96))
            _text.set_va("top")
            _text.set_fontsize(6.0)
        elif _text.get_text().startswith("– – within-regulator permutation"):
            _text.set_text(
                "– – within-regulator permutation null, per bin:\n"
                "2.33%, 2.26%, 2.54%, 4.56%"
            )
            _text.set_transform(ax_g.transAxes)
            _text.set_position((0.02, 0.89))
            _text.set_va("top")
            _text.set_fontsize(6.0)
    
    ax_h = fig.add_subplot(bottom[0, 2])
    panel_guide(ax_h, "h")
    
    for _text in ax_g.texts:
        if _text.get_text() == "g":
            _text.set_position((-0.18, 1.10))
    
    for _text in ax_h.texts:
        if _text.get_text().startswith("cluster-bootstrap 95% CI spans"):
            _text.set_text("CI: 0–100%\n4 clusters;\ndirectional only")
            _text.set_position((0.98, 0.62))
            _text.set_ha("right")
            _text.set_fontsize(6.0)
        elif _text.get_text().startswith("– – stable ceiling"):
            _text.set_text(f"ceiling {G['by_threshold']['>0']['stable_repro']*100:.1f}%")
            _text.set_transform(ax_h.transAxes)
            _text.set_position((0.02, 0.88))
            _text.set_va("top")
            _text.set_bbox(dict(boxstyle="square,pad=0.10", facecolor="white",
                                edgecolor="none", alpha=0.94))
            _text.set_fontsize(6.0)
        elif _text.get_text().startswith("– – permutation null"):
            _text.set_text(f"null {G['by_threshold']['>0']['perm_null_mean']*100:.2f}%")
            _text.set_transform(ax_h.transAxes)
            _text.set_position((0.02, 0.06))
            _text.set_va("bottom")
            _text.set_bbox(dict(boxstyle="square,pad=0.10", facecolor="white",
                                edgecolor="none", alpha=0.94))
            _text.set_fontsize(6.0)
        elif _text.get_text().startswith("single-guide DE"):
            _text.set_text("single-guide\nestimates noisy;\ndirectional only")
            _text.set_position((0.98, 0.17))
            _text.set_ha("right")
            _text.set_fontsize(6.0)
        elif _text.get_text() == "75.0%":
            _text.set_fontsize(7.0)
            _text.set_position((0.0, 107.0))
    
    # Preserve the predecessor typography outside the explicitly changed b-e panels.
    # Mathtext exponents need larger parent tick sizes to remain at the same effective floor.
    fig.canvas.draw()
    for _text in fig.findobj(match=plt.Text):
        if _text.get_text().strip() and _text.get_fontsize() < 7.0:
            _text.set_fontsize(7.0)
    # The 5-6 pt tier is restricted to true legends and in-panel notes.  Axis labels,
    # ticks, titles and panel tags remain at the 7 pt ordinary source floor.
    for _text in ax_b.texts:
        if _text.get_text().startswith("reversals track"):
            _text.set_fontsize(5.0)
    for _text in ax_c_top.texts + ax_c_bottom.texts:
        _text.set_fontsize(6.0)
    for _legend in (_leg, _dleg, ax_f.get_legend()):
        if _legend is not None:
            for _legend_text in _legend.get_texts():
                if _legend is _leg:
                    _legend_text.set_fontsize(5.0)
                else:
                    _legend_text.set_fontsize(6.0 if _legend is not ax_f.get_legend() else 5.6)
    if _dleg is not None:
        for _dtext in _dleg.get_texts():
            _dtext.set_fontsize(5.5)
    for _text in ax_d.texts:
        if _text.get_text().startswith("still on the original sign"):
            _text.set_fontsize(6.0)
        elif "KAT7 reversals lack" in _text.get_text():
            _text.set_fontsize(5.5)
    for _text in ax_g.texts:
        _glabel = _text.get_text()
        if "/" in _glabel and _glabel.replace("/", "").isdigit():
            _gx, _ = _text.get_position()
            _text.set_position((_gx, 15.5))
        elif (_glabel.startswith("– – stable-edge ceiling") or
                _text.get_text().startswith("– – within-regulator")):
            _text.set_fontsize(6.0)
        elif _glabel.endswith(" reg."):
            # In-bar sample-size labels are compact annotations, not ordinary axes.
            # Restore their intended small tier so they remain separate from k/n.
            _text.set_fontsize(6.0)
    # The four compact explanatory annotations in panel h are inset labels, not ordinary axis text.
    for _text in ax_h.texts:
        if (_text.get_text().startswith("CI:") or
                _text.get_text().startswith("ceiling ") or
                _text.get_text().startswith("null ") or
                _text.get_text().startswith("single-guide")):
            _text.set_fontsize(6.0)
    
    # Requested compact typography applies only to panels b-e.  Titles and panel tags
    # are prominent 7 pt text; every other visible text artist inside these axes is 5 pt.
    RIGHT_PANEL_AXES = (ax_b, ax_c_top, ax_c_bottom, ax_d, ax_e)
    RIGHT_PANEL_TAGS = ((ax_b, "b"), (ax_c_top, "c"), (ax_d, "d"), (ax_e, "e"))
    for _ax in RIGHT_PANEL_AXES:
        for _title_artist in (_ax.title, _ax._left_title, _ax._right_title):
            if _title_artist.get_text().strip():
                _title_artist.set_fontsize(7.0)
                _title_artist.set_fontweight("bold")
        _ax.xaxis.label.set_fontsize(5.0)
        _ax.yaxis.label.set_fontsize(5.0)
        _ax.xaxis.offsetText.set_fontsize(5.0)
        _ax.yaxis.offsetText.set_fontsize(5.0)
        _ax.tick_params(axis="both", which="both", labelsize=5.0)
        for _text in _ax.texts:
            _text.set_fontsize(5.0)
        _axis_legend = _ax.get_legend()
        if _axis_legend is not None:
            for _legend_text in _axis_legend.get_texts():
                _legend_text.set_fontsize(5.0)
    for _ax, _tag in RIGHT_PANEL_TAGS:
        _tags = [_text for _text in _ax.texts if _text.get_text() == _tag]
        assert len(_tags) == 1, f"panel {_tag} tag count changed: {len(_tags)}"
        _tags[0].set_fontsize(7.0)
        _tags[0].set_fontweight("bold")
    
    for _text in ax_f.get_xticklabels() + ax_f.get_yticklabels():
        _text.set_fontsize(10.5)
    fig.canvas.draw()
    
    EXPECTED_RIGHT_TITLES = {
        ax_b: ("center", "Edges that reverse sign\nbetween states"),
        ax_c_top: ("center", "A rare minority\nof shared edges"),
        ax_c_bottom: ("center", "Balanced in direction"),
        ax_d: ("left", "Median effects at 8 h\nnot reversed"),
        ax_e: ("left", "Burden is concentrated"),
    }
    for _ax, (_title_loc, _expected_title) in EXPECTED_RIGHT_TITLES.items():
        assert _ax.get_title(loc=_title_loc) == _expected_title, (
            f"right-panel title changed: {_ax.get_title(loc=_title_loc)!r} "
            f"!= {_expected_title!r}"
        )
        _title_artist = {
            "center": _ax.title,
            "left": _ax._left_title,
            "right": _ax._right_title,
        }[_title_loc]
        assert abs(_title_artist.get_fontsize() - 7.0) < 1e-9
        for _label in (_ax.xaxis.label, _ax.yaxis.label,
                       _ax.xaxis.offsetText, _ax.yaxis.offsetText):
            assert abs(_label.get_fontsize() - 5.0) < 1e-9
        for _tick_text in (
                _ax.get_xticklabels() + _ax.get_yticklabels()
                + _ax.get_xticklabels(minor=True) + _ax.get_yticklabels(minor=True)):
            if _tick_text.get_text().strip():
                assert abs(_tick_text.get_fontsize() - 5.0) < 1e-9
        _tag_values = {_tag for _tag_ax, _tag in RIGHT_PANEL_TAGS if _tag_ax is _ax}
        for _text in _ax.texts:
            _expected_size = 7.0 if _text.get_text() in _tag_values else 5.0
            assert abs(_text.get_fontsize() - _expected_size) < 1e-9, (
                f"right-panel text size drift: {_text.get_text()!r} "
                f"is {_text.get_fontsize():g} pt, expected {_expected_size:g}"
            )
        _axis_legend = _ax.get_legend()
        if _axis_legend is not None:
            for _legend_text in _axis_legend.get_texts():
                assert abs(_legend_text.get_fontsize() - 5.0) < 1e-9
    
    _fw, _fh = fig.get_size_inches()
    assert tuple(round(float(_v), 2) for _v in (_fw, _fh)) == PAGE_SIZE
    panel_axes = {
        "a": ax_a, "b": ax_b, "c_top": ax_c_top, "c_bottom": ax_c_bottom,
        "d": ax_d, "e": ax_e, "f": ax_f, "g": ax_g, "h": ax_h,
    }
    panel_geometry = {}
    for _name, _ax in panel_axes.items():
        _p = _ax.get_position()
        panel_geometry[_name] = {
            "width_in": round(_p.width * _fw, 5),
            "height_in": round(_p.height * _fh, 5),
            "bbox_fraction": [round(x, 8) for x in (_p.x0, _p.y0, _p.width, _p.height)],
        }
        print(f"  panel {_name}: {_p.width*_fw:.2f}w x {_p.height*_fh:.2f}h in")
    
    EXPECTED_PANEL_SIZE_IN = {
        "a": (4.18705, 3.68603),
        "b": (1.29140, 1.29140),
        "c_top": (1.19574, 0.64346),
        "c_bottom": (1.19574, 0.50190),
        "d": (1.29140, 1.72949),
        "e": (1.19574, 1.72949),
        "f": (2.09342, 2.23013),
        "g": (2.30998, 2.23013),
        "h": (1.44374, 2.23013),
    }
    for _name, (_expected_w, _expected_h) in EXPECTED_PANEL_SIZE_IN.items():
        _observed = panel_geometry[_name]
        assert abs(_observed["width_in"] - _expected_w) <= 0.00002, (
            f"panel {_name} width changed: {_observed['width_in']} != {_expected_w}"
        )
        assert abs(_observed["height_in"] - _expected_h) <= 0.00002, (
            f"panel {_name} height changed: {_observed['height_in']} != {_expected_h}"
        )
    assert panel_geometry["a"]["width_in"] > 3.89535
    assert panel_geometry["b"]["width_in"] < 1.41657
    assert panel_geometry["c_top"]["width_in"] < 1.31164
    
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{STEM}.{ext}", dpi=300, facecolor="white")
    print(f"wrote {STEM}")

if __name__ == '__main__':
    build()
