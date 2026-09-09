from matplotlib.collections import PathCollection
import textwrap
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import LinearSegmentedColormap
from fractions import Fraction
import os, sys, json, platform, argparse, hashlib, math
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as patheffects
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, FancyBboxPatch
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy import stats
import numpy as np
import pandas as pd

class figstyle: pass
class fig4_architecture: pass
class fig_target_activation_trajectory_v1: pass
class fig2g_trajectory_panel: pass
class fig2_biology_candidate: pass
class fig2h_activation_movement_similarity: pass
class fig2_eight_panel_candidate: pass
class fig2h_8hr_sign_states: pass
class fig2_eight_panel_8h_candidate: pass
class fig2h_8hr_sign_states_scatter: pass
class fig2_eight_panel_original_line_a4_candidate: pass
class fig2_panel_b_smallbar_visibility_candidate: pass
class fig2_panel_b_main_scale_visibility_candidate: pass
class fig_signed_kegg_transition_heatmap_1d_exploratory: pass
class fig2_signed_kegg_side_by_side_candidate_v1: pass
class fig2_signed_kegg_side_by_side_candidate_v2: pass
class fig2_signed_kegg_side_by_side_candidate_v3: pass
class fig2_signed_kegg_side_by_side_candidate_v4: pass
class fig2_signed_kegg_side_by_side_candidate_v5: pass
class fig2_signed_kegg_side_by_side_candidate_v7: pass
class fig2_signed_kegg_side_by_side_candidate_v8: pass
class fig2_signed_kegg_side_by_side_candidate_v9: pass
class fig2_signed_kegg_side_by_side_candidate_v10: pass
class fig2_signed_kegg_side_by_side_candidate_v11: pass

def init_figstyle():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
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
    
    
    # A4 portrait canvas in inches. Figure 1 is saved at exactly this size, so it must NOT be
    # written with bbox_inches="tight" (that crops to content and silently discards the page size).
    A4_PORTRAIT = (8.27, 11.69)
    
    OUT = ROOT / "figures"
    os.makedirs(OUT, exist_ok=True)
    
    def panel_tag(ax, s, dx=-0.04, dy=1.06):
        ax.text(dx, dy, s, transform=ax.transAxes, fontsize=12 * _SCALE[0], fontweight="bold",
                va="top", ha="right", color="#22262b")
    
    def save(fig, name):
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight",
                        facecolor="white")
        print("  saved", name + ".pdf/.png")
    
    for k, v in list(locals().items()):
        if not k.startswith('_') or k == '_SCALE': setattr(figstyle, k, v)
init_figstyle()

def init_fig4_architecture():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python
    """Regulator architecture panels used by the main and standalone figures.
    
    The standalone Fig4_architecture retains the original top-20 bar plot plus out-degree panel.
    The all-324 regulator views are exposed separately and written to separately named files.
    
    Refactored into panel functions (panel_count / panel_rate, each taking an axes + panel tag) so it
    renders standalone AND composes into the combined architecture+network figure."""
    S = figstyle
    S.setup()
    
    
    D = json.load(open(ROOT / "data/census/figure_data.json"))
    # per_regulator_rates.csv is SUPERSEDED: it drew n_both from B_masked (61,090 edges) but
    # n_inv from the full census, inflating every rate and omitting TADA2B and SUPT7L entirely.
    # perreg_rates_and_enrichment.py rebuilds both counts from one edge table.
    rr = pd.read_csv(ROOT / "data/census/per_regulator_rates.csv")
    fig2a_test = pd.read_csv(ROOT / "data/census/per_regulator_enrichment_all.csv")
    enr = dict(zip(rr.gene, rr.enrichment)); nboth = dict(zip(rr.gene, rr.n_both)); rate = dict(zip(rr.gene, rr.rate))
    top = D["top_regulators"]
    base = rr.n_inv.sum() / rr.n_both.sum()             # trans-edge baseline, self edges masked
    
    CLASS_COL = {
        "SAGA/co-activator (HAT)": "#6a51a3", "Mediator": "#2c6fbb",
        "CBP/p300 (HAT)": "#8d6cc4", "PAF1C (elongation)": "#3f8f7d",
        "signalling / activation TF": "#e07b39", "bHLH-PAS TF": "#a83e83",
        "cytoskeletal / immune-synapse": "#4f9130", "other (replication/translation/misc)": "#9aa0a6",
        "chromatin regulation": "#8c6d31", "RNA processing / surveillance": "#4d788e",
        "metabolic / mitochondrial": "#d39b2a",
        "other / unclassified": "#c7ccd1",
    }
    def col(k): return CLASS_COL.get(k, "#c7ccd1")
    
    # Complete annotation for the genes eligible to be labelled by the two-sided Figure 2 test.
    # This is descriptive annotation, not a tested class-level enrichment. PAF1C remains in the
    # shared palette for the standalone architecture plot but no PAF1C member enters the current
    # two-sided top 20, so an empty PAF1C item is not drawn in Figure 2.
    FIG2_CLASS = {
        "CCNC": "Mediator", "MED24": "Mediator",
        "SUPT7L": "SAGA/co-activator (HAT)",
        "AHR": "bHLH-PAS TF", "ARNT": "bHLH-PAS TF",
        "NCKAP1L": "cytoskeletal / immune-synapse",
        "DOCK2": "cytoskeletal / immune-synapse",
        "ARHGAP30": "cytoskeletal / immune-synapse",
        "STAT5B": "signalling / activation TF", "STAT3": "signalling / activation TF",
        "CBFB": "signalling / activation TF",
        "NFRKB": "chromatin regulation", "SETDB1": "chromatin regulation",
        "SIN3B": "chromatin regulation",
        "SMG1": "RNA processing / surveillance", "RBM10": "RNA processing / surveillance",
        "RPRD2": "RNA processing / surveillance",
        "DOLPP1": "metabolic / mitochondrial", "ATP5ME": "metabolic / mitochondrial",
        "STIP1": "other (replication/translation/misc)",
    }
    FIG2_CLASS_ORDER = [
        "Mediator", "SAGA/co-activator (HAT)", "bHLH-PAS TF",
        "signalling / activation TF", "cytoskeletal / immune-synapse",
        "chromatin regulation", "RNA processing / surveillance",
        "metabolic / mitochondrial", "other (replication/translation/misc)",
        "other / unclassified",
    ]
    
    # Direction-resolved counts come from the same authoritative Rest--Stim48 trans-edge census as
    # the regulator test. Loading them here prevents the bar panel from silently falling back to the
    # superseded 8 h subset.
    _edge_direction = pd.read_csv(
        ROOT / "data/census/census_master_edges.csv.gz",
        usecols=["regulator", "target", "kind", "direction", "is_self"])
    _edge_direction = _edge_direction[~_edge_direction.is_self]
    _edge_direction = _edge_direction[_edge_direction.kind == "rev"]
    direction_counts = (pd.crosstab(_edge_direction.regulator, _edge_direction.direction)
                        .reindex(columns=["lost_on_activation", "gained_on_activation"],
                                 fill_value=0))
    assert len(_edge_direction) == 4_379 and direction_counts.sum().sum() == 4_379
    
    
    def fig2_class_handles(genes):
        """Legend handles for only the functional annotations present in `genes`."""
        present = {FIG2_CLASS.get(g, "other / unclassified") for g in genes}
        return [Line2D([0], [0], marker="s", ls="", ms=6, color=col(category),
                       label=category)
                for category in FIG2_CLASS_ORDER if category in present]
    
    
    def fig2_direction_test_handles():
        """Keep edge direction and regulator-level inferential direction visually distinct."""
        return [
            Line2D([0], [0], color=S.PAL["LOST"], lw=4,
                   label="bar: lost on stimulation (+ to -)"),
            Line2D([0], [0], color=S.PAL["GAINED"], lw=4,
                   label="bar: gained on stimulation (- to +)"),
            Line2D([0], [0], marker="^", ls="", ms=5, color="black",
                   label="test: TOTAL reversal burden above expectation"),
            Line2D([0], [0], marker="v", ls="", ms=5, color="black",
                   label="test: TOTAL reversal burden below expectation"),
        ]
    
    
    def panel_top20_bar(ax, tag="a"):
        """Original top-20 raw-contributor bar panel, shaded by reversal rate."""
        sub = rr.nlargest(20, "n_inv").iloc[::-1]
        yy = np.arange(len(sub))
        pct = sub.pct.to_numpy()
    
        tested = pd.read_csv(ROOT / "data/census/per_regulator_tested_all.csv")
        tested = tested.set_index("gene")
        marks = {}
        for gene in sub.gene:
            q = tested.loc[gene, "q_degree_matched"]
            fold = tested.loc[gene, "enrichment_vs_decile"]
            if q < 0.05:
                stars = "***" if q < 0.001 else ("**" if q < 0.01 else "*")
                marks[gene] = (stars, "up" if fold > 1 else "down")
    
        cmap = plt.get_cmap("YlGnBu")
        norm = mpl.colors.Normalize(vmin=0, vmax=float(pct.max()))
        bar_cols = [cmap(0.18 + 0.82 * norm(v)) for v in pct]
        ax.barh(yy, sub.n_inv.to_numpy(), color=bar_cols, height=0.72,
                edgecolor="white", linewidth=0.4)
        for i, r in enumerate(sub.itertuples()):
            txt = f"{r.n_inv}   {r.pct:.1f}%  ({r.enrichment:.1f}\u00d7)"
            ax.text(r.n_inv + 3, i, txt, va="center", fontsize=6.4, color="#3d4348")
            if r.gene in marks:
                stars, direction = marks[r.gene]
                ax.text(0.995, i, ("\u25b2" if direction == "up" else "\u25bc") + stars,
                        transform=ax.get_yaxis_transform(), ha="right", va="center",
                        fontsize=6.0, fontweight="bold", color="black")
    
        ax.set_yticks(yy)
        ax.set_yticklabels(sub.gene, fontsize=7.5)
        cls_of = {r["gene"]: r["class"] for r in top}
        for lab in ax.get_yticklabels():
            lab.set_color(col(cls_of.get(lab.get_text(), "other / unclassified")))
        ax.set_xlabel("reversing target edges (count)")
        ax.set_title("Raw contributors, shaded by per-edge reversal rate")
        ax.set_xlim(0, sub.n_inv.max() * 1.52)
    
        cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                          fraction=0.030, pad=0.015)
        cb.set_label("per-edge reversal rate (%)", fontsize=6.4)
        cb.ax.tick_params(labelsize=6)
        cb.ax.axhline(base * 100, color="#c1272d", lw=1.2)
        cb.ax.text(1.9, base * 100, f"baseline {base*100:.1f}%", fontsize=5.6,
                   color="#c1272d", va="center", transform=cb.ax.get_yaxis_transform())
    
        classes_present = []
        for gene in sub.gene:
            category = cls_of.get(gene, "other / unclassified")
            if category not in classes_present:
                classes_present.append(category)
        handles = [Line2D([0], [0], marker="s", ls="", ms=6, color=col(category),
                          label=category) for category in classes_present]
        handles += [
            Line2D([0], [0], ls="", marker="^", ms=5, color="black",
                   label="reverses MORE than out-degree predicts"),
            Line2D([0], [0], ls="", marker="v", ms=5, color="black",
                   label="reverses LESS than out-degree predicts"),
        ]
        ax.legend(handles=handles, loc="lower right", fontsize=5.5, handletextpad=0.3,
                  labelspacing=0.28,
                  title="gene-label colour   |   * q<0.05  ** q<0.01  *** q<0.001",
                  title_fontsize=5.5)
        S.panel_tag(ax, tag, dx=-0.30, dy=1.04)
    
    
    def panel_count(ax, tag="a", detail_scale=1.0):
        """All reversal-bearing regulators: two-sided deviation from matched expectation."""
        sub = fig2a_test[fig2a_test.shown_in_fig2a].copy()
        sub["expected_n_inv"] = sub.n_both * sub.p_expected
        sig_up = sub.sig_two_sided & (sub.direction == "up")
        sig_down = sub.sig_two_sided & (sub.direction == "down")
        nonsig = ~sub.sig_two_sided
    
        ax.scatter(sub.loc[nonsig, "expected_n_inv"], sub.loc[nonsig, "n_inv"],
                   s=15 * detail_scale ** 2, color="#c7ccd1", edgecolor="#9aa0a6",
                   lw=0.35, alpha=0.82, zorder=2)
        ax.scatter(sub.loc[sig_up, "expected_n_inv"], sub.loc[sig_up, "n_inv"],
                   s=31 * detail_scale ** 2, marker="^", color="black", edgecolor="white",
                   lw=0.55, alpha=0.96, zorder=3)
        ax.scatter(sub.loc[sig_down, "expected_n_inv"], sub.loc[sig_down, "n_inv"],
                   s=31 * detail_scale ** 2, marker="v", color="black", edgecolor="white",
                   lw=0.55, alpha=0.96, zorder=3)
    
        lo = min(sub.expected_n_inv.min(), sub.n_inv.min()) * 0.78
        hi = max(sub.expected_n_inv.max(), sub.n_inv.max()) * 1.28
        ax.plot([lo, hi], [lo, hi], color="#565b61", lw=0.9, ls="--", zorder=1)
    
        # Fixed offsets keep the deterministic two-sided top-20 labels legible in the narrow
        # three-panel architecture row. Labels are coloured by functional annotation; marker shape
        # and colour, not label colour, carry the regulator-level test direction.
        offsets = {
            "CCNC": (-8, 9, "right"), "NCKAP1L": (-7, 7, "right"),
            "AHR": (-8, 11, "right"), "ARNT": (7, 8, "left"),
            "DOCK2": (-7, 8, "right"), "DOLPP1": (-7, -8, "right"),
            "SMG1": (-7, -8, "right"), "ARHGAP30": (6, 4, "left"),
            "STAT5B": (-8, 12, "right"), "STIP1": (8, 1, "left"),
            "STAT3": (7, -9, "left"), "CBFB": (-7, -15, "right"),
            "NFRKB": (7, 8, "left"), "SUPT7L": (7, -8, "left"),
            "ATP5ME": (-7, 8, "right"), "RBM10": (-7, 3, "right"),
            "MED24": (-8, 10, "right"), "RPRD2": (7, -8, "left"),
            "SETDB1": (-7, -10, "right"), "SIN3B": (7, -11, "left"),
        }
        labelled = sub[sub.label_rank_fig2a_two_sided > 0].sort_values(
            "label_rank_fig2a_two_sided")
        assert len(sub) == 324 and len(labelled) == 20 and labelled.sig_two_sided.all()
        for r in labelled.itertuples():
            dx, dy, ha = offsets[r.gene]
            ax.annotate(r.gene, (r.expected_n_inv, r.n_inv),
                        xytext=(dx * detail_scale, dy * detail_scale),
                        textcoords="offset points", ha=ha, va="center",
                        fontsize=5.9 * detail_scale, fontweight="bold",
                        color=col(FIG2_CLASS.get(r.gene, "other / unclassified")), zorder=4)
    
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("expected reversing edges (out-degree matched, log scale)")
        ax.set_ylabel("observed reversing edges (log scale)")
        n_up = int(sig_up.sum())
        n_down = int(sig_down.sum())
        n_nonsig = int(nonsig.sum())
        ax.set_title(f"All reversal-bearing regulators\n(n = {len(sub)})")
        handles = [
            Line2D([0], [0], ls="", marker="^", ms=5 * detail_scale,
                   markerfacecolor="black", markeredgecolor="white",
                   label=f"more than expected ({n_up})"),
            Line2D([0], [0], ls="", marker="v", ms=5 * detail_scale,
                   markerfacecolor="black", markeredgecolor="white",
                   label=f"fewer than expected ({n_down})"),
            Line2D([0], [0], ls="", marker="o", ms=4.5 * detail_scale, markerfacecolor="#c7ccd1",
                   markeredgecolor="#9aa0a6", label=f"not significant ({n_nonsig})"),
            Line2D([0], [0], color="#565b61", lw=0.9, ls="--",
                   label="degree-matched expectation"),
        ]
        ax.legend(handles=handles, loc="upper left", fontsize=5.7 * detail_scale,
                  handletextpad=0.35,
                  labelspacing=0.28, title="two-sided Fisher; BH over 617",
                  title_fontsize=5.7 * detail_scale)
        if tag:
            S.panel_tag(ax, tag, dx=-0.22, dy=1.04)
    
    
    def panel_fisher_labelled_bar(ax, tag="b", detail_scale=1.0):
        """Count/rate/direction view of the two-sided Fisher labels from panel a."""
        sub = (fig2a_test[fig2a_test.label_rank_fig2a_two_sided > 0]
               .sort_values(["n_inv", "p_two_sided", "gene"],
                            ascending=[False, True, True]).copy())
        assert len(sub) == 20 and sub.sig_two_sided.all()
        sub = sub.join(direction_counts, on="gene")
        assert (sub.lost_on_activation + sub.gained_on_activation == sub.n_inv).all()
    
        yy = np.arange(len(sub))
        lost = sub.lost_on_activation.to_numpy()
        gained = sub.gained_on_activation.to_numpy()
        ax.barh(yy, lost, color=S.PAL["LOST"], height=0.72,
                edgecolor="white", linewidth=0.35)
        ax.barh(yy, gained, left=lost, color=S.PAL["GAINED"], height=0.72,
                edgecolor="white", linewidth=0.35)
        for i, r in enumerate(sub.itertuples()):
            txt = (f"{int(r.n_inv)}/{int(r.n_both)}  {100*r.rate:.1f}%  "
                   f"{r.enrichment_vs_decile:.1f}\N{MULTIPLICATION SIGN}")
            # The 0--350 axis leaves only 62 units after CCNC's 288-edge bar. Put that single
            # annotation inside its teal segment; all other rows retain the common outside position.
            inside = r.gene == "CCNC"
            ax.text(r.n_inv - 3 if inside else r.n_inv + 3, i, txt,
                    ha="right" if inside else "left", va="center",
                    fontsize=5.85 * detail_scale,
                    color="white" if inside else "#3d4348",
                    fontweight="bold" if inside else "normal")
            stars = "***" if r.q_two_sided < 0.001 else (
                "**" if r.q_two_sided < 0.01 else "*")
            is_up = r.enrichment_vs_decile >= 1
            ax.text(0.995, i, ("\N{BLACK UP-POINTING TRIANGLE}" if is_up
                               else "\N{BLACK DOWN-POINTING TRIANGLE}") + stars,
                    transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=5.8 * detail_scale, fontweight="bold", color="black")
    
        ax.set_yticks(yy)
        ax.set_yticklabels(sub.gene, fontsize=7.0 * detail_scale, fontweight="bold")
        for lab in ax.get_yticklabels():
            lab.set_color(col(FIG2_CLASS.get(lab.get_text(), "other / unclassified")))
        ax.set_ylim(len(sub) - 0.35, -0.65)  # descending reversing-edge count at the top
        # Fixed to 350 so the 288-edge maximum retains annotation headroom while the three
        # top-row panels have comparable, near-square visual footprints.
        assert sub.n_inv.max() < 350
        ax.set_xlim(0, 350)
        ax.set_xlabel("reversing target edges (count)")
        ax.set_title("Two-sided Fisher priorities\n(20 labels from a; sorted by count)")
        if tag:
            S.panel_tag(ax, tag, dx=-0.18, dy=1.04)
    
    
    def panel_rate(ax, tag="b", show_note=True, detail_scale=1.0, compact=False):
        """Inversion RATE vs out-degree (the confound defused)."""
        ax.scatter(rr.n_both, rr.rate * 100, s=14 * detail_scale ** 2,
                   color="#c7ccd1", edgecolor="#9aa0a6", lw=0.3, zorder=2)
        ax.axhline(base * 100, ls="--", lw=1.1, color="#565b61", zorder=1)
        ax.text(rr.n_both.max()*0.92, base*100 + 0.4, f"census baseline {base*100:.1f}%",
                ha="right", fontsize=7 * detail_scale, color="#565b61")
        ENR = {"CCNC":"Mediator", "AHR":"bHLH-PAS TF", "ARNT":"bHLH-PAS TF",
               "NCKAP1L":"cytoskeletal / immune-synapse", "DOCK2":"cytoskeletal / immune-synapse",
               "CBFB":"signalling / activation TF", "STAT5B":"signalling / activation TF",
               "CREBBP":"CBP/p300 (HAT)", "ARHGAP30":"cytoskeletal / immune-synapse"}
        CCOL = {"CCNC":"#2c6fbb"}
        OFF = {"CCNC":(7,8,"left"), "AHR":(-7,8,"right"), "ARNT":(9,6,"left"), "ARHGAP30":(-9,6,"right"),
               "NCKAP1L":(10,4,"left"), "DOCK2":(-10,4,"right"), "CBFB":(-9,10,"right"),
               "STAT5B":(11,-12,"left"), "CREBBP":(10,-1,"left")}
        if compact:
            OFF.update({"STAT5B": (-8, -16, "right"), "CREBBP": (10, 7, "left")})
        for g, cls in ENR.items():
            if g not in nboth: continue
            x, y = nboth[g], rate[g]*100; c = CCOL.get(g, col(cls))
            dx, dy, ha = OFF[g]
            ax.scatter(x, y, s=42 * detail_scale ** 2,
                       color=c, edgecolor="white", lw=0.8, zorder=4)
            ax.annotate(g, (x, y), xytext=(dx * detail_scale, dy * detail_scale),
                        textcoords="offset points", ha=ha, fontsize=7 * detail_scale,
                        fontweight="bold", color=c, zorder=5)
        BASE_HD = {"MED12":"Mediator", "ATAD5":"other (replication/translation/misc)",
                   "SGF29":"SAGA/co-activator (HAT)", "SUPT20H":"SAGA/co-activator (HAT)"}
        for g, cls in BASE_HD.items():
            if g not in nboth: continue
            ax.scatter(nboth[g], rate[g]*100, s=42 * detail_scale ** 2,
                       color=col(cls), edgecolor="white", lw=0.8, zorder=4)
        if show_note:
            ax.annotate("high out-degree,\nbaseline rate:\nMED12, ATAD5,\nSGF29, SUPT20H",
                        xy=(1650, 5.6), xycoords="data",
                        xytext=((0.68, 0.86) if compact else (0.79, 0.86)),
                        textcoords="axes fraction", fontsize=6.8 * detail_scale,
                        color="#565b61", ha="left", va="center",
                        arrowprops=dict(arrowstyle="->", lw=0.9, color="#8a8f98"))
        ax.set_xscale("log")
        ax.set_xlabel("out-degree  (both-significant edges, log scale)")
        ax.set_ylabel("sign-reversal rate (%)")
        ax.set_title("Reversal rate is not set\nby out-degree" if compact
                     else "Reversal rate is not set by out-degree")
        ax.set_ylim(0, max(rr.rate)*100*1.08)
        if tag:
            S.panel_tag(ax, tag, dx=(-0.10 if compact else -0.13), dy=1.04)
    
    
    for k, v in list(locals().items()):
        setattr(fig4_architecture, k, v)
init_fig4_architecture()

def init_fig_target_activation_trajectory_v1():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Reversal scores and independent CD4 activation trajectories.
    
    Version 1.0, 2026-08-14.
    
    This is a standalone biological figure. It deliberately does not modify or compose
    Figures 1-3. It reads only the sealed GSE140244 target-gating outputs and writes:
    
        figures/Fig_target_activation_trajectory.pdf
        figures/Fig_target_activation_trajectory.png
    
    Panel a defines the reversal-specific question and constrained null. Panel b anchors
    inference with the all-target signed-score association and its two separately reported
    null distributions. Panel c shows the absolute 0-to-48 h donor medians (both groups rise)
    and raw/stratified donor AUCs in the unanimous subset. Panel d shows its time course.
    
    The >=5-reversal stratified attenuation is printed on the figure rather than hidden.
    The target-conditioned randomisation is the inferential anchor. The donor AUC panels are
    external reproducibility summaries, not the reversal-specific test.
    """
    
    
    
    # Headless/vector-safe rendering on macOS, CI and HPC nodes.
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    RESULT_DIR = PAPER / "data" / "target_gating"
    OUT_DIR = PAPER / "figures_final"
    
    S = figstyle
    
    
    RESULT_JSON = RESULT_DIR / "gse140244_target_gating_results.json"
    DONOR_TSV = RESULT_DIR / "gse140244_donor_target_gating.tsv"
    TIME_TSV = RESULT_DIR / "gse140244_timecourse_donor_auc.tsv"
    SPECIFICITY_JSON = RESULT_DIR / "gse140244_reversal_specificity_results.json"
    SPECIFICITY_TARGETS = RESULT_DIR / "gse140244_reversal_specificity_targets.tsv"
    SPECIFICITY_HIST = RESULT_DIR / "gse140244_reversal_specificity_null_histogram.tsv"
    
    OUT_PDF = OUT_DIR / "Fig_target_activation_trajectory.pdf"
    OUT_PNG = OUT_DIR / "Fig_target_activation_trajectory.png"
    
    # Existing paper palette: colour denotes reversal class, never activation direction.
    C_LOST = S.PAL["LOST"]
    C_GAINED = S.PAL["GAINED"]
    C_RAW = "#33475b"
    C_STRAT = "#d39b35"
    C_INK = S.PAL["INK"]
    C_MUTED = "#687078"
    C_LIGHT = "#dfe2e6"
    C_PANEL = "#f6f7f8"
    
    
    def sha256(path: Path) -> str:
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(block)
        return h.hexdigest()
    
    
    def load_and_validate():
        with RESULT_JSON.open() as handle:
            result = json.load(handle)
        donors = pd.read_csv(DONOR_TSV, sep="\t")
        timecourse = pd.read_csv(TIME_TSV, sep="\t")
        with SPECIFICITY_JSON.open() as handle:
            specificity = json.load(handle)
        specificity_targets = pd.read_csv(SPECIFICITY_TARGETS, sep="\t")
        specificity_hist = pd.read_csv(SPECIFICITY_HIST, sep="\t")
    
        primary = result["primary_and_sensitivities"]["minimum_2"]
        d2 = donors.loc[donors["minimum_reversals"] == 2].copy()
    
        # Fail loudly if a later rerun changes the sealed quantities this figure reports.
        assert result["mapping"]["frozen_targets"] == 633
        assert result["mapping"]["mapped_targets"] == 629
        assert len(d2) == 22 and d2["donor"].nunique() == 22
        assert primary["n_lost_predicted_induced"] == 351
        assert primary["n_gained_predicted_repressed"] == 278
        assert primary["donor_auc"]["positive"] == 22
        assert primary["donor_stratified_auc"]["positive"] == 22
        assert np.isclose(d2["auc"].median(), primary["donor_auc"]["median"])
        assert np.isclose(
            d2["stratified_auc"].median(),
            primary["donor_stratified_auc"]["median"],
        )
        assert (d2["median_delta_lost"] > 0).all()
        assert (d2["median_delta_gained"] > 0).all()
        assert sorted(timecourse["time"].unique().tolist()) == [0, 2, 4, 8, 12, 24, 48, 72]
        sp = specificity["primary_target_conditioned"]["result"]
        assert specificity["universe"]["mapped_targets"] == 903
        assert len(specificity_targets) == 903
        assert np.isclose(
            specificity_targets[["reversal_signed_score", "median_delta_0_to_48"]]
            .corr(method="spearman").iloc[0, 1],
            sp["observed_spearman_rho"],
        )
        expected_hist = {
            ("target_conditioned", "raw_tmm"): 100000,
            ("target_conditioned", "author_normalized"): 100000,
            ("regulator_signed_margin", "raw_tmm"): 10000,
            ("regulator_signed_margin", "author_normalized"): 10000,
        }
        for key, n_draws in expected_hist.items():
            part = specificity_hist.loc[
                (specificity_hist["null_model"] == key[0])
                & (specificity_hist["normalization"] == key[1])
            ]
            assert len(part) == 80 and int(part["count"].sum()) == n_draws
            assert part["total_draws"].eq(n_draws).all()
        return (
            result,
            d2.sort_values("donor"),
            timecourse.sort_values(["donor", "time"]),
            specificity,
            specificity_targets,
            specificity_hist,
        )
    
    
    def rounded_box(ax, xy, width, height, text, face, edge, fontsize=7.2,
                    text_color=C_INK, linewidth=1.0):
        box = FancyBboxPatch(
            xy,
            width,
            height,
            boxstyle="round,pad=0.014,rounding_size=0.025",
            transform=ax.transAxes,
            facecolor=face,
            edgecolor=edge,
            linewidth=linewidth,
        )
        ax.add_patch(box)
        ax.text(
            xy[0] + width / 2,
            xy[1] + height / 2,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=fontsize,
            color=text_color,
            linespacing=1.25,
        )
    
    
    def panel_hypothesis(ax):
        """Define the all-target score and the reversal-specific conditional null."""
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        S.panel_tag(ax, "a", dx=-0.035, dy=1.045)
        ax.text(
            0.0, 1.01, "Does reversal selection track target activation?",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=9.0,
            fontweight="bold", color=C_INK,
        )
    
        # A stylised target with a mixture of reversal and stable edges.
        ax.text(0.03, 0.91, "Zhu perturbation atlas", transform=ax.transAxes,
                fontsize=7.0, fontweight="bold", color=C_INK)
        target_xy = (0.31, 0.69)
        ax.scatter([target_xy[0]], [target_xy[1]], s=190, facecolor="#ffffff",
                   edgecolor=C_INK, lw=1.2, transform=ax.transAxes, zorder=5)
        ax.text(target_xy[0], target_xy[1] - 0.075, "one target", transform=ax.transAxes,
                ha="center", va="top", fontsize=6.2, fontweight="bold", color=C_INK)
        edges = [
            ((0.03, 0.83), C_LOST, "lost"),
            ((0.03, 0.70), C_LOST, "lost"),
            ((0.03, 0.57), C_GAINED, "gained"),
            ((0.57, 0.82), "#b9bec3", "stable"),
            ((0.59, 0.61), "#b9bec3", "stable"),
        ]
        for (x0, y0), colour, label in edges:
            ax.add_patch(FancyArrowPatch(
                (x0, y0), target_xy, transform=ax.transAxes,
                arrowstyle="-|>", mutation_scale=7, linewidth=1.4, color=colour,
                alpha=0.95,
            ))
            ax.text(x0, y0 + 0.027, label, transform=ax.transAxes, ha="center",
                    fontsize=5.9, color=colour)
    
        rounded_box(
            ax, (0.04, 0.37), 0.55, 0.14,
            "signed reversal score =\n(lost - gained) / all reversals",
            face="#ffffff", edge="#aeb4ba", fontsize=6.8,
        )
        rounded_box(
            ax, (0.66, 0.62), 0.31, 0.24,
            "Independent GSE140244\nactivation trajectory\n22 donors; 0-72 h\nmedian delta at 48 h",
            face=C_PANEL, edge=C_RAW, fontsize=6.5,
        )
        ax.add_patch(FancyArrowPatch(
            (0.60, 0.69), (0.65, 0.69), transform=ax.transAxes,
            arrowstyle="-|>", mutation_scale=8, linewidth=1.0, color=C_MUTED,
        ))
    
        rounded_box(
            ax, (0.04, 0.08), 0.92, 0.20,
            "Primary conditional null: within each target, reassign which edges reverse\n"
            "while preserving target identity, in-degree, reversal count, full Rest-sign\n"
            "composition and the independent trajectory (903/910 targets mapped)",
            face="#f6f7f8", edge="#aeb4ba", fontsize=6.3,
        )
        ax.add_patch(FancyArrowPatch(
            (0.50, 0.36), (0.50, 0.29), transform=ax.transAxes,
            arrowstyle="-|>", mutation_scale=8, linewidth=1.0, color=C_MUTED,
        ))
        ax.text(
            0.50, 0.015,
            "The test asks which edges reverse, beyond a target's ordinary Rest-sign bias",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=6.2,
            color=C_MUTED, fontstyle="italic",
        )
    
    
    def panel_specificity_relation(ax, targets, specificity):
        """Observed all-target relation; dots are displayed, not treated as IID replicates."""
        x = targets["reversal_signed_score"].to_numpy()
        y = targets["median_delta_0_to_48"].to_numpy()
        colours = np.where(x > 0, C_LOST, np.where(x < 0, C_GAINED, "#9da3a8"))
        ax.scatter(x, y, s=10, c=colours, alpha=0.25, edgecolor="none", rasterized=True)
        ax.axvline(0, color="#8a9095", lw=0.8, ls="--", zorder=0)
        ax.axhline(0, color="#8a9095", lw=0.8, ls="--", zorder=0)
    
        # Five fixed display bins summarise the monotone relation without redefining the test.
        bins = [(-1.01, -0.999, "-1"), (-0.999, -1e-12, "mixed -"),
                (-1e-12, 1e-12, "0"), (1e-12, 0.999, "mixed +"),
                (0.999, 1.01, "+1")]
        bx, med, lo, hi = [], [], [], []
        for left, right, _ in bins:
            keep = (x >= left) & (x <= right)
            values = y[keep]
            bx.append(np.median(x[keep]))
            med.append(np.median(values))
            lo.append(np.quantile(values, 0.25))
            hi.append(np.quantile(values, 0.75))
        ax.plot(bx, med, color=C_INK, lw=1.5, zorder=4)
        ax.errorbar(bx, med, yerr=[np.array(med) - np.array(lo),
                                   np.array(hi) - np.array(med)],
                    fmt="o", ms=4.2, color=C_INK, ecolor=C_INK, elinewidth=1.2,
                    capsize=2.0, zorder=5)
    
        primary = specificity["primary_target_conditioned"]["result"]
        regulator = specificity["regulator_margin_sensitivity"]["raw_tmm"]
        ax.set_xlim(-1.08, 1.08)
        ax.set_ylim(-6.2, 10.3)
        ax.set_xlabel("target signed reversal score\nmore gained  <  0  >  more lost")
        ax.set_ylabel("target activation response\nmedian log2 CPM (48 h - 0 h)")
        ax.set_title("Observed all-target relation", fontsize=8.3, pad=7)
        # rho never appears without the null centres it is measured against. Read against zero it
        # invites "0.099 is nothing"; the nulls sit at -0.214 and -0.107, so the displacement is
        # the result, not the coefficient. Both permutation p-values are given here as well as on
        # the null panel, because a reader can quote this panel on its own.
        ax.text(
            0.03, 0.97,
            f"n = {primary['n_targets']} targets,  Spearman rho = "
            f"{primary['observed_spearman_rho']:.3f}\n"
            f"target-conditioned null {primary['null_mean']:.3f},  "
            f"{_format_p_mc(primary['plus_one_one_sided_p'])}\n"
            f"regulator-margin null {regulator['null_mean']:.3f},  "
            f"{_format_p_mc(regulator['plus_one_one_sided_p'])}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2,
            color=C_MUTED, linespacing=1.5,
        )
        S.panel_tag(ax, "b", dx=-0.13, dy=1.055)
    
    
    def _hist_series(hist, null_model):
        part = hist.loc[
            (hist["null_model"] == null_model)
            & (hist["normalization"] == "raw_tmm")
        ].sort_values("bin_left")
        edges = np.r_[part["bin_left"].to_numpy(), part["bin_right"].iloc[-1]]
        width = part["bin_right"].to_numpy() - part["bin_left"].to_numpy()
        density = part["count"].to_numpy() / part["total_draws"].iloc[0] / width
        return edges, density
    
    
    def _format_p_mc(value):
        """Compact mathtext for Monte Carlo plus-one p-values near powers of ten."""
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10 ** exponent)
        if np.isclose(coefficient, 10.0, rtol=0.01):
            exponent += 1
            coefficient_text = "1"
        elif np.isclose(coefficient, 1.0, rtol=0.01):
            coefficient_text = "1"
        else:
            coefficient_text = f"{coefficient:.2g}"
        return rf"$p_{{MC}} = {coefficient_text} \times 10^{{{exponent}}}$"
    
    
    def panel_specificity_null(ax, hist, specificity):
        """Empirical conditional-null distributions and the observed rho."""
        tc_edges, tc_density = _hist_series(hist, "target_conditioned")
        rg_edges, rg_density = _hist_series(hist, "regulator_signed_margin")
        ax.stairs(tc_density, tc_edges, fill=True, color=C_RAW, alpha=0.45, lw=1.1,
                  label="target-conditioned null\n(100,000 draws)")
        ax.stairs(rg_density, rg_edges, fill=True, color=C_STRAT, alpha=0.38, lw=1.1,
                  label="regulator-margin null\n(10,000 draws; sensitivity)")
    
        primary = specificity["primary_target_conditioned"]["result"]
        regulator = specificity["regulator_margin_sensitivity"]["raw_tmm"]
        observed = primary["observed_spearman_rho"]
        ax.axvline(observed, color=C_LOST, lw=2.2, zorder=6)
        ax.text(observed - 0.008, ax.get_ylim()[1] * 0.92,
                f"observed\nrho = {observed:.3f}", ha="right", va="top",
                fontsize=6.7, fontweight="bold", color=C_LOST)
        ax.axvline(primary["null_mean"], color=C_RAW, lw=1.0, ls="--")
        ax.axvline(regulator["null_mean"], color="#a36e1e", lw=1.0, ls="--")
        ax.set_xlim(-0.33, 0.14)
        ax.set_xlabel("Spearman rho under constrained null")
        ax.set_ylabel("null density")
        ax.set_title("Weak association exceeds both nulls", fontsize=8.3, pad=7)
        ax.legend(loc="upper left", fontsize=6.1, handlelength=1.7,
                  borderaxespad=0.3, labelspacing=0.5)
        ax.text(
            0.98, 0.04,
            f"target null: {_format_p_mc(primary['plus_one_one_sided_p'])}\n"
            f"regulator-margin: {_format_p_mc(regulator['plus_one_one_sided_p'])}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.2,
            color=C_MUTED,
        )
    
    
    def add_median_iqr(ax, x, values, colour):
        q25, med, q75 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.plot([x, x], [q25, q75], color=C_INK, lw=2.2, zorder=5)
        ax.plot([x - 0.09, x + 0.09], [q25, q25], color=C_INK, lw=1.0, zorder=5)
        ax.plot([x - 0.09, x + 0.09], [q75, q75], color=C_INK, lw=1.0, zorder=5)
        ax.scatter([x], [med], s=55, marker="D", facecolor=colour, edgecolor="white",
                   lw=0.8, zorder=6)
        return q25, med, q75
    
    
    def panel_absolute(ax, donors):
        """Show donor-level absolute response so the relative AUC is not misread."""
        x_lost, x_gained = 0, 1
        for row in donors.itertuples(index=False):
            ax.plot(
                [x_lost, x_gained],
                [row.median_delta_lost, row.median_delta_gained],
                color="#b8bdc2",
                lw=0.65,
                alpha=0.62,
                zorder=1,
            )
        ax.scatter(
            np.full(len(donors), x_lost), donors["median_delta_lost"],
            s=21, facecolor=C_LOST, edgecolor="white", lw=0.45, alpha=0.83, zorder=3,
        )
        ax.scatter(
            np.full(len(donors), x_gained), donors["median_delta_gained"],
            s=21, facecolor=C_GAINED, edgecolor="white", lw=0.45, alpha=0.83, zorder=3,
        )
        _, med_lost, _ = add_median_iqr(ax, x_lost, donors["median_delta_lost"], C_LOST)
        _, med_gained, _ = add_median_iqr(ax, x_gained, donors["median_delta_gained"], C_GAINED)
    
        ax.axhline(0, color="#777d82", lw=0.9, ls="--", zorder=0)
        ax.set_xticks([x_lost, x_gained])
        ax.set_xticklabels(["lost-edge targets\n(n = 351)", "gained-edge targets\n(n = 278)"],
                           fontsize=7.0)
        ax.set_xlim(-0.38, 1.38)
        ax.set_ylim(-0.08, 1.06)
        ax.set_ylabel("median target response per donor\nlog2 CPM (48 h - 0 h)", labelpad=2)
        ax.set_title("Both unanimous groups rise at 48 h", fontsize=8.2, pad=8)
        ax.text(
            0.5,
            0.98,
            "both group medians > 0\nin all 22 donors",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.4,
            color=C_MUTED,
        )
        ax.text(x_lost - 0.05, med_lost + 0.10, f"{med_lost:.3f}", ha="right",
                va="center", fontsize=6.5, fontweight="bold", color=C_LOST)
        ax.text(x_gained + 0.05, med_gained - 0.09, f"{med_gained:.3f}", ha="left",
                va="center", fontsize=6.5, fontweight="bold", color=C_GAINED)
        S.panel_tag(ax, "c", dx=-0.20, dy=1.055)
    
    
    def panel_primary_auc(ax, donors, result):
        """Display donor consistency in the unanimous subset."""
        raw = donors["auc"].to_numpy()
        strat = donors["stratified_auc"].to_numpy()
        for a, b in zip(raw, strat):
            ax.plot([0, 1], [a, b], color="#b8bdc2", lw=0.65, alpha=0.62, zorder=1)
        ax.scatter(np.zeros(len(raw)), raw, s=21, facecolor=C_RAW, edgecolor="white",
                   lw=0.45, alpha=0.88, zorder=3)
        ax.scatter(np.ones(len(strat)), strat, s=22, marker="s", facecolor=C_STRAT,
                   edgecolor="white", lw=0.45, alpha=0.90, zorder=3)
        _, med_raw, _ = add_median_iqr(ax, 0, raw, C_RAW)
        _, med_strat, _ = add_median_iqr(ax, 1, strat, C_STRAT)
    
        ax.axhline(0.5, color="#777d82", lw=0.95, ls="--", zorder=0)
        ax.text(1.38, 0.5015, "null 0.5", fontsize=6.1, color=C_MUTED,
                ha="right", va="bottom")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["raw AUC", "stratified\nAUC"], fontsize=7.0)
        ax.set_xlim(-0.38, 1.42)
        ax.set_ylim(0.492, 0.648)
        ax.set_ylabel("donor AUC\nP(delta_lost > delta_gained)")
        ax.set_title("Donor rank separation", fontsize=8.2, pad=8)
        ax.text(0, med_raw + 0.014, f"median {med_raw:.3f}", ha="center", va="bottom",
                fontsize=6.3, fontweight="bold", color=C_RAW)
        ax.text(1, med_strat + 0.014, f"{med_strat:.3f}", ha="center", va="bottom",
                fontsize=6.3, fontweight="bold", color="#9a641d")
        p = result["primary_and_sensitivities"]["minimum_2"]["donor_auc"][
            "exact_one_sided_sign_p"
        ]
        ax.text(
            0.5,
            0.98,
            f"unanimous subset: 629/633 mapped\n22/22 donors above 0.5; sign p = {p:.2g}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.4,
            color=C_MUTED,
        )
    
    
    def panel_timecourse(ax, timecourse, result):
        """Show all donor AUC trajectories and the median/IQR time course."""
        hours = [0, 2, 4, 8, 12, 24, 48, 72]
        xpos = {hour: i for i, hour in enumerate(hours)}
        pivot = timecourse.pivot(index="donor", columns="time", values="auc")
        for _, row in pivot.iterrows():
            yy = np.array([row.get(hour, np.nan) for hour in hours], dtype=float)
            ax.plot(range(len(hours)), yy, color="#b6bbc0", lw=0.55, alpha=0.42, zorder=1)
            ax.scatter(range(len(hours)), yy, s=7, color="#b6bbc0", alpha=0.42,
                       edgecolor="none", zorder=1)
    
        summaries = []
        for hour in hours:
            values = timecourse.loc[timecourse["time"] == hour, "auc"].dropna().to_numpy()
            summaries.append((np.median(values), np.quantile(values, 0.25),
                              np.quantile(values, 0.75), len(values)))
        med = np.array([x[0] for x in summaries])
        q25 = np.array([x[1] for x in summaries])
        q75 = np.array([x[2] for x in summaries])
        xx = np.arange(len(hours))
    
        ax.axvspan(xpos[4] - 0.28, xpos[8] + 0.28, color="#eff3f5", zorder=0)
        ax.axvspan(xpos[48] - 0.18, xpos[48] + 0.18, color="#f7ecd2", zorder=0)
        ax.fill_between(xx, q25, q75, color=C_RAW, alpha=0.16, linewidth=0, zorder=2)
        ax.plot(xx, med, color=C_RAW, lw=2.1, marker="o", ms=4.4, zorder=4)
        ax.scatter([xpos[48]], [med[xpos[48]]], s=65, facecolor=C_STRAT,
                   edgecolor="white", lw=0.9, zorder=6)
        ax.axhline(0.5, color="#777d82", lw=0.95, ls="--", zorder=0)
        ax.text(7.12, 0.5015, "null 0.5", fontsize=6.1, color=C_MUTED,
                ha="right", va="bottom")
    
        ax.set_xticks(xx)
        ax.set_xticklabels(["0", "2", "4", "8", "12\n(n=19)", "24", "48", "72"])
        ax.set_xlim(-0.2, 7.2)
        ax.set_ylim(0.485, 0.735)
        ax.set_xlabel("hours after stimulation")
        ax.set_ylabel("donor AUC (median and IQR)")
        ax.set_title("Unanimous-subset separation peaks at 4-8 h", fontsize=8.5, pad=8)
        ax.text(
            xpos[8],
            med[xpos[8]] + 0.018,
            f"peak median {med[xpos[8]]:.3f}",
            ha="center",
            va="bottom",
            fontsize=6.4,
            fontweight="bold",
            color=C_RAW,
        )
        ax.annotate(
            f"48 h subset AUC\nmedian {med[xpos[48]]:.3f}",
            xy=(xpos[48], med[xpos[48]]),
            xytext=(xpos[48] + 0.45, med[xpos[48]] + 0.055),
            fontsize=6.3,
            color="#8a5b16",
            arrowprops=dict(arrowstyle="-", color="#b37a24", lw=0.8),
            ha="left",
            va="center",
        )
        S.panel_tag(ax, "d", dx=-0.08, dy=1.055)
    
    
    def main():
        S.setup(scale=1.0)
        (result, donors, timecourse, specificity,
         specificity_targets, specificity_hist) = load_and_validate()
        OUT_DIR.mkdir(parents=True, exist_ok=True)
    
        fig = plt.figure(figsize=(13.2, 8.35), facecolor="white")
        outer = fig.add_gridspec(
            2,
            2,
            width_ratios=[1.02, 1.28],
            height_ratios=[1.03, 0.97],
            left=0.065,
            right=0.985,
            top=0.84,
            bottom=0.20,
            wspace=0.24,
            hspace=0.38,
        )
        panel_hypothesis(fig.add_subplot(outer[0, 0]))
        anchor = outer[0, 1].subgridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.34)
        panel_specificity_relation(fig.add_subplot(anchor[0, 0]), specificity_targets,
                                   specificity)
        panel_specificity_null(fig.add_subplot(anchor[0, 1]), specificity_hist, specificity)
        donor = outer[1, 0].subgridspec(1, 2, wspace=0.40)
        panel_absolute(fig.add_subplot(donor[0, 0]), donors)
        panel_primary_auc(fig.add_subplot(donor[0, 1]), donors, result)
        panel_timecourse(fig.add_subplot(outer[1, 1]), timecourse, result)
    
        fig.text(
            0.5,
            0.965,
            "Reversal selection tracks an independent CD4 activation trajectory",
            ha="center",
            va="top",
            fontsize=14,
            fontweight="bold",
            color=C_INK,
        )
        fig.text(
            0.5,
            0.918,
            "The all-target conditional test anchors inference; unanimous-target panels show donor reproducibility and timing",
            ha="center",
            va="top",
            fontsize=8.0,
            color=C_MUTED,
        )
    
        sens5 = result["primary_and_sensitivities"]["minimum_5"]
        fig.text(
            0.5,
            0.103,
            "Interpretation boundary: PASS a weak relative-trajectory association under both separately reported nulls; "
            "FAIL a literal induced-versus-repressed split\n"
            "because both unanimous target-group medians increase at 48 h; association/organisation only, "
            "not causation or molecular mechanism.",
            ha="center",
            va="center",
            fontsize=7.2,
            color=C_INK,
            linespacing=1.35,
            wrap=True,
            bbox=dict(boxstyle="round,pad=0.48", facecolor=C_PANEL,
                      edgecolor="#d2d6da", linewidth=0.8),
        )
        fig.text(
            0.5,
            0.043,
            "The target null fixes each target's reversal count and full Rest-sign composition but not regulator margins; "
            "the complementary regulator-margin null does not fix target reversal counts.\n"
            "Neither controls arbitrary dyad-specific QC. "
            "For targets with >=5 reversals, raw AUC remains above 0.5 in "
            f"{sens5['donor_auc']['positive']}/22 donors, but the stratified estimate attenuates "
            f"(median {sens5['donor_stratified_auc']['median']:.3f}; "
            f"{sens5['donor_stratified_auc']['positive']}/22 above 0.5; exact sign p = "
            f"{sens5['donor_stratified_auc']['exact_one_sided_sign_p']:.3f}).",
            ha="center",
            va="center",
            fontsize=6.2,
            color=C_MUTED,
            linespacing=1.25,
            wrap=True,
        )
        fig.text(
            0.985,
            0.009,
            "GSE140244; TMM log2 CPM; donor summaries show every donor; diamonds and band show median and IQR",
            ha="right",
            va="bottom",
            fontsize=5.8,
            color="#7b8288",
        )
    
        fig.savefig(OUT_PDF, dpi=300, facecolor="white")
        fig.savefig(OUT_PNG, dpi=300, facecolor="white")
        plt.close(fig)
    
        print(f"wrote {OUT_PDF}")
        print(f"wrote {OUT_PNG}")
        print(f"input result JSON sha256: {sha256(RESULT_JSON)}")
        print(f"input specificity JSON sha256: {sha256(SPECIFICITY_JSON)}")
        print(f"input specificity targets sha256: {sha256(SPECIFICITY_TARGETS)}")
        print(f"input specificity null histogram sha256: {sha256(SPECIFICITY_HIST)}")
        print(f"input donor table sha256: {sha256(DONOR_TSV)}")
        print(f"input time-course table sha256: {sha256(TIME_TSV)}")
        print(f"output PDF sha256: {sha256(OUT_PDF)}")
        print(f"output PNG sha256: {sha256(OUT_PNG)}")
    
    
    for k, v in list(locals().items()):
        setattr(fig_target_activation_trajectory_v1, k, v)
init_fig_target_activation_trajectory_v1()

def init_fig2g_trajectory_panel():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Figure 2g: the one-direction regulators' reversal targets across Rest, 8 h and 48 h.
    
    Compact re-render of the diagnostic in
    data/figure_provenance/fig2_r49_direction_pure_trajectory_heatmap, sized for Figure 2's
    bottom row. The standalone renderer builds its own page with its own guard log and is left
    untouched; this module reads the same census and draws the same content into a supplied axes.
    
    CONTENT. The six regulators inside the 49-regulator enriched family whose reversals all run one
    way, and every reversal target they carry: 39 edges on 39 targets, 23 lost and 16 gained. Two
    blocks, lost on the left and gained on the right, three columns each for Rest, 8 h and 48 h.
    
    WHAT IS FIXED AND WHAT IS NOT. Rest and 48 h define the reversal and its direction, so the
    colour flip across those two columns is the selection rule. The 8 h column took no part in
    selection and is the only free one, which is why the marks sit there: a filled dot for the two
    edges reaching FDR below 0.10 and a cross for the one whose adjusted p is unavailable.
    
    THE 8 H COLLAPSE IS REAL BUT NARROW. 37 of 39 effects are smaller at 8 h than at both endpoints.
    That alone proves nothing, because a straight line from a positive Rest value to a negative 48 h
    value already passes near zero: interpolating these same edges predicts 23 of 39. Observed 8 h
    magnitudes are nonetheless closer to zero than that line (median 0.092 against 0.274, Wilcoxon
    p = 2.4e-08), so these edges do collapse sub-linearly. Census-wide there is no such effect
    (observed 0.753 against 0.732 predicted), so this belongs to these six regulators and rests on
    39 edges from six post-hoc rows. It is not evidence that reversals in general are non-monotonic.
    
    Writes data/census/fig2g_trajectory_panel.json
    """
    
    
    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    
    
    MASTER = ROOT / "data/census/census_master_edges.csv.gz"
    OUTJSON = ROOT / "data/census/fig2g_trajectory_panel.json"
    
    SIX_LOST = ["MCAT", "TIMM8B", "RB1CC1"]
    SIX_GAINED = ["CARD11", "NDUFA8", "SUGT1"]
    STATES = [("R", "lfc_Rest"), ("8", "lfc_Stim8hr"), ("48", "lfc_Stim48")]
    ANNOT = {"MCAT": "mito fatty-acid synth.", "TIMM8B": "mito import", "RB1CC1": "autophagy",
             "CARD11": "NF-kB scaffold", "NDUFA8": "complex I", "SUGT1": "HSP90 co-chaperone"}
    FDR = 0.10
    LOST, GAINED, INK, MUTED = "#c1272d", "#2a9d8f", "#22262B", "#666B70"
    
    
    def sha256(path: Path) -> str:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    
    
    def load() -> pd.DataFrame:
        m = pd.read_csv(MASTER)
        m = m.loc[~m["is_self"].astype(bool)]
        six = SIX_LOST + SIX_GAINED
        e = m.loc[m["regulator"].isin(six) & m["kind"].eq("rev")].copy()
        assert len(e) == 39 and e["target"].nunique() == 39
        assert (e.groupby("regulator")["direction"].nunique() == 1).all()
        e["q8"] = pd.to_numeric(e["padj_Stim8hr"], errors="coerce")
        e["block"] = np.where(e["regulator"].isin(SIX_LOST), "lost", "gained")
        order = {r: i for i, r in enumerate(six)}
        return e.sort_values(["block", "regulator", "target"],
                             key=lambda s: s.map(order) if s.name == "regulator" else s,
                             ascending=[False, True, True]).reset_index(drop=True)
    
    
    def panel(ax, edges: pd.DataFrame, fs: float = 7.0, vlim: float = 2.0) -> dict:
        """Draw the two blocks into ax, in the Figure 3d/3e idiom.
    
        Tile geometry is set by the x-unit scale rather than left to the axes aspect. One tile is one
        x-unit, and the label gutters are quoted in the same units, so widening the gutters narrows
        the tiles instead of stretching them. The earlier version used a 15-unit range across the full
        page, which gave tiles about five times wider than tall; here 27 units across two thirds of the
        page give the roughly 1.8:1 tiles Figure 3d and 3e use.
        """
        ax.set_axis_off()
        cmap = plt.get_cmap("RdBu_r")
        # x-units. GUTTER is the space each block reserves for its regulator bracket, the regulator
        # label and the target names, all right-aligned into it.
        GUTTER, NCOL = 10.2, len(STATES)
        # Both blocks are top-aligned on the taller one so their headings sit on one line; the
        # gained block simply stops seven rows earlier.
        tallest = max(int(edges["block"].eq(b).sum()) for b in ("lost", "gained"))
        counts = {}
        for side, regs, x0 in (("lost", SIX_LOST, GUTTER),
                               ("gained", SIX_GAINED, 2 * GUTTER + NCOL)):
            rows = []
            for reg in regs:
                sub = edges.loc[edges["regulator"].eq(reg)]
                rows.extend((reg, t) for t in sub["target"])
            counts[side] = len(rows)
            n = len(rows)
            for i, (reg, target) in enumerate(rows):
                y = tallest - 1 - i
                rec = edges.loc[edges["regulator"].eq(reg) & edges["target"].eq(target)].iloc[0]
                for j, (_, col) in enumerate(STATES):
                    v = float(rec[col])
                    ax.add_patch(plt.Rectangle((x0 + j, y), 1.0, 1.0,
                                               facecolor=cmap(0.5 + 0.5 * np.clip(v / vlim, -1, 1)),
                                               edgecolor="white", lw=0.35))
                ax.text(x0 - 0.35, y + 0.5, target, fontsize=fs - 1.2, ha="right", va="center",
                        color=INK)
                if pd.isna(rec["q8"]):
                    ax.plot([x0 + 1.5], [y + 0.5], marker="x", ms=3.0, color=INK, mew=0.9, zorder=4)
                elif rec["q8"] < FDR:
                    ax.plot([x0 + 1.5], [y + 0.5], marker="o", ms=2.6, color=INK, zorder=4)
            # white column dividers between the three states, as in Figure 3d
            for j in range(1, NCOL):
                ax.plot([x0 + j] * 2, [tallest - n, tallest], color="white", lw=0.9, zorder=3)
            # regulator brackets down the left of each block
            top = tallest
            colour = LOST if side == "lost" else GAINED
            for reg in regs:
                k = int(edges["regulator"].eq(reg).sum())
                lo, hi = top - k, top
                ax.plot([x0 - 4.35] * 2, [lo + 0.12, hi - 0.12], color=colour, lw=2.2,
                        solid_capstyle="butt")
                ax.text(x0 - 4.62, (lo + hi) / 2, f"{reg}\n{ANNOT[reg]}", fontsize=fs - 1.0,
                        ha="right", va="center", color=colour, linespacing=1.2)
                top = lo
            for j, (name, _) in enumerate(STATES):
                ax.text(x0 + j + 0.5, tallest + 0.28, name, fontsize=fs - 0.8, ha="center",
                        va="bottom", color=MUTED)
            # Two lines, not one: at this tile width a single-line heading is wider than its block
            # and the right-hand one ran off the axes.
            ax.text(x0 + NCOL / 2, tallest + 1.25, f"all {side}\non stimulation",
                    fontsize=fs, ha="center", va="bottom", color=colour, fontweight="bold",
                    linespacing=1.15)
    
        assert tallest == max(counts.values())
        ax.set_xlim(0.0, 2 * GUTTER + 2 * NCOL + 0.4)
        ax.set_ylim(-0.30, tallest + 3.05)
        return counts
    
    
    def add_colourbar(ax, fs: float = 7.0, vlim: float = 2.0) -> None:
        """Shared scale and marker key, drawn in the block of free rows under the gained column.
    
        Both blocks are top-aligned and the gained one is seven rows shorter, so the rectangle beneath
        it is empty by construction. Keeping the key inside the axes rather than below it means the
        whole slot height goes to the 23 rows, which is what decides whether the target names collide.
        """
        x0, x1 = 16.2, 23.4
        ybar0, ybar1 = 4.35, 5.05
        grad = np.linspace(0, 1, 256).reshape(1, -1)
        ax.imshow(grad, cmap="RdBu_r", aspect="auto", origin="lower",
                  extent=(x0, x1, ybar0, ybar1), zorder=3)
        ax.add_patch(plt.Rectangle((x0, ybar0), x1 - x0, ybar1 - ybar0, fill=False,
                                   edgecolor="#b8bfc5", lw=0.5, zorder=4))
        ax.text((x0 + x1) / 2, ybar1 + 0.22, "net knockdown coefficient (log2 fold change)",
                fontsize=fs - 1.0, ha="center", va="bottom", color=MUTED)
        for frac, label in ((0.0, f"$-${vlim:.0f}"), (0.5, "0"), (1.0, f"+{vlim:.0f}")):
            ax.text(x0 + frac * (x1 - x0), ybar0 - 0.22, label, fontsize=fs - 1.0,
                    ha="center", va="top", color=MUTED)
        ax.plot([x0 + 0.25], [2.75], marker="o", ms=2.6, color=INK, zorder=4)
        ax.text(x0 + 0.75, 2.75, "8 h FDR < 0.10", fontsize=fs - 1.0, ha="left", va="center",
                color=MUTED)
        ax.plot([x0 + 0.25], [1.65], marker="x", ms=3.0, color=INK, mew=0.9, zorder=4)
        ax.text(x0 + 0.75, 1.65, "8 h q unavailable", fontsize=fs - 1.0, ha="left", va="center",
                color=MUTED)
        ax.text(x0, 0.45, "columns R, 8 and 48: Rest, 8 h, 48 h",
                fontsize=fs - 1.0, ha="left", va="center", color=MUTED)
    
    
    def legend_note(edges: pd.DataFrame) -> str:
        a, b = edges["lfc_Rest"].abs(), edges["lfc_Stim48"].abs()
        c = edges["lfc_Stim8hr"].abs()
        lin = (edges["lfc_Rest"] + (edges["lfc_Stim48"] - edges["lfc_Rest"]) * (8 / 48)).abs()
        return (f"{int(((c < a) & (c < b)).sum())} of {len(edges)} effects are smaller at 8 h than at "
                f"both endpoints; straight-line interpolation predicts "
                f"{int(((lin < a) & (lin < b)).sum())}")
    
    
    def main() -> None:
        edges = load()
        # The standalone sets its own canvas, unlike the Figure 2 slot, which carries gridspec
        # margins. At 4.6 in with right=0.99 the gained heading ran past the page edge and was
        # clipped, so the canvas is wider and the right margin real.
        fig, ax = plt.subplots(figsize=(5.1, 3.0), facecolor="white")
        counts = panel(ax, edges, fs=7.0)
        add_colourbar(ax, fs=7.0)
        fig.subplots_adjust(left=0.015, right=0.955, top=0.95, bottom=0.03)
        out = ROOT / "data/figure_provenance/fig2g_trajectory_panel"
        out.mkdir(parents=True, exist_ok=True)
        for e in ("pdf", "png"):
            fig.savefig(out / f"Fig2g_trajectory_panel.{e}", dpi=300, facecolor="white")
        plt.close(fig)
        OUTJSON.write_text(json.dumps({
            "_what": "compact Figure 2g trajectory panel for the six one-direction regulators",
            "_script": "scripts/fig2g_trajectory_panel.py",
            "_claim_boundary": ("Rest and 48 h define the reversal and its direction, so the colour "
                                "flip is the selection rule. Only the 8 h column is free."),
            "inputs": {"census_master_edges.csv.gz": sha256(MASTER)},
            "n_edges": int(len(edges)), "block_sizes": counts,
            "n_8h_fdr_below_0_10": int((edges["q8"] < FDR).sum()),
            "n_8h_q_unavailable": int(edges["q8"].isna().sum()),
            "eight_hour_note": legend_note(edges),
        }, indent=2) + "\n")
        print(counts, legend_note(edges))
        print(f"wrote {out}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2g_trajectory_panel, k, v)
init_fig2g_trajectory_panel()

def init_fig2_biology_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Biology-first main Figure 2 on a full A4 portrait page.
    
    The established regulator analyses are retained in an A4 upper block:
    
      a  all 324 reversal-bearing regulators under the two-sided, degree-matched Fisher test;
      b  all 49 regulators with significantly elevated reversal burden, as horizontal rows;
      c  reversal rate versus out-degree, positioned below a.
    
    The lower block separates three distinct pieces of biological evidence:
    
      d  the 0--72 h donor trajectory for the frozen unanimous target groups;
      e  the absolute disclosure that both target groups rise after stimulation;
      f  every Rest/8 h/48 h reversal target of the six one-direction R49 regulators;
      g  the all-target reversal-orientation/activation relation against two frozen nulls.
    
    The target-conditioned null is a post-primary specificity analysis; the regulator-margin null
    is a separate sensitivity with different preserved margins. Neither establishes mechanism.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    OUT = ROOT / "figures"
    
    S = figstyle
    FA = fig4_architecture
    TG = fig_target_activation_trajectory_v1
    TRJ = fig2g_trajectory_panel
    
    
    SOURCE_WIDTH_IN, SOURCE_HEIGHT_IN = S.A4_PORTRAIT
    MIN_SOURCE_FONT_PT = 7.0
    FULL_PAGE_WIDTH_MM = 210.0
    TWO_COLUMN_WIDTH_MM = 165.1
    ROLE_SOURCE = (PAPER / "data/figure_provenance/fig2_r49_direction_pure_regulator_cards/"
                   "source_regulator_cards.tsv")
    ROLE_SOURCE_SHA256 = "bb40c74527f2e72f5690e08d0f17060929d0e212d932cd82c3f7d8e2fa31bd7b"
    
    ENRICHED_CLASS = {row["gene"]: row["class"] for row in FA.top}
    ENRICHED_CLASS.update(FA.FIG2_CLASS)
    ENRICHED_CLASS_ORDER = [
        "Mediator",
        "SAGA/co-activator (HAT)",
        "CBP/p300 (HAT)",
        "bHLH-PAS TF",
        "signalling / activation TF",
        "cytoskeletal / immune-synapse",
        "RNA processing / surveillance",
        "metabolic / mitochondrial",
        "other (replication/translation/misc)",
        "other / unclassified",
    ]
    ENRICHED_CLASS_SHORT = {
        "Mediator": "Mediator",
        "SAGA/co-activator (HAT)": "SAGA",
        "CBP/p300 (HAT)": "CBP/p300",
        "bHLH-PAS TF": "bHLH-PAS",
        "signalling / activation TF": "signalling",
        "cytoskeletal / immune-synapse": "immune-synapse",
        "RNA processing / surveillance": "RNA processing",
        "metabolic / mitochondrial": "metabolic",
        "other (replication/translation/misc)": "other annotated",
        "other / unclassified": "unclassified",
    }
    
    
    def _enriched_label_colour(category: str) -> str:
        """Use the annotation palette, with a readable grey for unclassified text."""
        return "#6f7780" if category == "other / unclassified" else FA.col(category)
    
    
    def _remove_panel_tag(ax, tag: str) -> None:
        """Remove a hard-coded panel tag from a reused standalone panel function."""
        for artist in list(ax.texts):
            if artist.get_text() == tag and artist.get_fontweight() == "bold":
                artist.remove()
    
    
    def _null_facet(ax, hist, specificity, null_model: str, *, show_xlabel: bool) -> None:
        """Draw one constrained null in its own facet so the p-values are never combined."""
        edges, density = TG._hist_series(hist, null_model)
        target_result = specificity["primary_target_conditioned"]["result"]
        observed = target_result["observed_spearman_rho"]
    
        if null_model == "target_conditioned":
            result = target_result
            colour = TG.C_RAW
            title = "Target-conditioned\nnull (specificity)"
            subtitle = "fixes target count\n+ Rest-sign mix"
            draws = 100_000
            p_text = r"$p_{\mathrm{MC}} = 1 \times 10^{-5}$"
        else:
            result = specificity["regulator_margin_sensitivity"]["raw_tmm"]
            colour = TG.C_STRAT
            title = "Regulator-margin\nnull (sensitivity)"
            subtitle = "fixes regulator\nsigned margins"
            draws = 10_000
            p_text = r"$p_{\mathrm{MC}} = 1 \times 10^{-4}$"
    
        ax.stairs(density, edges, fill=True, color=colour, alpha=0.48, lw=1.0)
        ax.axvline(result["null_mean"], color=colour, lw=1.0, ls="--")
        ax.axvspan(result["null_q025"], result["null_q975"], color=colour,
                   alpha=0.10, lw=0)
        ax.axvline(observed, color=TG.C_LOST, lw=2.0, zorder=5)
        ax.set_xlim(-0.33, 0.14)
        ax.set_yticks([])
        ax.set_ylabel("null density", fontsize=7.0, labelpad=1)
        ax.set_title(f"{title}\n{p_text}", fontsize=7.2, pad=2.0, loc="left",
                     linespacing=1.0)
        ax.text(0.02, 0.88, subtitle, transform=ax.transAxes, ha="left", va="top",
                fontsize=7.0, color=TG.C_MUTED, linespacing=1.02)
        ax.text(
            0.02,
            0.06,
            f"{draws:,} draws",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=7.0,
            color=TG.C_MUTED,
        )
        ax.text(observed - 0.006, 0.94, "observed", transform=ax.get_xaxis_transform(),
                ha="right", va="top", rotation=90, fontsize=7.0,
                fontweight="bold", color=TG.C_LOST)
        if show_xlabel:
            ax.set_xlabel("Spearman rho under null", fontsize=7.0, labelpad=1)
        else:
            ax.tick_params(axis="x", labelbottom=False)
    
    
    def _top_legend(fig) -> None:
        """Separate inferential encodings from descriptive gene-label colours."""
        test_handles = [
            Line2D([0], [0], marker="^", ls="", ms=5.5, color="black",
                   label="above (49)"),
            Line2D([0], [0], marker="v", ls="", ms=5.5, color="black",
                   label="below (57; 42 shown)"),
            Line2D([0], [0], marker="o", ls="", ms=5.0, markerfacecolor="#c7ccd1",
                   markeredgecolor="#9aa0a6", label="not significant"),
            Line2D([0], [0], color="#565b61", lw=1.0, ls="--",
                   label="expectation"),
            Line2D([0], [0], color=S.PAL["LOST"], lw=5,
                   label="lost (+ to -)"),
            Line2D([0], [0], color=S.PAL["GAINED"], lw=5,
                   label="gained (- to +)"),
        ]
        fig.legend(
            handles=test_handles,
            loc="upper center",
            bbox_to_anchor=(0.50, 0.996),
            ncol=6,
            fontsize=7.0,
            handletextpad=0.35,
            columnspacing=0.85,
            frameon=False,
        )
    
        class_handles = [
            Line2D([0], [0], marker="s", ls="", ms=5.2,
                   color=_enriched_label_colour(category),
                   label=ENRICHED_CLASS_SHORT[category])
            for category in ENRICHED_CLASS_ORDER
        ]
        fig.legend(
            handles=class_handles,
            loc="upper center",
            bbox_to_anchor=(0.50, 0.968),
            ncol=5,
            fontsize=7.0,
            handletextpad=0.35,
            columnspacing=0.90,
            labelspacing=0.28,
            frameon=False,
            title="Gene-label colour = functional annotation (descriptive only)",
            title_fontsize=7.0,
        )
    
    
    def _enforce_minimum_font(fig, minimum: float = MIN_SOURCE_FONT_PT) -> tuple[float, list[str]]:
        """Raise every non-empty visible Matplotlib Text object to the publication floor."""
        changed = []
        for artist in fig.findobj(match=Text):
            if not artist.get_visible() or not artist.get_text().strip():
                continue
            if artist.get_fontsize() < minimum:
                changed.append(artist.get_text().replace("\n", " ")[:70])
                artist.set_fontsize(minimum)
        fig.canvas.draw()
        sizes = [
            artist.get_fontsize()
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        return min(sizes), changed
    
    
    def _retune_panel_a(ax) -> None:
        """Remove the former top-20 labels because panel b now names all 49 enriched genes."""
        genes = set(
            FA.fig2a_test.loc[
                FA.fig2a_test.label_rank_fig2a_two_sided > 0, "gene"
            ]
        )
        for artist in list(ax.texts):
            if artist.get_text() in genes:
                artist.remove()
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        ax.set_title("Observed versus out-degree-matched expectation", fontsize=7.6)
    
    
    def _panel_all_enriched_bars(ax, ax_meta, tag: str = "b") -> None:
        """Show every significantly enriched regulator on one horizontal row.
    
        Bar length is the raw reversal burden on a fixed 0--300 count scale. Red and teal segments
        retain the lost-versus-gained split, while the tested denominator, rate and fold remain
        explicit in the adjacent text.
        """
        sub = FA.fig2a_test.loc[
            FA.fig2a_test.sig_two_sided & (FA.fig2a_test.direction == "up")
        ].copy()
        sub = sub.sort_values(
            ["n_inv", "p_two_sided", "gene"], ascending=[False, True, True]
        )
        assert len(sub) == 49 and sub.gene.nunique() == 49
        assert sub.q_two_sided.lt(0.05).all()
        sub = sub.join(FA.direction_counts, on="gene")
        assert (sub.lost_on_activation + sub.gained_on_activation == sub.n_inv).all()
        assert int(sub.n_inv.sum()) == 2_075
        assert int(sub.n_both.sum()) == 14_434
        assert int(sub.lost_on_activation.sum()) == 1_154
        assert int(sub.gained_on_activation.sum()) == 921
    
        yy = np.arange(len(sub))
        lost = sub.lost_on_activation.to_numpy()
        gained = sub.gained_on_activation.to_numpy()
        ax.barh(yy, lost, color=S.PAL["LOST"], height=0.66,
                edgecolor="white", linewidth=0.25)
        ax.barh(yy, gained, left=lost, color=S.PAL["GAINED"], height=0.66,
                edgecolor="white", linewidth=0.25)
        ax.set_yticks(yy)
        ax.set_yticklabels(sub.gene, fontsize=7.0, fontweight="bold", color=S.PAL["INK"])
        for label in ax.get_yticklabels():
            category = ENRICHED_CLASS.get(label.get_text(), "other / unclassified")
            label.set_color(_enriched_label_colour(category))
        ax.set_ylim(len(sub) - 0.35, -0.65)
        assert int(sub.n_inv.max()) == 288
        ax.set_xlim(0, 300)
        ax.set_xticks([0, 100, 200, 300])
        ax.set_xlabel("reversing target edges (count)", fontsize=7.0)
        ax_meta.set_xlim(0, 1)
        ax_meta.set_ylim(len(sub) - 0.35, -0.65)
        ax_meta.axis("off")
        ax_meta.text(0.01, 1.002, "reversals/tested; rate; fold",
                     transform=ax_meta.transAxes, ha="left", va="bottom",
                     fontsize=7.0, fontweight="bold", color="#565b61")
        # Everything that used to sit in the two figure-level legend rows now lives in the panel
        # it explains. Below the eighth row the bars are short, so the block right of x = 100 is
        # empty for four fifths of the panel height: the annotation key and the direction-purity
        # donut go there, and the page reclaims that whole strip at the top.
        n_lost_only = int(((lost > 0) & (gained == 0)).sum())
        n_gained_only = int(((lost == 0) & (gained > 0)).sum())
        n_mixed = len(sub) - n_lost_only - n_gained_only
        assert (n_mixed, n_lost_only, n_gained_only) == (43, 3, 3)
    
        ax.legend(
            handles=[Line2D([0], [0], marker="s", ls="", ms=5.2,
                            color=_enriched_label_colour(category),
                            label=ENRICHED_CLASS_SHORT[category])
                     for category in ENRICHED_CLASS_ORDER],
            # one column, not two: panel b's axes is 1.94 in wide and a two-column key runs
            # past its right edge into the reversals/tested text column beside it.
            loc="lower left", bbox_to_anchor=(0.145, 0.02), ncol=1, fontsize=7.0,
            handletextpad=0.35, labelspacing=0.26, frameon=False,
            title="label colour = annotation\n(descriptive only)",
            title_fontsize=7.0, alignment="left",
        )
    
        donut = ax.inset_axes([0.30, 0.355, 0.26, 0.26])
        counts = [n_mixed, n_lost_only, n_gained_only]
        donut.pie(
            counts, startangle=90, counterclock=False,
            colors=["#c7ccd1", S.PAL["LOST"], S.PAL["GAINED"]],
            wedgeprops=dict(width=0.44, edgecolor="white", linewidth=0.8),
        )
        donut.set_aspect("equal")
        donut.text(0, 0.06, f"{len(sub)}", ha="center", va="center", fontsize=8.5,
                   fontweight="bold", color=S.PAL["INK"])
        donut.text(0, -0.30, "regulators", ha="center", va="center", fontsize=7.0,
                   color="#565b61")
        for text, colour, y in (
                (f"both  {n_mixed} ({100 * n_mixed / len(sub):.0f}%)", "#7b8085", 1.32),
                (f"only lost  {n_lost_only}", S.PAL["LOST"], -1.36),
                (f"only gained  {n_gained_only}", S.PAL["GAINED"], -1.66)):
            donut.text(0, y, text, ha="center", va="center", fontsize=7.0,
                       fontweight="bold", color=colour)
        donut.set_title("bar colour: red lost (+ to -),\nteal gained (- to +)",
                        fontsize=7.0, color="#565b61", fontweight="normal", pad=15)
    
        for i, row in enumerate(sub.itertuples(index=False)):
            ax_meta.text(
                0.01,
                i,
                f"{int(row.n_inv)}/{int(row.n_both)}; {100 * row.rate:.1f}%; "
                f"{row.enrichment_vs_decile:.2f}x",
                ha="left",
                va="center",
                fontsize=7.0,
                color="#3d4348",
            )
    
    
    def _retune_panel_c(ax) -> None:
        """Keep the high-degree note inside the narrow publication-scale axes."""
        for artist in ax.texts:
            if artist.get_text().startswith("high out-degree"):
                artist.set_text("high out-degree, baseline rate:\nMED12, ATAD5,\nSGF29, SUPT20H")
                artist.set_position((0.97, 0.87))
                artist.set_ha("right")
                artist.set_va("center")
                artist.set_fontsize(8.0)
            elif artist.get_text() == "CCNC":
                # CCNC is the rightmost labelled point; the shared panel's outward label
                # is clipped in the three-column composition once fonts reach print size.
                artist.set_position((-8, -10))
                artist.set_ha("right")
            elif artist.get_text() == "AHR":
                # Separate the two nearly level high-degree labels after moving CCNC inward.
                artist.set_position((-10, 11))
                artist.set_ha("right")
            elif artist.get_text() == "ARNT":
                artist.set_position((6, 6))
            elif artist.get_text() == "CREBBP":
                # Three of these labels are placed in offset points, so shortening the axes moves
                # them together in data terms and they collide. CREBBP drops below CCNC, CBFB below
                # ARHGAP30 and STAT5B rises off the census-baseline caption.
                artist.set_position((13, 1))
                artist.set_ha("left")
            elif artist.get_text() == "CBFB":
                artist.set_position((-9, 9))
                artist.set_ha("right")
            elif artist.get_text() == "ARHGAP30":
                artist.set_position((-11, 10))
                artist.set_ha("right")
            elif artist.get_text() == "STAT5B":
                artist.set_position((-9, -10))
                artist.set_ha("right")
            elif artist.get_text().startswith("census baseline"):
                # the long caption reaches back under four of the labelled points
                artist.set_text("baseline 5.2%")
    
    
    def _merged_null_facet(ax, hist, specificity) -> None:
        """Both constrained nulls on one axes.
    
        They share an x axis, they barely overlap, and each previously carried its own copy of the
        same red observed line. One facet says the same thing and returns a square to the figure.
        """
        observed = specificity["primary_target_conditioned"]["result"]["observed_spearman_rho"]
        for model, key, colour, label in (
                ("target_conditioned", "primary_target_conditioned", TG.C_RAW,
                 "target-conditioned\n100,000 draws"),
                ("regulator_signed_margin", "regulator_margin_sensitivity", TG.C_STRAT,
                 "regulator-margin\n10,000 draws")):
            edges, density = TG._hist_series(hist, model)
            result = (specificity[key]["result"] if key == "primary_target_conditioned"
                      else specificity[key]["raw_tmm"])
            ax.stairs(density, edges, fill=True, color=colour, alpha=0.45, lw=1.0, label=label)
            ax.axvline(result["null_mean"], color=colour, lw=1.0, ls="--")
        ax.axvline(observed, color=TG.C_LOST, lw=2.0, zorder=5)
        ax.text(observed - 0.006, 0.97, "observed", transform=ax.get_xaxis_transform(),
                ha="right", va="top", rotation=90, fontsize=7.0, fontweight="bold",
                color=TG.C_LOST)
        ax.set_xlim(-0.33, 0.14)
        ax.set_yticks([])
        ax.set_ylabel("null density", fontsize=7.0, labelpad=1)
        ax.set_xlabel("Spearman rho under null", fontsize=7.0, labelpad=1)
        ax.set_title("Both nulls sit below zero", fontsize=7.2, pad=3, loc="left")
        ax.legend(loc="upper left", fontsize=7.0, frameon=False, handlelength=1.2,
                  labelspacing=0.45, borderaxespad=0.25)
        ax.text(0.02, 0.06,
                r"$p_{\mathrm{MC}} = 1 \times 10^{-5}$ and $1 \times 10^{-4}$",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=7.0, color=TG.C_MUTED)
    
    
    # The eight direction-pure regulators carrying at least six reversals. Purity below six is
    # arithmetic: 143 of the 155 pure regulators in the census carry five or fewer, and 97 carry
    # one. The functional labels are fixed here rather than derived, because they are read from
    # gene descriptions and are display text, not a computed quantity.
    # Short tags, not sentences: this is a square panel and the full descriptions ran off it.
    # The wording in full is in the legend.
    PURE_BIOLOGY = {
        "MCAT": "mito fatty-acid synth.",
        "TIMM8B": "mito import",
        "OXA1L": "mito insertase",
        "RB1CC1": "autophagy",
        "CXXC1": "H3K4me3 targeting",
        "NDUFA8": "complex I",
        "LCK": "TCR kinase",
        "YY1": "architectural TF",
    }
    # Read from the frozen result, never hardcoded. The first version of this panel carried
    # literals (4.95 and 0.009) produced by a throwaway Monte Carlo that was never saved as a
    # script, so the numbers could not be re-derived and the p was an MC estimate printed as if
    # exact. The committed analysis computes the null expectation and tail exactly, by rational
    # factorial moments and binomial-moment inversion, with no draws and no seed.
    PURITY_RESULT = json.loads(
        (PAPER / "data/census/regulator_direction_purity_posthoc.json").read_text())["primary"]
    PURITY_BUNDLE = json.loads(
        (PAPER / "data/census/regulator_direction_purity_posthoc.json").read_text())
    PURE_MIN_REVERSALS = int(PURITY_RESULT["cutoff_min_reversals"])
    PURE_NULL_MEAN = float(PURITY_RESULT["expected_direction_pure"]["decimal"])
    PURE_NULL_P = float(PURITY_RESULT["one_sided_exact_p"]["decimal"])
    PURE_N_ELIGIBLE = int(PURITY_RESULT["n_eligible_regulators"])
    
    
    def _assert_direction_purity_unchanged() -> None:
        """The purity result no longer has a panel, but the Results text still quotes it.
    
        Panel g used to draw these eight regulators. It was retired because what it showed was fixed
        by the rule that selects them, so it could not fail. The numbers survive in the Results
        (8 of 131 regulators with at least six reversals, 4.959 expected, exact one-sided p = 0.0106),
        and this check keeps the figure build failing if the underlying set ever moves.
        """
        edges = pd.read_csv(PAPER / "data/census/sign_inverting_edges.csv")
        edges["lost"] = edges["lfc_Rest"] > 0
        grouped = edges.groupby("regulator")["lost"].agg(["size", "sum"])
        grouped["gained"] = grouped["size"] - grouped["sum"]
        pure = grouped.loc[(grouped["size"] >= PURE_MIN_REVERSALS)
                           & ((grouped["sum"] == 0) | (grouped["gained"] == 0))]
        assert set(pure.index) == set(PURE_BIOLOGY), "the direction-pure set changed"
        assert int((grouped["size"] >= PURE_MIN_REVERSALS).sum()) == PURE_N_ELIGIBLE
        assert len(pure) == int(PURITY_RESULT["n_observed_direction_pure"])
        assert PURITY_BUNDLE["analysis_status"] == "POST_HOC_EXPLORATORY"
        assert np.isclose(PURE_NULL_MEAN, 4.958774250440917)
        assert np.isclose(PURE_NULL_P, 0.010558535828597558)
    
    
    def build(output_stem: str = "Fig2_biology_candidate") -> None:
        global _TRJ_EDGES
        _assert_direction_purity_unchanged()
        _TRJ_EDGES = TRJ.load()
        enriched = FA.fig2a_test.loc[
            FA.fig2a_test.sig_two_sided & FA.fig2a_test.direction.eq("up")
        ].join(FA.direction_counts, on="gene")
        enriched_pure = enriched.loc[
            enriched.lost_on_activation.eq(0) | enriched.gained_on_activation.eq(0)
        ]
        assert set(enriched_pure.gene) == set(TRJ.SIX_LOST + TRJ.SIX_GAINED)
        assert len(enriched_pure) == 6
        assert hashlib.sha256(ROLE_SOURCE.read_bytes()).hexdigest() == ROLE_SOURCE_SHA256
        roles = pd.read_csv(ROLE_SOURCE, sep="\t")
        assert set(roles["regulator"]) == set(TRJ.SIX_LOST + TRJ.SIX_GAINED)
        assert roles["canonical_role"].notna().all()
        S.setup(scale=0.92)
        (result, donors, timecourse, specificity,
         specificity_targets, specificity_hist) = TG.load_and_validate()
        OUT.mkdir(parents=True, exist_ok=True)
    
        # A single 49-row regulator list is legible only as a full A4-page figure. Saving without
        # bbox_inches="tight" preserves the exact A4 MediaBox.
        fig = plt.figure(figsize=(SOURCE_WIDTH_IN, SOURCE_HEIGHT_IN), facecolor="white")
        outer = fig.add_gridspec(
            2,
            1,
            height_ratios=[0.55, 0.45],
            hspace=0.175,
            left=0.080,
            right=0.980,
            top=0.962,
            bottom=0.065,
        )
    
        # Regulator block: b spans the full right column; removing the separate purity-null panel
        # gives a and c the height their one-line labels require. The exact 8/131 purity result is
        # still asserted above and remains available as a standalone diagnostic.
        top = outer[0].subgridspec(
            2, 2, width_ratios=[0.82, 1.58], height_ratios=[0.94, 1.06],
            wspace=0.29, hspace=0.38,
        )
        ax_a = fig.add_subplot(top[0, 0])
        FA.panel_count(ax_a, "a", detail_scale=0.95)
        _remove_panel_tag(ax_a, "a")
        S.panel_tag(ax_a, "a", dx=-0.19, dy=1.06)
        _retune_panel_a(ax_a)
        ax_a.set_xlabel("expected reversing edges (out-degree matched, log scale)")
        ax_a.set_ylabel("observed reversals (log scale)")
    
        bgrid = top[:, 1].subgridspec(1, 2, width_ratios=[1.0, 0.86], wspace=0.04)
        ax_b = fig.add_subplot(bgrid[0, 0])
        ax_b_meta = fig.add_subplot(bgrid[0, 1])
        _panel_all_enriched_bars(ax_b, ax_b_meta, "b")
        b_left = ax_b.get_position()
        b_right = ax_b_meta.get_position()
        fig.text(
            b_left.x0 - 0.019,
            b_left.y1 + 0.027,
            "b",
            ha="right",
            va="top",
            fontsize=11.0,
            fontweight="bold",
            color=S.PAL["INK"],
        )
        fig.text(
            (b_left.x0 + b_right.x1) / 2,
            b_left.y1 + 0.026,
            "All 49 regulators above degree-matched expectation\n"
            "2,075 reversals: 1,154 lost + 921 gained; Fisher/BH q < 0.05",
            ha="center",
            va="top",
            fontsize=7.6,
            fontweight="bold",
            color=S.PAL["INK"],
            linespacing=1.05,
        )
    
        ax_c = fig.add_subplot(top[1, 0])
        FA.panel_rate(
            ax_c,
            "c",
            show_note=True,
            detail_scale=0.92,
            compact=True,
        )
        _retune_panel_c(ax_c)
        ax_c.set_title("Reversal rate is not set by out-degree", fontsize=7.6, pad=4)
        ax_c.set_xlabel("out-degree (both-significant edges, log scale)")
    
        # panel a owns the encodings it uses; the lost/gained key moved into panel b's donut
        ax_a.legend(
            handles=[
                Line2D([0], [0], marker="^", ls="", ms=5.5, color="black", label="above (49)"),
                Line2D([0], [0], marker="v", ls="", ms=5.5, color="black",
                       label="below (57; 42 shown)"),
                Line2D([0], [0], marker="o", ls="", ms=5.0, markerfacecolor="#c7ccd1",
                       markeredgecolor="#9aa0a6", label="not significant"),
                Line2D([0], [0], color="#565b61", lw=1.0, ls="--", label="expectation"),
            ],
            loc="upper left", fontsize=7.0, handletextpad=0.35, labelspacing=0.30,
            borderaxespad=0.4, frameon=False,
        )
    
        # Biology-first target side. The exploratory coefficient-magnitude association is kept as a
        # standalone diagnostic rather than repeated beside the stronger, frozen orientation test.
        # The left column is deliberately wide enough for the time course and the 39-target heatmap;
        # the right column holds their compact interpretation/falsification panels.
        bottom = outer[1].subgridspec(
            2, 2, width_ratios=[1.94, 1.0], height_ratios=[0.86, 1.94],
            wspace=0.43, hspace=0.42,
        )
    
        # d: the time course, with the nuisance-stratified 48 h values overlaid.
        ax_d = fig.add_subplot(bottom[0, 0])
        TG.panel_timecourse(ax_d, timecourse, result)
        _remove_panel_tag(ax_d, "d")
        S.panel_tag(ax_d, "d", dx=-0.04, dy=1.10)
        ax_d.set_title("Relative separation peaks at 4-8 h", fontsize=7.5, pad=4)
        ax_d.set_ylabel("donor AUC (median and IQR)", fontsize=7.0)
        hours = [0, 2, 4, 8, 12, 24, 48, 72]
        donor_n = timecourse.groupby("time")["donor"].nunique().to_dict()
        assert donor_n == {0: 22, 2: 22, 4: 22, 8: 22, 12: 19, 24: 22, 48: 22, 72: 22}
        # the retired AUC scatter was exactly this panel's 48 h column, so only its stratified
        # counterpart carried information the time course did not already show
        strat = donors["stratified_auc"].to_numpy(float)
        raw48 = donors["auc"].to_numpy(float)
        assert np.allclose(np.sort(raw48),
                           np.sort(timecourse.loc[timecourse["time"].eq(48), "auc"].to_numpy(float)))
        x48 = hours.index(48) + 0.28
        ax_d.scatter(np.full(strat.size, x48), strat, s=7, marker="s", color="#6c7a89",
                     alpha=0.75, edgecolors="none", zorder=4)
        ax_d.plot([x48 - 0.12, x48 + 0.12], [np.median(strat)] * 2, color="#3f4a55", lw=1.4, zorder=5)
        # above the stratified points: below them it collided with the "null 0.5" label
        ax_d.text(x48 + 0.10, float(np.max(strat)) + 0.004, "nuisance stratified",
                  fontsize=6.4, color=TG.C_MUTED, ha="left", va="bottom", linespacing=1.2)
        ax_d.set_xticks(range(len(hours)))
        ax_d.set_xticklabels([str(hour) for hour in hours], fontsize=7.0)
        ax_d.set_xlabel("hours after stimulation (n=22; 12 h: n=19)", fontsize=7.0)
        ax_d.set_xlim(-0.25, 8.15)
        for artist in ax_d.texts:
            if artist.get_text().startswith("peak median"):
                # it was drawn just above the curve, which is above the axes top in this slot
                x, y = artist.get_position()
                artist.set_position((x, y - 0.012))
            if artist.get_text().startswith("48 h subset AUC"):
                artist.set_text("48 h median\n0.571")
                artist.set_position((5.55, 0.638))
                artist.set_ha("right")
                artist.set_va("center")
        # e: the absolute comparison, which is the falsification the relative AUC cannot make. It
        # failed as an inset on e at this figure's 7 pt floor; in its own slot it is legible.
        ax_e = fig.add_subplot(bottom[0, 1])
        TG.panel_absolute(ax_e, donors)
        _remove_panel_tag(ax_e, "c")
        S.panel_tag(ax_e, "e", dx=-0.19, dy=1.07)
        ax_e.set_title("Both unanimous groups rise at 48 h", fontsize=7.5, pad=4)
        ax_e.set_ylabel("donor median log2 CPM change", fontsize=7.0)
        ax_e.set_xticklabels(["lost targets\n(n = 351)", "gained targets\n(n = 278)"], fontsize=7.0)
        # headroom for the unanimity caption, which otherwise sits on the lost-group points
        ax_e.set_ylim(ax_e.get_ylim()[0], 1.28)
        for artist in ax_e.texts:
            # the unanimity caption and the lost-group median label were written for a wider box
            if artist.get_text().startswith("both group medians"):
                artist.set_position((0.50, 0.985))
                artist.set_ha("center")
                artist.set_va("top")
            elif artist.get_text() == "0.759":
                x, y = artist.get_position()
                artist.set_position((x, y - 0.235))
                artist.set_va("top")
    
        # f: the biology this figure is for. Redrawn in the Figure 3d/3e idiom, so the tiles are no
        # longer stretched across the page and the block only needs two of the three columns.
        ax_f = fig.add_subplot(bottom[1, 0])
        TRJ.panel(ax_f, _TRJ_EDGES, fs=MIN_SOURCE_FONT_PT)
        TRJ.add_colourbar(ax_f, fs=MIN_SOURCE_FONT_PT)
        for artist in ax_f.texts:
            if artist.get_text() == "net knockdown coefficient (log2 fold change)":
                artist.set_text("KD coefficient (log2 fold change)")
        S.panel_tag(ax_f, "f", dx=-0.03, dy=1.02)
        ax_f.set_title("All 39 targets of the six one-direction members of panel b",
                       fontsize=7.5, pad=6, loc="left")
    
        # g: the null-calibrated all-target orientation analysis. This is not a pooled form of the
        # removed magnitude diagnostic: it tests lost-versus-gained composition, not effect swing.
        ax_g = fig.add_subplot(bottom[1, 1])
        TG.panel_specificity_relation(ax_g, specificity_targets, specificity)
        _remove_panel_tag(ax_g, "b")
        S.panel_tag(ax_g, "g", dx=-0.19, dy=1.07)
        primary = specificity["primary_target_conditioned"]["result"]
        margin = specificity["regulator_margin_sensitivity"]["raw_tmm"]
        n_specificity_targets = int(primary["n_targets"])
        assert n_specificity_targets == len(specificity_targets)
        ax_g.set_title(f"Orientation across {n_specificity_targets} eligible targets",
                       fontsize=7.5, pad=4)
        for artist in ax_g.texts:
            if artist.get_text().startswith(f"n = {n_specificity_targets} targets"):
                artist.set_text(
                    f"rho = {primary['observed_spearman_rho']:.3f}\n"
                    f"target null {primary['null_mean']:.3f}, "
                    f"pMC={primary['plus_one_one_sided_p']:.0e}\n"
                    f"margin null {margin['null_mean']:.3f}, "
                    f"pMC={margin['plus_one_one_sided_p']:.0e}")
                artist.set_fontsize(7.0)
                artist.set_linespacing(1.35)
                artist.set_transform(ax_g.transAxes)
                artist.set_position((0.03, 0.97))
                artist.set_ha("left")
                artist.set_va("top")
                artist.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.2))
        ax_g.set_xlabel("signed reversal score (gained < 0 < lost)", fontsize=7.0)
        ax_g.set_ylabel("independent target response (48 h - 0 h)", fontsize=7.0)
        ax_g.tick_params(labelsize=7.0)
    
        minimum_font, raised = _enforce_minimum_font(fig)
        full_page_scale = FULL_PAGE_WIDTH_MM / (SOURCE_WIDTH_IN * 25.4)
        two_column_scale = TWO_COLUMN_WIDTH_MM / (SOURCE_WIDTH_IN * 25.4)
        print(
            f"source page: {SOURCE_WIDTH_IN:.2f} x {SOURCE_HEIGHT_IN:.2f} in; "
            f"{FULL_PAGE_WIDTH_MM:.1f}-mm scale: {full_page_scale:.4f}; "
            f"{TWO_COLUMN_WIDTH_MM:.1f}-mm scale: {two_column_scale:.4f}"
        )
        print(f"minimum visible source font: {minimum_font:.2f} pt")
        print(
            f"minimum font at {FULL_PAGE_WIDTH_MM:.1f} mm: "
            f"{minimum_font * full_page_scale:.2f} pt"
        )
        print(
            f"minimum font if reduced to {TWO_COLUMN_WIDTH_MM:.1f} mm: "
            f"{minimum_font * two_column_scale:.2f} pt"
        )
        print(f"text objects raised to {MIN_SOURCE_FONT_PT:.1f} pt: {len(raised)}")
    
        # Reflowing three panels into the space that held two makes label collisions the real risk,
        # and a downscaled PNG hides them. Measure on the canvas instead of looking.
        from figure_overlap_check import find_overlaps
        problems = find_overlaps(fig)
        for line in problems:
            print(f"OVERLAP {line}")
        print(f"measured text collisions: {len(problems)}")
    
        out_pdf = OUT / f"{output_stem}.pdf"
        out_png = OUT / f"{output_stem}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
    
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_biology_candidate, k, v)
init_fig2_biology_candidate()

def init_fig2h_activation_movement_similarity():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Corrected Figure 2h analysis: activation-movement correlations are not reversal-specific.
    
    The external response is the frozen, donor-aligned 22-donor GSE140244 TMM log2-CPM
    48-minus-0 h median.  The primary estimand is the raw within-regulator Spearman correlation
    between that external response and the absolute Rest-to-48 h coefficient movement, computed
    separately for reversing and stable edges and compared in paired regulators.
    
    The former se_Stim48-only residualisation is deliberately not used: it manufactured the
    positive separation that this corrected analysis falsifies.  A partial-Spearman sensitivity
    adjusting for both endpoint standard errors is retained in the JSON and caption only.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    
    
    PAPER = ROOT
    MASTER = PAPER / "data/census/census_master_edges.csv.gz"
    DELTA = (PAPER / "data/census/twopath_motif_exploratory/"
             "gse140244_all_measured_gene_delta.tsv.gz")
    DELTA_MANIFEST = (PAPER / "data/census/twopath_motif_exploratory/"
                      "gse140244_all_measured_gene_delta.manifest.json")
    OUT_TSV = PAPER / "data/census/fig2h_activation_movement_similarity.tsv"
    OUT_JSON = PAPER / "data/census/fig2h_activation_movement_similarity.json"
    OUTDIR = PAPER / "data/figure_provenance/fig2h_activation_movement_similarity"
    
    MIN_TARGETS_PRIMARY = 4
    MIN_TARGETS_SENSITIVITY = 8
    RAW_FLOOR_SWEEP = (4, 5, 6, 8, 10, 12)
    REV = "#c1272d"
    STA = "#7f8a94"
    INK = "#22262b"
    MUTED = "#687078"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _load_edges() -> pd.DataFrame:
        manifest = json.loads(DELTA_MANIFEST.read_text(encoding="utf-8"))
        assert manifest["paired_donors"] == 22
        assert manifest["processing"] == (
            "edgeR TMM log2 CPM prior.count=0.5; median paired 48-minus-0 delta"
        )
        assert manifest["sha256"] == sha256(DELTA)
    
        delta = pd.read_csv(DELTA, sep="\t")
        assert delta["target_ensg"].is_unique
        assert delta["n_paired_donors"].eq(22).all()
        response = delta.set_index("target_ensg")["median_delta_0_to_48"]
    
        edges = pd.read_csv(MASTER)
        edges = edges.loc[~edges["is_self"].astype(bool)].copy()
        assert len(edges) == 83_489
        edges["external_delta_48_minus_0"] = edges["target_ensg"].map(response)
        edges["abs_coefficient_movement"] = (
            edges["lfc_Stim48"] - edges["lfc_Rest"]
        ).abs()
        keep = [
            "regulator", "regulator_ensg", "target", "target_ensg", "kind",
            "external_delta_48_minus_0", "abs_coefficient_movement",
            "se_Rest", "se_Stim48",
        ]
        # The raw primary estimand does not condition on either endpoint SE.  Keep
        # those columns nullable here and apply their complete-case requirement only
        # inside the explicitly labelled both-endpoint-SE sensitivity analysis.
        return edges[keep].dropna(
            subset=["external_delta_48_minus_0", "abs_coefficient_movement"]
        ).reset_index(drop=True)
    
    
    def _raw_rho(edges: pd.DataFrame, kind: str, floor: int) -> pd.DataFrame:
        records: list[dict[str, object]] = []
        for (regulator, regulator_ensg), group in edges.loc[
                edges["kind"].eq(kind)].groupby(["regulator", "regulator_ensg"], sort=True):
            if len(group) < floor:
                continue
            rho = spearmanr(
                group["external_delta_48_minus_0"],
                group["abs_coefficient_movement"],
            ).statistic
            if np.isfinite(rho):
                records.append({
                    "regulator": regulator,
                    "regulator_ensg": regulator_ensg,
                    f"n_{kind}_edges": int(len(group)),
                    f"rho_{kind}_raw": float(rho),
                })
        return pd.DataFrame.from_records(records)
    
    
    def _rank_residual(values: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
        response = rankdata(values)
        design = np.column_stack(
            [np.ones(len(response))] + [rankdata(covariate) for covariate in covariates]
        )
        return response - design @ np.linalg.lstsq(design, response, rcond=None)[0]
    
    
    def _both_endpoint_se_rho(edges: pd.DataFrame, kind: str, floor: int) -> pd.DataFrame:
        records: list[dict[str, object]] = []
        for (regulator, regulator_ensg), group in edges.loc[
                edges["kind"].eq(kind)].groupby(["regulator", "regulator_ensg"], sort=True):
            group = group.dropna(
                subset=[
                    "external_delta_48_minus_0",
                    "abs_coefficient_movement",
                    "se_Rest",
                    "se_Stim48",
                ]
            )
            if len(group) < floor:
                continue
            covariates = [group["se_Rest"].to_numpy(), group["se_Stim48"].to_numpy()]
            x = _rank_residual(group["external_delta_48_minus_0"].to_numpy(), covariates)
            y = _rank_residual(group["abs_coefficient_movement"].to_numpy(), covariates)
            rho = pearsonr(x, y).statistic
            if np.isfinite(rho):
                records.append({
                    "regulator": regulator,
                    "regulator_ensg": regulator_ensg,
                    f"rho_{kind}_both_endpoint_se": float(rho),
                })
        return pd.DataFrame.from_records(records)
    
    
    def _paired_summary(table: pd.DataFrame, rev_col: str, stable_col: str) -> dict:
        paired = table[[rev_col, stable_col]].dropna()
        test = wilcoxon(paired[rev_col], paired[stable_col], alternative="two-sided")
        return {
            "n_paired_regulators": int(len(paired)),
            "median_reversing": float(paired[rev_col].median()),
            "median_stable": float(paired[stable_col].median()),
            "median_paired_difference": float((paired[rev_col] - paired[stable_col]).median()),
            "n_reversing_greater": int((paired[rev_col] > paired[stable_col]).sum()),
            "wilcoxon_two_sided_p": float(test.pvalue),
        }
    
    
    def compute() -> tuple[pd.DataFrame, dict]:
        edges = _load_edges()
        edge_counts = edges["kind"].value_counts().to_dict()
        assert edge_counts["rev"] == 4_344
        assert edge_counts["stable"] == 78_520
        reversing = _raw_rho(edges, "rev", MIN_TARGETS_PRIMARY)
        stable = _raw_rho(edges, "stable", MIN_TARGETS_PRIMARY)
        assert len(reversing) == 155
        source = reversing.merge(stable, on=["regulator", "regulator_ensg"], how="inner")
        source["rho_raw_difference_rev_minus_stable"] = (
            source["rho_rev_raw"] - source["rho_stable_raw"]
        )
        source = source.sort_values(["regulator_ensg", "regulator"], kind="stable").reset_index(drop=True)
    
        sensitivity_rev = _both_endpoint_se_rho(edges, "rev", MIN_TARGETS_SENSITIVITY)
        sensitivity_stable = _both_endpoint_se_rho(edges, "stable", MIN_TARGETS_SENSITIVITY)
        sensitivity = sensitivity_rev.merge(
            sensitivity_stable, on=["regulator", "regulator_ensg"], how="inner"
        )
        sensitivity["rho_both_endpoint_se_difference_rev_minus_stable"] = (
            sensitivity["rho_rev_both_endpoint_se"] - sensitivity["rho_stable_both_endpoint_se"]
        )
    
        primary = _paired_summary(source, "rho_rev_raw", "rho_stable_raw")
        both_se = _paired_summary(
            sensitivity, "rho_rev_both_endpoint_se", "rho_stable_both_endpoint_se"
        )
        floor_sweep: list[dict[str, object]] = []
        for floor in RAW_FLOOR_SWEEP:
            floor_rev = _raw_rho(edges, "rev", floor)
            floor_stable = _raw_rho(edges, "stable", floor)
            floor_table = floor_rev.merge(
                floor_stable, on=["regulator", "regulator_ensg"], how="inner"
            )
            floor_sweep.append({
                "minimum_edges_per_kind_per_regulator": floor,
                **_paired_summary(floor_table, "rho_rev_raw", "rho_stable_raw"),
            })
        # These are broad drift checks, not presentation rounding.
        assert primary["n_paired_regulators"] == 154
        assert np.isclose(primary["median_reversing"], 0.077, atol=0.002)
        assert np.isclose(primary["median_stable"], 0.055, atol=0.002)
        assert np.isclose(primary["wilcoxon_two_sided_p"], 0.572, atol=0.005)
    
        result = {
            "_what": (
                "paired within-regulator comparison of activation-movement correlations for "
                "reversing versus stable target edges"
            ),
            "_script": "scripts/fig2h_activation_movement_similarity.py",
            "analysis_status": "CORRECTED_NEGATIVE_RESULT",
            "design_provenance": "post-hoc paired-regulator comparison",
            "claim_boundary": (
                "No detectable difference was found between the reversing-edge and stable-edge "
                "correlations with the paired-22-donor TMM activation response. This is not an "
                "equivalence test, and the measured association does not establish mechanism."
            ),
            "invalidated_predecessor": {
                "script": "scripts/fig2g_census_swing_vs_induction.py",
                "reason": (
                    "se_Stim48-only residualisation induced the reported separation; it is not the "
                    "primary analysis and is not used in Figure 2"
                ),
            },
            "inputs": {
                "census_master_edges.csv.gz": sha256(MASTER),
                "gse140244_all_measured_gene_delta.tsv.gz": sha256(DELTA),
                "gse140244_all_measured_gene_delta.manifest.json": sha256(DELTA_MANIFEST),
            },
            "external_response": {
                "paired_donors": 22,
                "normalization": "edgeR TMM log2 CPM with prior.count=0.5",
                "summary": "per-gene median paired 48-minus-0 h response",
            },
            "mapped_edge_coverage": {
                "reversing_edges": int(edge_counts["rev"]),
                "stable_edges": int(edge_counts["stable"]),
                "reversing_regulators_with_raw_rho_at_floor_4": int(len(reversing)),
            },
            "primary": {
                "minimum_edges_per_kind_per_regulator": MIN_TARGETS_PRIMARY,
                "estimand": (
                    "raw within-regulator Spearman rho: external activation response versus absolute "
                    "Rest-to-48 h coefficient movement"
                ),
                **primary,
            },
            "raw_minimum_edge_floor_sweep": floor_sweep,
            "both_endpoint_se_sensitivity": {
                "display_status": "persisted sensitivity; not printed in panel h",
                "minimum_edges_per_kind_per_regulator": MIN_TARGETS_SENSITIVITY,
                "estimand": (
                    "partial Spearman via rank residuals adjusted jointly for se_Rest and se_Stim48"
                ),
                **both_se,
            },
            "caption_guidance": (
                "Raw within-regulator correlations were compared at a minimum of four edges of "
                "each kind; two-sided paired Wilcoxon p=0.572. Raw minimum-edge floors from 4 to "
                "12 gave no detectable separation. A secondary partial-Spearman sensitivity "
                "rank-residualised both Rest and 48 h standard errors and is retained in the source "
                "JSON, not panel h. Failure to detect a difference is not an equivalence test."
            ),
        }
        return source, result
    
    
    def persist(source: pd.DataFrame, result: dict) -> None:
        OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
        source.to_csv(OUT_TSV, sep="\t", index=False)
        result = dict(result)
        result["source_table"] = {
            "path": str(OUT_TSV.relative_to(PAPER)),
            "rows": int(len(source)),
            "sha256": sha256(OUT_TSV),
        }
        OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    
    
    def load_and_validate() -> tuple[pd.DataFrame, dict]:
        source = pd.read_csv(OUT_TSV, sep="\t")
        result = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        assert result["analysis_status"] == "CORRECTED_NEGATIVE_RESULT"
        assert result["design_provenance"] == "post-hoc paired-regulator comparison"
        assert result["source_table"]["sha256"] == sha256(OUT_TSV)
        assert result["source_table"]["rows"] == len(source) == 154
        assert result["inputs"]["census_master_edges.csv.gz"] == sha256(MASTER)
        assert result["inputs"]["gse140244_all_measured_gene_delta.tsv.gz"] == sha256(DELTA)
        return source, result
    
    
    def panel(ax, source: pd.DataFrame, result: dict, fs: float = 7.0, compact: bool = False) -> None:
        """Paired raw correlations; no highlighted post-hoc regulator subset."""
        x = np.array([0.0, 1.0])
        for row in source.itertuples(index=False):
            ax.plot(x, [row.rho_stable_raw, row.rho_rev_raw], color="#b8bec3",
                    lw=0.42, alpha=0.20, zorder=1)
        ax.scatter(np.zeros(len(source)), source["rho_stable_raw"], s=7 if compact else 10,
                   color=STA, alpha=0.34, edgecolors="none", zorder=2)
        ax.scatter(np.ones(len(source)), source["rho_rev_raw"], s=7 if compact else 10,
                   color=REV, alpha=0.34, edgecolors="none", zorder=2)
        stable_median = float(result["primary"]["median_stable"])
        reversal_median = float(result["primary"]["median_reversing"])
        ax.plot(x, [stable_median, reversal_median], color=INK, lw=1.2, zorder=4)
        ax.scatter(x, [stable_median, reversal_median], marker="D", s=32 if compact else 42,
                   color=[STA, REV], edgecolors="white", linewidths=0.7, zorder=5)
        ax.axhline(0, color="#7b838a", lw=0.75, ls=(0, (4, 2.5)), zorder=0)
        ax.set_xlim(-0.32, 1.32)
        ax.set_ylim(-1.02, 1.02)
        ax.set_xticks(x, ["stable", "reversing"])
        ax.set_yticks([-1, 0, 1] if compact else [-1, -0.5, 0, 0.5, 1])
        ax.tick_params(labelsize=fs)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        p = result["primary"]
        ax.text(
            0.5,
            0.97,
            "median rho\n"
            f"reversing {p['median_reversing']:+.3f}; stable {p['median_stable']:+.3f}\n"
            f"n={p['n_paired_regulators']}; two-sided p={p['wilcoxon_two_sided_p']:.3f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=fs,
            color=INK,
            linespacing=1.20,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.0),
            zorder=6,
        )
    
    
    def figure(source: pd.DataFrame, result: dict) -> None:
        OUTDIR.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4.4, 3.4), facecolor="white")
        panel(ax, source, result, fs=7.0, compact=False)
        ax.set_title(
            "No detected excess in activation–movement\n"
            "correlation for reversing edges",
            fontsize=8.2,
            fontweight="bold",
            pad=5,
            linespacing=1.0,
        )
        ax.set_ylabel("within-regulator Spearman rho", fontsize=7.0)
        fig.subplots_adjust(left=0.16, right=0.98, top=0.84, bottom=0.14)
        fig.savefig(OUTDIR / "Fig2h_activation_movement_similarity.pdf", dpi=300, facecolor="white")
        fig.savefig(OUTDIR / "Fig2h_activation_movement_similarity.png", dpi=300, facecolor="white")
        plt.close(fig)
    
    
    def main() -> None:
        source, result = compute()
        persist(source, result)
        source, result = load_and_validate()
        figure(source, result)
        p = result["primary"]
        s = result["both_endpoint_se_sensitivity"]
        print(
            f"primary: n={p['n_paired_regulators']}, median reversing "
            f"{p['median_reversing']:+.4f} vs stable {p['median_stable']:+.4f}, "
            f"two-sided Wilcoxon p={p['wilcoxon_two_sided_p']:.4g}"
        )
        print(
            f"both-SE sensitivity: n={s['n_paired_regulators']}, median reversing "
            f"{s['median_reversing']:+.4f} vs stable {s['median_stable']:+.4f}, "
            f"two-sided Wilcoxon p={s['wilcoxon_two_sided_p']:.4g}"
        )
        print(f"wrote {OUT_TSV}")
        print(f"wrote {OUT_JSON}")
        print(f"wrote {OUTDIR}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2h_activation_movement_similarity, k, v)
init_fig2h_activation_movement_similarity()

def init_fig2_eight_panel_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render a noncanonical eight-panel Figure 2 candidate.
    
    This keeps the promoted seven-panel Figure 2 unchanged and adds the corrected paired-donor
    activation-movement comparison as panel h.  Every panel is redrawn from its source data;
    no diagnostic PDF is embedded or cropped.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    BASE = fig2_biology_candidate
    MOVEMENT = fig2h_activation_movement_similarity
    
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_eight_panel_candidate"
    STEM = "Fig2_eight_panel_candidate"
    CANONICAL_PDF = PAPER / "figures/Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures/Fig2.png"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _validate_movement_inputs(source: pd.DataFrame, result: dict) -> None:
        """Fail closed if the corrected source table and report no longer agree."""
        primary = result["primary"]
        assert result["analysis_status"] == "CORRECTED_NEGATIVE_RESULT"
        assert len(source) == int(primary["n_paired_regulators"]) == 154
        assert np.isclose(source["rho_rev_raw"].median(),
                          float(primary["median_reversing"]), atol=5e-10)
        assert np.isclose(source["rho_stable_raw"].median(),
                          float(primary["median_stable"]), atol=5e-10)
        assert np.isclose(float(primary["wilcoxon_two_sided_p"]), 0.5717292019676893,
                          atol=5e-12)
    
    
    def _build() -> tuple[plt.Figure, dict]:
        BASE._assert_direction_purity_unchanged()
        trajectory_edges = BASE.TRJ.load()
        enriched = BASE.FA.fig2a_test.loc[
            BASE.FA.fig2a_test.sig_two_sided & BASE.FA.fig2a_test.direction.eq("up")
        ].join(BASE.FA.direction_counts, on="gene")
        enriched_pure = enriched.loc[
            enriched.lost_on_activation.eq(0) | enriched.gained_on_activation.eq(0)
        ]
        assert set(enriched_pure.gene) == set(BASE.TRJ.SIX_LOST + BASE.TRJ.SIX_GAINED)
        assert len(enriched_pure) == 6

    
        BASE.S.setup(scale=0.92)
        (result, donors, timecourse, specificity,
         specificity_targets, _) = BASE.TG.load_and_validate()
        movement_source, movement_result = MOVEMENT.load_and_validate()
        _validate_movement_inputs(movement_source, movement_result)
    
        fig = plt.figure(
            figsize=(BASE.SOURCE_WIDTH_IN, BASE.SOURCE_HEIGHT_IN),
            facecolor="white",
        )
        outer = fig.add_gridspec(
            2,
            1,
            height_ratios=[0.55, 0.45],
            hspace=0.175,
            left=0.080,
            right=0.980,
            top=0.962,
            bottom=0.065,
        )
    
        # Panels a-c are identical to the promoted seven-panel figure.
        top = outer[0].subgridspec(
            2,
            2,
            width_ratios=[0.82, 1.58],
            height_ratios=[0.94, 1.06],
            wspace=0.29,
            hspace=0.38,
        )
        ax_a = fig.add_subplot(top[0, 0])
        BASE.FA.panel_count(ax_a, "a", detail_scale=0.95)
        BASE._remove_panel_tag(ax_a, "a")
        BASE.S.panel_tag(ax_a, "a", dx=-0.19, dy=1.06)
        BASE._retune_panel_a(ax_a)
        ax_a.set_xlabel("expected reversing edges (out-degree matched, log scale)")
        ax_a.set_ylabel("observed reversals (log scale)")
        ax_a.legend(
            handles=[
                Line2D([0], [0], marker="^", ls="", ms=5.5, color="black", label="above (49)"),
                Line2D([0], [0], marker="v", ls="", ms=5.5, color="black",
                       label="below (57; 42 shown)"),
                Line2D([0], [0], marker="o", ls="", ms=5.0, markerfacecolor="#c7ccd1",
                       markeredgecolor="#9aa0a6", label="not significant"),
                Line2D([0], [0], color="#565b61", lw=1.0, ls="--", label="expectation"),
            ],
            loc="upper left",
            fontsize=7.0,
            handletextpad=0.35,
            labelspacing=0.30,
            borderaxespad=0.4,
            frameon=False,
        )
    
        bgrid = top[:, 1].subgridspec(1, 2, width_ratios=[1.0, 0.86], wspace=0.04)
        ax_b = fig.add_subplot(bgrid[0, 0])
        ax_b_meta = fig.add_subplot(bgrid[0, 1])
        BASE._panel_all_enriched_bars(ax_b, ax_b_meta, "b")
        b_left = ax_b.get_position()
        b_right = ax_b_meta.get_position()
        fig.text(
            b_left.x0 - 0.019,
            b_left.y1 + 0.027,
            "b",
            ha="right",
            va="top",
            fontsize=11.0,
            fontweight="bold",
            color=BASE.S.PAL["INK"],
        )
        fig.text(
            (b_left.x0 + b_right.x1) / 2,
            b_left.y1 + 0.026,
            "All 49 regulators above degree-matched expectation\n"
            "2,075 reversals: 1,154 lost + 921 gained; Fisher/BH q < 0.05",
            ha="center",
            va="top",
            fontsize=7.6,
            fontweight="bold",
            color=BASE.S.PAL["INK"],
            linespacing=1.05,
        )
    
        ax_c = fig.add_subplot(top[1, 0])
        BASE.FA.panel_rate(ax_c, "c", show_note=True, detail_scale=0.92, compact=True)
        BASE._retune_panel_c(ax_c)
        ax_c.set_title("Reversal rate is not set by out-degree", fontsize=7.6, pad=4)
        ax_c.set_xlabel("out-degree (both-significant edges, log scale)")
    
        # Preserve the promoted heatmap/orientation geometry in the tall row; split the short row
        # among the two expression summaries and the corrected paired negative result.
        lower = outer[1].subgridspec(2, 1, height_ratios=[1.94, 0.86], hspace=0.42)
        lower_top = lower[0].subgridspec(1, 2, width_ratios=[1.94, 1.0], wspace=0.43)
        lower_bottom = lower[1].subgridspec(
            1, 3, width_ratios=[1.60, 0.95, 0.95], wspace=0.28
        )
    
        ax_d = fig.add_subplot(lower_top[0, 0])
        BASE.TRJ.panel(ax_d, trajectory_edges, fs=BASE.MIN_SOURCE_FONT_PT)
        BASE.TRJ.add_colourbar(ax_d, fs=BASE.MIN_SOURCE_FONT_PT)
        for artist in ax_d.texts:
            if artist.get_text() == "net knockdown coefficient (log2 fold change)":
                artist.set_text("KD effect (log2 FC)")
            elif artist.get_text() == "columns R, 8 and 48: Rest, 8 h, 48 h":
                artist.set_text("R/8/48 = Rest/8 h/48 h")
            elif artist.get_text() == "CARD11\nNF-kB scaffold":
                artist.set_text("CARD11\nNF-κB\nscaffold")
                artist.set_linespacing(0.95)
            elif artist.get_text() == "SUGT1\nHSP90 co-chaperone":
                x, y = artist.get_position()
                artist.set_text("SUGT1\nHSP90\nco-chaperone")
                artist.set_position((x - 0.90, y))
                artist.set_linespacing(0.95)
        BASE.S.panel_tag(ax_d, "d", dx=-0.04, dy=1.02)
        ax_d.set_title("Six one-direction panel-b regulators: all 39 targets",
                       fontsize=7.5, pad=6, loc="left")
    
        ax_e = fig.add_subplot(lower_top[0, 1])
        BASE.TG.panel_specificity_relation(ax_e, specificity_targets, specificity)
        BASE._remove_panel_tag(ax_e, "b")
        BASE.S.panel_tag(ax_e, "e", dx=-0.25, dy=1.07)
        primary = specificity["primary_target_conditioned"]["result"]
        margin = specificity["regulator_margin_sensitivity"]["raw_tmm"]
        n_specificity_targets = int(primary["n_targets"])
        assert n_specificity_targets == len(specificity_targets) == 903
        ax_e.set_title(f"Target orientation (n = {n_specificity_targets})", fontsize=7.3, pad=4)
        for artist in ax_e.texts:
            if artist.get_text().startswith(f"n = {n_specificity_targets} targets"):
                artist.set_text(
                    f"rho = {primary['observed_spearman_rho']:.3f}\n"
                    f"target null {primary['null_mean']:.3f}; pMC={primary['plus_one_one_sided_p']:.0e}\n"
                    f"margin null {margin['null_mean']:.3f}; pMC={margin['plus_one_one_sided_p']:.0e}"
                )
                artist.set_fontsize(7.0)
                artist.set_linespacing(1.30)
                artist.set_transform(ax_e.transAxes)
                artist.set_position((0.03, 0.97))
                artist.set_ha("left")
                artist.set_va("top")
                artist.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.0))
        ax_e.set_xlabel("signed reversal score", fontsize=7.0)
        ax_e.set_ylabel("independent response (48 h - 0 h)", fontsize=7.0)
        ax_e.tick_params(labelsize=7.0)
    
        ax_f = fig.add_subplot(lower_bottom[0, 0])
        BASE.TG.panel_timecourse(ax_f, timecourse, result)
        BASE._remove_panel_tag(ax_f, "d")
        BASE.S.panel_tag(ax_f, "f", dx=-0.04, dy=1.10)
        ax_f.set_title("Relative separation peaks at 4-8 h", fontsize=7.5, pad=4)
        ax_f.set_ylabel("donor AUC (median and IQR)", fontsize=7.0)
        hours = [0, 2, 4, 8, 12, 24, 48, 72]
        donor_n = timecourse.groupby("time")["donor"].nunique().to_dict()
        assert donor_n == {0: 22, 2: 22, 4: 22, 8: 22, 12: 19, 24: 22, 48: 22, 72: 22}
        strat = donors["stratified_auc"].to_numpy(float)
        raw48 = donors["auc"].to_numpy(float)
        assert np.allclose(
            np.sort(raw48),
            np.sort(timecourse.loc[timecourse["time"].eq(48), "auc"].to_numpy(float)),
        )
        x48 = hours.index(48) + 0.28
        ax_f.scatter(np.full(strat.size, x48), strat, s=7, marker="s", color="#6c7a89",
                     alpha=0.75, edgecolors="none", zorder=4)
        ax_f.plot([x48 - 0.12, x48 + 0.12], [np.median(strat)] * 2,
                  color="#3f4a55", lw=1.4, zorder=5)
        ax_f.text(x48 + 0.10, float(np.max(strat)) + 0.004, "stratified",
                  fontsize=7.0, color=BASE.TG.C_MUTED, ha="left", va="bottom")
        ax_f.set_xticks(range(len(hours)))
        ax_f.set_xticklabels([str(hour) for hour in hours], fontsize=7.0)
        ax_f.set_xlabel("hours after stimulation (n=22; 12 h: n=19)", fontsize=7.0)
        ax_f.set_xlim(-0.25, 8.15)
        for artist in ax_f.texts:
            if artist.get_text().startswith("peak median"):
                x, y = artist.get_position()
                artist.set_position((x, y - 0.012))
            elif artist.get_text().startswith("48 h subset AUC"):
                artist.set_text("48 h median\n0.571")
                artist.set_position((5.55, 0.638))
                artist.set_ha("right")
                artist.set_va("center")
    
        ax_g = fig.add_subplot(lower_bottom[0, 1])
        BASE.TG.panel_absolute(ax_g, donors)
        BASE._remove_panel_tag(ax_g, "c")
        BASE.S.panel_tag(ax_g, "g", dx=-0.04, dy=1.10)
        ax_g.set_title("Both groups rise at 48 h", fontsize=7.5, pad=4)
        ax_g.set_ylabel("donor median log2 CPM change", fontsize=7.0)
        ax_g.set_xticklabels(["lost targets\n(n = 351)", "gained targets\n(n = 278)"],
                             fontsize=7.0)
        ax_g.set_ylim(ax_g.get_ylim()[0], 1.28)
        for artist in ax_g.texts:
            if artist.get_text().startswith("both group medians"):
                artist.set_position((0.50, 0.985))
                artist.set_ha("center")
                artist.set_va("top")
            elif artist.get_text() == "0.759":
                x, y = artist.get_position()
                artist.set_position((x, y - 0.235))
                artist.set_va("top")
    
        ax_h = fig.add_subplot(lower_bottom[0, 2])
        MOVEMENT.panel(ax_h, movement_source, movement_result, fs=7.0, compact=True)
        BASE.S.panel_tag(ax_h, "h", dx=-0.25, dy=1.07)
        ax_h.set_title(
            "No detected excess in\n"
            "activation–movement correlation",
            fontsize=7.3,
            pad=3,
            linespacing=1.0,
        )
        ax_h.set_ylabel("within-regulator rho", fontsize=7.0)
        ax_h.tick_params(labelsize=7.0)
    
        minimum_font, raised = BASE._enforce_minimum_font(fig)
        assert minimum_font >= BASE.MIN_SOURCE_FONT_PT
        return fig, {
            "minimum_source_font_pt": float(minimum_font),
            "text_objects_raised_to_floor": len(raised),
            "activation_movement_result": {
                "analysis_status": movement_result["analysis_status"],
                "paired_regulators": int(movement_result["primary"]["n_paired_regulators"]),
                "median_reversing": float(movement_result["primary"]["median_reversing"]),
                "median_stable": float(movement_result["primary"]["median_stable"]),
                "wilcoxon_two_sided_p": float(
                    movement_result["primary"]["wilcoxon_two_sided_p"]
                ),
            },
        }
    
    
    def build() -> None:
        canonical_before = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        fig, audit = _build()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out_pdf = OUTDIR / f"{STEM}.pdf"
        out_png = OUTDIR / f"{STEM}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
    
        canonical_after = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        assert canonical_before == canonical_after, "canonical Figure 2 changed during candidate build"
        manifest = {
            "schema_version": "fig2-eight-panel-candidate-v2-corrected-negative",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "script": "scripts/fig2_eight_panel_candidate.py",
            "page": {"width_mm": 210.0, "height_mm": 297.0, "format": "A4 portrait"},
            "panel_order": ["a", "b", "c", "d", "e", "f", "g", "h"],
            "renderer_closure": {
                path: sha256(PAPER / path)
                for path in (
                    "scripts/fig2_eight_panel_candidate.py",
                    "scripts/fig2h_activation_movement_similarity.py",
                    "scripts/fig2_biology_candidate.py",
                    "scripts/figstyle.py",
                    "scripts/fig4_architecture.py",
                    "scripts/fig_target_activation_trajectory_v1.py",
                    "scripts/fig2g_trajectory_panel.py",
                )
            },
            "added_panel": {
                "panel": "h",
                "source_script": "scripts/fig2h_activation_movement_similarity.py",
                "source_table": "data/census/fig2h_activation_movement_similarity.tsv",
                "source_result": "data/census/fig2h_activation_movement_similarity.json",
                "source_diagnostic_pdf": (
                    "data/figure_provenance/fig2h_activation_movement_similarity/"
                    "Fig2h_activation_movement_similarity.pdf"
                ),
                "source_hashes": {
                    "table_sha256": sha256(MOVEMENT.OUT_TSV),
                    "result_json_sha256": sha256(MOVEMENT.OUT_JSON),
                },
                "rendering": "direct panel-function and data reuse; no PDF embedding",
                "predecessor_excluded": "scripts/fig2g_census_swing_vs_induction.py",
                "claim": "no detected difference in activation-movement correlation",
            },
            "canonical_figure_2_unchanged": canonical_after,
            "audit": audit,
            "outputs": {
                "pdf": {"path": str(out_pdf.relative_to(PAPER)), "sha256": sha256(out_pdf)},
                "png": {"path": str(out_png.relative_to(PAPER)), "sha256": sha256(out_png)},
            },
        }
        manifest_path = OUTDIR / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"minimum visible source font: {audit['minimum_source_font_pt']:.2f} pt")
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
        print(f"wrote {manifest_path}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_eight_panel_candidate, k, v)
init_fig2_eight_panel_candidate()

def init_fig2h_8hr_sign_states():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the audited 8 h sign-state result as a compact Figure 2 panel.
    
    The reversal set is defined only by Rest and Stim48.  This renderer reads the
    already-frozen post-hoc 8 h sign-state artifacts; it performs no edge selection,
    testing, or resampling.  The unresolved majority is displayed explicitly.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    SIGN_STATES = PAPER / "data/kinetics/kinetics_8hr_sign_states.csv.gz"
    SUMMARY = PAPER / "data/kinetics/kinetics_8hr_sign_states.summary.json"
    PAIRED = PAPER / "data/kinetics/kinetics_8hr_sign_states.paired_regulators.tsv"
    SOURCE_MANIFEST = PAPER / "data/kinetics/kinetics_8hr_sign_states.manifest.json"
    OUTDIR = PAPER / "data/figure_provenance/fig2h_8hr_sign_states"
    
    EXPECTED_HASHES = {
        SIGN_STATES: "9d34bdb23acdc1a1556854bea0e148a8f1945fb5dfe14c41ad5162a4cd11c222",
        SUMMARY: "af034de06b5b0de8a414e8f8bd3c692377de5faa8a3d7899b6db83d0f2a70a9a",
        PAIRED: "a8bc41f6745335396fdd51726011533b558f3a7ced452b67fd878174c97ec3ba",
        SOURCE_MANIFEST: "2b26eb8dd424eff79a2912cd99c045ed562628f3d12848e29360c178e46ac212",
    }
    
    LOST = "lost_on_activation"
    GAINED = "gained_on_activation"
    REST = "Rest-sign-retaining"
    STIM48 = "Stim48-sign-matching"
    UNRESOLVED = "Unresolved"
    
    INK = "#24292e"
    MUTED = "#66717a"
    RESOLVED = "#59636d"
    UNRESOLVED_COLOUR = "#d6dadd"
    # Paper-wide state and reversal-direction conventions.
    REST_COLOUR = "#2c6fbb"
    STIM48_COLOUR = "#e07b39"
    LOST_COLOUR = "#c1272d"
    GAINED_COLOUR = "#2a9d8f"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def load_and_validate() -> tuple[pd.DataFrame, dict]:
        for path, expected in EXPECTED_HASHES.items():
            observed = sha256(path)
            if observed != expected:
                raise AssertionError(f"source hash changed for {path}: {observed}")
    
        states = pd.read_csv(SIGN_STATES)
        summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
        manifest = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
        assert summary["analysis_status"] == "POST_HOC_EXPLORATORY"
        assert manifest["analysis_status"] == "POST_HOC_EXPLORATORY"
        assert len(states) == 4_377
        assert states["eight_hour_resolved"].sum() == 1_279
    
        observed = pd.crosstab(states["direction"], states["eight_hour_sign_state"])
        expected = {
            GAINED: {REST: 295, STIM48: 280, UNRESOLVED: 1_506},
            LOST: {REST: 529, STIM48: 175, UNRESOLVED: 1_592},
        }
        for direction, counts in expected.items():
            for state, count in counts.items():
                assert int(observed.loc[direction, state]) == count
    
        coverage = summary["coverage"]
        assert coverage == {
            "n_Rest_sign_retaining": 824,
            "n_Stim48_sign_matching": 455,
            "n_endpoint_defined_reversals_in_8h_input": 4_377,
            "n_resolved": 1_279,
            "n_unresolved": 3_098,
        }
        paired = summary["paired_regulator_result"]
        assert paired["n_regulators_with_both_directions"] == 75
        assert np.isclose(paired["p_value"], 0.0049970904, atol=5e-11)
        return states, summary
    
    
    def _label_segment(ax, left: float, width: float, y: float, text: str, colour: str) -> None:
        ax.text(
            left + width / 2,
            y,
            text,
            ha="center",
            va="center",
            fontsize=7.0,
            fontweight="bold",
            color=colour,
            linespacing=0.85,
            clip_on=False,
        )
    
    
    def panel(ax, states: pd.DataFrame, summary: dict, fs: float = 7.0) -> None:
        """Coverage plus resolved-edge sign states; no switch-time claim."""
        del states  # validated upstream; displayed values are read from the sealed summary.
        coverage = summary["coverage"]
        direction = summary["direction_description"]
        paired = summary["paired_regulator_result"]
    
        resolved_pct = 100 * coverage["n_resolved"] / coverage[
            "n_endpoint_defined_reversals_in_8h_input"
        ]
        unresolved_pct = 100 - resolved_pct
        bar_height = 0.28
    
        # All endpoint-defined reversals: the unresolved majority remains visible.
        ax.barh(1.45, resolved_pct, left=0, height=bar_height, color=RESOLVED, linewidth=0)
        ax.barh(
            1.45,
            unresolved_pct,
            left=resolved_pct,
            height=bar_height,
            color=UNRESOLVED_COLOUR,
            linewidth=0,
        )
        _label_segment(
            ax,
            0,
            resolved_pct,
            1.45,
            f"resolved\n{coverage['n_resolved']:,} ({resolved_pct:.1f}%)",
            "white",
        )
        _label_segment(
            ax,
            resolved_pct,
            unresolved_pct,
            1.45,
            f"unresolved\n{coverage['n_unresolved']:,} ({unresolved_pct:.1f}%)",
            INK,
        )
        ax.text(
            -3.0,
            1.45,
            "8 h-covered\n(n=4,377)",
            ha="right",
            va="center",
            fontsize=fs,
            color=MUTED,
            linespacing=0.95,
        )
    
        # Conditional composition among resolvable edges, split by reversal orientation.
        rows = [
            (0.76, LOST, "Resolved lost", 704, LOST_COLOUR),
            (0.24, GAINED, "Resolved gained", 575, GAINED_COLOUR),
        ]
        for y, key, label, n_resolved, row_colour in rows:
            rest_pct = 100 * float(direction[key]["fraction_Rest_sign_retaining"])
            stim_pct = 100 - rest_pct
            ax.barh(y, rest_pct, left=0, height=bar_height, color=REST_COLOUR, linewidth=0)
            ax.barh(
                y,
                stim_pct,
                left=rest_pct,
                height=bar_height,
                color=STIM48_COLOUR,
                linewidth=0,
            )
            rest_count = int(direction[key]["n_Rest_sign_retaining"])
            stim_count = int(direction[key]["n_resolved"]) - rest_count
            _label_segment(
                ax,
                0,
                rest_pct,
                y,
                f"Rest sign\n{rest_count:,} ({rest_pct:.1f}%)",
                "white",
            )
            _label_segment(
                ax,
                rest_pct,
                stim_pct,
                y,
                f"48 h sign\n{stim_count:,} ({stim_pct:.1f}%)",
                "white",
            )
            ax.text(
                -3.0,
                y,
                f"{label}\n(n={n_resolved})",
                ha="right",
                va="center",
                fontsize=fs,
                color=row_colour,
                fontweight="bold",
                linespacing=0.95,
            )
    
        ax.set_xlim(-3.0, 100.0)
        ax.set_ylim(-0.05, 1.75)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.text(
            0.5,
            0.00,
            f"paired regulators n={paired['n_regulators_with_both_directions']}; "
            f"p={paired['p_value']:.4f}",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=fs,
            color=MUTED,
            clip_on=False,
        )
    
    
    def build_standalone() -> None:
        states, summary = load_and_validate()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4.0, 2.75), facecolor="white")
        panel(ax, states, summary, fs=7.0)
        ax.set_title(
            "Lost reversals retain the Rest sign more often at 8 h",
            fontsize=8.2,
            fontweight="bold",
            color=INK,
            pad=5,
        )
        fig.subplots_adjust(left=0.19, right=0.98, top=0.86, bottom=0.19)
        out_pdf = OUTDIR / "Fig2h_8hr_sign_states.pdf"
        out_png = OUTDIR / "Fig2h_8hr_sign_states.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
    
        manifest = {
            "schema_version": "fig2h-8hr-sign-states-render-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "analysis_status": summary["analysis_status"],
            "renderer": {
                "path": "scripts/fig2h_8hr_sign_states.py",
                "sha256": sha256(Path(__file__).resolve()),
            },
            "sources": {
                str(path.relative_to(PAPER)): sha256(path) for path in EXPECTED_HASHES
            },
            "displayed_contract": {
                "covered": 4_377,
                "resolved": 1_279,
                "unresolved": 3_098,
                "lost_resolved": 704,
                "lost_rest_sign_retaining": 529,
                "gained_resolved": 575,
                "gained_rest_sign_retaining": 295,
                "paired_regulators": 75,
                "wilcoxon_two_sided_p": 0.0049970904,
            },
            "claim_boundary": summary["claim_boundary"],
            "outputs": {
                "pdf": {
                    "path": str(out_pdf.relative_to(PAPER)),
                    "sha256": sha256(out_pdf),
                },
                "png": {
                    "path": str(out_png.relative_to(PAPER)),
                    "sha256": sha256(out_png),
                },
            },
        }
        (OUTDIR / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2h_8hr_sign_states, k, v)
init_fig2h_8hr_sign_states()

def init_fig2_eight_panel_8h_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render Figure 2 with the audited 8 h sign-state result in panel h.
    
    This is a noncanonical successor candidate.  It reuses the promoted a--g
    renderer byte-for-byte, clears only its final panel axis in memory, and redraws
    that axis from the frozen 8 h sign-state artifacts.  No existing file is moved,
    deleted, or overwritten.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    PREDECESSOR = fig2_eight_panel_candidate
    SIGN8 = fig2h_8hr_sign_states
    
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_eight_panel_8h_candidate"
    STEM = "Fig2_eight_panel_8h_candidate"
    CANONICAL_PDF = PAPER / "figures/Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures/Fig2.png"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _build() -> tuple[plt.Figure, dict]:
        # The predecessor reconstructs the exact promoted a--g geometry.  Its h axis
        # is then cleared before any output is written.
        fig, _ = PREDECESSOR._build()
        assert len(fig.axes) == 9
        ax_h = fig.axes[-1]
        assert "activation–movement" in ax_h.get_title()
        ax_h.clear()
    
        states, summary = SIGN8.load_and_validate()
        SIGN8.panel(ax_h, states, summary, fs=7.0)
        PREDECESSOR.BASE.S.panel_tag(ax_h, "h", dx=-0.25, dy=1.07)
        ax_h.set_title(
            "Lost reversals more often retain\n"
            "the Rest sign at 8 h",
            fontsize=7.3,
            pad=3,
            linespacing=1.0,
        )
    
        minimum_font, raised = PREDECESSOR.BASE._enforce_minimum_font(fig)
        assert minimum_font >= PREDECESSOR.BASE.MIN_SOURCE_FONT_PT
        all_text = "\n".join(text.get_text() for text in fig.findobj(match=Text))
        assert "p=0.572" not in all_text
        assert "75.1%" in all_text and "51.3%" in all_text
        assert "unresolved" in all_text and "70.8%" in all_text
        return fig, {
            "minimum_source_font_pt": float(minimum_font),
            "text_objects_raised_to_floor": len(raised),
            "panel_h": {
                "analysis_status": summary["analysis_status"],
                "covered_reversals": 4_377,
                "resolved": 1_279,
                "unresolved": 3_098,
                "lost_resolved": 704,
                "lost_rest_sign_retaining": 529,
                "lost_rest_sign_retaining_pct": 75.1,
                "gained_resolved": 575,
                "gained_rest_sign_retaining": 295,
                "gained_rest_sign_retaining_pct": 51.3,
                "paired_regulators": 75,
                "wilcoxon_two_sided_p": 0.0049970904,
            },
        }
    
    
    def build() -> None:
        canonical_before = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        fig, audit = _build()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out_pdf = OUTDIR / f"{STEM}.pdf"
        out_png = OUTDIR / f"{STEM}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
    
        canonical_after = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        assert canonical_before == canonical_after
        sources = {
            "data/kinetics/kinetics_8hr_sign_states.csv.gz": SIGN8.SIGN_STATES,
            "data/kinetics/kinetics_8hr_sign_states.summary.json": SIGN8.SUMMARY,
            "data/kinetics/kinetics_8hr_sign_states.paired_regulators.tsv": SIGN8.PAIRED,
            "data/kinetics/kinetics_8hr_sign_states.manifest.json": (
                SIGN8.SOURCE_MANIFEST
            ),
        }
        renderer_paths = (
            "scripts/fig2_eight_panel_8h_candidate.py",
            "scripts/fig2h_8hr_sign_states.py",
            "scripts/fig2_eight_panel_candidate.py",
            "scripts/fig2h_activation_movement_similarity.py",
            "scripts/fig2_biology_candidate.py",
            "scripts/figstyle.py",
            "scripts/fig4_architecture.py",
            "scripts/fig_target_activation_trajectory_v1.py",
            "scripts/fig2g_trajectory_panel.py",
        )
        manifest = {
            "schema_version": "fig2-eight-panel-8h-sign-state-candidate-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "panel_order": ["a", "b", "c", "d", "e", "f", "g", "h"],
            "replacement": {
                "panel": "h",
                "superseded_visual_in_candidate_only": (
                    "post-hoc activation-movement comparison; its files remain untouched"
                ),
                "new_content": "8 h sign states among endpoint-defined reversals",
                "rendering": "direct from sealed tables; no PDF embedding",
            },
            "renderer_closure": {
                path: sha256(PAPER / path) for path in renderer_paths
            },
            "panel_h_sources": {
                path: sha256(source) for path, source in sources.items()
            },
            "canonical_figure_2_unchanged": canonical_after,
            "audit": audit,
            "claim_boundary": (
                "Panel h is post-hoc and conditional on endpoint-defined reversals. "
                "It compares detectable 8 h sign states, not exact switch times or mechanism."
            ),
            "outputs": {
                "pdf": {
                    "path": str(out_pdf.relative_to(PAPER)),
                    "sha256": sha256(out_pdf),
                },
                "png": {
                    "path": str(out_png.relative_to(PAPER)),
                    "sha256": sha256(out_png),
                },
            },
        }
        manifest_path = OUTDIR / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"minimum visible source font: {audit['minimum_source_font_pt']:.2f} pt")
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
        print(f"wrote {manifest_path}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_eight_panel_8h_candidate, k, v)
init_fig2_eight_panel_8h_candidate()

def init_fig2h_8hr_sign_states_scatter():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the Figure 2e 8 h result as a retention plot plus composition key.
    
    Two lines show the share of resolved lost and gained reversals retaining the
    Rest sign across the three sampled conditions.  Rest and 48 h are open markers
    because their values are fixed by the endpoint reversal definition; only the
    8 h values discriminate timing.  Adjacent 100-percent bars show the complete
    8 h compositions and the resolved/unresolved classification coverage.
    
    The historical ``scatter`` filename is retained so this visual-only revision
    overwrites its own noncanonical diagnostic rather than creating another copy.
    """
    
    
    
    matplotlib.use("Agg")
    
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    
    SOURCE = fig2h_8hr_sign_states
    
    
    PAPER = ROOT
    OUTDIR = PAPER / "data/figure_provenance/fig2h_8hr_sign_states_scatter"
    
    INK = "#24292e"
    MUTED = "#66717a"
    # Paper-wide reversal-direction convention (scripts/figstyle.py).
    LOST = "#c1272d"
    GAINED = "#2a9d8f"
    TIED = "#8f989f"
    GRID = "#dde1e4"
    REST_SIGN = "#2c6fbb"
    STIM48_SIGN = "#e07b39"
    RESOLVED = "#59636d"
    UNRESOLVED_COLOUR = "#d6dadd"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def load_and_validate() -> tuple[pd.DataFrame, dict, pd.DataFrame]:
        states, summary = SOURCE.load_and_validate()
        paired = pd.read_csv(SOURCE.PAIRED, sep="\t")
        required = {
            "regulator_ensg",
            "regulator",
            "lost_n_resolved",
            "lost_n_rest_sign_retaining",
            "lost_rest_sign_retaining_fraction",
            "gained_n_resolved",
            "gained_n_rest_sign_retaining",
            "gained_rest_sign_retaining_fraction",
            "lost_minus_gained_exact_fraction",
        }
        assert required.issubset(paired.columns)
        assert len(paired) == 75
        assert paired["regulator_ensg"].is_unique
    
        exact_difference = paired["lost_minus_gained_exact_fraction"].map(Fraction)
        paired = paired.assign(
            relation=np.select(
                [exact_difference > 0, exact_difference < 0],
                ["lost_higher", "gained_higher"],
                default="tied",
            ),
            gained_pct=100 * paired["gained_rest_sign_retaining_fraction"],
            lost_pct=100 * paired["lost_rest_sign_retaining_fraction"],
        )
        assert paired["relation"].value_counts().to_dict() == {
            "lost_higher": 40,
            "tied": 21,
            "gained_higher": 14,
        }
        # These are the paired-regulator analysis rows, not all resolved edges.
        assert int(paired["lost_n_resolved"].sum()) == 651
        assert int(paired["lost_n_rest_sign_retaining"].sum()) == 485
        assert int(paired["gained_n_resolved"].sum()) == 506
        assert int(paired["gained_n_rest_sign_retaining"].sum()) == 257
        return states, summary, paired
    
    
    def panel(
        ax: plt.Axes,
        states: pd.DataFrame,
        summary: dict,
        paired: pd.DataFrame,
        fs: float = 7.0,
        compact: bool = False,
        show_side_summary: bool = True,
    ) -> None:
        """Rest-sign retention across sampled states, with fixed endpoints marked."""
        del states  # Its validated coverage is reported explicitly below.
        direction = summary["direction_description"]
        assert len(paired) == 75
        assert paired["relation"].value_counts().to_dict() == {
            "lost_higher": 40,
            "tied": 21,
            "gained_higher": 14,
        }
        lost = direction["lost_on_activation"]
        gained = direction["gained_on_activation"]
        assert (
            lost["n_Rest_sign_retaining"],
            lost["n_resolved"],
            float(lost["percentage_Rest_sign_retaining"]),
        ) == (
            529, 704, 75.1)
        assert (
            gained["n_Rest_sign_retaining"],
            gained["n_resolved"],
            float(gained["percentage_Rest_sign_retaining"]),
        ) == (
            295, 575, 51.3)
        lost_pct = 100.0 * lost["n_Rest_sign_retaining"] / lost["n_resolved"]
        gained_pct = 100.0 * gained["n_Rest_sign_retaining"] / gained["n_resolved"]
        x = np.asarray([0.0, 1.0, 2.0])
        lost_values = np.asarray([100.0, lost_pct, 0.0])
        gained_values = np.asarray([100.0, gained_pct, 0.0])
        ax.plot(x, lost_values, color=LOST, lw=1.7, zorder=2)
        ax.plot(
            x, gained_values, color=GAINED, lw=1.7, ls="--", zorder=2,
        )
        # Open endpoint markers disclose values fixed by the reversal definition.
        ax.scatter(
            [0, 2], [100, 0], s=27, facecolor="white", edgecolor=LOST,
            linewidth=1.0, marker="o", zorder=3,
        )
        ax.scatter(
            [0, 2], [100, 0], s=13, facecolor="white", edgecolor=GAINED,
            linewidth=0.9, marker="s", zorder=4,
        )
        ax.scatter([1], [lost_pct], s=24, color=LOST, marker="o", linewidth=0, zorder=4)
        ax.scatter(
            [1], [gained_pct], s=22, color=GAINED, marker="s", linewidth=0, zorder=4,
        )
        ax.text(
            1, lost_pct + 3.0, f"Lost {lost_pct:.1f}%", ha="center", va="bottom",
            fontsize=fs, fontweight="bold", color=LOST,
        )
        ax.text(
            1, gained_pct - 3.0, f"Gained {gained_pct:.1f}%", ha="center", va="top",
            fontsize=fs, fontweight="bold", color=GAINED,
        )
        ax.text(
            1.0, 96.0, "open endpoints fixed", ha="center", va="top",
            fontsize=fs, color=MUTED,
        )
        ax.set_xlim(-0.16, 2.16)
        ax.set_ylim(-6, 106)
        ax.set_xticks([0, 1, 2], ["Rest", "8 h", "48 h"])
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel("resolved edges retaining Rest sign (%)", fontsize=fs, labelpad=2)
        ax.tick_params(labelsize=fs, colors=MUTED, length=2.2, width=0.65)
        ax.grid(axis="y", color=GRID, linewidth=0.45, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(MUTED)
        ax.spines[["left", "bottom"]].set_linewidth(0.65)
    
    
    def mini_key_panel(ax: plt.Axes, summary: dict, fs: float = 7.0) -> None:
        """Compact rendering of the prior stacked sign-state panel as a visual key."""
        coverage = summary["coverage"]
        resolved_pct = 100 * coverage["n_resolved"] / coverage[
            "n_endpoint_defined_reversals_in_8h_input"
        ]
        direction = summary["direction_description"]
        lost = direction["lost_on_activation"]
        gained = direction["gained_on_activation"]
        lost_pct = 100 * lost["n_Rest_sign_retaining"] / lost["n_resolved"]
        gained_pct = 100 * gained["n_Rest_sign_retaining"] / gained["n_resolved"]
        rows = (
            (2.0, "8 h-covered (n=4,377)", resolved_pct, 100 - resolved_pct,
             RESOLVED, UNRESOLVED_COLOUR,
             "resolved\n1,279",
             "unresolved\n3,098"),
            (1.0, "resolved lost (n=704)", lost_pct, 100 - lost_pct,
             REST_SIGN, STIM48_SIGN,
             f"Rest\n{lost_pct:.1f}%", f"48 h\n{100 - lost_pct:.1f}%"),
            (0.0, "resolved gained (n=575)", gained_pct, 100 - gained_pct,
             REST_SIGN, STIM48_SIGN,
             f"Rest\n{gained_pct:.1f}%", f"48 h\n{100 - gained_pct:.1f}%"),
        )
        for y, label, left, right, colour_left, colour_right, text_left, text_right in rows:
            ax.barh(y, left, height=0.42, color=colour_left, linewidth=0)
            ax.barh(y, right, left=left, height=0.42, color=colour_right, linewidth=0)
            row_label = ax.text(
                0,
                y + 0.29,
                label,
                ha="left",
                va="bottom",
                fontsize=fs,
                fontweight="bold" if label.startswith("resolved") else "normal",
                color=(LOST if label.startswith("resolved lost") else
                       GAINED if label.startswith("resolved gained") else MUTED),
            )
            row_label.set_gid("mini_key_row_label")
            left_text = ax.text(
                left / 2,
                y,
                text_left,
                ha="center",
                va="center",
                fontsize=fs,
                fontweight="bold",
                color="white",
                linespacing=0.92,
            )
            left_text.set_gid("mini_key_segment_label")
            right_text = ax.text(
                left + right / 2,
                y,
                text_right,
                ha="center",
                va="center",
                fontsize=fs,
                fontweight="bold",
                color=INK if colour_right == UNRESOLVED_COLOUR else "white",
                linespacing=0.92,
            )
            right_text.set_gid("mini_key_segment_label")
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.92, 2.70)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title("8 h classification", fontsize=fs, fontweight="bold", pad=3)
        ax.text(
            50,
            -0.66,
            "p=0.0050\n75 paired regulators",
            ha="center",
            va="center",
            fontsize=fs,
            color=MUTED,
            linespacing=1.02,
        )
    
    
    def build_standalone() -> None:
        states, summary, paired = load_and_validate()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        fig = plt.figure(figsize=(5.8, 3.25), facecolor="white")
        ax = fig.add_axes([0.09, 0.17, 0.34, 0.68])
        key_ax = fig.add_axes([0.52, 0.20, 0.44, 0.60])
        panel(ax, states, summary, paired, fs=7.0)
        mini_key_panel(key_ax, summary, fs=7.0)
        ax.set_title("Rest-sign retention", fontsize=8.0, fontweight="bold", pad=4)
        fig.suptitle(
            "Resolved lost and gained reversals differ at 8 h",
            fontsize=8.2,
            fontweight="bold",
            color=INK,
            y=0.95,
        )
        out_pdf = OUTDIR / "Fig2h_8hr_sign_states_scatter.pdf"
        out_png = OUTDIR / "Fig2h_8hr_sign_states_scatter.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
    
        manifest = {
            "schema_version": "fig2h-8hr-sign-states-retention-plus-composition-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "analysis_status": summary["analysis_status"],
            "renderer": {
                "path": "scripts/fig2h_8hr_sign_states_scatter.py",
                "sha256": sha256(Path(__file__).resolve()),
            },
            "source_renderer": {
                "path": "scripts/fig2h_8hr_sign_states.py",
                "sha256": sha256(Path(SOURCE.__file__).resolve()),
            },
            "sources": {
                str(path.relative_to(PAPER)): sha256(path)
                for path in SOURCE.EXPECTED_HASHES
            },
            "display_contract": {
                "retention_values": {
                    "8h_lost_Rest_sign": "529/704 resolved lost reversals (75.1%)",
                    "8h_lost_48h_sign": "175/704 resolved lost reversals (24.9%)",
                    "8h_gained_Rest_sign": "295/575 resolved gained reversals (51.3%)",
                    "8h_gained_48h_sign": "280/575 resolved gained reversals (48.7%)",
                },
                "display_unit": "resolved reversal edges",
                "p_value_unit": "75 paired regulator-level fractions",
                "colour_convention": {
                    "lost": LOST,
                    "gained": GAINED,
                    "Rest_sign": REST_SIGN,
                    "48h_sign": STIM48_SIGN,
                },
                "retention_lines": {
                    "lost": [100.0, 100.0 * 529 / 704, 0.0],
                    "gained": [100.0, 100.0 * 295 / 575, 0.0],
                },
                "endpoint_marker_contract": (
                    "Rest and 48 h are open because their sign-state shares are "
                    "fixed by the endpoint-defined reversal criterion"
                ),
                "paired_regulators": 75,
                "lost_higher": 40,
                "gained_higher": 14,
                "tied": 21,
                "caption_edge_counts_not_plotted_as_the_inferential_unit": {
                    "gained": "295/575 (51.3%)",
                    "lost": "529/704 (75.1%)",
                },
                "resolved": 1_279,
                "unresolved": 3_098,
                "wilcoxon_two_sided_p": 0.0049970904,
            },
            "claim_boundary": summary["claim_boundary"],
            "outputs": {
                "pdf": {
                    "path": str(out_pdf.relative_to(PAPER)),
                    "sha256": sha256(out_pdf),
                },
                "png": {
                    "path": str(out_png.relative_to(PAPER)),
                    "sha256": sha256(out_png),
                },
            },
        }
        (OUTDIR / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2h_8hr_sign_states_scatter, k, v)
init_fig2h_8hr_sign_states_scatter()

def init_fig2_eight_panel_original_line_a4_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the original Figure 2 arrangement with the first simple 8 h line panel.
    
    This candidate deliberately keeps panels a--d and the single ordered 49-row
    panel b in their established A4 geometry.  It swaps the 903-target orientation
    analysis to h, restores the first two-line Rest--8 h--48 h display in e, and
    uses only the existing lower whitespace to give f--h more height.  It never
    writes canonical Figure 2 paths and never removes prior diagnostics.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    PREDECESSOR = fig2_eight_panel_8h_candidate
    PROFILE8 = fig2h_8hr_sign_states_scatter
    
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_eight_panel_original_line_a4_candidate"
    STEM = "Fig2_eight_panel_original_line_a4_candidate"
    CANONICAL_PDF = PAPER / "figures/Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures/Fig2.png"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _build() -> tuple[plt.Figure, dict]:
        fig, _ = PREDECESSOR._build()
        assert tuple(fig.get_size_inches()) == (8.27, 11.69)
        assert len(fig.axes) == 9
        ax_a, ax_b, ax_b_meta, ax_c, ax_d = fig.axes[:5]
        ax_orientation = next(
            ax for ax in fig.axes if ax.get_title().startswith("Target orientation")
        )
        ax_f, ax_g, ax_line = fig.axes[-3], fig.axes[-2], fig.axes[-1]
        assert "Rest sign at 8 h" in ax_line.get_title()
        ax_line.clear()
    
        base = PREDECESSOR.PREDECESSOR.BASE
    
        # Preserve a--d and the complete one-column panel b exactly.  Only enlarge
        # the existing lower row into whitespace and move the unchanged orientation
        # analysis into its h slot.
        left, right = 0.080, 0.980
        row_y, row_h = 0.045, 0.145
        gap_fg, gap_gh = 0.035, 0.050
        usable = right - left - gap_fg - gap_gh
        ratios = (1.55, 0.80, 1.15)
        unit = usable / sum(ratios)
        width_f, width_g, width_h = (ratio * unit for ratio in ratios)
        x_f = left
        x_g = x_f + width_f + gap_fg
        x_h = x_g + width_g + gap_gh
        ax_f.set_position([x_f, row_y, width_f, row_h])
        ax_g.set_position([x_g, row_y, width_g, row_h])
        ax_orientation.set_position([x_h, row_y, width_h, row_h])
    
        base._remove_panel_tag(ax_f, "f")
        base._remove_panel_tag(ax_g, "g")
        base._remove_panel_tag(ax_orientation, "e")
        base.S.panel_tag(ax_f, "f", dx=-0.04, dy=1.05)
        base.S.panel_tag(ax_g, "g", dx=-0.04, dy=1.05)
        base.S.panel_tag(ax_orientation, "h", dx=0.03, dy=1.14)
        ax_f.set_title("Relative separation peaks at 4–8 h", fontsize=7.5, pad=3)
        ax_g.set_title("Both groups rise at 48 h", fontsize=7.5, pad=3)
        ax_orientation.set_title("Target orientation (n=903)", fontsize=7.3, pad=3)
        ax_orientation.set_xlabel("signed reversal score", fontsize=7.0, labelpad=2)
        ax_orientation.set_ylabel("48−0 h", fontsize=7.0, labelpad=1)
        ax_orientation.tick_params(labelsize=7.0)
        for artist in ax_orientation.texts:
            if artist.get_text().startswith("rho ="):
                artist.set_text("rho=.099\nnull −.214/−.107; pMC 1e−5/1e−4")
                artist.set_position((0.03, 0.96))
                artist.set_fontsize(7.0)
                artist.set_linespacing(1.20)
    
        middle = ax_d.get_position()
        # Keep the 100% tick clear of panel d's right-hand heading while retaining
        # the established d/e arrangement.
        ax_line.set_position([0.620, middle.y0, 0.1472, middle.height])
        key_ax = fig.add_axes([
            0.7805,
            middle.y0 + 0.0170,
            0.1995,
            middle.height - 0.0383,
        ])
        states, summary, paired = PROFILE8.load_and_validate()
        PROFILE8.panel(ax_line, states, summary, paired, fs=7.0)
        PROFILE8.mini_key_panel(key_ax, summary, fs=7.0)
        ax_line.set_ylabel("Rest-sign retention (%)", fontsize=7.0, labelpad=1)
        base.S.panel_tag(ax_line, "e", dx=-0.02, dy=1.08)
        ax_line.set_title("Rest-sign retention", fontsize=7.3, pad=3)
    
        # Only bar-internal and row-header text is compacted; all ordinary labels
        # stay at the established 7 pt source size.
        for artist in key_ax.texts:
            if artist.get_gid() in {"mini_key_segment_label", "mini_key_row_label"}:
                artist.set_fontsize(5.5)
        visible_sizes = [
            artist.get_fontsize()
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        assert min(visible_sizes) >= 5.5
        all_text = "\n".join(
            artist.get_text() for artist in fig.findobj(match=Text)
            if artist.get_visible()
        )
        assert "75.1%" in all_text and "51.3%" in all_text
        assert "24.9%" in all_text and "48.7%" in all_text
        assert "75 paired regulators" in all_text
        assert "p=0.572" not in all_text
        return fig, {
            "canvas": {"format": "A4 portrait", "inches": [8.27, 11.69]},
            "layout": "established Figure 2 arrangement; panel b unchanged",
            "minimum_source_font_pt": float(min(visible_sizes)),
            "ordinary_text_floor_pt": 7.0,
            "compact_key_text_pt": 5.5,
            "panel_e": {
                "encoding": (
                    "first simple Lost/Gained Rest--8 h--48 h retention lines; "
                    "selection-fixed endpoints open; exact classification bars adjacent"
                ),
                "resolved": 1_279,
                "unresolved": 3_098,
                "paired_regulators": 75,
                "wilcoxon_two_sided_p": 0.0049970904,
            },
            "panel_h": "unchanged 903-target orientation analysis; geometry only",
            "lower_row_height_inches": row_h * 11.69,
        }
    
    
    def build() -> None:
        canonical_before = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        fig, audit = _build()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out_pdf = OUTDIR / f"{STEM}.pdf"
        out_png = OUTDIR / f"{STEM}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
        canonical_after = {"pdf": sha256(CANONICAL_PDF), "png": sha256(CANONICAL_PNG)}
        assert canonical_before == canonical_after
    
        renderer_paths = (
            "scripts/fig2_eight_panel_original_line_a4_candidate.py",
            "scripts/fig2_eight_panel_8h_candidate.py",
            "scripts/fig2h_8hr_sign_states_scatter.py",
            "scripts/fig2h_8hr_sign_states.py",
            "scripts/fig2_eight_panel_candidate.py",
            "scripts/fig2h_activation_movement_similarity.py",
            "scripts/fig2_biology_candidate.py",
            "scripts/figstyle.py",
            "scripts/fig4_architecture.py",
            "scripts/fig_target_activation_trajectory_v1.py",
            "scripts/fig2g_trajectory_panel.py",
        )
        source_paths = (
            PROFILE8.SOURCE.SIGN_STATES,
            PROFILE8.SOURCE.SUMMARY,
            PROFILE8.SOURCE.PAIRED,
            PROFILE8.SOURCE.SOURCE_MANIFEST,
        )
        manifest = {
            "schema_version": "fig2-original-arrangement-first-line-a4-candidate-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "renderer_closure": {
                str(path): sha256(PAPER / path) for path in renderer_paths
            },
            "eight_hour_sources": {
                str(path.relative_to(PAPER)): sha256(path) for path in source_paths
            },
            "canonical_figure_2_unchanged": canonical_after,
            "audit": audit,
            "claim_boundary": (
                "Panel e is post-hoc and conditional on endpoint-defined reversals; "
                "it compares detectable 8 h sign retention, not exact switch time or mechanism."
            ),
            "outputs": {
                "pdf": {"path": str(out_pdf.relative_to(PAPER)), "sha256": sha256(out_pdf)},
                "png": {"path": str(out_png.relative_to(PAPER)), "sha256": sha256(out_png)},
            },
        }
        manifest_path = OUTDIR / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"minimum visible source font: {audit['minimum_source_font_pt']:.2f} pt")
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
        print(f"wrote {manifest_path}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_eight_panel_original_line_a4_candidate, k, v)
init_fig2_eight_panel_original_line_a4_candidate()

def init_fig2_panel_b_smallbar_visibility_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render Figure 2 with a visibility-only adjustment to panel b.
    
    The full 0--300 linear count panel is unchanged.  A second aligned 0--15
    linear inset enlarges every row with at most 10 reversals without changing the
    main bars.
    Canonical Figure 2 files are never written by this candidate renderer.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    ORIGINAL = fig2_eight_panel_original_line_a4_candidate
    
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_panel_b_smallbar_visibility_candidate"
    STEM = "Fig2_panel_b_smallbar_visibility_candidate"
    CANONICAL_PDF = PAPER / "figures/Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures/Fig2.png"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _build() -> tuple[plt.Figure, dict]:
        fig, predecessor_audit = ORIGINAL._build()
        assert len(fig.axes) == 10
        ax_b = fig.axes[1]
        ax_meta = fig.axes[2]
        assert ax_b.get_xlabel() == "reversing target edges (count)"
    
        base = ORIGINAL.PREDECESSOR.PREDECESSOR.BASE
        sub = base.FA.fig2a_test.loc[
            base.FA.fig2a_test.sig_two_sided
            & base.FA.fig2a_test.direction.eq("up")
        ].copy()
        sub = sub.sort_values(
            ["n_inv", "p_two_sided", "gene"], ascending=[False, True, True]
        ).join(base.FA.direction_counts, on="gene")
        assert len(sub) == 49
        assert int(sub.n_inv.sum()) == 2_075
        assert int(sub.n_inv.min()) == 3 and int(sub.n_inv.max()) == 288
    
        lost_only = sub.gained_on_activation.eq(0)
        gained_only = sub.lost_on_activation.eq(0)
        assert int(lost_only.sum()) == 3 and int(gained_only.sum()) == 3
        assert int(sub.loc[lost_only, "n_inv"].sum()) == 23
        assert int(sub.loc[gained_only, "n_inv"].sum()) == 16
    
        header = next(
            artist for artist in ax_meta.texts
            if artist.get_text() == "reversals/tested; rate; fold"
        )
        header.set_text("rev/test; rate; fold; L/G")
        for row in sub.itertuples(index=False):
            old = (
                f"{int(row.n_inv)}/{int(row.n_both)}; {100 * row.rate:.1f}%; "
                f"{row.enrichment_vs_decile:.2f}x"
            )
            artist = next(text for text in ax_meta.texts if text.get_text() == old)
            artist.set_text(
                f"{old}; {int(row.lost_on_activation)}/{int(row.gained_on_activation)}"
            )
    
        # Reclaim the unused space between the smallest bars and the zoom without
        # dropping the descriptive gene-label key.
        annotation_key = ax_b.get_legend()
        assert annotation_key is not None
        annotation_key.set_bbox_to_anchor((0.07, 0.02), transform=ax_b.transAxes)
        annotation_key.get_title().set_text("label colour (descriptive)")
        annotation_key.get_title().set_fontsize(6.0)
        for text in annotation_key.get_texts():
            text.set_fontsize(6.0)
    
        # Duplicate every row with at most 10 reversals in an aligned linear zoom.
        # Using the threshold rather than an arbitrary row count retains all rows
        # tied at 10, including MCAT.
        # inset occupies their own lower-right whitespace, so it does not cover any
        # original bar.  Both axes start at zero and retain additive stacked widths.
        zoom = sub.loc[sub.n_inv.le(10)].copy()
        assert len(zoom) == 22
        assert int(zoom.n_inv.max()) == 10 and int(zoom.n_inv.min()) == 3
        zoom_y = np.arange(len(sub) - len(zoom), len(sub))
        ax_zoom = ax_b.inset_axes([0.64, 0.00, 0.34, 0.45], zorder=8)
        ax_zoom.set_facecolor("#fbfbfb")
        ax_zoom.barh(
            zoom_y,
            zoom.lost_on_activation.to_numpy(),
            color=base.S.PAL["LOST"],
            height=0.66,
            edgecolor="#3f454b",
            linewidth=0.55,
            zorder=3,
        )
        ax_zoom.barh(
            zoom_y,
            zoom.gained_on_activation.to_numpy(),
            left=zoom.lost_on_activation.to_numpy(),
            color=base.S.PAL["GAINED"],
            height=0.66,
            edgecolor="#3f454b",
            linewidth=0.55,
            zorder=3,
        )
        ax_zoom.set_xlim(0, 15)
        ax_zoom.set_xticks([0, 5, 10, 15])
        ax_zoom.set_ylim(len(sub) - 0.35, len(sub) - len(zoom) - 0.65)
        ax_zoom.set_yticks([])
        ax_zoom.xaxis.tick_top()
        ax_zoom.tick_params(axis="x", labelsize=5.5, pad=1, length=2)
        ax_zoom.set_title("Totals ≤10 (0–15)", fontsize=5.5, pad=9)
        for spine in ax_zoom.spines.values():
            spine.set_color("#9aa0a6")
            spine.set_linewidth(0.5)
    
        return fig, {
            "predecessor": predecessor_audit,
            "panel_b": {
                "regulators": 49,
                "both_orientations": 43,
                "lost_only_regulators": 3,
                "lost_only_edges": 23,
                "gained_only_regulators": 3,
                "gained_only_edges": 16,
                "one_direction_edges": 39,
                "total_edges": 2_075,
                "metadata_fields": (
                    "rev/test; rate; fold; L/G"
                ),
                "main_scale": "linear 0–300 count scale, unchanged",
                "zoom_scale": "linear 0–15 count scale",
                "zoom_rows": 22,
                "zoom_total_range": [3, 10],
                "visible_data_baseline": 0,
                "smallest_bar": 3,
                "change": (
                    "aligned 0–15 linear inset for all 22 rows with totals at most "
                    "10, thin neutral outlines around each nonzero segment, and exact "
                    "lost/gained metadata"
                ),
            },
        }
    
    
    def build() -> None:
        canonical_before = {
            "pdf": sha256(CANONICAL_PDF),
            "png": sha256(CANONICAL_PNG),
        }
        fig, audit = _build()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out_pdf = OUTDIR / f"{STEM}.pdf"
        out_png = OUTDIR / f"{STEM}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
        canonical_after = {
            "pdf": sha256(CANONICAL_PDF),
            "png": sha256(CANONICAL_PNG),
        }
        assert canonical_before == canonical_after
    
        renderer_paths = (
            "scripts/fig2_panel_b_smallbar_visibility_candidate.py",
            "scripts/fig2_eight_panel_original_line_a4_candidate.py",
            "scripts/fig2_eight_panel_8h_candidate.py",
            "scripts/fig2h_8hr_sign_states_scatter.py",
            "scripts/fig2h_8hr_sign_states.py",
            "scripts/fig2_eight_panel_candidate.py",
            "scripts/fig2h_activation_movement_similarity.py",
            "scripts/fig2_biology_candidate.py",
            "scripts/figstyle.py",
            "scripts/fig4_architecture.py",
            "scripts/fig_target_activation_trajectory_v1.py",
            "scripts/fig2g_trajectory_panel.py",
        )
        source_paths = (
            "data/census/per_regulator_enrichment_all.csv",
            "data/census/census_master_edges.csv.gz",
        )
        manifest = {
            "schema_version": "fig2-panel-b-smallbar-visibility-candidate-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "renderer_closure": {
                path: sha256(PAPER / path) for path in renderer_paths
            },
            "source_closure": {
                path: sha256(PAPER / path) for path in source_paths
            },
            "canonical_figure_2_unchanged": canonical_after,
            "audit": audit,
            "claim_boundary": (
                "Panel b retains its original linear 0–300 count axis and adds a "
                "separate linear 0–15 inset for all 22 rows with totals at most 10. "
                "Both axes "
                "preserve exact stacked endpoints, ordering, membership, and inference."
            ),
            "outputs": {
                "pdf": {
                    "path": str(out_pdf.relative_to(PAPER)),
                    "sha256": sha256(out_pdf),
                },
                "png": {
                    "path": str(out_png.relative_to(PAPER)),
                    "sha256": sha256(out_png),
                },
            },
        }
        manifest_path = OUTDIR / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
        print(f"wrote {manifest_path}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_panel_b_smallbar_visibility_candidate, k, v)
init_fig2_panel_b_smallbar_visibility_candidate()

def init_fig2_panel_b_main_scale_visibility_candidate():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render Figure 2 with honest one-edge visibility marks on panel b's main scale.
    
    The predecessor's complete 0--300 linear bars and 0--15 inset are retained.
    On the main scale only, every one-edge segment receives a slim colour tick at
    its exact segment midpoint.
    The tick indicates nonzero presence; the underlying bar endpoints remain the
    count encoding. Canonical Figure 2 files are never written by this renderer.
    """
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    PREDECESSOR = fig2_panel_b_smallbar_visibility_candidate
    
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_panel_b_main_scale_visibility_candidate"
    STEM = "Fig2_panel_b_main_scale_visibility_candidate"
    CANONICAL_PDF = PAPER / "figures/Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures/Fig2.png"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()
    
    
    def _build() -> tuple[plt.Figure, dict]:
        fig, predecessor_audit = PREDECESSOR._build()
        assert len(fig.axes) == 10
        ax_b = fig.axes[1]
        assert ax_b.get_xlabel() == "reversing target edges (count)"
        assert np.allclose(ax_b.get_xlim(), (0.0, 300.0))
    
        base = PREDECESSOR.ORIGINAL.PREDECESSOR.PREDECESSOR.BASE
        sub = base.FA.fig2a_test.loc[
            base.FA.fig2a_test["sig_two_sided"]
            & base.FA.fig2a_test["direction"].eq("up")
        ].copy()
        sub = sub.sort_values(
            ["n_inv", "p_two_sided", "gene"], ascending=[False, True, True]
        ).join(base.FA.direction_counts, on="gene")
        assert len(sub) == 49
    
        rasal3 = sub.loc[sub["gene"].eq("RASAL3")]
        assert len(rasal3) == 1
        assert int(rasal3.iloc[0]["lost_on_activation"]) == 1
        assert int(rasal3.iloc[0]["gained_on_activation"]) == 3
    
        one_edge_marks: list[tuple[str, str, float, int, str]] = []
        for y, (_, row) in enumerate(sub.iterrows()):
            gene = str(row["gene"])
            lost = int(row["lost_on_activation"])
            gained = int(row["gained_on_activation"])
            if lost == 1:
                one_edge_marks.append((gene, "lost", 0.5, y, base.S.PAL["LOST"]))
            if gained == 1:
                one_edge_marks.append(
                    (gene, "gained", lost + 0.5, y, base.S.PAL["GAINED"])
                )
    
        expected = {
            ("ACTR8", "lost"),
            ("GSS", "lost"),
            ("ATP5ME", "gained"),
            ("ATP5PO", "gained"),
            ("ATP5F1A", "gained"),
            ("CCNF", "lost"),
            ("NDUFS8", "lost"),
            ("ZAP70", "lost"),
            ("RASAL3", "lost"),
            ("WAPL", "gained"),
        }
        assert {(gene, direction) for gene, direction, *_ in one_edge_marks} == expected
    
        # The coloured core sits at the exact midpoint of each one-edge segment;
        # the neutral halo separates it from the adjacent segment at native size.
        # The axis geometry and all category-coloured gene labels remain untouched.
        ax_b.spines["left"].set_zorder(1)
        for _, _, x, y, colour in one_edge_marks:
            ax_b.vlines(
                x, y - 0.29, y + 0.29, colors="#3f454b", linewidth=1.55,
                zorder=8, clip_on=False,
            )
            ax_b.vlines(
                x, y - 0.29, y + 0.29, colors=colour, linewidth=0.90,
                zorder=9, clip_on=False,
            )
    
        assert np.allclose(ax_b.get_xlim(), (0.0, 300.0))
        return fig, {
            "predecessor": predecessor_audit,
            "panel_b_main_scale_visibility": {
                "main_scale": "linear 0–300 count scale, unchanged",
                "bar_endpoints_unchanged": True,
                "left_spine_position_unchanged": True,
                "one_edge_segment_count": 10,
                "one_edge_rule": "slim colour tick at the exact segment midpoint",
                "tick_outer_colour": "#3f454b",
                "tick_outer_width_pt": 1.55,
                "tick_inner_width_pt": 0.90,
                "tick_role": "nonzero visibility only; bar length remains the count encoding",
                "one_edge_segments": [
                    {"gene": gene, "direction": direction, "x": x}
                    for gene, direction, x, _, _ in one_edge_marks
                ],
                "rasal3": {
                    "lost": 1,
                    "gained": 3,
                    "lost_tick_x": 0.5,
                },
            },
        }
    
    
    def build() -> None:
        canonical_before = {
            "pdf": sha256(CANONICAL_PDF),
            "png": sha256(CANONICAL_PNG),
        }
        fig, audit = _build()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out_pdf = OUTDIR / f"{STEM}.pdf"
        out_png = OUTDIR / f"{STEM}.png"
        fig.savefig(out_pdf, dpi=300, facecolor="white")
        fig.savefig(out_png, dpi=300, facecolor="white")
        plt.close(fig)
        canonical_after = {
            "pdf": sha256(CANONICAL_PDF),
            "png": sha256(CANONICAL_PNG),
        }
        assert canonical_before == canonical_after
    
        renderer_paths = (
            "scripts/fig2_panel_b_main_scale_visibility_candidate.py",
            "scripts/fig2_panel_b_smallbar_visibility_candidate.py",
            "scripts/fig2_eight_panel_original_line_a4_candidate.py",
            "scripts/fig2_eight_panel_8h_candidate.py",
            "scripts/fig2h_8hr_sign_states_scatter.py",
            "scripts/fig2h_8hr_sign_states.py",
            "scripts/fig2_eight_panel_candidate.py",
            "scripts/fig2h_activation_movement_similarity.py",
            "scripts/fig2_biology_candidate.py",
            "scripts/figstyle.py",
            "scripts/fig4_architecture.py",
            "scripts/fig_target_activation_trajectory_v1.py",
            "scripts/fig2g_trajectory_panel.py",
        )
        source_paths = (
            "data/census/per_regulator_enrichment_all.csv",
            "data/census/census_master_edges.csv.gz",
        )
        manifest = {
            "schema_version": "fig2-panel-b-main-scale-visibility-candidate-v1",
            "status": "NONCANONICAL_VISUAL_CANDIDATE",
            "renderer_closure": {path: sha256(PAPER / path) for path in renderer_paths},
            "source_closure": {path: sha256(PAPER / path) for path in source_paths},
            "canonical_figure_2_unchanged": canonical_after,
            "audit": audit,
            "claim_boundary": (
                "Panel b retains exact 0–300 bar endpoints. Main-scale colour ticks "
                "mark the exact midpoint of each one-edge segment because these segments "
                "are otherwise subpixel at page scale; ticks encode presence, not length."
            ),
            "outputs": {
                "pdf": {"path": str(out_pdf.relative_to(PAPER)), "sha256": sha256(out_pdf)},
                "png": {"path": str(out_png.relative_to(PAPER)), "sha256": sha256(out_png)},
            },
        }
        manifest_path = OUTDIR / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {out_pdf}")
        print(f"wrote {out_png}")
        print(f"wrote {manifest_path}")
    
    
    for k, v in list(locals().items()):
        setattr(fig2_panel_b_main_scale_visibility_candidate, k, v)
init_fig2_panel_b_main_scale_visibility_candidate()

def init_fig_signed_kegg_transition_heatmap_1d_exploratory():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Build a reference-inspired 1D heatmap for descriptive KEGG ORA candidates.
    
    This is a separate exploratory diagnostic. It reads the frozen filtered-ORA and
    signed-orientation results, but it does not recompute either analysis and never
    modifies the existing signed KEGG transition map or manuscript artifacts.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    
    SIGNED_JSON = ROOT / "data" / "census" / "signed_kegg_transition_map_exploratory.json"
    SIGNED_TERMS = ROOT / "data" / "census" / "signed_kegg_transition_map_exploratory.terms.tsv"
    FILTERED_ORA = ROOT / "data" / "census" / "unique_reversal_target_kegg_ora_filtered_exploratory.tsv"
    FILTERED_ORA_MANIFEST = ROOT / "data" / "census" / "unique_reversal_target_kegg_ora_filtered_figure_manifest.json"
    # Layout reference was an unrelated figure from a separate project on the author's
    # machine. It is not an input to this figure and is not distributed.
    REFERENCE_PDF = Path("external/Figure_Functional_Shift_1D.pdf")
    
    OUT_PLOT = ROOT / "data" / "census" / "signed_kegg_transition_heatmap_1d_exploratory.plot.tsv"
    OUT_MANIFEST = ROOT / "data" / "census" / "signed_kegg_transition_heatmap_1d_figure_manifest.json"
    OUT_PDF = ROOT / "figures_final" / "diagnostics" / "Fig_signed_kegg_transition_heatmap_1D_exploratory.pdf"
    OUT_PNG = ROOT / "figures_final" / "diagnostics" / "Fig_signed_kegg_transition_heatmap_1D_exploratory.png"
    
    EXPECTED_SHA256 = {
        SIGNED_JSON: "1c1aa3fc6916114e90cbca0594111dc0419e38971ed10e57232ae12fcb390f90",
        SIGNED_TERMS: "33bcd8fc5949b48be940c45525515e9d5e62376e503bdb8e79f6442f5bd18f5d",
        FILTERED_ORA: "6904b877e3f04cfa874f08aee1a70403bd546b1c8d8df48a8e62732a1b1667f4",
        FILTERED_ORA_MANIFEST: "5485572369e2e28ca21c4df34942899f2884fdeef589aa02c8d57173efe6386a",
        REFERENCE_PDF: "df3fff67ee2fb33020b220d5b992b453b462be9c6e89d77c7af8527323c1df7d",
    }
    
    N_CANDIDATES = 14
    SIGNED_SCALE = 5.276379894164125
    ORA_GLOBAL_P = 0.18439078046097696
    PAGE_INCHES = (11.69, 8.27)
    PNG_DPI = 300
    LOST = "#c1272d"
    GAINED = "#2a9d8f"
    INK = "#22262b"
    GREY = "#667078"
    HAIR = "#dfe3e6"
    DISAGREE = "#a65a00"
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    
    
    def display_path(path: Path) -> str:
        try:
            return str(path.relative_to(ROOT))
        except ValueError:
            return str(path)
    
    
    def file_record(path: Path) -> dict[str, object]:
        return {
            "path": display_path(path),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
    
    
    def atomic_text(path: Path, payload: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(payload)
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
        os.close(fd)
        try:
            frame.to_csv(temporary, sep="\t", index=False, float_format="%.12g")
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def atomic_figure(fig: plt.Figure, path: Path, **kwargs: object) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{path.stem}.", suffix=path.suffix, dir=path.parent)
        os.close(fd)
        try:
            fig.savefig(temporary, **kwargs)
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def verify_sources() -> dict[str, dict[str, object]]:
        observed: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_SHA256.items():
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = sha256(path)
            if actual != expected:
                raise AssertionError(f"source drift: {path}: {actual} != {expected}")
            observed[display_path(path)] = file_record(path)
    
        ora_manifest = json.loads(FILTERED_ORA_MANIFEST.read_text(encoding="utf-8"))
        recorded_tsv = ora_manifest["source_data"]["tsv"]
        if recorded_tsv["sha256"] != EXPECTED_SHA256[FILTERED_ORA]:
            raise AssertionError("filtered-ORA manifest does not identify the pinned TSV")
        if ora_manifest["guardrails"]["terms_displayable"] != 113:
            raise AssertionError("filtered-ORA family drifted from 113 terms")
        return observed
    
    
    def load_plot_data() -> tuple[pd.DataFrame, dict[str, object]]:
        signed_summary = json.loads(SIGNED_JSON.read_text(encoding="utf-8"))
        signed = pd.read_csv(SIGNED_TERMS, sep="\t")
        ora = pd.read_csv(FILTERED_ORA, sep="\t")
    
        assert signed_summary["status"] == "EXPLORATORY_LOCALISATION_PASSES_BOTH_SINGLE_MARGIN_GATES"
        assert len(signed) == 113 and int(signed["signed_test_eligible"].sum()) == 60
        assert int(signed["both_single_margin_maxT_supported"].sum()) == 1
        assert len(ora) == 113 and ora["term"].is_unique and signed["term"].is_unique
        assert int(ora["conventional_BH_descriptive"].sum()) == N_CANDIDATES
        assert int(ora["regulator_maxT_term"].sum()) == 0
        assert int(ora["joint_sensitivity_supported"].sum()) == 0
        assert not bool(ora["regulator_global_gate"].any())
        assert math.isclose(float(ora["regulator_maxT_p"].min()), ORA_GLOBAL_P, rel_tol=1e-12)
        assert math.isclose(
            float(signed.loc[signed["signed_test_eligible"], "target_observed_T"].abs().max()),
            SIGNED_SCALE,
            rel_tol=1e-9,
        )
    
        ora_columns = [
            "term",
            "kegg_id",
            "kegg_display_name",
            "category",
            "subcategory",
            "term_background_size",
            "query_overlap",
            "odds_ratio",
            "fisher_p_greater",
            "fisher_bh_q",
            "conventional_BH_descriptive",
            "regulator_maxT_term",
            "joint_sensitivity_supported",
        ]
        signed_columns = [
            "term",
            "n_reversal_dyads",
            "n_lost_dyads",
            "n_gained_dyads",
            "pathway_Z",
            "signed_test_eligible",
            "target_observed_T",
            "target_maxT_p",
            "target_excess_sign",
            "regulator_observed_T",
            "regulator_maxT_p",
            "regulator_excess_sign",
            "concordant_excess_sign",
            "both_single_margin_maxT_supported",
        ]
        selected = ora.loc[ora["conventional_BH_descriptive"], ora_columns].merge(
            signed[signed_columns],
            on="term",
            how="left",
            validate="one_to_one",
        )
        if len(selected) != N_CANDIDATES or selected["target_observed_T"].isna().any():
            raise AssertionError("exact-term merge did not retain all 14 candidates")
        if not bool(selected["signed_test_eligible"].all()):
            raise AssertionError("a descriptive ORA candidate is signed-test-ineligible")
        if not np.isfinite(selected["target_observed_T"]).all():
            raise AssertionError("non-finite tile value")
        if float(selected["target_observed_T"].abs().max()) > SIGNED_SCALE + 1e-12:
            raise AssertionError("a tile would be clipped by the frozen scale")
    
        selected["regulator_sign_relation"] = np.where(
            selected["target_excess_sign"].eq(selected["regulator_excess_sign"]),
            "agrees",
            "disagrees",
        )
        selected = selected.sort_values(
            ["target_observed_T", "kegg_id"],
            ascending=[False, True],
            kind="stable",
        ).reset_index(drop=True)
        selected.insert(0, "display_order", np.arange(1, len(selected) + 1, dtype=int))
        selected["tile_value"] = selected["target_observed_T"]
        selected["tile_encoding"] = "target-conditioned reversal-orientation T"
        selected["ora_interpretation"] = "directionless descriptive target representation"
        selected["support_marker"] = np.where(
            selected["both_single_margin_maxT_supported"],
            "passes both orientation maxT gates",
            "",
        )
    
        antigen = selected.loc[selected["kegg_id"].eq("hsa04612")]
        assert len(antigen) == 1
        antigen = antigen.iloc[0]
        assert bool(antigen["both_single_margin_maxT_supported"])
        assert int(antigen["n_lost_dyads"]) == 12 and int(antigen["n_gained_dyads"]) == 43
        assert math.isclose(float(antigen["target_maxT_p"]), 0.0034498275086245687, rel_tol=1e-9)
        assert math.isclose(float(antigen["regulator_maxT_p"]), 1 / 20_001, rel_tol=1e-9)
        assert int(selected["both_single_margin_maxT_supported"].sum()) == 1
        assert int(selected["regulator_sign_relation"].eq("disagrees").sum()) == 1
    
        coverage = signed_summary["pathway_union_coverage"]
        exact_coverage = {
            "descriptive_family_reversal_targets_covered": 796,
            "descriptive_family_reversal_targets_total": 2_577,
            "descriptive_family_reversal_dyads_covered": 1_443,
            "descriptive_family_reversal_dyads_total": 4_379,
            "signed_test_informative_targets_covered": 604,
            "signed_test_informative_targets_total": 2_259,
            "signed_test_informative_reversal_dyads_covered": 1_138,
            "signed_test_informative_reversal_dyads_total": 3_943,
        }
        for key, value in exact_coverage.items():
            assert coverage[key] == value, f"coverage drift: {key}"
        return selected, exact_coverage
    
    
    def make_figure(plot: pd.DataFrame) -> plt.Figure:
        matplotlib.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "font.size": 7.5,
                "axes.titlesize": 9.0,
                "axes.labelsize": 8.0,
                "xtick.labelsize": 7.0,
                "ytick.labelsize": 7.0,
                "legend.fontsize": 7.0,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "figure.facecolor": "white",
                "axes.facecolor": "white",
            }
        )
        fig = plt.figure(figsize=PAGE_INCHES, facecolor="white", constrained_layout=False)
        label_ax = fig.add_axes([0.035, 0.205, 0.405, 0.59])
        tile_ax = fig.add_axes([0.455, 0.205, 0.030, 0.59])
        note_ax = fig.add_axes([0.505, 0.205, 0.465, 0.59])
        n_rows = len(plot)
        y_positions = np.arange(n_rows)
        for axis in (label_ax, note_ax):
            axis.set_xlim(0, 1)
            axis.set_ylim(n_rows - 0.5, -0.5)
            axis.axis("off")
    
        for row in plot.itertuples(index=False):
            y = row.display_order - 1
            wrapped = "\n".join(textwrap.wrap(row.kegg_display_name, width=43))
            label_ax.text(
                0.99,
                y,
                wrapped,
                ha="right",
                va="center",
                fontsize=7.5,
                linespacing=1.05,
                fontweight="bold" if row.kegg_id == "hsa04612" else "normal",
                color=INK,
            )
            note_ax.axhline(y + 0.5, color=HAIR, linewidth=0.45, zorder=0)
            note_ax.text(0.00, y, f"lost {row.n_lost_dyads:,} · gained {row.n_gained_dyads:,}", ha="left", va="center", fontsize=7.1, color=INK)
            note_ax.text(0.225, y, f"Fisher q={row.fisher_bh_q:.3g}", ha="left", va="center", fontsize=7.1, color=GREY)
            relation_color = GREY if row.regulator_sign_relation == "agrees" else DISAGREE
            note_ax.text(0.405, y, f"regulator sign {row.regulator_sign_relation}", ha="left", va="center", fontsize=7.1, color=relation_color)
            if row.kegg_id == "hsa04612":
                note_ax.scatter([0.665], [y], marker="D", s=28, color=GAINED, edgecolor="#176c63", linewidth=0.5, zorder=3)
                note_ax.text(0.69, y, "passes both orientation\nmaxT gates", ha="left", va="center", fontsize=7.0, linespacing=1.05, fontweight="bold", color="#176c63")
    
        cmap = LinearSegmentedColormap.from_list("gained_white_lost", [GAINED, "#ffffff", LOST])
        norm = TwoSlopeNorm(vmin=-SIGNED_SCALE, vcenter=0.0, vmax=SIGNED_SCALE)
        values = plot["tile_value"].to_numpy(float)[:, None]
        tile_ax.pcolormesh(
            np.array([-0.5, 0.5]),
            np.arange(n_rows + 1) - 0.5,
            values,
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors="white",
            linewidth=0.8,
        )
        tile_ax.set_xlim(-0.5, 0.5)
        tile_ax.set_ylim(n_rows - 0.5, -0.5)
        tile_ax.axis("off")
    
        fig.text(0.035, 0.965, "Signed orientation within descriptive KEGG ORA candidates", ha="left", va="top", fontsize=15.0, fontweight="bold", color=INK)
        fig.text(0.035, 0.925, "All 14 Fisher-BH q < 0.10 candidates (descriptive list); tile color encodes target-conditioned reversal orientation.", ha="left", va="top", fontsize=9.0, color=GREY)
        fig.text(0.035, 0.895, "ORA is directionless; calibrated regulator-opportunity global p = 0.18439; 0 ORA maxT-supported terms.", ha="left", va="top", fontsize=8.2, color=GREY)
    
        fig.text(0.435, 0.815, "KEGG pathway", ha="right", va="bottom", fontsize=7.2, fontweight="bold", color=GREY)
        fig.text(0.470, 0.815, "orientation", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=GREY)
        fig.text(0.505, 0.815, "reversal dyads", ha="left", va="bottom", fontsize=7.2, fontweight="bold", color=GREY)
        fig.text(0.610, 0.815, "descriptive ORA", ha="left", va="bottom", fontsize=7.2, fontweight="bold", color=GREY)
        fig.text(0.695, 0.815, "regulator sensitivity", ha="left", va="bottom", fontsize=7.2, fontweight="bold", color=GREY)
    
        color_ax = fig.add_axes([0.455, 0.125, 0.245, 0.018])
        colorbar = ColorbarBase(color_ax, cmap=cmap, norm=norm, orientation="horizontal")
        colorbar.set_ticks([-SIGNED_SCALE, 0, SIGNED_SCALE])
        colorbar.set_ticklabels([f"−{SIGNED_SCALE:.2f}", "0", f"+{SIGNED_SCALE:.2f}"])
        colorbar.outline.set_edgecolor("#90979d")
        colorbar.outline.set_linewidth(0.6)
        color_ax.tick_params(length=2.5, width=0.6, pad=2, labelsize=7.0)
        fig.text(0.455, 0.154, "Target-conditioned orientation T", ha="left", va="bottom", fontsize=7.5, fontweight="bold", color=INK)
        fig.text(0.455, 0.103, "teal: gained-oriented     white: null centre     red: lost-oriented", ha="left", va="top", fontsize=7.0, color=GREY)
        fig.text(0.735, 0.145, "Antigen processing: PSME1 + TAP1\ncontribute 79.8% of its signed score.", ha="left", va="center", fontsize=7.2, linespacing=1.08, color="#176c63", fontweight="bold")
    
        fig.text(0.035, 0.058, "Coverage: 113-term family, 796/2,577 reversal targets and 1,443/4,379 reversal dyads; 60 signed terms, 604/2,259 informative targets and 1,138/3,943 informative dyads. Memberships overlap.", ha="left", va="bottom", fontsize=7.0, color=GREY)
        fig.text(0.035, 0.030, "Claim boundary: ORA describes target representation; tile color describes orientation. Neither establishes pathway up/down, activation, suppression, movement, or rewiring; the two gates are not a combined p-value.", ha="left", va="bottom", fontsize=7.0, color=GREY)
        return fig
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true", help="replace only this diagnostic's outputs")
        args = parser.parse_args()
        outputs = [OUT_PLOT, OUT_MANIFEST, OUT_PDF, OUT_PNG]
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to overwrite existing outputs without --force: {existing}")
    
        sources = verify_sources()
        plot, coverage = load_plot_data()
        atomic_tsv(plot, OUT_PLOT)
        figure = make_figure(plot)
        atomic_figure(
            figure,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Signed orientation within descriptive KEGG ORA candidates",
                "Subject": "Exploratory repository diagnostic; not manuscript evidence",
            },
        )
        atomic_figure(figure, OUT_PNG, format="png", dpi=PNG_DPI, facecolor="white")
        plt.close(figure)
    
        manifest = {
            "schema_version": "signed-kegg-transition-heatmap-1d-figure-v1",
            "status": "EXPLORATORY_REPOSITORY_DIAGNOSTIC_NOT_MANUSCRIPT_EVIDENCE",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "selection": {
                "rule": "all and only conventional_BH_descriptive == True in the corrected 113-term filtered ORA",
                "n_rows": N_CANDIDATES,
                "all_rows_signed_test_eligible": True,
                "presentation_order": "target_observed_T descending, then KEGG ID",
                "outcome_sorted_for_presentation": True,
            },
            "encoding": {
                "tile": "target_observed_T only",
                "fixed_symmetric_scale": [-SIGNED_SCALE, SIGNED_SCALE],
                "clipping": False,
                "negative": "teal, gained-on-activation orientation",
                "zero": "white, null centre",
                "positive": "red, lost-on-activation orientation",
                "ora_directionality": "directionless descriptive target representation",
                "regulator_sensitivity": "separate sign-agreement text; not averaged with or combined into the tile",
                "support_marker": "filled diamond only for antigen processing and presentation, the sole term passing both complementary orientation maxT gates",
            },
            "ora_calibration": {
                "regulator_opportunity_global_p": ORA_GLOBAL_P,
                "regulator_opportunity_global_gate": False,
                "ora_maxT_supported_terms": 0,
            },
            "coverage": coverage,
            "claim_boundary": "ORA describes target representation and the tile describes reversal orientation. Neither establishes pathway up/down, activation, suppression, movement, rewiring, mechanism, or joint conditioning.",
            "page": {
                "inches": list(PAGE_INCHES),
                "orientation": "landscape A4",
                "png_dpi": PNG_DPI,
                "bbox_inches": None,
                "opaque_background": True,
                "minimum_intended_font_pt": 7.0,
            },
            "sources": {
                **sources,
                "script": file_record(Path(__file__).resolve()),
            },
            "outputs": {
                "plot_tsv": file_record(OUT_PLOT),
                "pdf": file_record(OUT_PDF),
                "png": file_record(OUT_PNG),
            },
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "matplotlib": matplotlib.__version__,
            },
            "alt_text": "Landscape one-dimensional heatmap of all 14 descriptive Fisher-BH KEGG candidates. Rows are sorted from positive lost-oriented to negative gained-oriented target-conditioned T. Right-hand text reports lost and gained reversal-dyad counts, descriptive Fisher q values, and whether the regulator-margin excess sign agrees. Only antigen processing and presentation carries a filled support diamond; it has gained-oriented excess and passes both complementary orientation maxT gates.",
        }
        atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "rows": len(plot),
                    "supported_rows": int(plot["both_single_margin_maxT_supported"].sum()),
                    "pdf": display_path(OUT_PDF),
                    "png": display_path(OUT_PNG),
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig_signed_kegg_transition_heatmap_1d_exploratory, k, v)
init_fig_signed_kegg_transition_heatmap_1d_exploratory()

def init_fig2_signed_kegg_side_by_side_candidate_v1():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render a diagnostic Figure 2 with compact regulator and signed-KEGG panels.
    
    The current promoted Figure 2 is reconstructed in memory from its active
    renderer.  Only the top-right panel-b rectangle is replaced:
    
    * b: one continuous ranks 1--49 horizontal lost/gained reversal-bar list;
    * d: a native compact redraw of the frozen 14-row full-census signed-KEGG list.
    
    Panels a and c retain their geometry.  The lower scientific panels retain
    their geometry and content, with their visible tags relabelled d--h to e--i.
    No PDF is embedded or cropped.  Canonical Figure 2 paths are never written.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    BASE = fig2_biology_candidate
    ACTIVE = fig2_panel_b_main_scale_visibility_candidate
    KEGG = fig_signed_kegg_transition_heatmap_1d_exploratory
    
    
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v1"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v1"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    
    CANONICAL_PDF = PAPER / "figures_final" / "Fig2.pdf"
    CANONICAL_PNG = PAPER / "figures_final" / "Fig2.png"
    EXPECTED_CANONICAL = {
        "pdf": "fdba720a016927f17d4d2ddea4bc01023102b71674d377264d207363cc6b8004",
        "png": "58deb555438b5640f4110c40786b27e21ad4ae222a03abdac3c94ebd3bacb098",
    }
    
    FONT_MIN_PT = 5.0
    FONT_MAX_PT = 7.0
    PNG_DPI = 300
    INK = "#22262b"
    MUTED = "#667078"
    HAIR = "#dfe3e6"
    LOST = BASE.S.PAL["LOST"]
    GAINED = BASE.S.PAL["GAINED"]
    DISAGREE = "#a65a00"
    
    RENDERER_PATHS = (
        "scripts/fig2_signed_kegg_side_by_side_candidate_v1.py",
        "scripts/fig2_panel_b_main_scale_visibility_candidate.py",
        "scripts/fig2_panel_b_smallbar_visibility_candidate.py",
        "scripts/fig2_eight_panel_original_line_a4_candidate.py",
        "scripts/fig2_eight_panel_8h_candidate.py",
        "scripts/fig2h_8hr_sign_states_scatter.py",
        "scripts/fig2h_8hr_sign_states.py",
        "scripts/fig2_eight_panel_candidate.py",
        "scripts/fig2h_activation_movement_similarity.py",
        "scripts/fig2_biology_candidate.py",
        "scripts/figstyle.py",
        "scripts/fig4_architecture.py",
        "scripts/fig_target_activation_trajectory_v1.py",
        "scripts/fig2g_trajectory_panel.py",
        "scripts/fig_signed_kegg_transition_heatmap_1d_exploratory.py",
    )
    
    SCIENTIFIC_SOURCE_PATHS = (
        PAPER / "data/census/per_regulator_enrichment_all.csv",
        PAPER / "data/census/census_master_edges.csv.gz",
        KEGG.SIGNED_JSON,
        KEGG.SIGNED_TERMS,
        KEGG.FILTERED_ORA,
        KEGG.FILTERED_ORA_MANIFEST,
    )
    
    
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    
    
    def relative(path: Path) -> str:
        try:
            return str(path.resolve().relative_to(PAPER.resolve()))
        except ValueError:
            return str(path.resolve())
    
    
    def file_record(path: Path) -> dict[str, object]:
        return {
            "path": relative(path),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
    
    
    def atomic_text(path: Path, payload: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                handle.write(payload)
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        os.close(descriptor)
        try:
            frame.to_csv(temporary, sep="\t", index=False, float_format="%.12g")
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def atomic_figure(fig: plt.Figure, path: Path, **kwargs: object) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{path.stem}.", suffix=path.suffix, dir=path.parent
        )
        os.close(descriptor)
        try:
            fig.savefig(temporary, **kwargs)
            os.replace(temporary, path)
            path.chmod(0o644)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    
    
    def canonical_hashes() -> dict[str, str]:
        observed = {
            "pdf": sha256(CANONICAL_PDF),
            "png": sha256(CANONICAL_PNG),
        }
        if observed != EXPECTED_CANONICAL:
            raise AssertionError(
                f"active canonical Figure 2 drifted: {observed} != {EXPECTED_CANONICAL}"
            )
        return observed
    
    
    def verify_kegg_scientific_sources() -> dict[str, dict[str, object]]:
        """Verify only scientific KEGG inputs, not the historical style PDF."""
        scientific = (
            KEGG.SIGNED_JSON,
            KEGG.SIGNED_TERMS,
            KEGG.FILTERED_ORA,
            KEGG.FILTERED_ORA_MANIFEST,
        )
        records: dict[str, dict[str, object]] = {}
        for path in scientific:
            expected = KEGG.EXPECTED_SHA256[path]
            actual = sha256(path)
            if actual != expected:
                raise AssertionError(f"signed-KEGG source drift: {path}: {actual} != {expected}")
            records[relative(path)] = file_record(path)
        return records
    
    
    def load_regulators() -> pd.DataFrame:
        tested = BASE.FA.fig2a_test
        rows = tested.loc[
            tested["sig_two_sided"] & tested["direction"].eq("up")
        ].copy()
        rows = rows.sort_values(
            ["n_inv", "p_two_sided", "gene"],
            ascending=[False, True, True],
            kind="stable",
        ).join(BASE.FA.direction_counts, on="gene")
        rows = rows.reset_index(drop=True)
        rows.insert(0, "display_rank", np.arange(1, len(rows) + 1, dtype=int))
    
        assert len(rows) == 49
        assert rows["gene"].nunique() == 49
        assert rows["q_two_sided"].lt(0.05).all()
        assert int(rows["n_inv"].sum()) == 2_075
        assert int(rows["n_both"].sum()) == 14_434
        assert int(rows["lost_on_activation"].sum()) == 1_154
        assert int(rows["gained_on_activation"].sum()) == 921
        assert (
            rows["lost_on_activation"] + rows["gained_on_activation"]
            == rows["n_inv"]
        ).all()
        assert int(rows["n_inv"].min()) == 3
        assert int(rows["n_inv"].max()) == 288
        return rows[
            [
                "display_rank",
                "gene",
                "n_inv",
                "n_both",
                "rate",
                "enrichment_vs_decile",
                "p_two_sided",
                "q_two_sided",
                "lost_on_activation",
                "gained_on_activation",
            ]
        ].copy()
    
    
    def reletter_lower_panels(fig: plt.Figure) -> dict[str, str]:
        mapping = {"d": "e", "e": "f", "f": "g", "g": "h", "h": "i"}
        located: dict[str, Text] = {}
        for old in mapping:
            matches = [
                artist
                for artist in fig.findobj(match=Text)
                if artist.get_visible()
                and artist.get_text() == old
                and artist.get_fontweight() == "bold"
                and artist.get_fontsize() >= 8.0
            ]
            if len(matches) != 1:
                raise AssertionError(f"expected one current panel tag {old!r}; found {len(matches)}")
            located[old] = matches[0]
        for old, new in mapping.items():
            located[old].set_text(new)
        return mapping
    
    
    def draw_regulator_panel(
        parent: mpl.axes.Axes,
        rows: pd.DataFrame,
    ) -> dict[str, object]:
        """Draw all 49 horizontal bars in one uninterrupted ranked list."""
        bar_ax = parent.inset_axes([0.015, 0.075, 0.335, 0.825])
        count_ax = parent.inset_axes([0.360, 0.075, 0.090, 0.825])
        y = np.arange(len(rows))
        lost = rows["lost_on_activation"].to_numpy(dtype=int)
        gained = rows["gained_on_activation"].to_numpy(dtype=int)
    
        bar_ax.barh(
            y,
            lost,
            color=LOST,
            height=0.62,
            edgecolor="#3f454b",
            linewidth=0.18,
            zorder=3,
        )
        bar_ax.barh(
            y,
            gained,
            left=lost,
            color=GAINED,
            height=0.62,
            edgecolor="#3f454b",
            linewidth=0.18,
            zorder=3,
        )
    
        one_edge_marks: list[dict[str, object]] = []
        for row_index, (_, row) in enumerate(rows.iterrows()):
            lost_count = int(row["lost_on_activation"])
            gained_count = int(row["gained_on_activation"])
            if lost_count == 1:
                one_edge_marks.append(
                    {
                        "gene": str(row["gene"]),
                        "direction": "lost",
                        "x": 0.5,
                    }
                )
                bar_ax.vlines(
                    0.5,
                    row_index - 0.27,
                    row_index + 0.27,
                    color="#3f454b",
                    linewidth=1.10,
                    zorder=7,
                )
                bar_ax.vlines(
                    0.5,
                    row_index - 0.27,
                    row_index + 0.27,
                    color=LOST,
                    linewidth=0.65,
                    zorder=8,
                )
            if gained_count == 1:
                midpoint = lost_count + 0.5
                one_edge_marks.append(
                    {
                        "gene": str(row["gene"]),
                        "direction": "gained",
                        "x": midpoint,
                    }
                )
                bar_ax.vlines(
                    midpoint,
                    row_index - 0.27,
                    row_index + 0.27,
                    color="#3f454b",
                    linewidth=1.10,
                    zorder=7,
                )
                bar_ax.vlines(
                    midpoint,
                    row_index - 0.27,
                    row_index + 0.27,
                    color=GAINED,
                    linewidth=0.65,
                    zorder=8,
                )
    
        expected_marks = {
            ("ACTR8", "lost"),
            ("GSS", "lost"),
            ("ATP5ME", "gained"),
            ("ATP5PO", "gained"),
            ("ATP5F1A", "gained"),
            ("CCNF", "lost"),
            ("NDUFS8", "lost"),
            ("ZAP70", "lost"),
            ("RASAL3", "lost"),
            ("WAPL", "gained"),
        }
        assert {
            (str(record["gene"]), str(record["direction"]))
            for record in one_edge_marks
        } == expected_marks
    
        bar_ax.set_xlim(0, 300)
        bar_ax.set_xticks([0, 100, 200, 300])
        bar_ax.set_ylim(len(rows) - 0.35, -0.65)
        bar_ax.set_yticks(y)
        bar_ax.set_yticklabels(
            rows["gene"].astype(str).tolist(),
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=INK,
        )
        bar_ax.tick_params(axis="y", length=0, pad=1.0)
        bar_ax.tick_params(axis="x", labelsize=FONT_MIN_PT, length=2.0, pad=1.0)
        bar_ax.set_xlabel("reversal edges (linear 0–300)", fontsize=FONT_MIN_PT, labelpad=1.0)
        bar_ax.grid(axis="x", color=HAIR, linewidth=0.35, zorder=0)
        bar_ax.set_axisbelow(True)
        bar_ax.spines["top"].set_visible(False)
        bar_ax.spines["right"].set_visible(False)
    
        count_ax.set_xlim(0, 1)
        count_ax.set_ylim(len(rows) - 0.35, -0.65)
        count_ax.axis("off")
        count_ax.text(
            0.50,
            -1.15,
            "L / G",
            ha="center",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
        for row_index, (_, row) in enumerate(rows.iterrows()):
            count_ax.text(
                0.46,
                row_index,
                f"{int(row['lost_on_activation'])}",
                ha="right",
                va="center",
                fontsize=FONT_MIN_PT,
                color=LOST,
            )
            count_ax.text(
                0.50,
                row_index,
                "/",
                ha="center",
                va="center",
                fontsize=FONT_MIN_PT,
                color=MUTED,
            )
            count_ax.text(
                0.54,
                row_index,
                f"{int(row['gained_on_activation'])}",
                ha="left",
                va="center",
                fontsize=FONT_MIN_PT,
                color=GAINED,
            )
    
        parent.text(
            0.000,
            0.995,
            "b",
            transform=parent.transAxes,
            ha="right",
            va="top",
            fontsize=FONT_MAX_PT,
            fontweight="bold",
            color=INK,
        )
        parent.text(
            0.225,
            0.995,
            "49 reversal-enriched regulators",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MAX_PT,
            fontweight="bold",
            color=INK,
        )
        parent.text(
            0.225,
            0.952,
            "2,075 reversals; Fisher/BH q<0.05; ranks 1–49",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        parent.text(
            0.225,
            0.018,
            "red lost; teal gained; exact counts shown as L/G",
            transform=parent.transAxes,
            ha="center",
            va="bottom",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        return {
            "regulators": 49,
            "display_order": "descending total reversals, then p, then gene",
            "bar_scale": [0, 300],
            "tested_edges": 14_434,
            "reversal_edges": 2_075,
            "lost_edges": 1_154,
            "gained_edges": 921,
            "metadata_in_artwork": "gene, total bar, lost/gained split, exact L/G counts",
            "metadata_in_source_only": "tested edges, rate, fold, p and q for every regulator",
            "one_edge_marks": one_edge_marks,
            "one_edge_marker_role": "visibility only; exact bar endpoints remain the count encoding",
        }
    
    
    def draw_kegg_panel(
        parent: mpl.axes.Axes,
        plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> dict[str, object]:
        """Redraw the signed-KEGG diagnostic natively in the right half."""
        ax = parent.inset_axes([0.525, 0.075, 0.475, 0.825])
        n_rows = len(plot)
        ax.set_xlim(0, 1)
        ax.set_ylim(n_rows - 0.5, -1.15)
        ax.axis("off")
    
        cmap = LinearSegmentedColormap.from_list(
            "gained_white_lost_compact", [GAINED, "#ffffff", LOST]
        )
        norm = TwoSlopeNorm(
            vmin=-KEGG.SIGNED_SCALE,
            vcenter=0.0,
            vmax=KEGG.SIGNED_SCALE,
        )
    
        columns = {
            "pathway": 0.475,
            "tile_left": 0.505,
            "tile_width": 0.055,
            "counts": 0.585,
            "q": 0.715,
            "relation": 0.895,
            "support": 0.970,
        }
        ax.text(
            columns["pathway"],
            -0.82,
            "KEGG pathway",
            ha="right",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
        ax.text(
            columns["tile_left"] + columns["tile_width"] / 2,
            -0.82,
            "T",
            ha="center",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
        ax.text(
            columns["counts"],
            -0.82,
            "L/G",
            ha="left",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
        ax.text(
            columns["q"],
            -0.82,
            "BH q",
            ha="left",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
        ax.text(
            (columns["relation"] + columns["support"]) / 2,
            -0.82,
            "markers",
            ha="center",
            va="center",
            fontsize=FONT_MIN_PT,
            fontweight="bold",
            color=MUTED,
        )
    
        for row_index, (_, row) in enumerate(plot.iterrows()):
            name = str(row["kegg_display_name"])
            wrapped = "\n".join(textwrap.wrap(name, width=25))
            supported = bool(row["both_single_margin_maxT_supported"])
            ax.text(
                columns["pathway"],
                row_index,
                wrapped,
                ha="right",
                va="center",
                fontsize=FONT_MIN_PT,
                linespacing=0.92,
                fontweight="bold" if supported else "normal",
                color=INK,
            )
            ax.add_patch(
                Rectangle(
                    (columns["tile_left"], row_index - 0.34),
                    columns["tile_width"],
                    0.68,
                    facecolor=cmap(norm(float(row["target_observed_T"]))),
                    edgecolor="#aeb4b9",
                    linewidth=0.28,
                )
            )
            ax.text(
                columns["counts"],
                row_index,
                f"{int(row['n_lost_dyads'])}/{int(row['n_gained_dyads'])}",
                ha="left",
                va="center",
                fontsize=FONT_MIN_PT,
                color=INK,
            )
            ax.text(
                columns["q"],
                row_index,
                f"{float(row['fisher_bh_q']):.3g}",
                ha="left",
                va="center",
                fontsize=FONT_MIN_PT,
                color=MUTED,
            )
            relation = str(row["regulator_sign_relation"])
            if relation == "agrees":
                ax.scatter(
                    [columns["relation"]],
                    [row_index],
                    marker="o",
                    s=8,
                    facecolor="white",
                    edgecolor=MUTED,
                    linewidth=0.55,
                    zorder=4,
                )
            elif relation == "disagrees":
                ax.scatter(
                    [columns["relation"]],
                    [row_index],
                    marker="x",
                    s=11,
                    color=DISAGREE,
                    linewidth=0.70,
                    zorder=4,
                )
            else:
                raise AssertionError(f"unexpected regulator sign relation: {relation}")
            if supported:
                ax.scatter(
                    [columns["support"]],
                    [row_index],
                    marker="D",
                    s=13,
                    color=GAINED,
                    edgecolor="#176c63",
                    linewidth=0.45,
                    zorder=5,
                )
            ax.axhline(row_index + 0.5, color=HAIR, linewidth=0.30, zorder=0)
    
        parent.text(
            0.500,
            0.995,
            "d",
            transform=parent.transAxes,
            ha="right",
            va="top",
            fontsize=FONT_MAX_PT,
            fontweight="bold",
            color=INK,
        )
        parent.text(
            0.762,
            0.995,
            "Full-census signed KEGG",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MAX_PT,
            fontweight="bold",
            color=INK,
        )
        parent.text(
            0.762,
            0.952,
            "Full census: 2,577 targets / 4,379 dyads — not T279",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        parent.text(
            0.762,
            0.925,
            "14 descriptive ORA q<0.10; global p=0.18439; ORA maxT 0/113",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        parent.text(
            0.762,
            0.044,
            "○ regulator sign agrees; × differs   ◆ both orientation gates",
            transform=parent.transAxes,
            ha="center",
            va="bottom",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        parent.text(
            0.762,
            0.022,
            "T: −5.28 gained-oriented   0 null centre   +5.28 lost-oriented",
            transform=parent.transAxes,
            ha="center",
            va="bottom",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
        parent.text(
            0.762,
            0.000,
            "ORA is directionless; tiles do not show pathway activation.",
            transform=parent.transAxes,
            ha="center",
            va="top",
            fontsize=FONT_MIN_PT,
            color=MUTED,
        )
    
        antigen = plot.loc[plot["kegg_id"].eq("hsa04612")]
        assert len(antigen) == 1
        assert bool(antigen.iloc[0]["both_single_margin_maxT_supported"])
        assert int(plot["both_single_margin_maxT_supported"].sum()) == 1
        assert int(plot["regulator_sign_relation"].eq("disagrees").sum()) == 1
        assert coverage["descriptive_family_reversal_targets_total"] == 2_577
        assert coverage["descriptive_family_reversal_dyads_total"] == 4_379
        return {
            "analysis_scope": "full census; not the 279-target result",
            "reversal_targets": 2_577,
            "reversal_dyads": 4_379,
            "displayed_terms": 14,
            "selection": "conventional Fisher-BH q<0.10 within the fixed 113-term family",
            "regulator_opportunity_global_p": KEGG.ORA_GLOBAL_P,
            "ora_maxT_supported_terms": 0,
            "orientation_both_single_margin_maxT_terms": 1,
            "orientation_supported_term": "Antigen processing and presentation",
            "tile": "target-conditioned reversal-orientation T",
            "regulator_relation_marker": "circle agrees; cross disagrees",
            "claim_boundary": (
                "ORA is directionless and tiles describe reversal orientation; neither "
                "establishes pathway activation, suppression, movement, or rewiring"
            ),
        }
    
    
    def clamp_visible_fonts(fig: plt.Figure) -> dict[str, object]:
        changed_low: list[str] = []
        changed_high: list[str] = []
        # A draw can lazily create tick-label Text objects, so clamp after drawing
        # and repeat until no new out-of-contract text appears.
        for _ in range(3):
            fig.canvas.draw()
            changed_this_pass = 0
            for artist in fig.findobj(match=Text):
                if not artist.get_visible() or not artist.get_text().strip():
                    continue
                size = float(artist.get_fontsize())
                if size < FONT_MIN_PT:
                    artist.set_fontsize(FONT_MIN_PT)
                    changed_low.append(artist.get_text().replace("\n", " ")[:80])
                    changed_this_pass += 1
                elif size > FONT_MAX_PT:
                    artist.set_fontsize(FONT_MAX_PT)
                    changed_high.append(artist.get_text().replace("\n", " ")[:80])
                    changed_this_pass += 1
            if changed_this_pass == 0:
                break
        fig.canvas.draw()
        sizes = [
            float(artist.get_fontsize())
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        if not sizes:
            raise AssertionError("candidate contains no visible text")
        if min(sizes) < FONT_MIN_PT or max(sizes) > FONT_MAX_PT:
            raise AssertionError(f"font range outside 5--7 pt: {min(sizes)}--{max(sizes)}")
        return {
            "minimum_visible_font_pt": min(sizes),
            "maximum_visible_font_pt": max(sizes),
            "raised_to_5_pt": len(changed_low),
            "lowered_to_7_pt": len(changed_high),
        }
    
    
    def axes_geometry(axes: list[mpl.axes.Axes]) -> list[dict[str, object]]:
        records = []
        for index, axis in enumerate(axes):
            bounds = tuple(float(value) for value in axis.get_position().bounds)
            records.append(
                {
                    "retained_index": index,
                    "title": axis.get_title(),
                    "xlabel": axis.get_xlabel(),
                    "ylabel": axis.get_ylabel(),
                    "bounds": list(bounds),
                }
            )
        return records
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, active_audit = ACTIVE._build()
        if not np.allclose(fig.get_size_inches(), (8.27, 11.69)):
            raise AssertionError(f"unexpected active canvas: {fig.get_size_inches()}")
        if len(fig.axes) != 10:
            raise AssertionError(f"unexpected active axes count: {len(fig.axes)}")
    
        old_bar = fig.axes[1]
        old_meta = fig.axes[2]
        if old_bar.get_xlabel() != "reversing target edges (count)":
            raise AssertionError("active panel-b bar axes not found")
        header_matches = [
            artist
            for artist in old_meta.texts
            if artist.get_text() == "rev/test; rate; fold; L/G"
        ]
        if len(header_matches) != 1:
            raise AssertionError("active panel-b metadata axes not found")
    
        left = float(old_bar.get_position().x0)
        bottom = float(old_bar.get_position().y0)
        right = float(old_meta.get_position().x1)
        top = float(old_bar.get_position().y1)
        panel_b_bounds = [left, bottom, right - left, top - bottom]
        expected_bounds = [
            0.46253275109170305,
            0.5083448275862068,
            0.5174672489082972,
            0.45365517241379316,
        ]
        if not np.allclose(panel_b_bounds, expected_bounds, atol=1e-12):
            raise AssertionError(f"active panel-b geometry drifted: {panel_b_bounds}")
    
        removed_axes = {old_bar, old_meta, *list(old_bar.child_axes)}
        retained_axes = [axis for axis in fig.axes if axis not in removed_axes]
        geometry_before = [tuple(axis.get_position().bounds) for axis in retained_axes]
        retained_audit = axes_geometry(retained_axes)
    
        relettering = reletter_lower_panels(fig)
        for child in list(old_bar.child_axes):
            child.remove()
        old_bar.remove()
        old_meta.remove()
    
        removed_figure_text: list[str] = []
        for artist in list(fig.texts):
            text = artist.get_text()
            if text == "b" or text.startswith("All 49 regulators above degree-matched expectation"):
                removed_figure_text.append(text)
                artist.remove()
        if len(removed_figure_text) != 2:
            raise AssertionError(
                f"expected to remove active b tag and heading; removed {removed_figure_text}"
            )
    
        parent = fig.add_axes(panel_b_bounds)
        parent.set_axis_off()
        panel_b = draw_regulator_panel(parent, regulators)
        panel_d = draw_kegg_panel(parent, kegg_plot, coverage)
    
        geometry_after = [tuple(axis.get_position().bounds) for axis in retained_axes]
        for before, after in zip(geometry_before, geometry_after):
            if not np.allclose(before, after, atol=1e-14):
                raise AssertionError(f"retained scientific-panel geometry changed: {before} -> {after}")
    
        font_audit = clamp_visible_fonts(fig)
        visible_text = "\n".join(
            artist.get_text()
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        )
        visible_text_normalized = " ".join(visible_text.split())
        required_text = (
            "49 reversal-enriched regulators",
            "Full census: 2,577 targets / 4,379 dyads — not T279",
            "14 descriptive ORA q<0.10; global p=0.18439; ORA maxT 0/113",
            "Antigen processing and presentation",
        )
        for text in required_text:
            if " ".join(text.split()) not in visible_text_normalized:
                raise AssertionError(f"required visible scope text missing: {text}")
    
        panel_tags = {
            letter: sum(
                1
                for artist in fig.findobj(match=Text)
                if artist.get_visible()
                and artist.get_text() == letter
                and artist.get_fontweight() == "bold"
            )
            for letter in "abcdefghi"
        }
        if panel_tags != {letter: 1 for letter in "abcdefghi"}:
            raise AssertionError(f"panel-tag audit failed: {panel_tags}")
    
        return fig, {
            "active_predecessor_audit": active_audit,
            "canvas_inches": [8.27, 11.69],
            "replaced_rectangle": panel_b_bounds,
            "retained_axes_geometry": retained_audit,
            "retained_axes_geometry_unchanged": True,
            "relettering": relettering,
            "panel_tags": panel_tags,
            "font_audit": font_audit,
            "panel_b": panel_b,
            "panel_d": panel_d,
        }
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "--force",
            action="store_true",
            help="replace only this diagnostic bundle's generated files",
        )
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(
                f"refusing to replace existing diagnostic outputs without --force: {existing}"
            )
    
        canonical_before = canonical_hashes()
        kegg_source_records = verify_kegg_scientific_sources()
        regulators = load_regulators()
        kegg_plot, coverage = KEGG.load_plot_data()
    
        assert len(kegg_plot) == 14
        assert int(kegg_plot["conventional_BH_descriptive"].sum()) == 14
        assert int(kegg_plot["regulator_maxT_term"].sum()) == 0
        # `regulator_maxT_p` in this merged display table is the signed-orientation
        # regulator-margin result.  The separate ORA opportunity gate is bound by
        # KEGG.ORA_GLOBAL_P and the zero `regulator_maxT_term` indicators above.
        assert np.isclose(KEGG.ORA_GLOBAL_P, 0.18439078046097696)
    
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        atomic_tsv(regulators, OUT_REGULATORS)
        atomic_tsv(kegg_plot, OUT_KEGG)
        atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG side-by-side diagnostic candidate",
                "Subject": "Noncanonical visual candidate; full-census signed KEGG is not T279",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        atomic_figure(
            fig,
            OUT_PNG,
            format="png",
            dpi=PNG_DPI,
            facecolor="white",
        )
        plt.close(fig)
    
        canonical_after = canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed during diagnostic render")
    
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v1",
            "status": "NONCANONICAL_DIAGNOSTIC_VISUAL_CANDIDATE",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "panel_order": ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
            "replacement": {
                "scope": "only the active top-right panel-b rectangle",
                "left": "b, one continuous ranks 1-49 lost/gained bar list",
                "right": "d, compact full-census 14-row signed-KEGG list",
                "unchanged_geometry": "a, c and all lower scientific axes",
                "relettered_only": "former d-h are e-i",
                "pdf_embedding_or_cropping": False,
            },
            "font_contract": {
                "minimum_pt": FONT_MIN_PT,
                "maximum_pt": FONT_MAX_PT,
                **audit["font_audit"],
            },
            "scientific_scope": {
                "signed_kegg_population": "all 2,577 full-census reversal targets",
                "signed_kegg_reversal_dyads": 4_379,
                "not_population": "the separate 279 recurrent direction-unanimous targets",
                "displayed_descriptive_terms": 14,
                "fixed_family_terms": 113,
                "regulator_opportunity_global_p": KEGG.ORA_GLOBAL_P,
                "ora_maxT_supported_terms": 0,
                "orientation_both_single_margin_maxT_terms": 1,
            },
            "audit": audit,
            "canonical_figure_2": {
                "before": canonical_before,
                "after": canonical_after,
                "written": False,
            },
            "renderer_closure": {
                path: sha256(PAPER / path) for path in RENDERER_PATHS
            },
            "scientific_source_closure": {
                relative(path): file_record(path) for path in SCIENTIFIC_SOURCE_PATHS
            },
            "verified_kegg_scientific_sources": kegg_source_records,
            "bundle_sources": {
                "panel_b_regulators": file_record(OUT_REGULATORS),
                "panel_d_signed_kegg": file_record(OUT_KEGG),
            },
            "outputs": {
                "pdf": file_record(OUT_PDF),
                "png": file_record(OUT_PNG),
            },
            "claim_boundary": (
                "Panel d is a full-census descriptive target-pathway summary. Conventional "
                "ORA is directionless; tile colour describes reversal orientation. The panel "
                "is not the T279 result and does not establish pathway activation, suppression, "
                "movement, coordinated regulator disruption, or rewiring."
            ),
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }
        atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "pdf": file_record(OUT_PDF),
                    "png": file_record(OUT_PNG),
                    "manifest": relative(OUT_MANIFEST),
                    "font_range_pt": [
                        audit["font_audit"]["minimum_visible_font_pt"],
                        audit["font_audit"]["maximum_visible_font_pt"],
                    ],
                    "canonical_unchanged": canonical_before == canonical_after,
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v1, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v1()

def init_fig2_signed_kegg_side_by_side_candidate_v2():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the frozen v2 semantic delta of the Figure 2 signed-KEGG candidate.
    
    Version 1 is reconstructed in memory and preserved byte-for-byte.  Version 2
    changes only the auditor-requested presentation semantics: the panel-b title,
    two panel-d headers, and a sealed panel-d note disclosing fixed-family
    coverage and the exact antigen-processing orientation evidence.  The panel-d
    table is shifted upward within its existing half-panel only to make room for
    the note.  No canonical, manuscript, legend, guard, or v1 path is written.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V1 = fig2_signed_kegg_side_by_side_candidate_v1
    
    
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v2"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v2"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    
    V1_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v1.py",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v1"
        / "Fig2_signed_kegg_side_by_side_candidate_v1.pdf",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v1"
        / "Fig2_signed_kegg_side_by_side_candidate_v1.png",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v1"
        / "source_panel_b_regulators.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v1"
        / "source_panel_d_signed_kegg.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v1"
        / "manifest.json",
    )
    EXPECTED_V1_SHA256 = {
        V1_PATHS[0]: "f43e58b9523aa619846f7d285f937732bf6560008b17f0c678ea66452a298d5c",
        V1_PATHS[1]: "be1d8acdbea16f2c6e77fbb398324c7adfcd13e75499d73e1015bb89794c12cf",
        V1_PATHS[2]: "80daa5cbfbdccac502f5b5e1c4a7300af836be2dcf3156a9c80108ded060b570",
        V1_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V1_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V1_PATHS[5]: "4f4c818d88a4235eac0c0a01043e09302169169ac8f8d02210741aec0d8ea244",
    }
    
    PARENT_BOUNDS = np.array(
        [
            0.46253275109170305,
            0.5083448275862068,
            0.5174672489082972,
            0.45365517241379316,
        ],
        dtype=float,
    )
    NOTE_TEXT = (
        "Fixed-family coverage: 796/2,577 targets;\n"
        "1,443/4,379 dyads; pathway memberships overlap.\n"
        "◆ Antigen processing: target-conditioned\n"
        "maxT p=.00344983; regulator-margin\n"
        "maxT p=.0000499975.\n"
        "PSME1+TAP1=79.8% of signed score;\n"
        "complementary single-margin tests,\n"
        "not a joint test.\n"
        "○ regulator sign agrees; × differs\n"
        "T: −5.28 gained · 0 · +5.28 lost\n"
        "ORA directionless; no pathway activation."
    )
    
    
    def verify_v1() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V1_SHA256.items():
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"frozen v1 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def find_unique_axis_with_text(fig: plt.Figure, value: str) -> mpl.axes.Axes:
        matches = [
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if any(artist.get_text() == value for artist in axis.texts)
        ]
        if len(matches) != 1:
            raise AssertionError(f"expected one axes containing {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def find_parent_axis(fig: plt.Figure) -> mpl.axes.Axes:
        matches = [
            axis
            for axis in fig.axes
            if np.allclose(axis.get_position().bounds, PARENT_BOUNDS, atol=1e-12)
        ]
        if len(matches) != 1:
            raise AssertionError(f"expected one v1 replacement parent; found {len(matches)}")
        return matches[0]
    
    
    def replace_unique_text(
        artists: list[Text],
        old: str,
        new: str,
    ) -> Text:
        matches = [artist for artist in artists if artist.get_text() == old]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {old!r}; found {len(matches)}")
        matches[0].set_text(new)
        return matches[0]
    
    
    def row_text_snapshot(
        axis: mpl.axes.Axes,
    ) -> list[tuple[str, tuple[float, float]]]:
        snapshot: list[tuple[str, tuple[float, float]]] = []
        for artist in axis.texts:
            x, y = artist.get_position()
            if float(y) >= -0.5:
                snapshot.append((artist.get_text(), (float(x), float(y))))
        return snapshot
    
    
    def apply_v2_delta(
        fig: plt.Figure,
        v1_audit: dict[str, object],
    ) -> dict[str, object]:
        parent = find_parent_axis(fig)
        kegg_axis = find_unique_axis_with_text(fig, "KEGG pathway")
        retained_axes = [axis for axis in fig.findobj(match=mpl.axes.Axes) if axis is not kegg_axis]
        retained_before = [tuple(axis.get_position().bounds) for axis in retained_axes]
        kegg_before = tuple(float(value) for value in kegg_axis.get_position().bounds)
        row_text_before = row_text_snapshot(kegg_axis)
    
        replace_unique_text(
            list(parent.texts),
            "49 reversal-enriched regulators",
            "49 above degree-matched expectation",
        )
        header_lg = replace_unique_text(list(kegg_axis.texts), "L/G", "L/G\ndyads")
        header_lg.set_linespacing(0.92)
        header_q = replace_unique_text(
            list(kegg_axis.texts),
            "BH q",
            "target ORA\nq",
        )
        header_q.set_linespacing(0.92)
    
        # The exact table values and ordering are untouched.  Only its internal
        # position changes, within panel d, to open a sealed explanatory strip.
        parent_left, parent_bottom, parent_width, parent_height = PARENT_BOUNDS
        kegg_after = (
            float(parent_left + 0.525 * parent_width),
            float(parent_bottom + 0.205 * parent_height),
            float(0.475 * parent_width),
            float(0.695 * parent_height),
        )
        kegg_axis.set_axes_locator(None)
        kegg_axis.set_position(kegg_after)
    
        old_footer = (
            "○ regulator sign agrees; × differs   ◆ both orientation gates",
            "T: −5.28 gained-oriented   0 null centre   +5.28 lost-oriented",
            "ORA is directionless; tiles do not show pathway activation.",
        )
        removed: list[str] = []
        for text in old_footer:
            matches = [artist for artist in list(parent.texts) if artist.get_text() == text]
            if len(matches) != 1:
                raise AssertionError(f"expected one v1 footer line {text!r}; found {len(matches)}")
            removed.append(text)
            matches[0].remove()
    
        note_box = Rectangle(
            (0.515, 0.004),
            0.480,
            0.188,
            transform=parent.transAxes,
            facecolor="#f7f8f9",
            edgecolor="#cdd2d6",
            linewidth=0.45,
            clip_on=False,
            zorder=1,
        )
        parent.add_patch(note_box)
        parent.text(
            0.530,
            0.184,
            NOTE_TEXT,
            transform=parent.transAxes,
            ha="left",
            va="top",
            fontsize=V1.FONT_MIN_PT,
            linespacing=1.05,
            color=V1.MUTED,
            zorder=2,
        )
    
        row_text_after = row_text_snapshot(kegg_axis)
        if row_text_after != row_text_before:
            raise AssertionError("panel-d row text or row coordinates changed")
        retained_after = [tuple(axis.get_position().bounds) for axis in retained_axes]
        for before, after in zip(retained_before, retained_after):
            if not np.allclose(before, after, atol=1e-14):
                raise AssertionError(f"non-table axes geometry changed: {before} -> {after}")
    
        font_audit = V1.clamp_visible_fonts(fig)
        visible = " ".join(
            "\n".join(
                artist.get_text()
                for artist in fig.findobj(match=Text)
                if artist.get_visible() and artist.get_text().strip()
            ).split()
        )
        required = (
            "49 above degree-matched expectation",
            "L/G dyads",
            "target ORA q",
            "796/2,577 targets",
            "1,443/4,379 dyads",
            "pathway memberships overlap",
            "target-conditioned maxT p=.00344983",
            "regulator-margin maxT p=.0000499975",
            "PSME1+TAP1=79.8% of signed score",
            "complementary single-margin tests",
            "not a joint test",
        )
        for value in required:
            if " ".join(value.split()) not in visible:
                raise AssertionError(f"required v2 disclosure missing: {value}")
        if "49 reversal-enriched regulators" in visible:
            raise AssertionError("v1 panel-b title survived in v2")
    
        if v1_audit["panel_tags"] != {letter: 1 for letter in "abcdefghi"}:
            raise AssertionError("v1 panel-tag audit drifted")
        panel_tags = {
            letter: sum(
                1
                for artist in fig.findobj(match=Text)
                if artist.get_visible()
                and artist.get_text() == letter
                and artist.get_fontweight() == "bold"
            )
            for letter in "abcdefghi"
        }
        if panel_tags != {letter: 1 for letter in "abcdefghi"}:
            raise AssertionError(f"v2 panel tags drifted: {panel_tags}")
    
        return {
            "predecessor_v1_audit": v1_audit,
            "delta_only": True,
            "panel_b_title": {
                "before": "49 reversal-enriched regulators",
                "after": "49 above degree-matched expectation",
            },
            "panel_d_headers": {
                "before": ["L/G", "BH q"],
                "after": ["L/G dyads", "target ORA q"],
            },
            "panel_d_table": {
                "row_text_and_coordinates_unchanged": True,
                "order_unchanged": True,
                "values_unchanged": True,
                "bounds_before": list(kegg_before),
                "bounds_after": list(kegg_after),
                "geometry_change_reason": "essential room for the sealed disclosure strip",
            },
            "sealed_note": {
                "text": NOTE_TEXT,
                "bounds_within_parent": [0.515, 0.004, 0.480, 0.188],
                "removed_v1_footer_lines": removed,
            },
            "all_other_axes_geometry_unchanged": True,
            "panel_tags": panel_tags,
            "font_audit": font_audit,
        }
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, v1_audit = V1.build_figure(regulators, kegg_plot, coverage)
        delta_audit = apply_v2_delta(fig, v1_audit)
        return fig, delta_audit
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "--force",
            action="store_true",
            help="replace only this v2 diagnostic bundle's generated files",
        )
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(
                f"refusing to replace existing v2 outputs without --force: {existing}"
            )
    
        canonical_before = V1.canonical_hashes()
        v1_before = verify_v1()
        kegg_source_records = V1.verify_kegg_scientific_sources()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        assert len(regulators) == 49
        assert len(kegg_plot) == 14
        assert int(kegg_plot["n_lost_dyads"].sum()) == 495
        assert int(kegg_plot["n_gained_dyads"].sum()) == 504
    
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG side-by-side diagnostic candidate v2",
                "Subject": "Frozen semantic delta; noncanonical and not T279",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        V1.atomic_figure(
            fig,
            OUT_PNG,
            format="png",
            dpi=V1.PNG_DPI,
            facecolor="white",
        )
        plt.close(fig)
    
        canonical_after = V1.canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed during v2 render")
        v1_after = verify_v1()
        if v1_after != v1_before:
            raise AssertionError("frozen v1 artifact changed during v2 render")
    
        renderer_paths = list(
            dict.fromkeys(
                [
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                ]
            )
        )
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v2",
            "status": "FROZEN_NONCANONICAL_DIAGNOSTIC_CANDIDATE",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "predecessor": {
                "version": "v1",
                "preserved": True,
                "before": v1_before,
                "after": v1_after,
            },
            "canonical_figure_2": {
                "before": canonical_before,
                "after": canonical_after,
                "written": False,
            },
            "delta": {
                "panel_b_title": "49 above degree-matched expectation",
                "panel_d_headers": ["L/G dyads", "target ORA q"],
                "fixed_family_coverage": {
                    "targets": "796/2,577",
                    "dyads": "1,443/4,379",
                    "pathway_memberships_overlap": True,
                },
                "antigen_processing_note": {
                    "target_conditioned_maxT_p": 0.00344983,
                    "regulator_margin_maxT_p": 0.0000499975,
                    "psme1_tap1_signed_score_fraction": 0.798,
                    "tests": "complementary single-margin tests, not a joint test",
                },
            },
            "scientific_scope": {
                "population": "full reversal census, 2,577 targets and 4,379 dyads",
                "not_population": "the separate T279 result",
                "fixed_kegg_family": 113,
                "displayed_descriptive_ora_q_lt_0_10": 14,
                "regulator_opportunity_global_p": V1.KEGG.ORA_GLOBAL_P,
                "ora_maxT_supported_terms": 0,
                "pathway_memberships_overlap": True,
            },
            "audit": audit,
            "font_contract": {
                "minimum_pt": V1.FONT_MIN_PT,
                "maximum_pt": V1.FONT_MAX_PT,
                **audit["font_audit"],
            },
            "renderer_closure": {
                path: V1.sha256(PAPER / path) for path in renderer_paths
            },
            "scientific_source_closure": {
                V1.relative(path): V1.file_record(path)
                for path in V1.SCIENTIFIC_SOURCE_PATHS
            },
            "verified_kegg_scientific_sources": kegg_source_records,
            "bundle_sources": {
                "panel_b_regulators": V1.file_record(OUT_REGULATORS),
                "panel_d_signed_kegg": V1.file_record(OUT_KEGG),
            },
            "outputs": {
                "pdf": V1.file_record(OUT_PDF),
                "png": V1.file_record(OUT_PNG),
            },
            "claim_boundary": (
                "Panel d is a partial fixed-family annotation of the full reversal "
                "census, with overlapping memberships. ORA is directionless. The two "
                "antigen-processing maxT values come from complementary single-margin "
                "tests and are not a joint test. The panel is not the T279 analysis and "
                "does not establish pathway activation or coordinated rewiring."
            ),
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "pdf": V1.file_record(OUT_PDF),
                    "png": V1.file_record(OUT_PNG),
                    "manifest": V1.relative(OUT_MANIFEST),
                    "font_range_pt": [
                        audit["font_audit"]["minimum_visible_font_pt"],
                        audit["font_audit"]["maximum_visible_font_pt"],
                    ],
                    "canonical_unchanged": canonical_before == canonical_after,
                    "v1_preserved": v1_before == v1_after,
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v2, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v2()

def init_fig2_signed_kegg_side_by_side_candidate_v3():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the frozen typography-only v3 delta of the Figure 2 candidate.
    
    Version 2 is reconstructed in memory and preserved byte-for-byte.  Version 3
    makes only the two corrections requested by the independent PDF-bbox audit:
    move panel tag b left of its longer title, and remove the redundant ``markers``
    header whose meaning is already stated in the sealed panel-d note.  No data,
    row order, scientific text, axes geometry, canonical, manuscript, legend,
    guard, v1, or v2 path is written.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V2 = fig2_signed_kegg_side_by_side_candidate_v2
    
    
    V1 = V2.V1
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v3"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v3"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    
    V2_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v2"
        / "Fig2_signed_kegg_side_by_side_candidate_v2.pdf",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v2"
        / "Fig2_signed_kegg_side_by_side_candidate_v2.png",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v2"
        / "source_panel_b_regulators.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v2"
        / "source_panel_d_signed_kegg.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v2"
        / "manifest.json",
    )
    EXPECTED_V2_SHA256 = {
        V2_PATHS[0]: "84b1cfdca5e22176f4f53313a96018342c265ab0441945e7e68ddf47bb032afe",
        V2_PATHS[1]: "315a6a8feaa595025b0571d9588d1d40d6f2c0f88951c2345dd935a731750093",
        V2_PATHS[2]: "b7479a08dfb43032fb8e933d78a3a1e89cf8feb71c0589c589eaa3e39c393978",
        V2_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V2_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V2_PATHS[5]: "2b4950c64a874d292c45ded4c13be67e9ced6d561cb6b8c10cec8ca7f58314f4",
    }
    
    
    def verify_v2() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V2_SHA256.items():
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"frozen v2 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def unique_text(artists: list[Text], value: str) -> Text:
        matches = [artist for artist in artists if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def axes_snapshot(fig: plt.Figure) -> list[tuple[float, float, float, float]]:
        return [
            tuple(float(value) for value in axis.get_position().bounds)
            for axis in fig.findobj(match=mpl.axes.Axes)
        ]
    
    
    def apply_v3_delta(
        fig: plt.Figure,
        v2_audit: dict[str, object],
    ) -> dict[str, object]:
        parent = V2.find_parent_axis(fig)
        kegg_axis = V2.find_unique_axis_with_text(fig, "KEGG pathway")
        axes_before = axes_snapshot(fig)
        rows_before = V2.row_text_snapshot(kegg_axis)
    
        tag_b = unique_text(
            [
                artist
                for artist in list(parent.texts)
                if artist.get_text() == "b" and artist.get_fontweight() == "bold"
            ],
            "b",
        )
        title_b = unique_text(
            list(parent.texts),
            "49 above degree-matched expectation",
        )
        tag_before = tuple(float(value) for value in tag_b.get_position())
        if not np.allclose(tag_before, (0.0, 0.995), atol=1e-14):
            raise AssertionError(f"unexpected v2 panel-b tag position: {tag_before}")
        tag_b.set_x(-0.05)
        tag_after = tuple(float(value) for value in tag_b.get_position())
    
        marker_header = unique_text(list(kegg_axis.texts), "markers")
        marker_header_position = tuple(
            float(value) for value in marker_header.get_position()
        )
        marker_header.remove()
    
        rows_after = V2.row_text_snapshot(kegg_axis)
        if rows_after != rows_before:
            raise AssertionError("panel-d row text or row coordinates changed")
        axes_after = axes_snapshot(fig)
        if len(axes_after) != len(axes_before):
            raise AssertionError("axes inventory changed")
        for before, after in zip(axes_before, axes_after):
            if not np.allclose(before, after, atol=1e-14):
                raise AssertionError(f"axes geometry changed: {before} -> {after}")
    
        font_audit = V1.clamp_visible_fonts(fig)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        tag_bbox = tag_b.get_window_extent(renderer=renderer)
        title_bbox = title_b.get_window_extent(renderer=renderer)
        clearance_pt = float((title_bbox.x0 - tag_bbox.x1) * 72.0 / fig.dpi)
        if clearance_pt < 4.0:
            raise AssertionError(
                f"panel-b tag/title clearance is only {clearance_pt:.3f} pt"
            )
    
        visible = " ".join(
            "\n".join(
                artist.get_text()
                for artist in fig.findobj(match=Text)
                if artist.get_visible() and artist.get_text().strip()
            ).split()
        )
        required = (
            "49 above degree-matched expectation",
            "L/G dyads",
            "target ORA q",
            "796/2,577 targets",
            "1,443/4,379 dyads",
            "pathway memberships overlap",
            "target-conditioned maxT p=.00344983",
            "regulator-margin maxT p=.0000499975",
            "PSME1+TAP1=79.8% of signed score",
            "complementary single-margin tests",
            "not a joint test",
            "regulator sign agrees; × differs",
        )
        for value in required:
            if " ".join(value.split()) not in visible:
                raise AssertionError(f"required v3 text missing: {value}")
        if any(
            artist.get_visible() and artist.get_text() == "markers"
            for artist in fig.findobj(match=Text)
        ):
            raise AssertionError("redundant markers header survived in v3")
    
        panel_tags = {
            letter: sum(
                1
                for artist in fig.findobj(match=Text)
                if artist.get_visible()
                and artist.get_text() == letter
                and artist.get_fontweight() == "bold"
            )
            for letter in "abcdefghi"
        }
        if panel_tags != {letter: 1 for letter in "abcdefghi"}:
            raise AssertionError(f"v3 panel tags drifted: {panel_tags}")
    
        return {
            "predecessor_v2_audit": v2_audit,
            "delta_only": True,
            "panel_b_tag": {
                "position_before_parent_axes": list(tag_before),
                "position_after_parent_axes": list(tag_after),
                "agg_bbox_clearance_pt": clearance_pt,
                "minimum_required_clearance_pt": 4.0,
            },
            "panel_d_header": {
                "removed": "markers",
                "removed_position_data": list(marker_header_position),
                "meaning_retained_in_sealed_note": True,
                "required_headers_retained": ["L/G dyads", "target ORA q"],
            },
            "panel_d_rows_unchanged": True,
            "all_axes_geometry_unchanged": True,
            "panel_tags": panel_tags,
            "font_audit": font_audit,
        }
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, v2_audit = V2.build_figure(regulators, kegg_plot, coverage)
        audit = apply_v3_delta(fig, v2_audit)
        return fig, audit
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "--force",
            action="store_true",
            help="replace only this v3 diagnostic bundle's generated files",
        )
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(
                f"refusing to replace existing v3 outputs without --force: {existing}"
            )
    
        canonical_before = V1.canonical_hashes()
        v1_before = V2.verify_v1()
        v2_before = verify_v2()
        kegg_source_records = V1.verify_kegg_scientific_sources()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        assert len(regulators) == 49
        assert len(kegg_plot) == 14
        assert int(kegg_plot["n_lost_dyads"].sum()) == 495
        assert int(kegg_plot["n_gained_dyads"].sum()) == 504
    
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG side-by-side diagnostic candidate v3",
                "Subject": "Typography-only delta; noncanonical and not T279",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        V1.atomic_figure(
            fig,
            OUT_PNG,
            format="png",
            dpi=V1.PNG_DPI,
            facecolor="white",
        )
        plt.close(fig)
    
        canonical_after = V1.canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed during v3 render")
        v1_after = V2.verify_v1()
        if v1_after != v1_before:
            raise AssertionError("frozen v1 artifact changed during v3 render")
        v2_after = verify_v2()
        if v2_after != v2_before:
            raise AssertionError("frozen v2 artifact changed during v3 render")
    
        renderer_paths = list(
            dict.fromkeys(
                [
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                ]
            )
        )
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v3",
            "status": "FROZEN_NONCANONICAL_DIAGNOSTIC_CANDIDATE",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "predecessor": {
                "version": "v2",
                "preserved": True,
                "before": v2_before,
                "after": v2_after,
            },
            "v1_preserved": {
                "before": v1_before,
                "after": v1_after,
            },
            "canonical_figure_2": {
                "before": canonical_before,
                "after": canonical_after,
                "written": False,
            },
            "delta": {
                "panel_b_tag_parent_x": {"before": 0.0, "after": -0.05},
                "panel_d_header_removed": "markers",
                "scientific_or_data_changes": False,
            },
            "scientific_scope": {
                "population": "full reversal census, 2,577 targets and 4,379 dyads",
                "not_population": "the separate T279 result",
                "fixed_kegg_family": 113,
                "displayed_descriptive_ora_q_lt_0_10": 14,
                "regulator_opportunity_global_p": V1.KEGG.ORA_GLOBAL_P,
                "ora_maxT_supported_terms": 0,
                "fixed_family_coverage_targets": "796/2,577",
                "fixed_family_coverage_dyads": "1,443/4,379",
                "pathway_memberships_overlap": True,
            },
            "audit": audit,
            "font_contract": {
                "minimum_pt": V1.FONT_MIN_PT,
                "maximum_pt": V1.FONT_MAX_PT,
                **audit["font_audit"],
            },
            "renderer_closure": {
                path: V1.sha256(PAPER / path) for path in renderer_paths
            },
            "scientific_source_closure": {
                V1.relative(path): V1.file_record(path)
                for path in V1.SCIENTIFIC_SOURCE_PATHS
            },
            "verified_kegg_scientific_sources": kegg_source_records,
            "bundle_sources": {
                "panel_b_regulators": V1.file_record(OUT_REGULATORS),
                "panel_d_signed_kegg": V1.file_record(OUT_KEGG),
            },
            "outputs": {
                "pdf": V1.file_record(OUT_PDF),
                "png": V1.file_record(OUT_PNG),
            },
            "claim_boundary": (
                "Typography-only v3 of the full-census, partial fixed-family "
                "annotation. ORA is directionless; pathway memberships overlap; "
                "the two antigen-processing tests are complementary single-margin "
                "tests, not a joint test; this is not the T279 analysis."
            ),
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "pdf": V1.file_record(OUT_PDF),
                    "png": V1.file_record(OUT_PNG),
                    "manifest": V1.relative(OUT_MANIFEST),
                    "font_range_pt": [
                        audit["font_audit"]["minimum_visible_font_pt"],
                        audit["font_audit"]["maximum_visible_font_pt"],
                    ],
                    "panel_b_tag_title_clearance_pt": audit["panel_b_tag"][
                        "agg_bbox_clearance_pt"
                    ],
                    "canonical_unchanged": canonical_before == canonical_after,
                    "v1_preserved": v1_before == v1_after,
                    "v2_preserved": v2_before == v2_after,
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v3, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v3()

def init_fig2_signed_kegg_side_by_side_candidate_v4():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render Figure 2 signed-KEGG candidate v4 with complete regulator metadata.
    
    This is a narrow delta on the audited v3 member of the existing
    ``fig2_signed_kegg_side_by_side_candidate_v2`` lineage.  The full-census KEGG
    panel and every other scientific panel are reconstructed unchanged.  Panel b
    retains all 49 lost/gained bars and adds, for every displayed regulator, the
    requested ``rev/test; rate; fold; L/G`` fields.  Here ``fold`` is the observed
    reversal rate divided by the regulator's out-degree-decile expectation.
    
    Only a noncanonical diagnostic bundle is written.  Canonical Figure 2,
    manuscript, legends, and all predecessor bundles are read-only guards.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V3 = fig2_signed_kegg_side_by_side_candidate_v3
    
    
    V2 = V3.V2
    V1 = V2.V1
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    V3_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v3"
        / "Fig2_signed_kegg_side_by_side_candidate_v3.pdf",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v3"
        / "Fig2_signed_kegg_side_by_side_candidate_v3.png",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v3"
        / "source_panel_b_regulators.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v3"
        / "source_panel_d_signed_kegg.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v3"
        / "manifest.json",
    )
    EXPECTED_V3_SHA256 = {
        V3_PATHS[0]: "76e0240840aaa0b512bb1c63c13943fdfcce96976803442de9b0a599459e2378",
        V3_PATHS[1]: "bc06e56ea7f81f4fa0e600bb2ecef17cedcd9c910800593bb7c5d95b280b4daf",
        V3_PATHS[2]: "3c9a2fc305049665178a0116f86d2b36fdaf6ba1d200ff568a3c07cf3e74bd8a",
        V3_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V3_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V3_PATHS[5]: "b395df7fa223279d2e3b0a718d6af6c28a33554e059fa4697c948d3f3b4f1aa4",
    }
    
    PANEL_B_BAR_RELATIVE_BOUNDS = (0.015, 0.075, 0.190, 0.825)
    PANEL_B_META_RELATIVE_BOUNDS = (0.215, 0.075, 0.370, 0.825)
    PANEL_B_TITLE_CENTER_X = 0.300
    PANEL_D_TABLE_RELATIVE_BOUNDS = (0.610, 0.205, 0.390, 0.695)
    PANEL_D_TITLE_CENTER_X = 0.805
    PANEL_D_NOTE_RELATIVE_BOUNDS = (0.605, 0.004, 0.390, 0.188)
    
    
    def current_canonical_hashes() -> dict[str, str]:
        """Snapshot the current canonical bytes without invoking stale hash pins."""
        return {
            "pdf": V1.sha256(V1.CANONICAL_PDF),
            "png": V1.sha256(V1.CANONICAL_PNG),
        }
    
    
    def verify_v3() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V3_SHA256.items():
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"frozen v3 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def unique_text(artists: list[Text], value: str) -> Text:
        matches = [artist for artist in artists if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def relative_to_figure(
        parent: mpl.axes.Axes,
        bounds: tuple[float, float, float, float],
    ) -> tuple[float, float, float, float]:
        parent_left, parent_bottom, parent_width, parent_height = (
            float(value) for value in parent.get_position().bounds
        )
        left, bottom, width, height = bounds
        return (
            parent_left + left * parent_width,
            parent_bottom + bottom * parent_height,
            width * parent_width,
            height * parent_height,
        )
    
    
    def all_axes_snapshot(
        fig: plt.Figure,
    ) -> dict[int, tuple[float, float, float, float]]:
        return {
            id(axis): tuple(float(value) for value in axis.get_position().bounds)
            for axis in fig.findobj(match=mpl.axes.Axes)
        }
    
    
    def metadata_bbox_audit(
        fig: plt.Figure,
        row_artists: list[list[Text]],
        metadata_axis: mpl.axes.Axes,
    ) -> dict[str, object]:
        """Confirm that every numeric field is visible and horizontally separated."""
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        axis_bbox = metadata_axis.get_window_extent(renderer=renderer)
        minimum_gap_pt = float("inf")
        minimum_left_margin_pt = float("inf")
        minimum_right_margin_pt = float("inf")
        for artists in row_artists:
            boxes = [artist.get_window_extent(renderer=renderer) for artist in artists]
            for left_box, right_box in zip(boxes[:-1], boxes[1:]):
                gap_pt = float((right_box.x0 - left_box.x1) * 72.0 / fig.dpi)
                minimum_gap_pt = min(minimum_gap_pt, gap_pt)
            minimum_left_margin_pt = min(
                minimum_left_margin_pt,
                float((boxes[0].x0 - axis_bbox.x0) * 72.0 / fig.dpi),
            )
            minimum_right_margin_pt = min(
                minimum_right_margin_pt,
                float((axis_bbox.x1 - boxes[-1].x1) * 72.0 / fig.dpi),
            )
        if minimum_gap_pt < 0.20:
            raise AssertionError(
                f"panel-b numeric fields overlap; minimum gap {minimum_gap_pt:.3f} pt"
            )
        if minimum_left_margin_pt < -0.20 or minimum_right_margin_pt < -0.20:
            raise AssertionError(
                "panel-b metadata extends outside its axes: "
                f"left {minimum_left_margin_pt:.3f} pt, "
                f"right {minimum_right_margin_pt:.3f} pt"
            )
        return {
            "rows_checked": len(row_artists),
            "minimum_horizontal_field_gap_pt": minimum_gap_pt,
            "minimum_left_margin_pt": minimum_left_margin_pt,
            "minimum_right_margin_pt": minimum_right_margin_pt,
        }
    
    
    def apply_v4_delta(
        fig: plt.Figure,
        regulators: pd.DataFrame,
        v3_audit: dict[str, object],
    ) -> dict[str, object]:
        parent = V2.find_parent_axis(fig)
        bar_matches = [
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if axis.get_xlabel() == "reversal edges (linear 0–300)"
        ]
        if len(bar_matches) != 1:
            raise AssertionError(f"expected one panel-b bar axis; found {len(bar_matches)}")
        bar_axis = bar_matches[0]
        metadata_matches = [
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if any(artist.get_text() == "L / G" for artist in axis.texts)
        ]
        if len(metadata_matches) != 1:
            raise AssertionError(
                f"expected one panel-b L/G axis; found {len(metadata_matches)}"
            )
        metadata_axis = metadata_matches[0]
        kegg_axis = V2.find_unique_axis_with_text(fig, "KEGG pathway")
    
        axes_before = all_axes_snapshot(fig)
        kegg_rows_before = V2.row_text_snapshot(kegg_axis)
        kegg_bounds_before = tuple(float(value) for value in kegg_axis.get_position().bounds)
        bar_bounds_before = tuple(float(value) for value in bar_axis.get_position().bounds)
        metadata_bounds_before = tuple(
            float(value) for value in metadata_axis.get_position().bounds
        )
    
        bar_axis.set_axes_locator(None)
        metadata_axis.set_axes_locator(None)
        kegg_axis.set_axes_locator(None)
        bar_bounds_after = relative_to_figure(parent, PANEL_B_BAR_RELATIVE_BOUNDS)
        metadata_bounds_after = relative_to_figure(parent, PANEL_B_META_RELATIVE_BOUNDS)
        kegg_bounds_after = relative_to_figure(parent, PANEL_D_TABLE_RELATIVE_BOUNDS)
        bar_axis.set_position(bar_bounds_after)
        metadata_axis.set_position(metadata_bounds_after)
        kegg_axis.set_position(kegg_bounds_after)
        bar_axis.set_xlabel(
            "reversal edges (0–300)",
            fontsize=V1.FONT_MIN_PT,
            labelpad=1.0,
        )
    
        for artist in list(metadata_axis.texts):
            artist.remove()
        metadata_axis.set_xlim(0.0, 1.0)
        metadata_axis.set_ylim(len(regulators) - 0.35, -0.65)
        metadata_axis.axis("off")
    
        header_y = -1.15
        metadata_axis.text(
            0.01,
            header_y,
            "rev/test",
            ha="left",
            va="center",
            fontsize=V1.FONT_MIN_PT,
            fontweight="bold",
            color=V1.MUTED,
        )
        metadata_axis.text(
            0.50,
            header_y,
            "rate",
            ha="right",
            va="center",
            fontsize=V1.FONT_MIN_PT,
            fontweight="bold",
            color=V1.MUTED,
        )
        metadata_axis.text(
            0.71,
            header_y,
            "fold",
            ha="right",
            va="center",
            fontsize=V1.FONT_MIN_PT,
            fontweight="bold",
            color=V1.MUTED,
        )
        metadata_axis.text(
            0.88,
            header_y,
            "L/G",
            ha="center",
            va="center",
            fontsize=V1.FONT_MIN_PT,
            fontweight="bold",
            color=V1.MUTED,
        )
    
        row_artists: list[list[Text]] = []
        row_strings: list[dict[str, str]] = []
        for row_index, (_, row) in enumerate(regulators.iterrows()):
            reversal_tested = f"{int(row['n_inv'])}/{int(row['n_both'])}"
            rate = f"{100.0 * float(row['rate']):.1f}%"
            fold = f"{float(row['enrichment_vs_decile']):.2f}×"
            lost = str(int(row["lost_on_activation"]))
            gained = str(int(row["gained_on_activation"]))
            rev_artist = metadata_axis.text(
                0.01,
                row_index,
                reversal_tested,
                ha="left",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.INK,
            )
            rate_artist = metadata_axis.text(
                0.50,
                row_index,
                rate,
                ha="right",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.INK,
            )
            fold_artist = metadata_axis.text(
                0.71,
                row_index,
                fold,
                ha="right",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.INK,
            )
            lost_artist = metadata_axis.text(
                0.865,
                row_index,
                lost,
                ha="right",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.LOST,
            )
            slash_artist = metadata_axis.text(
                0.88,
                row_index,
                "/",
                ha="center",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.MUTED,
            )
            gained_artist = metadata_axis.text(
                0.895,
                row_index,
                gained,
                ha="left",
                va="center",
                fontsize=V1.FONT_MIN_PT,
                color=V1.GAINED,
            )
            row_artists.append(
                [
                    rev_artist,
                    rate_artist,
                    fold_artist,
                    lost_artist,
                    slash_artist,
                    gained_artist,
                ]
            )
            row_strings.append(
                {
                    "gene": str(row["gene"]),
                    "rev_test": reversal_tested,
                    "rate": rate,
                    "fold": fold,
                    "lost_gained": f"{lost}/{gained}",
                }
            )
    
        title = unique_text(list(parent.texts), "49 above degree-matched expectation")
        subtitle = unique_text(
            list(parent.texts),
            "2,075 reversals; Fisher/BH q<0.05; ranks 1–49",
        )
        footer = unique_text(
            list(parent.texts),
            "red lost; teal gained; exact counts shown as L/G",
        )
        title.set_x(PANEL_B_TITLE_CENTER_X)
        subtitle.set_x(PANEL_B_TITLE_CENTER_X)
        footer.set_x(PANEL_B_TITLE_CENTER_X)
        subtitle.set_text("2,075/14,434 reversals/tested; BH q<0.05; ranks 1–49")
        footer.set_text("red lost; teal gained; fold = rate / out-degree-decile expectation")
    
        tag_d = unique_text(
            [
                artist
                for artist in list(parent.texts)
                if artist.get_text() == "d" and artist.get_fontweight() == "bold"
            ],
            "d",
        )
        tag_d.set_x(0.595)
        for value in (
            "Full-census signed KEGG",
            "Full census: 2,577 targets / 4,379 dyads — not T279",
            "14 descriptive ORA q<0.10; global p=0.18439; ORA maxT 0/113",
        ):
            unique_text(list(parent.texts), value).set_x(PANEL_D_TITLE_CENTER_X)
        note_boxes = [
            patch
            for patch in parent.patches
            if isinstance(patch, mpl.patches.Rectangle)
            and np.allclose(patch.get_xy(), (0.515, 0.004), atol=1e-14)
            and np.isclose(float(patch.get_width()), 0.480, atol=1e-14)
            and np.isclose(float(patch.get_height()), 0.188, atol=1e-14)
        ]
        if len(note_boxes) != 1:
            raise AssertionError(f"expected one panel-d note box; found {len(note_boxes)}")
        note_box = note_boxes[0]
        note_box.set_x(PANEL_D_NOTE_RELATIVE_BOUNDS[0])
        note_box.set_y(PANEL_D_NOTE_RELATIVE_BOUNDS[1])
        note_box.set_width(PANEL_D_NOTE_RELATIVE_BOUNDS[2])
        note_box.set_height(PANEL_D_NOTE_RELATIVE_BOUNDS[3])
        note_text = unique_text(list(parent.texts), V2.NOTE_TEXT)
        note_text.set_x(0.615)
    
        font_audit = V1.clamp_visible_fonts(fig)
        bbox_audit = metadata_bbox_audit(fig, row_artists, metadata_axis)
        kegg_rows_after = V2.row_text_snapshot(kegg_axis)
        if kegg_rows_after != kegg_rows_before:
            raise AssertionError("full-census KEGG row text or row coordinates changed")
    
        axes_after = all_axes_snapshot(fig)
        changed_axes = {id(bar_axis), id(metadata_axis), id(kegg_axis)}
        for axis_id, before in axes_before.items():
            if axis_id in changed_axes:
                continue
            after = axes_after[axis_id]
            if not np.allclose(before, after, atol=1e-14):
                raise AssertionError(f"non-panel-b axes geometry changed: {before} -> {after}")
    
        visible = " ".join(
            "\n".join(
                artist.get_text()
                for artist in fig.findobj(match=Text)
                if artist.get_visible() and artist.get_text().strip()
            ).split()
        )
        required = (
            "rev/test rate fold L/G",
            "2,075/14,434 reversals/tested",
            "fold = rate / out-degree-decile expectation",
            "Full census: 2,577 targets / 4,379 dyads — not T279",
            "14 descriptive ORA q<0.10; global p=0.18439; ORA maxT 0/113",
            "target-conditioned maxT p=.00344983",
            "regulator-margin maxT p=.0000499975",
        )
        for value in required:
            if " ".join(value.split()) not in visible:
                raise AssertionError(f"required v4 disclosure missing: {value}")
        for row in row_strings:
            for field in ("rev_test", "rate", "fold"):
                if row[field] not in visible:
                    raise AssertionError(
                        f"panel-b row field absent for {row['gene']}: {field}={row[field]}"
                    )
            lost, gained = row["lost_gained"].split("/")
            displayed_lost_gained = f"{lost} / {gained}"
            if displayed_lost_gained not in visible:
                raise AssertionError(
                    "panel-b L/G field absent for "
                    f"{row['gene']}: {displayed_lost_gained}"
                )
    
        panel_tags = {
            letter: sum(
                1
                for artist in fig.findobj(match=Text)
                if artist.get_visible()
                and artist.get_text() == letter
                and artist.get_fontweight() == "bold"
            )
            for letter in "abcdefghi"
        }
        if panel_tags != {letter: 1 for letter in "abcdefghi"}:
            raise AssertionError(f"v4 panel tags drifted: {panel_tags}")
    
        return {
            "predecessor_v3_audit": v3_audit,
            "delta_only": True,
            "panel_b": {
                "regulators": 49,
                "fields_visible_for_every_row": ["rev/test", "rate", "fold", "L/G"],
                "fold_definition": "observed rate divided by out-degree-decile expectation",
                "row_values": row_strings,
                "bar_bounds_before": list(bar_bounds_before),
                "bar_bounds_after": list(bar_bounds_after),
                "metadata_bounds_before": list(metadata_bounds_before),
                "metadata_bounds_after": list(metadata_bounds_after),
                "bbox_audit": bbox_audit,
            },
            "panel_d_full_census_kegg": {
                "row_text_and_coordinates_unchanged": True,
                "table_bounds_before": list(kegg_bounds_before),
                "table_bounds_after": list(kegg_bounds_after),
                "note_bounds_after_parent_axes": list(PANEL_D_NOTE_RELATIVE_BOUNDS),
                "displayed_terms": 14,
                "fixed_family": 113,
                "conventional_ora_q_lt_0_10": 14,
                "regulator_opportunity_global_p": 0.18439,
                "ora_maxT_supported_terms": 0,
            },
            "all_other_axes_geometry_unchanged": True,
            "panel_tags": panel_tags,
            "font_audit": font_audit,
        }
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, v3_audit = V3.build_figure(regulators, kegg_plot, coverage)
        audit = apply_v4_delta(fig, regulators, v3_audit)
        return fig, audit
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "--force",
            action="store_true",
            help="replace only this v4 diagnostic bundle's generated files",
        )
        args = parser.parse_args()
        outputs = (
            OUT_PDF,
            OUT_PNG,
            OUT_REGULATORS,
            OUT_KEGG,
            OUT_MANIFEST,
            OUT_PROVENANCE,
        )
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(
                f"refusing to replace existing v4 outputs without --force: {existing}"
            )
    
        canonical_before = current_canonical_hashes()
        v1_before = V2.verify_v1()
        v2_before = V3.verify_v2()
        v3_before = verify_v3()
        kegg_source_records = V1.verify_kegg_scientific_sources()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        assert len(regulators) == 49
        assert len(kegg_plot) == 14
        assert int(regulators["n_inv"].sum()) == 2_075
        assert int(regulators["n_both"].sum()) == 14_434
        assert int(regulators["lost_on_activation"].sum()) == 1_154
        assert int(regulators["gained_on_activation"].sum()) == 921
        assert int(kegg_plot["n_lost_dyads"].sum()) == 495
        assert int(kegg_plot["n_gained_dyads"].sum()) == 504
    
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG candidate v4 with full regulator metadata",
                "Subject": "Noncanonical full-census KEGG and 49-regulator diagnostic",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        V1.atomic_figure(
            fig,
            OUT_PNG,
            format="png",
            dpi=V1.PNG_DPI,
            facecolor="white",
        )
        plt.close(fig)
    
        canonical_after = current_canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed during v4 render")
        v1_after = V2.verify_v1()
        v2_after = V3.verify_v2()
        v3_after = verify_v3()
        if v1_after != v1_before or v2_after != v2_before or v3_after != v3_before:
            raise AssertionError("a frozen predecessor changed during v4 render")
    
        renderer_paths = list(
            dict.fromkeys(
                [
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                ]
            )
        )
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v4-full-metadata",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "predecessors": {
                "v1": {"preserved": True, "before": v1_before, "after": v1_after},
                "v2": {"preserved": True, "before": v2_before, "after": v2_after},
                "v3": {"preserved": True, "before": v3_before, "after": v3_after},
            },
            "canonical_figure_2": {
                "before": canonical_before,
                "after": canonical_after,
                "written": False,
            },
            "delta": {
                "panel_b_fields_added_to_artwork_for_all_49": [
                    "rev/test",
                    "rate",
                    "fold",
                    "L/G",
                ],
                "panel_b_fold_definition": (
                    "observed reversal rate divided by out-degree-decile expectation"
                ),
                "panel_d_scientific_content_changed": False,
                "other_scientific_panels_changed": False,
            },
            "scientific_scope": {
                "regulator_panel": {
                    "degree_matched_bh_q_lt_0_05": 49,
                    "tested_edges": 14_434,
                    "reversal_edges": 2_075,
                    "lost_edges": 1_154,
                    "gained_edges": 921,
                },
                "full_census_kegg_panel": {
                    "population": "2,577 reversal-bearing targets among 9,554 eligible census targets",
                    "fixed_kegg_legacy_family": 113,
                    "displayed_conventional_ora_bh_q_lt_0_10": 14,
                    "regulator_opportunity_global_p": V1.KEGG.ORA_GLOBAL_P,
                    "ora_maxT_supported_terms": 0,
                    "interpretation": "descriptive conventional ORA sensitivity analysis",
                },
            },
            "audit": audit,
            "font_contract": {
                "minimum_pt": V1.FONT_MIN_PT,
                "maximum_pt": V1.FONT_MAX_PT,
                **audit["font_audit"],
            },
            "renderer_closure": {
                path: V1.sha256(PAPER / path) for path in renderer_paths
            },
            "scientific_source_closure": {
                V1.relative(path): V1.file_record(path)
                for path in V1.SCIENTIFIC_SOURCE_PATHS
            },
            "verified_kegg_scientific_sources": kegg_source_records,
            "bundle_sources": {
                "panel_b_regulators": V1.file_record(OUT_REGULATORS),
                "panel_d_signed_kegg": V1.file_record(OUT_KEGG),
            },
            "outputs": {
                "pdf": V1.file_record(OUT_PDF),
                "png": V1.file_record(OUT_PNG),
            },
            "claim_boundary": (
                "The 14 displayed KEGG Legacy terms pass only conventional Fisher/BH "
                "q<0.10. The full-family opportunity-aware test is not significant "
                "(global p=0.18439), and 0/113 terms survive ORA maxT. Panel d is a "
                "descriptive sensitivity analysis and does not establish pathway "
                "enrichment, activity, activation, suppression, or rewiring."
            ),
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        provenance = f"""# Figure 2 signed-KEGG candidate v4 provenance
    
    - Status: noncanonical diagnostic candidate for review; no manuscript or canonical file was written.
    - Renderer: `scripts/fig2_signed_kegg_side_by_side_candidate_v4.py`.
    - Lineage: v4 is a narrow presentation delta on v3, which is the audited continuation of `fig2_signed_kegg_side_by_side_candidate_v2.py`.
    - Panel b: all 49 degree-matched enriched regulators; every row shows `rev/test`, reversal `rate`, degree-decile `fold`, and lost/gained (`L/G`) counts.
    - Panel d: full reversal census, 2,577 reversal-bearing targets among 9,554 eligible targets; 14 conventional ORA rows shown from a fixed 113-term KEGG Legacy family.
    - Claim boundary: global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms. The KEGG list is descriptive, not supported pathway enrichment or pathway activity.
    - Canonical Figure 2 remained byte-identical during rendering: PDF `{canonical_after['pdf']}`, PNG `{canonical_after['png']}`.
    - Full hashes, source closure, font audit, and row-level geometry checks are recorded in `manifest.json`.
    """
        V1.atomic_text(OUT_PROVENANCE, provenance)
        print(
            json.dumps(
                {
                    "pdf": V1.file_record(OUT_PDF),
                    "png": V1.file_record(OUT_PNG),
                    "manifest": V1.relative(OUT_MANIFEST),
                    "provenance": V1.relative(OUT_PROVENANCE),
                    "font_range_pt": [
                        audit["font_audit"]["minimum_visible_font_pt"],
                        audit["font_audit"]["maximum_visible_font_pt"],
                    ],
                    "panel_b_rows_with_all_fields": audit["panel_b"]["bbox_audit"][
                        "rows_checked"
                    ],
                    "canonical_unchanged": canonical_before == canonical_after,
                    "v1_v2_v3_preserved": True,
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v4, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v4()

def init_fig2_signed_kegg_side_by_side_candidate_v5():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the thin, DOCX-sized v5 successor of signed-KEGG Figure 2.
    
    Version 5 inherits every data artist and statistic from v4/v3.  It redistributes
    the top row across available page width, removes the long in-artwork KEGG note,
    keeps the essential ORA correction statement and support symbols visible, and
    enforces the established 5--7 pt contract at the actual 5.833333-inch DOCX
    placement.  Only a noncanonical diagnostic bundle is written.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V4 = fig2_signed_kegg_side_by_side_candidate_v4
    
    
    V3 = V4.V3
    V2 = V4.V2
    V1 = V4.V1
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v5_docx_ready"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v5_docx_ready"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    V4_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "Fig2_signed_kegg_side_by_side_candidate_v4_full_metadata.pdf",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "Fig2_signed_kegg_side_by_side_candidate_v4_full_metadata.png",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "source_panel_b_regulators.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "source_panel_d_signed_kegg.tsv",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "manifest.json",
        PAPER
        / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v4_full_metadata"
        / "PROVENANCE.md",
    )
    EXPECTED_V4_SHA256 = {
        V4_PATHS[0]: "d1a7389d3ee8790e9be4c805671b2d852118f086742f135bed8db4150e29f6c8",
        V4_PATHS[1]: "59f594736e367d0ea9d28f902412cc7692b5f64ac1ce3c92013f1e2d84fa5c68",
        V4_PATHS[2]: "f98e79f38ec72dd8bba4a274c63b6fd92d8cd27eaaad6c4e435a695c136b78b5",
        V4_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V4_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V4_PATHS[5]: "cea21c762fe713f31df7973630407386366b3af662e7c42b82f7b326df56c687",
        V4_PATHS[6]: "b8d2ee6d2480945b2d05c0405f7bb8498c024a2521bc4157adbd12c3898b1228",
    }
    
    SOURCE_PAGE_WIDTH_IN = 8.27
    DELIVERED_WIDTH_IN = 5.833333333333333
    DELIVERED_SCALE = DELIVERED_WIDTH_IN / SOURCE_PAGE_WIDTH_IN
    SOURCE_FONT_MIN_PT = 7.09
    SOURCE_FONT_MAX_PT = 9.92
    
    BAR_BOUNDS = (0.412, 0.5424, 0.135, 0.3742)
    META_BOUNDS = (0.552, 0.5424, 0.198, 0.3742)
    KEGG_BOUNDS = (0.770, 0.5424, 0.210, 0.3420)
    
    
    def current_canonical_hashes() -> dict[str, str]:
        return {
            "pdf": V1.sha256(V1.CANONICAL_PDF),
            "png": V1.sha256(V1.CANONICAL_PNG),
        }
    
    
    def verify_v4() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V4_SHA256.items():
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"frozen v4 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def page_to_parent_x(parent: mpl.axes.Axes, page_x: float) -> float:
        left, _, width, _ = (float(value) for value in parent.get_position().bounds)
        return (page_x - left) / width
    
    
    def unique_text(artists: list[Text], value: str) -> Text:
        matches = [artist for artist in artists if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def enforce_delivered_typography(fig: plt.Figure) -> dict[str, object]:
        fig.canvas.draw()
        visible = [
            artist
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        before = [float(artist.get_fontsize()) for artist in visible]
        raised = 0
        lowered = 0
        for artist in visible:
            size = float(artist.get_fontsize())
            if size < SOURCE_FONT_MIN_PT:
                artist.set_fontsize(SOURCE_FONT_MIN_PT)
                raised += 1
            elif size > SOURCE_FONT_MAX_PT:
                artist.set_fontsize(SOURCE_FONT_MAX_PT)
                lowered += 1
        fig.canvas.draw()
        after = [
            float(artist.get_fontsize())
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        delivered = [size * DELIVERED_SCALE for size in after]
        if min(after) < SOURCE_FONT_MIN_PT or max(after) > SOURCE_FONT_MAX_PT:
            raise AssertionError(f"source font contract failed: {min(after)}--{max(after)}")
        if min(delivered) < 5.0 or max(delivered) > 7.0:
            raise AssertionError(
                f"delivered font contract failed: {min(delivered)}--{max(delivered)}"
            )
        return {
            "source_page_width_in": SOURCE_PAGE_WIDTH_IN,
            "delivered_width_in": DELIVERED_WIDTH_IN,
            "placement_scale": DELIVERED_SCALE,
            "visible_text_artists": len(after),
            "source_font_range_before_pt": [min(before), max(before)],
            "source_font_range_after_pt": [min(after), max(after)],
            "delivered_font_range_pt": [min(delivered), max(delivered)],
            "artists_raised_to_floor": raised,
            "artists_lowered_to_ceiling": lowered,
        }
    
    
    def kegg_bbox_audit(
        fig: plt.Figure,
        kegg_axis: mpl.axes.Axes,
    ) -> dict[str, object]:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        q_texts = [
            artist
            for artist in kegg_axis.texts
            if float(artist.get_position()[0]) == 0.685
            and float(artist.get_position()[1]) >= 0
        ]
        marker_x_pixels: dict[int, float] = {}
        for collection in kegg_axis.collections:
            if not isinstance(collection, PathCollection):
                continue
            offsets = np.asarray(collection.get_offsets(), dtype=float)
            for x_value, y_value in offsets:
                if np.isclose(x_value, 0.945) or np.isclose(x_value, 0.990):
                    x_pixel = float(kegg_axis.transData.transform((x_value, y_value))[0])
                    row = int(round(float(y_value)))
                    marker_x_pixels[row] = min(marker_x_pixels.get(row, x_pixel), x_pixel)
        gaps: list[float] = []
        for artist in q_texts:
            row = int(round(float(artist.get_position()[1])))
            if row in marker_x_pixels:
                q_right = float(artist.get_window_extent(renderer=renderer).x1)
                gaps.append((marker_x_pixels[row] - q_right) * 72.0 / fig.dpi)
        if len(q_texts) != 14 or len(gaps) != 14:
            raise AssertionError(
                f"unexpected KEGG q/marker audit inventory: q={len(q_texts)}, gaps={len(gaps)}"
            )
        if min(gaps) < 1.0:
            raise AssertionError(f"KEGG q/marker gap only {min(gaps):.3f} pt")
        figure_bbox = fig.bbox
        all_visible = []
        for artist in fig.findobj(match=Text):
            if not artist.get_visible() or not artist.get_text().strip():
                continue
            bbox_artist = artist.get_window_extent(renderer=renderer)
            if (bbox_artist.y0 + bbox_artist.y1) / 2.0 >= 0.50 * figure_bbox.y1:
                all_visible.append(artist)
        right_edges = [
            (
                float(artist.get_window_extent(renderer=renderer).x1 - figure_bbox.x1),
                artist.get_text(),
            )
            for artist in all_visible
        ]
        right_overflow = max(value for value, _ in right_edges)
        if right_overflow > 0.5:
            offenders = sorted(right_edges, reverse=True)[:5]
            raise AssertionError(
                f"visible text extends {right_overflow:.3f}px past page right: {offenders}"
            )
        return {
            "q_rows_checked": len(q_texts),
            "minimum_q_to_marker_center_gap_pt": min(gaps),
            "maximum_text_right_overflow_px": max(0.0, right_overflow),
        }
    
    
    def apply_v5_delta(
        fig: plt.Figure,
        kegg_plot: pd.DataFrame,
        v4_audit: dict[str, object],
    ) -> dict[str, object]:
        parent = V2.find_parent_axis(fig)
        bar_axis = next(
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if axis.get_xlabel() == "reversal edges (0–300)"
        )
        metadata_axis = next(
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if any(artist.get_text() == "rev/test" for artist in axis.texts)
        )
        kegg_axis = V2.find_unique_axis_with_text(fig, "KEGG pathway")
    
        before = {
            "bar": list(float(value) for value in bar_axis.get_position().bounds),
            "metadata": list(float(value) for value in metadata_axis.get_position().bounds),
            "kegg": list(float(value) for value in kegg_axis.get_position().bounds),
        }
        for axis, bounds in (
            (bar_axis, BAR_BOUNDS),
            (metadata_axis, META_BOUNDS),
            (kegg_axis, KEGG_BOUNDS),
        ):
            axis.set_axes_locator(None)
            axis.set_position(bounds)
    
        note_text = unique_text(list(parent.texts), V2.NOTE_TEXT)
        note_text.remove()
        note_boxes = [
            patch
            for patch in list(parent.patches)
            if isinstance(patch, mpl.patches.Rectangle)
            and np.isclose(float(patch.get_height()), 0.188, atol=1e-14)
        ]
        if len(note_boxes) != 1:
            raise AssertionError(f"expected one long KEGG note box; found {len(note_boxes)}")
        note_boxes[0].remove()
    
        title_b = unique_text(list(parent.texts), "49 above degree-matched expectation")
        subtitle_b = unique_text(
            list(parent.texts),
            "2,075/14,434 reversals/tested; BH q<0.05; ranks 1–49",
        )
        footer_b = unique_text(
            list(parent.texts),
            "red lost; teal gained; fold = rate / out-degree-decile expectation",
        )
        title_b.set_x(page_to_parent_x(parent, (0.412 + 0.750) / 2.0))
        subtitle_b.set_x(title_b.get_position()[0])
        footer_b.set_x(title_b.get_position()[0])
        subtitle_b.set_text("2,075/14,434 rev/test; BH q<.05; ranks 1–49")
        footer_b.set_text("red lost; teal gained; fold vs degree-decile expectation")
        tag_b = unique_text(
            [
                artist
                for artist in parent.texts
                if artist.get_text() == "b" and artist.get_fontweight() == "bold"
            ],
            "b",
        )
        tag_b.set_x(page_to_parent_x(parent, 0.397))
    
        title_d = unique_text(list(parent.texts), "Full-census signed KEGG")
        subtitle_d = unique_text(
            list(parent.texts),
            "Full census: 2,577 targets / 4,379 dyads — not T279",
        )
        statistic_d = unique_text(
            list(parent.texts),
            "14 descriptive ORA q<0.10; global p=0.18439; ORA maxT 0/113",
        )
        panel_d_center = page_to_parent_x(parent, (0.770 + 0.980) / 2.0)
        title_d.set_x(panel_d_center)
        title_d.set_text("Full-census KEGG ORA")
        subtitle_d.set_x(panel_d_center)
        subtitle_d.set_text("2,577 reversal targets / 9,554 eligible")
        statistic_d.set_x(panel_d_center)
        statistic_d.set_text("14 Fisher/BH q<.10\nglobal p=.184; maxT 0/113")
        statistic_d.set_linespacing(0.94)
        tag_d = unique_text(
            [
                artist
                for artist in parent.texts
                if artist.get_text() == "d" and artist.get_fontweight() == "bold"
            ],
            "d",
        )
        tag_d.set_x(page_to_parent_x(parent, 0.758))
        parent.text(
            panel_d_center,
            0.018,
            "○ agrees   × differs   ◆ both tests",
            transform=parent.transAxes,
            ha="center",
            va="bottom",
            fontsize=V1.FONT_MIN_PT,
            color=V1.MUTED,
        )
    
        header_pathway = unique_text(list(kegg_axis.texts), "KEGG pathway")
        header_pathway.set_text("pathway")
        header_lg = unique_text(list(kegg_axis.texts), "L/G\ndyads")
        header_lg.set_text("L/G")
        header_q = unique_text(list(kegg_axis.texts), "target ORA\nq")
        header_q.set_text("ORA\nq")
        header_q.set_x(0.685)
        header_lg.set_x(0.505)
    
        for row_index, (_, row) in enumerate(kegg_plot.iterrows()):
            pathway_matches = [
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.475)
            ]
            if len(pathway_matches) != 1:
                raise AssertionError(
                    f"expected one KEGG pathway label at row {row_index}; "
                    f"found {len(pathway_matches)}"
                )
            pathway_matches[0].set_text(
                "\n".join(textwrap.wrap(str(row["kegg_display_name"]), width=18))
            )
            count_matches = [
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.585)
            ]
            q_matches = [
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.715)
            ]
            if len(count_matches) != 1 or len(q_matches) != 1:
                raise AssertionError(f"KEGG numeric row inventory failed at {row_index}")
            count_matches[0].set_x(0.505)
            q_matches[0].set_x(0.685)
    
        for collection in kegg_axis.collections:
            if not isinstance(collection, PathCollection):
                continue
            offsets = np.asarray(collection.get_offsets(), dtype=float)
            if offsets.size == 0:
                continue
            changed = offsets.copy()
            changed[np.isclose(changed[:, 0], 0.895), 0] = 0.945
            changed[np.isclose(changed[:, 0], 0.970), 0] = 0.990
            collection.set_offsets(changed)
    
        # Reuse the previously audited one-sided panel-h gutter correction before
        # increasing source typography for the actual DOCX placement.
        orientation_axis = next(
            axis
            for axis in fig.axes
            if axis.get_title().startswith("Target orientation")
        )
        panel_h_before = orientation_axis.get_position().frozen()
        panel_h_shift = 0.012
        orientation_axis.set_position(
            [
                panel_h_before.x0 + panel_h_shift,
                panel_h_before.y0,
                panel_h_before.width - panel_h_shift,
                panel_h_before.height,
            ]
        )
    
        typography = enforce_delivered_typography(fig)
        bbox = kegg_bbox_audit(fig, kegg_axis)
        visible = " ".join(
            "\n".join(
                artist.get_text()
                for artist in fig.findobj(match=Text)
                if artist.get_visible() and artist.get_text().strip()
            ).split()
        )
        for required in (
            "rev/test rate fold L/G",
            "2,075/14,434 rev/test",
            "14 Fisher/BH q<.10 global p=.184; maxT 0/113",
            "2,577 reversal targets / 9,554 eligible",
            "Antigen processing and presentation",
            "agrees × differs ◆ both tests",
        ):
            if " ".join(required.split()) not in visible:
                raise AssertionError(f"required v5 text missing: {required}")
    
        return {
            "predecessor_v4_audit": v4_audit,
            "delta_only": True,
            "top_row_page_bounds": {
                "panel_b_bar": list(BAR_BOUNDS),
                "panel_b_numeric_table": list(META_BOUNDS),
                "panel_d_kegg": list(KEGG_BOUNDS),
                "before": before,
            },
            "panel_b": {
                "regulators": 49,
                "fields_visible_for_all_rows": ["rev/test", "rate", "fold", "L/G"],
            },
            "panel_d": {
                "long_note_box_removed": True,
                "support_symbols_retained": True,
                "visible_scope": "14 Fisher/BH q<.10; global p=.184; maxT 0/113",
                "displayed_terms": 14,
                "fixed_family": 113,
                "bbox_audit": bbox,
            },
            "panel_h_left_shift_fraction_page": panel_h_shift,
            "typography": typography,
        }
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, v4_audit = V4.build_figure(regulators, kegg_plot, coverage)
        audit = apply_v5_delta(fig, kegg_plot, v4_audit)
        return fig, audit
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument(
            "--force",
            action="store_true",
            help="replace only this v5 diagnostic bundle's generated files",
        )
        args = parser.parse_args()
        outputs = (
            OUT_PDF,
            OUT_PNG,
            OUT_REGULATORS,
            OUT_KEGG,
            OUT_MANIFEST,
            OUT_PROVENANCE,
        )
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(
                f"refusing to replace existing v5 outputs without --force: {existing}"
            )
    
        canonical_before = current_canonical_hashes()
        v4_before = verify_v4()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG candidate v5 DOCX-ready",
                "Subject": "Noncanonical full-census KEGG and 49-regulator diagnostic",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        V1.atomic_figure(
            fig,
            OUT_PNG,
            format="png",
            dpi=V1.PNG_DPI,
            facecolor="white",
        )
        plt.close(fig)
    
        canonical_after = current_canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed during v5 render")
        v4_after = verify_v4()
        if v4_after != v4_before:
            raise AssertionError("frozen v4 predecessor changed during v5 render")
    
        renderer_paths = list(
            dict.fromkeys(
                [
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                ]
            )
        )
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v5-docx-ready",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "predecessor_v4": {
                "preserved": True,
                "before": v4_before,
                "after": v4_after,
            },
            "canonical_figure_2": {
                "before": canonical_before,
                "after": canonical_after,
                "written": False,
            },
            "audit": audit,
            "scientific_scope": {
                "regulators": {
                    "displayed": 49,
                    "reversal_edges": 2_075,
                    "tested_edges": 14_434,
                    "lost": 1_154,
                    "gained": 921,
                },
                "kegg": {
                    "reversal_targets": 2_577,
                    "eligible_census_targets": 9_554,
                    "fixed_legacy_family": 113,
                    "conventional_fisher_bh_q_lt_0_10": 14,
                    "opportunity_aware_global_p": V1.KEGG.ORA_GLOBAL_P,
                    "ora_maxT_supported_terms": 0,
                },
            },
            "font_contract": audit["typography"],
            "renderer_closure": {
                path: V1.sha256(PAPER / path) for path in renderer_paths
            },
            "scientific_source_closure": {
                V1.relative(path): V1.file_record(path)
                for path in V1.SCIENTIFIC_SOURCE_PATHS
            },
            "bundle_sources": {
                "panel_b_regulators": V1.file_record(OUT_REGULATORS),
                "panel_d_signed_kegg": V1.file_record(OUT_KEGG),
            },
            "outputs": {
                "pdf": V1.file_record(OUT_PDF),
                "png": V1.file_record(OUT_PNG),
            },
            "claim_boundary": (
                "Panel d shows a conventional Fisher/BH KEGG sensitivity analysis: "
                "14 terms have q<0.10, but the opportunity-aware global test is null "
                "(p=0.18439) and 0/113 survive ORA maxT. It is not supported pathway "
                "enrichment or pathway activity."
            ),
            "software": {
                "python": platform.python_version(),
                "matplotlib": matplotlib.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        provenance = f"""# Figure 2 signed-KEGG candidate v5 provenance
    
    - Status: noncanonical diagnostic candidate; canonical Figure 2 and manuscript were not written.
    - Lineage: thin successor to v4/v3/v2; no analysis was recomputed.
    - Panel b: all 49 rows show `rev/test`, `rate`, degree-decile `fold`, and `L/G`.
    - Panel d: 14 conventional Fisher/BH KEGG rows from 2,577 reversal targets against 9,554 eligible census targets.
    - Visible correction: global opportunity-aware p=.184 and 0/113 ORA maxT.
    - Typography: source {SOURCE_FONT_MIN_PT:.2f}–{SOURCE_FONT_MAX_PT:.2f} pt, delivering 5–7 pt at {DELIVERED_WIDTH_IN:.6f} inches.
    - Canonical hashes remained PDF `{canonical_after['pdf']}` and PNG `{canonical_after['png']}`.
    - Full sources, geometry, hashes and claim boundary are in `manifest.json`.
    """
        V1.atomic_text(OUT_PROVENANCE, provenance)
        print(
            json.dumps(
                {
                    "pdf": V1.file_record(OUT_PDF),
                    "png": V1.file_record(OUT_PNG),
                    "manifest": V1.relative(OUT_MANIFEST),
                    "provenance": V1.relative(OUT_PROVENANCE),
                    "delivered_font_range_pt": audit["typography"][
                        "delivered_font_range_pt"
                    ],
                    "canonical_unchanged": True,
                },
                sort_keys=True,
            )
        )
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v5, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v5()

def init_fig2_signed_kegg_side_by_side_candidate_v7():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the clean v7 signed-KEGG Figure 2 review candidate.
    
    This thin successor reuses the v4/v3/v2 scientific lineage.  It applies v5's
    clean top-row redistribution without its whole-figure typography pass, so all
    unrelated panels retain their predecessor typography and geometry.  Only the
    new b/d material is sized for 5--7 pt at the actual DOCX placement.  Panel b
    uses safe vertical bounds and aligned numeric columns; nonessential footers are
    removed.  Canonical Figure 2 and manuscript files are never written.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V5 = fig2_signed_kegg_side_by_side_candidate_v5
    
    
    V4 = V5.V4
    V3 = V5.V3
    V2 = V5.V2
    V1 = V5.V1
    OUTDIR = (
        PAPER
        / "figures_final"
        / "diagnostics"
        / "fig2_signed_kegg_side_by_side_candidate_v7_layout_clean"
    )
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v7_layout_clean"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    BAR_BOUNDS = (0.412, 0.495, 0.135, 0.422)
    META_BOUNDS = (0.552, 0.495, 0.198, 0.422)
    ROW_FONT_SOURCE_PT = V5.SOURCE_FONT_MIN_PT
    HEADING_FONT_SOURCE_PT = V5.SOURCE_FONT_MAX_PT
    
    
    def current_canonical_hashes() -> dict[str, str]:
        return {
            "pdf": V1.sha256(V1.CANONICAL_PDF),
            "png": V1.sha256(V1.CANONICAL_PNG),
        }
    
    
    def unique_text(artists: list[Text], value: str) -> Text:
        matches = [artist for artist in artists if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def no_global_typography(fig: plt.Figure) -> dict[str, object]:
        visible = [
            artist
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
        sizes = [float(artist.get_fontsize()) for artist in visible]
        return {
            "whole_figure_typography_changed": False,
            "visible_text_artists": len(visible),
            "source_font_range_unchanged_pt": [min(sizes), max(sizes)],
        }
    
    
    def metadata_bbox_audit(
        fig: plt.Figure,
        metadata_axis: mpl.axes.Axes,
    ) -> dict[str, object]:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        horizontal_gaps: list[float] = []
        for row_index in range(49):
            artists = [
                artist
                for artist in metadata_axis.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            if len(artists) != 6:
                raise AssertionError(
                    f"row {row_index} has {len(artists)} metadata fields, expected 6"
                )
            boxes = [artist.get_window_extent(renderer=renderer) for artist in artists]
            horizontal_gaps.extend(
                (right.x0 - left.x1) * 72.0 / fig.dpi
                for left, right in zip(boxes[:-1], boxes[1:])
            )
        y_pixels = [
            float(metadata_axis.transData.transform((0.0, float(row_index)))[1])
            for row_index in range(49)
        ]
        row_pitches = [
            abs(right - left) * 72.0 / fig.dpi
            for left, right in zip(y_pixels[:-1], y_pixels[1:])
        ]
        if min(horizontal_gaps) < 0.35:
            raise AssertionError(
                f"panel-b fields overlap; minimum gap {min(horizontal_gaps):.3f} pt"
            )
        if min(row_pitches) < ROW_FONT_SOURCE_PT:
            raise AssertionError(
                f"panel-b row pitch {min(row_pitches):.3f} pt is below font size"
            )
        return {
            "rows_checked": 49,
            "minimum_horizontal_gap_pt": min(horizontal_gaps),
            "minimum_row_pitch_pt": min(row_pitches),
        }
    
    
    def apply_v7_delta(
        fig: plt.Figure,
        kegg_plot: pd.DataFrame,
        v4_audit: dict[str, object],
    ) -> dict[str, object]:
        parent = V2.find_parent_axis(fig)
        outside_axes = [
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if axis is not parent
            and axis.get_xlabel() != "reversal edges (0–300)"
            and not any(artist.get_text() == "rev/test" for artist in axis.texts)
            and not any(artist.get_text() == "KEGG pathway" for artist in axis.texts)
        ]
        outside_positions = {
            id(axis): tuple(float(value) for value in axis.get_position().bounds)
            for axis in outside_axes
        }
        outside_fonts = {
            id(artist): float(artist.get_fontsize())
            for axis in outside_axes
            for artist in axis.findobj(match=Text)
        }
    
        original_typography = V5.enforce_delivered_typography
        V5.enforce_delivered_typography = no_global_typography
        try:
            v5_layout_audit = V5.apply_v5_delta(fig, kegg_plot, v4_audit)
        finally:
            V5.enforce_delivered_typography = original_typography
    
        for axis in outside_axes:
            axis.set_position(outside_positions[id(axis)])
            for artist in axis.findobj(match=Text):
                if id(artist) in outside_fonts:
                    artist.set_fontsize(outside_fonts[id(artist)])
    
        bar_axis = next(
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if axis.get_xlabel() == "reversal edges (0–300)"
        )
        metadata_axis = next(
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if any(artist.get_text() == "rev/test" for artist in axis.texts)
        )
        kegg_axis = V2.find_unique_axis_with_text(fig, "pathway")
        for axis, bounds in ((bar_axis, BAR_BOUNDS), (metadata_axis, META_BOUNDS)):
            axis.set_axes_locator(None)
            axis.set_position(bounds)
        bar_axis.set_xlabel("")
    
        footer_b = unique_text(
            list(parent.texts),
            "red lost; teal gained; fold vs degree-decile expectation",
        )
        footer_b.remove()
        symbol_key = unique_text(
            list(parent.texts),
            "○ agrees   × differs   ◆ both tests",
        )
        symbol_key.remove()
    
        header_rev = unique_text(list(metadata_axis.texts), "rev/test")
        header_rate = unique_text(list(metadata_axis.texts), "rate")
        header_fold = unique_text(list(metadata_axis.texts), "fold")
        header_lg = unique_text(list(metadata_axis.texts), "L/G")
        header_rev.set_x(0.00)
        header_rate.set_x(0.55)
        header_rate.set_text("rate %")
        header_fold.set_x(0.77)
        header_fold.set_text("fold ×")
        header_lg.set_x(0.92)
        for row_index in range(49):
            artists = [
                artist
                for artist in metadata_axis.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            if len(artists) != 6:
                raise AssertionError(f"unexpected metadata inventory at row {row_index}")
            rev, rate, fold, lost, slash, gained = artists
            rev.set_x(0.00)
            rate.set_x(0.55)
            rate.set_text(rate.get_text().removesuffix("%"))
            fold.set_x(0.77)
            fold.set_text(fold.get_text().removesuffix("×"))
            lost.set_x(0.905)
            slash.set_x(0.920)
            gained.set_x(0.935)
    
        for axis in (bar_axis, metadata_axis, kegg_axis):
            for artist in axis.findobj(match=Text):
                if artist.get_visible() and artist.get_text().strip():
                    artist.set_fontsize(ROW_FONT_SOURCE_PT)
        for artist in (header_rev, header_rate, header_fold, header_lg):
            artist.set_fontsize(ROW_FONT_SOURCE_PT)
            artist.set_fontweight("bold")
    
        top_titles = (
            unique_text(list(parent.texts), "49 above degree-matched expectation"),
            unique_text(list(parent.texts), "Full-census KEGG ORA"),
        )
        top_subtitles = (
            unique_text(list(parent.texts), "2,075/14,434 rev/test; BH q<.05; ranks 1–49"),
            unique_text(list(parent.texts), "2,577 reversal targets / 9,554 eligible"),
            unique_text(list(parent.texts), "14 Fisher/BH q<.10\nglobal p=.184; maxT 0/113"),
        )
        for artist in top_titles:
            artist.set_fontsize(HEADING_FONT_SOURCE_PT)
        for artist in top_subtitles:
            artist.set_fontsize(ROW_FONT_SOURCE_PT)
        for tag in ("b", "d"):
            matches = [
                artist
                for artist in parent.texts
                if artist.get_text() == tag and artist.get_fontweight() == "bold"
            ]
            if len(matches) != 1:
                raise AssertionError(f"expected one top panel tag {tag}")
            matches[0].set_fontsize(HEADING_FONT_SOURCE_PT)
    
        metadata_audit = metadata_bbox_audit(fig, metadata_axis)
        kegg_audit = V5.kegg_bbox_audit(fig, kegg_axis)
        delivered = [
            float(artist.get_fontsize()) * V5.DELIVERED_SCALE
            for axis in (bar_axis, metadata_axis, kegg_axis)
            for artist in axis.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ] + [
            float(artist.get_fontsize()) * V5.DELIVERED_SCALE
            for artist in (*top_titles, *top_subtitles)
        ]
        if min(delivered) < 5.0 or max(delivered) > 7.0:
            raise AssertionError(f"new-panel delivered font range failed: {delivered}")
    
        return {
            "predecessor_v4_audit": v4_audit,
            "v5_layout_audit": v5_layout_audit,
            "outside_panels_geometry_restored": True,
            "outside_panels_typography_restored": True,
            "panel_b": {
                "bar_bounds": list(BAR_BOUNDS),
                "metadata_bounds": list(META_BOUNDS),
                "fields": ["rev/test", "rate %", "fold ×", "L/G"],
                "bbox_audit": metadata_audit,
                "xlabel_removed": True,
                "footer_removed": True,
            },
            "panel_d": {
                "long_note_removed": True,
                "symbol_key_removed_from_artwork": True,
                "visible_correction": "14 Fisher/BH q<.10; global p=.184; maxT 0/113",
                "bbox_audit": kegg_audit,
            },
            "new_panel_delivered_font_range_pt": [min(delivered), max(delivered)],
        }
    
    
    def build_figure(
        regulators: pd.DataFrame,
        kegg_plot: pd.DataFrame,
        coverage: dict[str, object],
    ) -> tuple[plt.Figure, dict[str, object]]:
        fig, v4_audit = V4.build_figure(regulators, kegg_plot, coverage)
        return fig, apply_v7_delta(fig, kegg_plot, v4_audit)
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true")
        args = parser.parse_args()
        outputs = (
            OUT_PDF,
            OUT_PNG,
            OUT_REGULATORS,
            OUT_KEGG,
            OUT_MANIFEST,
            OUT_PROVENANCE,
        )
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to replace v7 outputs: {existing}")
    
        canonical_before = current_canonical_hashes()
        v4_before = V5.verify_v4()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(
            fig,
            OUT_PDF,
            format="pdf",
            facecolor="white",
            metadata={
                "Title": "Figure 2 signed-KEGG candidate v7 layout clean",
                "Subject": "Noncanonical full-census KEGG and 49-regulator diagnostic",
                "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        V1.atomic_figure(fig, OUT_PNG, format="png", dpi=V1.PNG_DPI, facecolor="white")
        plt.close(fig)
        canonical_after = current_canonical_hashes()
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 2 changed")
        v4_after = V5.verify_v4()
        if v4_after != v4_before:
            raise AssertionError("frozen v4 changed")
    
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v7-layout-clean",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "canonical_figure_2": {"before": canonical_before, "after": canonical_after, "written": False},
            "predecessor_v4_preserved": {"before": v4_before, "after": v4_after},
            "audit": audit,
            "scientific_scope": {
                "panel_b": {"regulators": 49, "reversals": 2_075, "tested": 14_434, "lost": 1_154, "gained": 921},
                "panel_d": {"reversal_targets": 2_577, "eligible_targets": 9_554, "kegg_family": 113, "fisher_bh_q_lt_0_10": 14, "global_p": V1.KEGG.ORA_GLOBAL_P, "maxT_supported": 0},
            },
            "renderer_closure": {
                path: V1.sha256(PAPER / path)
                for path in (
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                )
            },
            "scientific_source_closure": {
                V1.relative(path): V1.file_record(path) for path in V1.SCIENTIFIC_SOURCE_PATHS
            },
            "bundle_sources": {"panel_b_regulators": V1.file_record(OUT_REGULATORS), "panel_d_signed_kegg": V1.file_record(OUT_KEGG)},
            "outputs": {"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG)},
            "claim_boundary": "The 14 KEGG rows are conventional Fisher/BH sensitivity results; global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms do not support pathway enrichment or activity.",
            "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        V1.atomic_text(
            OUT_PROVENANCE,
            "# Figure 2 signed-KEGG candidate v7 provenance\n\n"
            "- Noncanonical diagnostic; no canonical or manuscript write.\n"
            "- Thin v4/v3/v2 lineage; v5 top-row geometry reused without global typography changes.\n"
            "- Panel b shows rev/test, rate, fold and L/G for all 49 regulators.\n"
            "- Panel d shows 14 conventional KEGG Fisher/BH rows with visible global p=.184 and 0/113 maxT.\n"
            "- Full hashes, sources and collision audits are in manifest.json.\n",
        )
        print(json.dumps({
            "pdf": V1.file_record(OUT_PDF),
            "png": V1.file_record(OUT_PNG),
            "manifest": V1.relative(OUT_MANIFEST),
            "metadata_min_gap_pt": audit["panel_b"]["bbox_audit"]["minimum_horizontal_gap_pt"],
            "metadata_row_pitch_pt": audit["panel_b"]["bbox_audit"]["minimum_row_pitch_pt"],
            "new_panel_delivered_font_range_pt": audit["new_panel_delivered_font_range_pt"],
            "canonical_unchanged": True,
        }, sort_keys=True))
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v7, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v7()

def init_fig2_signed_kegg_side_by_side_candidate_v8():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render v8: the title-clearance-only successor to Figure 2 candidate v7."""
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V7 = fig2_signed_kegg_side_by_side_candidate_v7
    
    
    V5 = V7.V5
    V2 = V7.V2
    V1 = V7.V1
    OUTDIR = PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review"
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v8_final_review"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    TITLE_SOURCE_PT = 8.5
    
    V7_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/Fig2_signed_kegg_side_by_side_candidate_v7_layout_clean.pdf",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/Fig2_signed_kegg_side_by_side_candidate_v7_layout_clean.png",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/source_panel_b_regulators.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/source_panel_d_signed_kegg.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/manifest.json",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v7_layout_clean/PROVENANCE.md",
    )
    EXPECTED_V7 = {
        V7_PATHS[0]: "a7e4949603d012ba8c909635030565b907916cc03b1b7778983c3ace55fab238",
        V7_PATHS[1]: "b02b83b323d38746054b07db7e0513d8899c4f3c8930bcc706f42f8e30a3894c",
        V7_PATHS[2]: "78c5bb7ead527d5654058c22226d6a249d429ed2cafbb0ccc65f93e4d391f326",
        V7_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V7_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V7_PATHS[5]: "46340a0d0bfdfe72a12eea16b4bec4c47ccc407177c2718fa80f428fedeefcaf",
        V7_PATHS[6]: "b40996f7de29d2df0645c121985bf82c1845f6e95b433f91c1b36d206ffdb640",
    }
    
    
    def verify_v7() -> dict[str, dict[str, object]]:
        records = {}
        for path, expected in EXPECTED_V7.items():
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"v7 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def current_canonical() -> dict[str, str]:
        return {"pdf": V1.sha256(V1.CANONICAL_PDF), "png": V1.sha256(V1.CANONICAL_PNG)}
    
    
    def unique(parent, value: str) -> Text:
        matches = [artist for artist in parent.texts if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def apply_v8(fig: plt.Figure, v7_audit: dict[str, object]) -> dict[str, object]:
        parent = V2.find_parent_axis(fig)
        title_b = unique(parent, "49 above degree-matched expectation")
        title_d = unique(parent, "Full-census KEGG ORA")
        subtitle_b = unique(parent, "2,075/14,434 rev/test; BH q<.05; ranks 1–49")
        subtitle_d = unique(parent, "2,577 reversal targets / 9,554 eligible")
        statistic_d = unique(parent, "14 Fisher/BH q<.10\nglobal p=.184; maxT 0/113")
        tag_b = unique(parent, "b")
        tag_d = unique(parent, "d")
        for artist in (title_b, title_d, tag_b, tag_d):
            artist.set_fontsize(TITLE_SOURCE_PT)
        for artist in (subtitle_b, subtitle_d, statistic_d):
            artist.set_fontsize(V5.SOURCE_FONT_MIN_PT)
    
        title_b.set_x(V5.page_to_parent_x(parent, 0.580))
        subtitle_b.set_x(title_b.get_position()[0])
        tag_b.set_x(V5.page_to_parent_x(parent, 0.397))
        title_d.set_x(V5.page_to_parent_x(parent, 0.875))
        subtitle_d.set_x(title_d.get_position()[0])
        statistic_d.set_x(title_d.get_position()[0])
        tag_d.set_x(V5.page_to_parent_x(parent, 0.758))
    
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ordered = (tag_b, title_b, tag_d, title_d)
        boxes = [artist.get_window_extent(renderer=renderer) for artist in ordered]
        gaps = [
            (right.x0 - left.x1) * 72.0 / fig.dpi
            for left, right in zip(boxes[:-1], boxes[1:])
        ]
        right_margin = (fig.bbox.x1 - boxes[-1].x1) * 72.0 / fig.dpi
        if min(gaps) < 2.0 or right_margin < 2.0:
            raise AssertionError(
                f"top heading clearance failed: gaps={gaps}, right={right_margin}"
            )
        delivered_heading = TITLE_SOURCE_PT * V5.DELIVERED_SCALE
        if not 5.0 <= delivered_heading <= 7.0:
            raise AssertionError(f"delivered heading {delivered_heading} pt")
        return {
            "predecessor_v7_audit": v7_audit,
            "delta_only": True,
            "source_title_and_tag_pt": TITLE_SOURCE_PT,
            "delivered_title_and_tag_pt": delivered_heading,
            "ordered_clearance_pairs": ["b-tag to b-title", "b-title to d-tag", "d-tag to d-title"],
            "clearance_pt": gaps,
            "right_margin_pt": right_margin,
            "all_scientific_panels_and_values_unchanged": True,
        }
    
    
    def build_figure(regulators: pd.DataFrame, kegg_plot: pd.DataFrame, coverage: dict[str, object]):
        fig, v7_audit = V7.build_figure(regulators, kegg_plot, coverage)
        return fig, apply_v8(fig, v7_audit)
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true")
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST, OUT_PROVENANCE)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to replace v8 outputs: {existing}")
        canonical_before = current_canonical()
        v7_before = verify_v7()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        V1.atomic_figure(fig, OUT_PDF, format="pdf", facecolor="white", metadata={"Title": "Figure 2 signed-KEGG candidate v8 final review", "Creator": "CoDEG_Tcell deterministic Matplotlib renderer", "CreationDate": None, "ModDate": None})
        V1.atomic_figure(fig, OUT_PNG, format="png", dpi=V1.PNG_DPI, facecolor="white")
        plt.close(fig)
        canonical_after = current_canonical()
        v7_after = verify_v7()
        if canonical_after != canonical_before or v7_after != v7_before:
            raise AssertionError("canonical or v7 changed during v8 render")
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v8-final-review",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "canonical_figure_2": {"before": canonical_before, "after": canonical_after, "written": False},
            "predecessor_v7": {"before": v7_before, "after": v7_after, "preserved": True},
            "audit": audit,
            "scientific_scope": {"panel_b": {"regulators": 49, "reversals": 2_075, "tested": 14_434, "fields": ["rev/test", "rate", "fold", "L/G"]}, "panel_d": {"reversal_targets": 2_577, "eligible_targets": 9_554, "fisher_bh_terms": 14, "global_p": 0.18439, "maxT": "0/113"}},
            "renderer_closure": {path: V1.sha256(PAPER / path) for path in ("scripts/fig2_signed_kegg_side_by_side_candidate_v8.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py", *V1.RENDERER_PATHS)},
            "scientific_source_closure": {V1.relative(path): V1.file_record(path) for path in V1.SCIENTIFIC_SOURCE_PATHS},
            "bundle_sources": {"panel_b": V1.file_record(OUT_REGULATORS), "panel_d": V1.file_record(OUT_KEGG)},
            "outputs": {"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG)},
            "claim_boundary": "The 14 KEGG rows are conventional Fisher/BH sensitivity results; global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms do not support pathway enrichment or activity.",
            "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        V1.atomic_text(OUT_PROVENANCE, "# Figure 2 signed-KEGG candidate v8 provenance\n\n- Noncanonical review candidate; no canonical/manuscript write.\n- Title-clearance-only successor to v7; all data and statistics unchanged.\n- See manifest.json for complete hashes, sources, collision audit, and claim boundary.\n")
        print(json.dumps({"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG), "manifest": V1.relative(OUT_MANIFEST), "title_clearance_pt": audit["clearance_pt"], "right_margin_pt": audit["right_margin_pt"], "canonical_unchanged": True}, sort_keys=True))
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v8, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v8()

def init_fig2_signed_kegg_side_by_side_candidate_v9():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the compact, exact-typography v9 Figure 2 review candidate.
    
    This is a presentation-only successor to the frozen v8 candidate.  It narrows
    panels a/c and panel-b metadata, gives the recovered width to the full-census
    KEGG panel, restores the canonical functional colours of all regulator labels,
    and enforces the delivered 5/6/7-point typography contract.  No scientific
    value, row, panel, canonical figure, manuscript, or predecessor is changed.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V8 = fig2_signed_kegg_side_by_side_candidate_v8
    
    
    V7 = V8.V7
    V5 = V8.V5
    V2 = V8.V2
    V1 = V8.V1
    BASE = V1.BASE
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main"
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v9_compact_main"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    V8_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v8.py",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/Fig2_signed_kegg_side_by_side_candidate_v8_final_review.pdf",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/Fig2_signed_kegg_side_by_side_candidate_v8_final_review.png",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/source_panel_b_regulators.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/source_panel_d_signed_kegg.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/manifest.json",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v8_final_review/PROVENANCE.md",
    )
    EXPECTED_V8 = {
        V8_PATHS[0]: "0e8c6d6e83a206725dd150fb0089fef3c32a993822a7eab3fffbeb5a8168f2b4",
        V8_PATHS[1]: "039038e5a4e8852a217f6b22e2a7761b199d5651714163050008b23ad9c18bb6",
        V8_PATHS[2]: "d5cb55ab48cc01332e7ed35f5760ded51079ffe14efbf9ab112892b8f1a14854",
        V8_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V8_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V8_PATHS[5]: "5724cd63675ac87b4d486522c00a400408d360e243ee19595c28db97c7e3c161",
        V8_PATHS[6]: "ec9f9bb4de9aab06c511a2f31ecd40f9fd6ad33eecdc31bab74cdffccb85c0e0",
    }
    
    SOURCE_5_PT = 5.0 / V5.DELIVERED_SCALE
    SOURCE_6_PT = 6.0 / V5.DELIVERED_SCALE
    SOURCE_7_PT = 7.0 / V5.DELIVERED_SCALE
    
    PANEL_A_WIDTH = 0.226
    PANEL_C_WIDTH = 0.226
    BAR_BOUNDS = (0.377, 0.495, 0.135, 0.422)
    META_BOUNDS = (0.517, 0.495, 0.185, 0.422)
    KEGG_BOUNDS = (0.718, 0.5424, 0.262, 0.3420)
    
    PANEL_TITLES = {
        "Observed versus out-degree-matched expectation",
        "Observed versus out-degree-matched\nexpectation",
        "Reversal rate is not set by out-degree",
        "Reversal rate is not set\nby out-degree",
        "49 above degree-matched expectation",
        "49 regulators above\ndegree-matched expectation",
        "Full-census KEGG ORA",
        "Six one-direction panel-b regulators: all 39 targets",
        "Rest-sign retention",
        "8 h classification",
        "Relative separation peaks at 4–8 h",
        "Both groups rise at 48 h",
        "Target orientation (n=903)",
    }
    
    
    def verify_v8() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V8.items():
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"v8 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def unique_text(fig: plt.Figure, value: str) -> Text:
        matches = [artist for artist in fig.findobj(match=Text) if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def visible_text(fig: plt.Figure) -> list[Text]:
        return [
            artist
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
    
    
    def metadata_axis(fig: plt.Figure) -> mpl.axes.Axes:
        return next(
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if any(artist.get_text() == "rev/test" for artist in axis.texts)
        )
    
    
    def bar_axis(fig: plt.Figure) -> mpl.axes.Axes:
        matches = [
            axis
            for axis in fig.findobj(match=mpl.axes.Axes)
            if len(axis.get_yticklabels()) == 49
            and np.allclose(axis.get_xlim(), (0.0, 300.0), atol=1e-12)
        ]
        if len(matches) != 1:
            raise AssertionError(f"expected one regulator bar axis; found {len(matches)}")
        return matches[0]
    
    
    def set_all_text_roles(fig: plt.Figure) -> dict[str, object]:
        texts = visible_text(fig)
        for artist in texts:
            artist.set_fontsize(SOURCE_5_PT)
    
        title_artists: list[Text] = []
        for artist in texts:
            if artist.get_text() in PANEL_TITLES:
                artist.set_fontsize(SOURCE_7_PT)
                title_artists.append(artist)
            elif artist.get_text() in set("abcdefghi") and artist.get_fontweight() == "bold":
                artist.set_fontsize(SOURCE_7_PT)
    
        axis_labels: list[Text] = []
        for axis in fig.findobj(match=mpl.axes.Axes):
            for artist in (axis.xaxis.label, axis.yaxis.label):
                if artist.get_visible() and artist.get_text().strip():
                    artist.set_fontsize(SOURCE_6_PT)
                    axis_labels.append(artist)
    
        delivered = [artist.get_fontsize() * V5.DELIVERED_SCALE for artist in visible_text(fig)]
        allowed = (5.0, 6.0, 7.0)
        off_contract = [value for value in delivered if not any(np.isclose(value, target, atol=1e-9) for target in allowed)]
        if off_contract:
            raise AssertionError(f"visible text outside exact 5/6/7 contract: {off_contract[:10]}")
        return {
            "visible_text_artists": len(delivered),
            "panel_title_artists": len(title_artists),
            "axis_title_artists": len(axis_labels),
            "delivered_values_pt": sorted(set(round(value, 9) for value in delivered)),
        }
    
    
    def tighten_metadata(fig: plt.Figure) -> dict[str, object]:
        axis = metadata_axis(fig)
        axis.set_axes_locator(None)
        axis.set_position(META_BOUNDS)
        header_rev = next(artist for artist in axis.texts if artist.get_text() == "rev/test")
        header_rate = next(artist for artist in axis.texts if artist.get_text() == "rate %")
        header_fold = next(artist for artist in axis.texts if artist.get_text() == "fold ×")
        header_lg = next(artist for artist in axis.texts if artist.get_text() == "L/G")
        header_rate.set_text("rate")
        header_fold.set_text("fold")
        for artist, x in (
            (header_rev, 0.00),
            (header_rate, 0.55),
            (header_fold, 0.77),
            (header_lg, 0.920),
        ):
            artist.set_x(x)
            artist.set_fontsize(SOURCE_5_PT)
            artist.set_fontstretch("condensed")
    
        for row_index in range(49):
            artists = [
                artist
                for artist in axis.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            if len(artists) != 6:
                raise AssertionError(f"metadata row {row_index}: expected 6 fields, found {len(artists)}")
            rev, rate, fold, lost, slash, gained = artists
            rev.set_x(0.00)
            rate.set_x(0.55)
            fold.set_x(0.77)
            lost.set_x(0.905)
            slash.set_x(0.920)
            gained.set_x(0.935)
            for artist in artists:
                artist.set_fontsize(SOURCE_5_PT)
                artist.set_fontstretch("condensed")
    
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        gaps: list[float] = []
        for row_index in range(49):
            artists = [
                artist
                for artist in axis.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            boxes = [artist.get_window_extent(renderer=renderer) for artist in artists]
            gaps.extend((right.x0 - left.x1) * 72.0 / fig.dpi for left, right in zip(boxes[:-1], boxes[1:]))
        delivered_sizes = {
            round(artist.get_fontsize() * V5.DELIVERED_SCALE, 9)
            for artist in axis.texts
            if artist.get_visible() and artist.get_text().strip()
        }
        if delivered_sizes != {5.0}:
            raise AssertionError(f"metadata not exactly 5 pt delivered: {delivered_sizes}")
        if min(gaps) < 0.35:
            raise AssertionError(f"metadata fields overlap: minimum gap {min(gaps):.3f} pt")
        return {
            "bounds": list(META_BOUNDS),
            "headers": [header_rev.get_text(), header_rate.get_text(), header_fold.get_text(), header_lg.get_text()],
            "delivered_font_pt": 5.0,
            "minimum_horizontal_gap_pt": min(gaps),
        }
    
    
    def restore_regulator_labels(fig: plt.Figure) -> dict[str, object]:
        axis = bar_axis(fig)
        axis.set_axes_locator(None)
        axis.set_position(BAR_BOUNDS)
        labels = axis.get_yticklabels()
        if len(labels) != 49:
            raise AssertionError(f"expected 49 regulator labels, found {len(labels)}")
        observed: dict[str, str] = {}
        for label in labels:
            gene = label.get_text()
            category = BASE.ENRICHED_CLASS.get(gene, "other / unclassified")
            colour = BASE._enriched_label_colour(category)
            label.set_color(colour)
            label.set_fontweight("normal")
            label.set_fontsize(SOURCE_5_PT)
            observed[gene] = colour
        forbidden = {mpl.colors.to_hex(V1.LOST), mpl.colors.to_hex(V1.GAINED)}
        collisions = sorted({mpl.colors.to_hex(value) for value in observed.values()} & forbidden)
        if collisions:
            raise AssertionError(f"functional labels collide with lost/gained colours: {collisions}")
    
        return {
            "labels": len(labels),
            "fontweight": "normal",
            "delivered_font_pt": 5.0,
            "functional_colours": observed,
            "lost_gained_colour_collisions": collisions,
            "compact_key_in_artwork": False,
            "key_owner": "figure legend; omitted from artwork to preserve panel-d width",
        }
    
    
    def widen_kegg(fig: plt.Figure, kegg_plot: pd.DataFrame) -> dict[str, object]:
        axis = V2.find_unique_axis_with_text(fig, "pathway")
        old = list(float(value) for value in axis.get_position().bounds)
        axis.set_axes_locator(None)
        axis.set_position(KEGG_BOUNDS)
        for row_index, (_, row) in enumerate(kegg_plot.iterrows()):
            matches = [
                artist
                for artist in axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.475)
            ]
            if len(matches) != 1:
                raise AssertionError(f"KEGG row {row_index}: pathway label inventory {len(matches)}")
            matches[0].set_text("\n".join(textwrap.wrap(str(row["kegg_display_name"]), width=23)))
        return {"old_bounds": old, "new_bounds": list(KEGG_BOUNDS), "width_gain_fraction": KEGG_BOUNDS[2] / old[2] - 1.0}
    
    
    def reflow_top(fig: plt.Figure) -> dict[str, object]:
        axis_a = next(axis for axis in fig.axes if axis.get_title() == "Observed versus out-degree-matched expectation")
        axis_c = next(axis for axis in fig.axes if axis.get_title() == "Reversal rate is not set by out-degree")
        old_a = tuple(float(value) for value in axis_a.get_position().bounds)
        old_c = tuple(float(value) for value in axis_c.get_position().bounds)
        axis_a.set_position((old_a[0], old_a[1], PANEL_A_WIDTH, old_a[3]))
        axis_c.set_position((old_c[0], old_c[1], PANEL_C_WIDTH, old_c[3]))
        axis_a.set_title("Observed versus out-degree-matched\nexpectation")
        axis_a.set_xlabel("expected reversing edges\n(out-degree matched, log scale)")
        axis_c.set_title("Reversal rate is not set\nby out-degree")
        axis_c.set_xlabel("out-degree (both-significant edges,\nlog scale)")
    
        parent = V2.find_parent_axis(fig)
        title_b = unique_text(fig, "49 above degree-matched expectation")
        subtitle_b = unique_text(fig, "2,075/14,434 rev/test; BH q<.05; ranks 1–49")
        title_d = unique_text(fig, "Full-census KEGG ORA")
        subtitle_d = unique_text(fig, "2,577 reversal targets / 9,554 eligible")
        statistic_d = unique_text(fig, "14 Fisher/BH q<.10\nglobal p=.184; maxT 0/113")
        title_b.set_text("49 regulators above\ndegree-matched expectation")
        title_b.set_linespacing(0.90)
        tags = {
            value: next(
                artist
                for artist in parent.texts
                if artist.get_text() == value and artist.get_fontweight() == "bold"
            )
            for value in ("b", "d")
        }
        b_center = (BAR_BOUNDS[0] + META_BOUNDS[0] + META_BOUNDS[2]) / 2.0
        d_center = KEGG_BOUNDS[0] + KEGG_BOUNDS[2] / 2.0
        for artist in (title_b, subtitle_b):
            artist.set_x(V5.page_to_parent_x(parent, b_center))
        for artist in (title_d, subtitle_d, statistic_d):
            artist.set_x(V5.page_to_parent_x(parent, d_center))
        tags["b"].set_x(V5.page_to_parent_x(parent, BAR_BOUNDS[0] - 0.017))
        tags["d"].set_x(V5.page_to_parent_x(parent, KEGG_BOUNDS[0] - 0.013))
        return {
            "panel_a": {"old": list(old_a), "new": list(axis_a.get_position().bounds)},
            "panel_c": {"old": list(old_c), "new": list(axis_c.get_position().bounds)},
            "panel_b_center_page_x": b_center,
            "panel_d_center_page_x": d_center,
        }
    
    
    def geometry_audit(fig: plt.Figure, functional: dict[str, object]) -> dict[str, object]:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        all_text = visible_text(fig)
        page = fig.bbox
        overflow = []
        tick_text_ids = {
            id(artist)
            for axis in fig.findobj(match=mpl.axes.Axes)
            for artist in (*axis.get_xticklabels(), *axis.get_yticklabels())
        }
        wholly_offpage_ticks = 0
        for artist in all_text:
            box = artist.get_window_extent(renderer=renderer)
            if not mpl.transforms.Bbox.overlaps(box, page):
                if id(artist) not in tick_text_ids:
                    overflow.append((artist.get_text(), tuple(float(value) for value in box.bounds)))
                else:
                    wholly_offpage_ticks += 1
                continue
            if box.x0 < page.x0 - 0.5 or box.x1 > page.x1 + 0.5 or box.y0 < page.y0 - 0.5 or box.y1 > page.y1 + 0.5:
                overflow.append((artist.get_text(), tuple(float(value) for value in box.bounds)))
        if overflow:
            raise AssertionError(f"visible text clips page: {overflow[:5]}")
    
        parent = V2.find_parent_axis(fig)
        ordered = [
            next(artist for artist in parent.texts if artist.get_text() == "b" and artist.get_fontweight() == "bold"),
            unique_text(fig, "49 regulators above\ndegree-matched expectation"),
            next(artist for artist in parent.texts if artist.get_text() == "d" and artist.get_fontweight() == "bold"),
            unique_text(fig, "Full-census KEGG ORA"),
        ]
        boxes = [artist.get_window_extent(renderer=renderer) for artist in ordered]
        heading_gaps = [(right.x0 - left.x1) * 72.0 / fig.dpi for left, right in zip(boxes[:-1], boxes[1:])]
        if min(heading_gaps) < 2.0:
            raise AssertionError(f"top headings collide: {heading_gaps}")
    
        axis_a = next(axis for axis in fig.axes if axis.get_title() == "Observed versus out-degree-matched\nexpectation")
        axis_c = next(axis for axis in fig.axes if axis.get_title() == "Reversal rate is not set\nby out-degree")
        labels = bar_axis(fig).get_yticklabels()
        label_left = min(label.get_window_extent(renderer=renderer).x0 for label in labels)
        ac_right = max(axis_a.bbox.x1, axis_c.bbox.x1)
        ac_to_b_gap = (label_left - ac_right) * 72.0 / fig.dpi
        if ac_to_b_gap < 2.0:
            raise AssertionError(f"a/c to panel-b label gap only {ac_to_b_gap:.3f} pt")
    
        return {
            "page_text_overflow_count": len(overflow),
            "wholly_offpage_clipped_tick_count": wholly_offpage_ticks,
            "top_heading_minimum_gap_pt": min(heading_gaps),
            "panel_ac_to_b_label_gap_pt": ac_to_b_gap,
            "functional_key_in_artwork": False,
            "functional_key_owner": functional["key_owner"],
        }
    
    
    def apply_v9(fig: plt.Figure, kegg_plot: pd.DataFrame, v8_audit: dict[str, object]) -> dict[str, object]:
        reflow = reflow_top(fig)
        functional = restore_regulator_labels(fig)
        metadata = tighten_metadata(fig)
        kegg = widen_kegg(fig, kegg_plot)
        typography = set_all_text_roles(fig)
        geometry = geometry_audit(fig, functional)
        return {
            "predecessor_v8_audit": v8_audit,
            "presentation_only": True,
            "reflow": reflow,
            "panel_b_bar_bounds": list(BAR_BOUNDS),
            "panel_b_metadata": metadata,
            "panel_b_regulator_labels": functional,
            "panel_d": kegg,
            "typography": typography,
            "geometry": geometry,
            "scientific_values_and_panels_preserved": True,
        }
    
    
    def build_figure(regulators: pd.DataFrame, kegg_plot: pd.DataFrame, coverage: dict[str, object]):
        fig, v8_audit = V8.build_figure(regulators, kegg_plot, coverage)
        return fig, apply_v9(fig, kegg_plot, v8_audit)
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true")
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST, OUT_PROVENANCE)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to replace v9 outputs: {existing}")
        canonical_before = V8.current_canonical()
        v8_before = verify_v8()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        metadata = {
            "Title": "Figure 2 signed-KEGG candidate v9 compact main",
            "Creator": "CoDEG_Tcell deterministic Matplotlib renderer",
            "CreationDate": None,
            "ModDate": None,
        }
        V1.atomic_figure(fig, OUT_PDF, format="pdf", facecolor="white", metadata=metadata)
        V1.atomic_figure(fig, OUT_PNG, format="png", dpi=V1.PNG_DPI, facecolor="white")
        plt.close(fig)
        canonical_after = V8.current_canonical()
        v8_after = verify_v8()
        if canonical_after != canonical_before or v8_after != v8_before:
            raise AssertionError("canonical or v8 changed during v9 render")
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v9-compact-main",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "canonical_figure_2": {"before": canonical_before, "after": canonical_after, "written": False},
            "predecessor_v8": {"before": v8_before, "after": v8_after, "preserved": True},
            "audit": audit,
            "scientific_scope": {
                "panel_b": {"regulators": 49, "reversals": 2_075, "tested": 14_434, "lost": 1_154, "gained": 921, "fields": ["rev/test", "rate", "fold", "L/G"]},
                "panel_d": {"reversal_targets": 2_577, "eligible_targets": 9_554, "fisher_bh_terms": 14, "global_p": 0.18439, "maxT": "0/113"},
            },
            "renderer_closure": {
                path: V1.sha256(PAPER / path)
                for path in (
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v9.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v8.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py",
                    "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py",
                    *V1.RENDERER_PATHS,
                )
            },
            "scientific_source_closure": {V1.relative(path): V1.file_record(path) for path in V1.SCIENTIFIC_SOURCE_PATHS},
            "bundle_sources": {"panel_b": V1.file_record(OUT_REGULATORS), "panel_d": V1.file_record(OUT_KEGG)},
            "outputs": {"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG)},
            "claim_boundary": "The 14 KEGG rows are conventional Fisher/BH sensitivity results; global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms do not support pathway enrichment or activity.",
            "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        V1.atomic_text(
            OUT_PROVENANCE,
            "# Figure 2 signed-KEGG candidate v9 provenance\n\n"
            "- Noncanonical review candidate; no canonical or manuscript write.\n"
            "- Presentation-only successor to v8; all panels, rows, and scientific values are preserved.\n"
            "- Panels a/c are modestly narrower; b shifts left, its numeric table is compact, and panel d uses the recovered width.\n"
            "- All regulator labels use the canonical functional colours at regular weight; the key is legend-owned because a 5-pt in-artwork key did not fit without taking width from panel d.\n"
            "- Delivered typography is exact: panel titles/letters 7 pt, axis titles 6 pt, and all other visible text 5 pt.\n"
            "- Panel d visibly retains the full ORA disclosure: 14 Fisher/BH q<.10; global p=.184; maxT 0/113.\n"
            "- See manifest.json for hashes, source closure, and geometry audits.\n",
        )
        print(json.dumps({"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG), "manifest": V1.relative(OUT_MANIFEST), "audit": {"metadata_min_gap_pt": audit["panel_b_metadata"]["minimum_horizontal_gap_pt"], "kegg_width_gain_fraction": audit["panel_d"]["width_gain_fraction"], "typography": audit["typography"]["delivered_values_pt"], "geometry": audit["geometry"]}, "canonical_unchanged": True}, sort_keys=True))
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v9, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v9()

def init_fig2_signed_kegg_side_by_side_candidate_v10():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render the collision-cleared v10 successor of Figure 2 candidate v9.
    
    Only visual geometry changes from v9: larger vertical separation between
    panels a/c, an opaque panel-a legend, a wider h/i gutter, vertical clearance
    for the panel-b heading, and separated columns in the widened KEGG panel.
    Scientific data, rows, colours, statistics, and the exact delivered 5/6/7 pt
    typography contract are unchanged.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V9 = fig2_signed_kegg_side_by_side_candidate_v9
    
    
    V8 = V9.V8
    V5 = V9.V5
    V2 = V9.V2
    V1 = V9.V1
    
    OUTDIR = PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v10_collision_clear"
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v10_collision_clear"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    V9_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v9.py",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/Fig2_signed_kegg_side_by_side_candidate_v9_compact_main.pdf",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/Fig2_signed_kegg_side_by_side_candidate_v9_compact_main.png",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/source_panel_b_regulators.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/source_panel_d_signed_kegg.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/manifest.json",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v9_compact_main/PROVENANCE.md",
    )
    EXPECTED_V9 = {
        V9_PATHS[0]: "98e5af803ac5931b180718144a69023da37ccd8b8a125dfb063245264db4952b",
        V9_PATHS[1]: "d9b6397164d730d82b239dad3cb5745ae48f671fc56860a3d0c8c1a9856e8cc8",
        V9_PATHS[2]: "4d5196184efd5c6f6ce2d7302aa8283572d586909ebfd631d776c0db171c719b",
        V9_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V9_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
        V9_PATHS[5]: "7112db4a27201be002e104a61bff2f64c2f43e7da792372b215de4baed160bf7",
        V9_PATHS[6]: "bff004ed5ab31c6fe874868a9d8079c44e6d0e92e2b2416fe59af7b243dfab52",
    }
    
    
    def verify_v9() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V9.items():
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"v9 drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def unique_text(fig: plt.Figure, value: str) -> Text:
        matches = [artist for artist in fig.findobj(match=Text) if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def find_axis(fig: plt.Figure, title: str) -> mpl.axes.Axes:
        matches = [axis for axis in fig.axes if axis.get_title() == title]
        if len(matches) != 1:
            raise AssertionError(f"expected one axis titled {title!r}; found {len(matches)}")
        return matches[0]
    
    
    def apply_v10(fig: plt.Figure, kegg_plot: pd.DataFrame, v9_audit: dict[str, object]) -> dict[str, object]:
        axis_a = find_axis(fig, "Observed versus out-degree-matched\nexpectation")
        axis_c = find_axis(fig, "Reversal rate is not set\nby out-degree")
        axis_h = find_axis(fig, "Both groups rise at 48 h")
        axis_i = find_axis(fig, "Target orientation (n=903)")
    
        old_positions = {
            "a": list(float(value) for value in axis_a.get_position().bounds),
            "c": list(float(value) for value in axis_c.get_position().bounds),
            "i": list(float(value) for value in axis_i.get_position().bounds),
        }
        # Keep panel-a top fixed while shortening it; move and shorten c so both
        # inter-panel gutters grow without changing their data limits.
        axis_a.set_position((0.0800, 0.7928, V9.PANEL_A_WIDTH, 0.1692))
        axis_c.set_position((0.0800, 0.4983, V9.PANEL_C_WIDTH, 0.2000))
        # Keep panel-i right edge fixed while moving its left edge right.
        old_i = old_positions["i"]
        axis_i.set_position((0.7322, old_i[1], 0.9800 - 0.7322, old_i[3]))
    
        legend_a = axis_a.get_legend()
        if legend_a is None:
            raise AssertionError("panel-a legend missing")
        frame = legend_a.get_frame()
        frame.set_visible(True)
        frame.set_facecolor("white")
        frame.set_edgecolor("none")
        frame.set_alpha(0.96)
    
        parent = V2.find_parent_axis(fig)
        title_b = unique_text(fig, "49 regulators above\ndegree-matched expectation")
        subtitle_b = unique_text(fig, "2,075/14,434 rev/test; BH q<.05; ranks 1–49")
        title_b.set_y(1.025)
        title_b.set_linespacing(0.90)
        tag_b = next(
            artist
            for artist in parent.texts
            if artist.get_text() == "b" and artist.get_fontweight() == "bold"
        )
        tag_b.set_y(1.008)
    
        kegg_axis = V2.find_unique_axis_with_text(fig, "pathway")
        header_pathway = next(artist for artist in kegg_axis.texts if artist.get_text() == "pathway" and float(artist.get_position()[1]) < 0)
        header_t = next(artist for artist in kegg_axis.texts if artist.get_text() == "T" and float(artist.get_position()[1]) < 0)
        header_lg = next(artist for artist in kegg_axis.texts if artist.get_text() == "L/G" and float(artist.get_position()[1]) < 0)
        header_q = next(artist for artist in kegg_axis.texts if artist.get_text() == "ORA\nq" and float(artist.get_position()[1]) < 0)
        header_pathway.set_x(0.490)
        header_t.set_x(0.5525)
        header_lg.set_x(0.594)
        header_lg.set_ha("left")
        header_q.set_x(0.746)
        header_q.set_ha("left")
    
        tiles = [
            patch
            for patch in kegg_axis.patches
            if isinstance(patch, Rectangle)
            and np.isclose(float(patch.get_x()), 0.505)
            and np.isclose(float(patch.get_width()), 0.055)
        ]
        if len(tiles) != 14:
            raise AssertionError(f"expected 14 KEGG tiles, found {len(tiles)}")
        for tile in tiles:
            tile.set_x(0.525)
    
        for row_index in range(14):
            pathway = next(
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.475)
            )
            count = next(
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.505)
            )
            q_value = next(
                artist
                for artist in kegg_axis.texts
                if int(round(float(artist.get_position()[1]))) == row_index
                and np.isclose(float(artist.get_position()[0]), 0.685)
            )
            pathway.set_x(0.490)
            count.set_x(0.594)
            q_value.set_x(0.746)
    
        for collection in kegg_axis.collections:
            if not isinstance(collection, PathCollection):
                continue
            offsets = np.asarray(collection.get_offsets(), dtype=float)
            if offsets.size == 0:
                continue
            moved = offsets.copy()
            moved[np.isclose(moved[:, 0], 0.945), 0] = 0.963
            moved[np.isclose(moved[:, 0], 0.990), 0] = 0.982
            collection.set_offsets(moved)
    
        typography = V9.set_all_text_roles(fig)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    
        def vertical_gap(upper: Text, lower: Text) -> float:
            upper_box = upper.get_window_extent(renderer=renderer)
            lower_box = lower.get_window_extent(renderer=renderer)
            return float((upper_box.y0 - lower_box.y1) * 72.0 / fig.dpi)
    
        a_c_gap = vertical_gap(axis_a.xaxis.label, axis_c.title)
        b_heading_gap = vertical_gap(title_b, subtitle_b)
        if a_c_gap < 2.0:
            raise AssertionError(f"panel-a xlabel / panel-c title gap {a_c_gap:.3f} pt")
        if b_heading_gap < 2.0:
            raise AssertionError(f"panel-b title / subtitle gap {b_heading_gap:.3f} pt")
    
        i_ylabel_box = axis_i.yaxis.label.get_window_extent(renderer=renderer)
        h_right_gap = float((i_ylabel_box.x0 - axis_h.bbox.x1) * 72.0 / fig.dpi)
        if h_right_gap < 2.0:
            raise AssertionError(f"panel-h / panel-i ylabel gutter {h_right_gap:.3f} pt")
    
        legend_box = legend_a.get_window_extent(renderer=renderer)
        legend_inside = (
            legend_box.x0 >= axis_a.bbox.x0
            and legend_box.x1 <= axis_a.bbox.x1
            and legend_box.y0 >= axis_a.bbox.y0
            and legend_box.y1 <= axis_a.bbox.y1
        )
        if not legend_inside:
            raise AssertionError("panel-a legend does not fit inside narrowed axis")
    
        header_boxes = [artist.get_window_extent(renderer=renderer) for artist in (header_t, header_lg, header_q)]
        header_gaps = [float((right.x0 - left.x1) * 72.0 / fig.dpi) for left, right in zip(header_boxes[:-1], header_boxes[1:])]
        if min(header_gaps) < 2.0:
            raise AssertionError(f"KEGG headers collide: {header_gaps}")
    
        count_q_gaps: list[float] = []
        q_marker_gaps: list[float] = []
        marker_x_by_row: dict[int, float] = {}
        for collection in kegg_axis.collections:
            if not isinstance(collection, PathCollection):
                continue
            for x_value, y_value in np.asarray(collection.get_offsets(), dtype=float):
                if np.isclose(x_value, 0.963) or np.isclose(x_value, 0.982):
                    row = int(round(float(y_value)))
                    x_pixel = float(kegg_axis.transData.transform((x_value, y_value))[0])
                    marker_x_by_row[row] = min(marker_x_by_row.get(row, x_pixel), x_pixel)
        for row_index in range(14):
            count = next(artist for artist in kegg_axis.texts if np.isclose(float(artist.get_position()[0]), 0.594) and int(round(float(artist.get_position()[1]))) == row_index)
            q_value = next(artist for artist in kegg_axis.texts if np.isclose(float(artist.get_position()[0]), 0.746) and int(round(float(artist.get_position()[1]))) == row_index)
            count_box = count.get_window_extent(renderer=renderer)
            q_box = q_value.get_window_extent(renderer=renderer)
            count_q_gaps.append(float((q_box.x0 - count_box.x1) * 72.0 / fig.dpi))
            if row_index not in marker_x_by_row:
                raise AssertionError(f"missing KEGG marker row {row_index}")
            q_marker_gaps.append(float((marker_x_by_row[row_index] - q_box.x1) * 72.0 / fig.dpi))
        if min(count_q_gaps) < 2.0 or min(q_marker_gaps) < 2.0:
            raise AssertionError(f"KEGG row columns collide: count/q={min(count_q_gaps):.3f}, q/marker={min(q_marker_gaps):.3f}")
    
        # Re-run v9's full-page, heading, and a/c-to-b label checks after moving axes.
        inherited_geometry = V9.geometry_audit(
            fig,
            {"key_owner": "figure legend; omitted from artwork to preserve panel-d width"},
        )
        return {
            "predecessor_v9_audit": v9_audit,
            "delta_only": True,
            "old_positions": old_positions,
            "new_positions": {
                "a": list(float(value) for value in axis_a.get_position().bounds),
                "c": list(float(value) for value in axis_c.get_position().bounds),
                "i": list(float(value) for value in axis_i.get_position().bounds),
            },
            "panel_a_legend_opaque_and_inside": legend_inside,
            "panel_a_xlabel_to_c_title_gap_pt": a_c_gap,
            "panel_b_title_to_subtitle_gap_pt": b_heading_gap,
            "panel_h_to_i_ylabel_gap_pt": h_right_gap,
            "kegg_header_minimum_gap_pt": min(header_gaps),
            "kegg_count_to_q_minimum_gap_pt": min(count_q_gaps),
            "kegg_q_to_marker_minimum_gap_pt": min(q_marker_gaps),
            "typography": typography,
            "inherited_geometry": inherited_geometry,
            "scientific_values_and_panels_preserved": True,
        }
    
    
    def build_figure(regulators: pd.DataFrame, kegg_plot: pd.DataFrame, coverage: dict[str, object]):
        fig, v9_audit = V9.build_figure(regulators, kegg_plot, coverage)
        return fig, apply_v10(fig, kegg_plot, v9_audit)
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true")
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST, OUT_PROVENANCE)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to replace v10 outputs: {existing}")
        canonical_before = V8.current_canonical()
        v9_before = verify_v9()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        metadata = {"Title": "Figure 2 signed-KEGG candidate v10 collision clear", "Creator": "CoDEG_Tcell deterministic Matplotlib renderer", "CreationDate": None, "ModDate": None}
        V1.atomic_figure(fig, OUT_PDF, format="pdf", facecolor="white", metadata=metadata)
        V1.atomic_figure(fig, OUT_PNG, format="png", dpi=V1.PNG_DPI, facecolor="white")
        plt.close(fig)
        canonical_after = V8.current_canonical()
        v9_after = verify_v9()
        if canonical_after != canonical_before or v9_after != v9_before:
            raise AssertionError("canonical or v9 changed during v10 render")
        manifest = {
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v10-collision-clear",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "canonical_figure_2": {"before": canonical_before, "after": canonical_after, "written": False},
            "predecessor_v9": {"before": v9_before, "after": v9_after, "preserved": True},
            "audit": audit,
            "scientific_scope": {"panel_b": {"regulators": 49, "reversals": 2_075, "tested": 14_434, "lost": 1_154, "gained": 921, "fields": ["rev/test", "rate", "fold", "L/G"]}, "panel_d": {"reversal_targets": 2_577, "eligible_targets": 9_554, "fisher_bh_terms": 14, "global_p": 0.18439, "maxT": "0/113"}},
            "renderer_closure": {path: V1.sha256(PAPER / path) for path in ("scripts/fig2_signed_kegg_side_by_side_candidate_v10.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v9.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v8.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py", *V1.RENDERER_PATHS)},
            "scientific_source_closure": {V1.relative(path): V1.file_record(path) for path in V1.SCIENTIFIC_SOURCE_PATHS},
            "bundle_sources": {"panel_b": V1.file_record(OUT_REGULATORS), "panel_d": V1.file_record(OUT_KEGG)},
            "outputs": {"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG)},
            "claim_boundary": "The 14 KEGG rows are conventional Fisher/BH sensitivity results; global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms do not support pathway enrichment or activity.",
            "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        }
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        V1.atomic_text(OUT_PROVENANCE, "# Figure 2 signed-KEGG candidate v10 provenance\n\n- Noncanonical review candidate; no canonical or manuscript write.\n- Thin collision-clearing successor to v9; all scientific panels, rows, values, functional colours, and exact delivered 5/6/7 pt typography are preserved.\n- The panel-a legend is opaque; a/c have a larger vertical gutter; i shifts right; b title/subtitle and KEGG columns have audited clearances.\n- The functional-colour key remains legend-owned and is omitted from the artwork to preserve panel-d width.\n- Panel d visibly retains 14 Fisher/BH q<.10; global p=.184; maxT 0/113.\n- See manifest.json for hashes and cross-axis bbox audits.\n")
        print(json.dumps({"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG), "manifest": V1.relative(OUT_MANIFEST), "clearance": {"a_to_c_pt": audit["panel_a_xlabel_to_c_title_gap_pt"], "b_title_to_subtitle_pt": audit["panel_b_title_to_subtitle_gap_pt"], "h_to_i_ylabel_pt": audit["panel_h_to_i_ylabel_gap_pt"], "kegg_header_pt": audit["kegg_header_minimum_gap_pt"], "kegg_count_q_pt": audit["kegg_count_to_q_minimum_gap_pt"], "kegg_q_marker_pt": audit["kegg_q_to_marker_minimum_gap_pt"]}, "canonical_unchanged": True}, sort_keys=True))
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v10, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v10()

def init_fig2_signed_kegg_side_by_side_candidate_v11():
    ROOT = Path(__file__).resolve().parents[1]
    OUT = ROOT / 'figures'
    #!/usr/bin/env python3
    """Render v11 with the final overlap-clearance presentation pass.
    
    The scientific panels and values are inherited unchanged from v10.  This pass
    only clears residual cross-panel collisions, removes the panel-c annotation
    leader from the labelled point cloud, and adds the shared functional-colour key
    requested for panels b and c.
    """
    
    
    
    
    matplotlib.use("Agg")
    
    
    PAPER = ROOT
    
    V10 = fig2_signed_kegg_side_by_side_candidate_v10
    
    
    V9 = V10.V9
    V8 = V10.V8
    V1 = V10.V1
    OUTDIR = PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v11_final_review"
    STEM = "Fig2_signed_kegg_side_by_side_candidate_v11_final_review"
    OUT_PDF = OUTDIR / f"{STEM}.pdf"
    OUT_PNG = OUTDIR / f"{STEM}.png"
    OUT_REGULATORS = OUTDIR / "source_panel_b_regulators.tsv"
    OUT_KEGG = OUTDIR / "source_panel_d_signed_kegg.tsv"
    OUT_MANIFEST = OUTDIR / "manifest.json"
    OUT_PROVENANCE = OUTDIR / "PROVENANCE.md"
    
    V10_PARTIAL_PATHS = (
        PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v10.py",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v10_collision_clear/Fig2_signed_kegg_side_by_side_candidate_v10_collision_clear.pdf",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v10_collision_clear/Fig2_signed_kegg_side_by_side_candidate_v10_collision_clear.png",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v10_collision_clear/source_panel_b_regulators.tsv",
        PAPER / "data/figure_provenance/fig2_signed_kegg_side_by_side_candidate_v10_collision_clear/source_panel_d_signed_kegg.tsv",
    )
    EXPECTED_V10_PARTIAL = {
        V10_PARTIAL_PATHS[0]: "323c974e72c77ca55215b7260c94c5511fba77deef1053320ade7a5c40a23778",
        V10_PARTIAL_PATHS[1]: "8cd70e0452cad60df7951f96e343cdac68ebd7c051888da233fda02e446d10a5",
        V10_PARTIAL_PATHS[2]: "279102e615389e73d34532a3154d7e14ed6b9f18e33faa3e8ee3ea4cd7a0a469",
        V10_PARTIAL_PATHS[3]: "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
        V10_PARTIAL_PATHS[4]: "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
    }
    
    
    def native(value):
        """Convert numpy scalar containers to JSON-native values without changing them."""
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {str(key): native(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [native(item) for item in value]
        return value
    
    
    def verify_v10_partial() -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for path, expected in EXPECTED_V10_PARTIAL.items():
            actual = V1.sha256(path)
            if actual != expected:
                raise AssertionError(f"v10 partial drift: {path}: {actual} != {expected}")
            records[V1.relative(path)] = V1.file_record(path)
        return records
    
    
    def build_figure(regulators: pd.DataFrame, kegg_plot: pd.DataFrame, coverage: dict[str, object]):
        fig, audit = V10.build_figure(regulators, kegg_plot, coverage)
        return fig, native(apply_v11_overlap_clearance(fig, audit))
    
    
    def _visible_text(fig: plt.Figure) -> list[Text]:
        return [
            artist
            for artist in fig.findobj(match=Text)
            if artist.get_visible() and artist.get_text().strip()
        ]
    
    
    def _unique_text(fig: plt.Figure, value: str) -> Text:
        matches = [artist for artist in _visible_text(fig) if artist.get_text() == value]
        if len(matches) != 1:
            raise AssertionError(f"expected one visible text {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def _panel_tag(axis: mpl.axes.Axes, value: str) -> Text:
        matches = [
            artist
            for artist in axis.texts
            if artist.get_visible()
            and artist.get_text() == value
            and artist.get_fontweight() == "bold"
        ]
        if len(matches) != 1:
            raise AssertionError(f"expected one panel tag {value!r}; found {len(matches)}")
        return matches[0]
    
    
    def _boxes_overlap(left: mpl.transforms.Bbox, right: mpl.transforms.Bbox) -> bool:
        return bool(
            min(left.x1, right.x1) > max(left.x0, right.x0)
            and min(left.y1, right.y1) > max(left.y0, right.y0)
        )
    
    
    def apply_v11_overlap_clearance(
        fig: plt.Figure,
        predecessor_audit: dict[str, object],
    ) -> dict[str, object]:
        """Apply and audit the final presentation-only collision fixes."""
        axis_b = V9.bar_axis(fig)
        axis_meta = V9.metadata_axis(fig)
        axis_c = V10.find_axis(fig, "Reversal rate is not set\nby out-degree")
        axis_d = V10.V2.find_unique_axis_with_text(fig, "pathway")
        axis_f = V10.find_axis(fig, "Rest-sign retention")
        axis_g = V10.find_axis(fig, "Relative separation peaks at 4–8 h")
        axis_h = V10.find_axis(fig, "Both groups rise at 48 h")
        axis_i = V10.find_axis(fig, "Target orientation (n=903)")
    
        old_positions = {
            key: list(float(value) for value in axis.get_position().bounds)
            for key, axis in (("meta", axis_meta), ("d", axis_d), ("c", axis_c), ("h", axis_h))
        }
    
        # Pack rev/test, rate, fold and L/G into one compact block and hand the
        # recovered width to panel d.  Text remains 5 pt and every row field stays
        # explicit.
        axis_meta.set_position((0.5170, 0.4950, 0.1740, 0.4220))
        axis_d.set_position((0.7000, 0.5424, 0.2800, 0.3420))
        header_by_text = {
            artist.get_text(): artist
            for artist in axis_meta.texts
            if artist.get_text() in {"rev/test", "rate", "fold", "L/G"}
        }
        for key, x_value in (("rev/test", 0.000), ("rate", 0.540), ("fold", 0.752), ("L/G", 0.910)):
            header_by_text[key].set_x(x_value)
        for row_index in range(49):
            row = [
                artist
                for artist in axis_meta.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            if len(row) != 6:
                raise AssertionError(f"panel-b metadata row {row_index}: found {len(row)} fields")
            for artist, x_value in zip(row, (0.000, 0.540, 0.752, 0.895, 0.910, 0.925)):
                artist.set_x(x_value)
    
        # Keep panel c's upper boundary fixed while opening the gutter to panel e.
        axis_c.set_position((0.0800, 0.5035, V9.PANEL_C_WIDTH, 0.1948))
        # Shift only the left boundary of h; its right edge and all data limits stay fixed.
        old_h = old_positions["h"]
        new_h_left = 0.4860
        axis_h.set_position((new_h_left, old_h[1], old_h[0] + old_h[2] - new_h_left, old_h[3]))
        # Matplotlib retained a 12.5 tick object above panel i's displayed y range;
        # remove only ticks outside the fixed limits so the panel tag can sit beside
        # its title instead of floating in the inter-row gutter.
        i_ymin, i_ymax = axis_i.get_ylim()
        axis_i.set_yticks([
            value for value in axis_i.get_yticks()
            if i_ymin <= float(value) <= i_ymax
        ])
    
        # The wording is unchanged.  Moving the line break leftward keeps long KEGG
        # labels clear of panel b without narrowing either panel.
        kegg_rewrap = {
            "Antigen processing and\npresentation": "Antigen processing\nand presentation",
            "MAPK signaling pathway": "MAPK signaling\npathway",
        }
        for before, after in kegg_rewrap.items():
            artist = _unique_text(fig, before)
            artist.set_text(after)
            artist.set_linespacing(0.92)
    
        # Preserve the named high-out-degree examples while removing the leader that
        # crossed the highlighted point labels.  The dashed baseline is labelled at
        # its quiet left end instead of inside the coloured cluster.
        high_degree = next(
            artist
            for artist in axis_c.texts
            if artist.get_text().startswith("high out-degree")
        )
        high_degree.set_text("high out-degree, baseline rate:\nMED12, ATAD5, SGF29, SUPT20H")
        high_degree.set_position((0.98, 0.97))
        high_degree.set_ha("right")
        high_degree.set_va("top")
        high_degree.set_linespacing(0.94)
        high_degree.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.90, pad=0.35))
        if getattr(high_degree, "arrow_patch", None) is not None:
            high_degree.arrow_patch.set_visible(False)
    
        baseline = _unique_text(fig, "baseline 5.2%")
        baseline.set_position((1.22, 5.55))
        baseline.set_ha("left")
        baseline.set_va("bottom")
        baseline.set_bbox(dict(facecolor="white", edgecolor="none", alpha=0.88, pad=0.20))
    
        # One compact legend explains the label colours in both b and c.  It sits in
        # the unused right side of the short-bar rows and therefore costs panel d no
        # width.
        class_handles = [
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                markersize=3.5,
                color=V9.BASE._enriched_label_colour(category),
                label=V9.BASE.ENRICHED_CLASS_SHORT[category],
            )
            for category in V9.BASE.ENRICHED_CLASS_ORDER
        ]
        functional_legend = axis_b.legend(
            handles=class_handles,
            loc="lower left",
            bbox_to_anchor=(0.120, 0.330),
            ncol=1,
            fontsize=V9.SOURCE_5_PT,
            handlelength=0.25,
            handletextpad=0.25,
            labelspacing=0.12,
            borderaxespad=0.0,
            borderpad=0.20,
            frameon=True,
            framealpha=0.94,
            facecolor="white",
            edgecolor="none",
            title="b,c label colour\nfunctional class",
            title_fontsize=V9.SOURCE_5_PT,
            alignment="left",
        )
    
        # Move lower-row panel letters into their gutters.  These are figure labels,
        # not data or headings.
        tag_positions = {
            "f": (-0.17, 1.08),
            "g": (-0.13, 1.13),
            "h": (-0.195, 1.13),
            "i": (-0.09, 1.14),
        }
        tag_axes = {"f": axis_f, "g": axis_g, "h": axis_h, "i": axis_i}
        for value, position in tag_positions.items():
            _panel_tag(tag_axes[value], value).set_position(position)
    
        typography = V9.set_all_text_roles(fig)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    
        def bbox(artist: Text) -> mpl.transforms.Bbox:
            return artist.get_window_extent(renderer=renderer)
    
        def vertical_gap(upper: Text, lower: Text) -> float:
            return float((bbox(upper).y0 - bbox(lower).y1) * 72.0 / fig.dpi)
    
        # Panel c / e gutter.
        title_e = _unique_text(fig, "Six one-direction panel-b regulators: all 39 targets")
        c_e_gap = vertical_gap(axis_c.xaxis.label, title_e)
        if c_e_gap < 2.0:
            raise AssertionError(f"panel-c xlabel / panel-e title gap {c_e_gap:.3f} pt")
    
        # Panel b numbers / panel d pathway names.  Only vertically overlapping rows
        # can collide, so compare those exact displayed boxes.
        meta_artists = [artist for artist in axis_meta.texts if artist.get_visible() and artist.get_text().strip()]
        pathway_artists = [
            artist
            for artist in axis_d.texts
            if artist.get_visible()
            and artist.get_text().strip()
            and np.isclose(float(artist.get_position()[0]), 0.490)
        ]
        b_d_gaps: list[float] = []
        for left_artist in meta_artists:
            left_box = bbox(left_artist)
            for right_artist in pathway_artists:
                right_box = bbox(right_artist)
                if min(left_box.y1, right_box.y1) > max(left_box.y0, right_box.y0):
                    b_d_gaps.append(float((right_box.x0 - left_box.x1) * 72.0 / fig.dpi))
        if not b_d_gaps or min(b_d_gaps) < 2.0:
            raise AssertionError(f"panel-b metadata / panel-d pathway gap {min(b_d_gaps):.3f} pt")
    
        metadata_row_gaps: list[float] = []
        for row_index in range(49):
            row = [
                artist
                for artist in axis_meta.texts
                if np.isclose(float(artist.get_position()[1]), float(row_index))
            ]
            row_boxes = [bbox(artist) for artist in row]
            metadata_row_gaps.extend(
                float((right.x0 - left.x1) * 72.0 / fig.dpi)
                for left, right in zip(row_boxes[:-1], row_boxes[1:])
            )
        if min(metadata_row_gaps) < 0.35:
            raise AssertionError(f"panel-b compact metadata gap {min(metadata_row_gaps):.3f} pt")
    
        # Shared legend containment and separation from all visible bar rectangles.
        legend_box = functional_legend.get_window_extent(renderer=renderer)
        legend_inside = bool(
            legend_box.x0 >= axis_b.bbox.x0
            and legend_box.x1 <= axis_b.bbox.x1
            and legend_box.y0 >= axis_b.bbox.y0
            and legend_box.y1 <= axis_b.bbox.y1
        )
        if not legend_inside:
            raise AssertionError("panel-b/c functional legend is not inside panel b")
        legend_bar_overlaps = 0
        for patch in axis_b.patches:
            patch_box = patch.get_window_extent(renderer=renderer)
            if patch.get_visible() and _boxes_overlap(legend_box, patch_box):
                legend_bar_overlaps += 1
        if legend_bar_overlaps:
            raise AssertionError(f"functional legend overlaps {legend_bar_overlaps} reversal bars")
    
        # Lower-panel tags must clear their own titles and visible tick labels.
        lower_tag_gaps: dict[str, float] = {}
        tag_tick_overlaps: dict[str, int] = {}
        for value, axis in tag_axes.items():
            tag = _panel_tag(axis, value)
            title = axis.title
            tag_box = bbox(tag)
            title_box = bbox(title)
            lower_tag_gaps[value] = float((title_box.x0 - tag_box.x1) * 72.0 / fig.dpi)
            if _boxes_overlap(tag_box, title_box) or lower_tag_gaps[value] < 2.0:
                raise AssertionError(f"panel-{value} tag / title gap {lower_tag_gaps[value]:.3f} pt")
            ticks = [*axis.get_xticklabels(), *axis.get_yticklabels()]
            tag_tick_overlaps[value] = sum(
                _boxes_overlap(tag_box, bbox(tick))
                for tick in ticks
                if tick.get_visible() and tick.get_text().strip()
            )
            if tag_tick_overlaps[value]:
                raise AssertionError(f"panel-{value} tag overlaps tick labels")
    
        g_h_gap = float((bbox(axis_h.yaxis.label).x0 - bbox(axis_g.title).x1) * 72.0 / fig.dpi)
        if g_h_gap < 2.0:
            raise AssertionError(f"panel-g title / panel-h ylabel gap {g_h_gap:.3f} pt")
    
        # The panel-c note now has no visible arrow and clears every highlighted gene label.
        highlighted = {
            "CCNC", "AHR", "ARNT", "NCKAP1L", "DOCK2", "CBFB", "STAT5B", "CREBBP", "ARHGAP30"
        }
        highlighted_artists = [artist for artist in axis_c.texts if artist.get_text() in highlighted]
        # Annotation.get_window_extent can retain the hidden arrow's historical
        # envelope; call Text's implementation to audit the visible note itself.
        note_text_box = Text.get_window_extent(high_degree, renderer=renderer)
        note_label_overlaps = sum(_boxes_overlap(note_text_box, bbox(artist)) for artist in highlighted_artists)
        if note_label_overlaps:
            raise AssertionError(f"panel-c high-degree note overlaps {note_label_overlaps} highlighted labels")
        arrow_visible = bool(
            getattr(high_degree, "arrow_patch", None) is not None
            and high_degree.arrow_patch.get_visible()
        )
        if arrow_visible:
            raise AssertionError("panel-c high-degree leader remains visible")
    
        inherited_geometry = V9.geometry_audit(
            fig,
            {
                "key_owner": "shared in-artwork panel-b legend for panels b and c",
            },
        )
        inherited_geometry["functional_key_in_artwork"] = True
        inherited_geometry["functional_key_owner"] = "shared panel-b legend for panels b and c"
        return {
            "predecessor_v10_audit": predecessor_audit,
            "predecessor_v9_audit": predecessor_audit["predecessor_v9_audit"],
            "presentation_only": True,
            "delta_only": True,
            "old_positions": old_positions,
            "new_positions": {
                "meta": list(float(value) for value in axis_meta.get_position().bounds),
                "d": list(float(value) for value in axis_d.get_position().bounds),
                "c": list(float(value) for value in axis_c.get_position().bounds),
                "h": list(float(value) for value in axis_h.get_position().bounds),
            },
            "panel_b_metadata_minimum_internal_gap_pt": min(metadata_row_gaps),
            "panel_d_width_gain_from_v10_fraction": (
                axis_d.get_position().width / old_positions["d"][2] - 1.0
            ),
            "kegg_rewrap": kegg_rewrap,
            "panel_b_to_d_pathway_minimum_gap_pt": min(b_d_gaps),
            "panel_c_xlabel_to_e_title_gap_pt": c_e_gap,
            "panel_g_title_to_h_ylabel_gap_pt": g_h_gap,
            "panel_a_xlabel_to_c_title_gap_pt": predecessor_audit["panel_a_xlabel_to_c_title_gap_pt"],
            "panel_b_title_to_subtitle_gap_pt": predecessor_audit["panel_b_title_to_subtitle_gap_pt"],
            "panel_h_to_i_ylabel_gap_pt": predecessor_audit["panel_h_to_i_ylabel_gap_pt"],
            "kegg_header_minimum_gap_pt": predecessor_audit["kegg_header_minimum_gap_pt"],
            "kegg_count_to_q_minimum_gap_pt": predecessor_audit["kegg_count_to_q_minimum_gap_pt"],
            "kegg_q_to_marker_minimum_gap_pt": predecessor_audit["kegg_q_to_marker_minimum_gap_pt"],
            "lower_panel_tag_to_title_gaps_pt": lower_tag_gaps,
            "lower_panel_tag_tick_overlap_counts": tag_tick_overlaps,
            "shared_functional_legend": {
                "present": True,
                "applies_to_panels": ["b", "c"],
                "categories": len(class_handles),
                "inside_panel_b": legend_inside,
                "bar_overlap_count": legend_bar_overlaps,
                "delivered_font_pt": 5.0,
            },
            "panel_c_annotation": {
                "high_degree_examples_retained": ["MED12", "ATAD5", "SGF29", "SUPT20H"],
                "arrow_visible": arrow_visible,
                "highlighted_label_overlap_count": note_label_overlaps,
                "baseline_label_retained": True,
            },
            "typography": typography,
            "inherited_geometry": inherited_geometry,
            "scientific_values_and_panels_preserved": True,
        }
    
    
    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--force", action="store_true")
        args = parser.parse_args()
        outputs = (OUT_PDF, OUT_PNG, OUT_REGULATORS, OUT_KEGG, OUT_MANIFEST, OUT_PROVENANCE)
        existing = [path for path in outputs if path.exists()]
        if existing and not args.force:
            raise FileExistsError(f"refusing to replace v11 outputs: {existing}")
        canonical_before = V8.current_canonical()
        v9_before = V10.verify_v9()
        v10_before = verify_v10_partial()
        regulators = V1.load_regulators()
        kegg_plot, coverage = V1.KEGG.load_plot_data()
        fig, audit = build_figure(regulators, kegg_plot, coverage)
        OUTDIR.mkdir(parents=True, exist_ok=True)
        V1.atomic_tsv(regulators, OUT_REGULATORS)
        V1.atomic_tsv(kegg_plot, OUT_KEGG)
        metadata = {"Title": "Figure 2 signed-KEGG candidate v11 final review", "Creator": "CoDEG_Tcell deterministic Matplotlib renderer", "CreationDate": None, "ModDate": None}
        V1.atomic_figure(fig, OUT_PDF, format="pdf", facecolor="white", metadata=metadata)
        V1.atomic_figure(fig, OUT_PNG, format="png", dpi=V1.PNG_DPI, facecolor="white")
        plt.close(fig)
        canonical_after = V8.current_canonical()
        v9_after = V10.verify_v9()
        v10_after = verify_v10_partial()
        if canonical_after != canonical_before or v9_after != v9_before or v10_after != v10_before:
            raise AssertionError("canonical or predecessor changed during v11 render")
        manifest = native({
            "schema_version": "fig2-signed-kegg-side-by-side-candidate-v11-final-review",
            "status": "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "delta_from_v10": "presentation-only overlap clearance, compact panel-b metadata, wider panel d, and a shared b/c functional-colour legend; no data, statistics, rows, panels, or typography change",
            "canonical_figure_2": {"before": canonical_before, "after": canonical_after, "written": False},
            "predecessor_v9": {"before": v9_before, "after": v9_after, "preserved": True},
            "predecessor_v10_partial": {"before": v10_before, "after": v10_after, "preserved": True},
            "audit": audit,
            "scientific_scope": {"panel_b": {"regulators": 49, "reversals": 2_075, "tested": 14_434, "lost": 1_154, "gained": 921, "fields": ["rev/test", "rate", "fold", "L/G"]}, "panel_d": {"reversal_targets": 2_577, "eligible_targets": 9_554, "fisher_bh_terms": 14, "global_p": 0.18439, "maxT": "0/113"}},
            "renderer_closure": {path: V1.sha256(PAPER / path) for path in ("scripts/fig2_signed_kegg_side_by_side_candidate_v11.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v10.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v9.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v8.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v7.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v5.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v4.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v3.py", "scripts/fig2_signed_kegg_side_by_side_candidate_v2.py", *V1.RENDERER_PATHS)},
            "scientific_source_closure": {V1.relative(path): V1.file_record(path) for path in V1.SCIENTIFIC_SOURCE_PATHS},
            "bundle_sources": {"panel_b": V1.file_record(OUT_REGULATORS), "panel_d": V1.file_record(OUT_KEGG)},
            "outputs": {"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG)},
            "claim_boundary": "The 14 KEGG rows are conventional Fisher/BH sensitivity results; global opportunity-aware p=0.18439 and 0/113 ORA maxT-supported terms do not support pathway enrichment or activity.",
            "software": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        })
        V1.atomic_text(OUT_MANIFEST, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        V1.atomic_text(OUT_PROVENANCE, "# Figure 2 signed-KEGG candidate v11 provenance\n\n- Noncanonical review candidate; no canonical or manuscript write.\n- Presentation-only successor to v10: residual overlaps are cleared, panel-b metadata is compacted, and the recovered width is assigned to panel d.\n- Exact delivered typography: panel titles/letters 7 pt, axis titles 6 pt, all other visible text 5 pt.\n- All 49 regular-weight regulator labels use canonical functional colours; one compact in-artwork legend explains the colours used in panels b and c.\n- Panel c retains the high-out-degree examples and baseline label without a leader crossing the highlighted point labels.\n- Panel d visibly reports 14 Fisher/BH q<.10; global p=.184; maxT 0/113.\n- See manifest.json for source closure, hashes, and explicit cross-panel bbox audits.\n")
        print(json.dumps({"pdf": V1.file_record(OUT_PDF), "png": V1.file_record(OUT_PNG), "manifest": V1.file_record(OUT_MANIFEST), "sources": {"panel_b": V1.file_record(OUT_REGULATORS), "panel_d": V1.file_record(OUT_KEGG)}, "renderer": V1.file_record(PAPER / "scripts/fig2_signed_kegg_side_by_side_candidate_v11.py"), "canonical_unchanged": True}, sort_keys=True))
    
    
    for k, v in list(locals().items()):
        setattr(fig2_signed_kegg_side_by_side_candidate_v11, k, v)
init_fig2_signed_kegg_side_by_side_candidate_v11()


if __name__ == '__main__':
    V11 = fig2_signed_kegg_side_by_side_candidate_v11
    V10 = fig2_signed_kegg_side_by_side_candidate_v10
    V9 = fig2_signed_kegg_side_by_side_candidate_v9
    V2 = fig2_signed_kegg_side_by_side_candidate_v2
    V1 = fig2_signed_kegg_side_by_side_candidate_v1
    KEGG = fig_signed_kegg_transition_heatmap_1d_exploratory
    
    regulators = V1.load_regulators()
    kegg_plot, coverage = KEGG.load_plot_data()
    fig, audit = V11.build_figure(regulators, kegg_plot, coverage)
    
    axis_b = V9.bar_axis(fig)
    axis_meta = V9.metadata_axis(fig)
    axis_c = V10.find_axis(fig, "Reversal rate is not set\nby out-degree")
    axis_d = V2.find_unique_axis_with_text(fig, "pathway")
    axis_f = V10.find_axis(fig, "Rest-sign retention")
    axis_g = V10.find_axis(fig, "Relative separation peaks at 4–8 h")
    axis_h = V10.find_axis(fig, "Both groups rise at 48 h")
    axis_i = V10.find_axis(fig, "Target orientation (n=903)")
    axis_c.set_position((0.0800, 0.5035, V9.PANEL_C_WIDTH, 0.1948))
    
    typography = V9.set_all_text_roles(fig)
    
    for ext in ("pdf", "png"):
        fig.savefig(V1.OUT / f"Fig2.{ext}", dpi=V1.PNG_DPI, facecolor="white")
    print("wrote Fig2")
