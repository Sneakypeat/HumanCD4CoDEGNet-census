#!/usr/bin/env python
"""Figure v2: reversal pattern overlaid on the CORRECTED, promoted excess-burden enrichment
result (Supplementary Results R7 / Figure S5), not the superseded exploratory v1.

Reads results_data/census/gsea_pathway_reversal_overlay.json and asserts its anchors; computes
nothing. Scope: the 9 KEGG Legacy pathways meeting the repeatable rule and the 10 GO Biological
Process terms additionally surviving the strict global maxT test in both permutation banks
(global_familywise_stable) -- the two sets R7 itself names as its headline results.

Panel a  each pathway's members split into lost / gained / mixed orientation, ranked within
         library by q, with members carrying no reversal shown as remainder.
Panel b  targets reversal-positive in 4+ of the 19 pathways: which pathways, and each target's
         own orientation.

Writes figures_final/diagnostics/gsea_pathway_reversal_overlay_v2/
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results_data/census/gsea_pathway_reversal_overlay.json"
OUT = ROOT / "figures_final/diagnostics/gsea_pathway_reversal_overlay_v2"
STEM = "gsea_pathway_reversal_overlay_v2"

LOST, GAINED, MIXED, ZERO = "#c1272d", "#2a9d8f", "#8a6fb0", "#e3e6e9"
INK, MUTED, GRID = "#22262B", "#666B70", "#c7ccd1"
HEAD_PT, LABEL_PT, SMALL_PT = 7.0, 7.0, 7.0
FIG_W, FIG_H = 7.20, 11.80

SHORT = {
    "Proteasome": "Proteasome*",
    "Cytokine Cytokine Receptor Interaction": "Cytokine-cytokine receptor interaction",
    "P53 Signaling Pathway": "p53 signalling",
    "Apoptosis": "Apoptosis",
    "Antigen Processing And Presentation": "Antigen processing and presentation",
    "Jak Stat Signaling Pathway": "JAK-STAT signalling",
    "Prion Diseases": "Prion diseases*",
    "Thyroid Cancer": "Thyroid cancer",
    "Chronic Myeloid Leukemia": "Chronic myeloid leukaemia",
    "Mitotic Cytokinesis": "Mitotic cytokinesis",
    "Apoptotic Process": "Apoptotic process",
    "Organelle Organization": "Organelle organisation",
    "Regulation Of Type II Interferon Production": "Regulation of type II IFN production",
    "Positive Regulation Of Apoptotic Process": "Positive reg. of apoptotic process",
    "Proteasome-Mediated Ubiquitin-Dependent Protein Catabolic Process": "Proteasome-mediated ubiquitin catabolism",
    "Intrinsic Apoptotic Signaling Pathway In Response To Endoplasmic Reticulum Stress": "Intrinsic apoptotic signalling, ER stress",
    "Cellular Response To Lipopolysaccharide": "Cellular response to LPS",
    "Cellular Response To Cytokine Stimulus": "Cellular response to cytokine",
    "Intrinsic Apoptotic Signaling Pathway": "Intrinsic apoptotic signalling",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pretty(term: str) -> str:
    return SHORT.get(term, term)


def load() -> dict:
    d = json.loads(SRC.read_text())
    assert len(d["pathways"]) == 19
    assert sum(1 for r in d["pathways"] if r["library"] == "KEGG_Legacy") == 9
    assert sum(1 for r in d["pathways"] if r["library"].startswith("GO_Biological")) == 10
    for r in d["pathways"]:
        assert r["n_lost"] + r["n_gained"] + r["n_mixed"] == r["n_members_with_reversal"]
    return d


def ordered_rows(d: dict) -> list[dict]:
    kegg = sorted((r for r in d["pathways"] if r["library"] == "KEGG_Legacy"), key=lambda r: r["q"])
    go = sorted((r for r in d["pathways"] if r["library"].startswith("GO_Biological")),
                key=lambda r: r["q"])
    return kegg + go


def panel_pathways(ax, d: dict) -> None:
    rows = ordered_rows(d)
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        n_mem = r["n_members_in_universe"]
        silent = n_mem - r["n_members_with_reversal"]
        left = 0.0
        for val, col in ((r["n_lost"], LOST), (r["n_gained"], GAINED), (r["n_mixed"], MIXED),
                          (silent, ZERO)):
            frac = val / n_mem
            ax.barh(yy, frac, left=left, height=0.62, color=col, edgecolor="white", lw=0.4)
            left += frac
        ax.text(1.05, yy, f"n={n_mem:<3d}  NES {r['NES']:+.2f}  q={r['q']:.2f}",
                fontsize=SMALL_PT, color=INK, va="center", ha="left", family="monospace")
    ax.set_yticks(y)
    labels = [f"{i + 1}. {pretty(r['term'])}" for i, r in enumerate(rows)]
    ax.set_yticklabels(labels, fontsize=LABEL_PT, color=INK)
    n_kegg = sum(1 for r in rows if r["library"] == "KEGG_Legacy")
    div_y = y[n_kegg - 1] - 0.5
    ax.axhline(div_y, color=INK, lw=1.0)
    ax.set_xlim(0, 1.62)
    ax.set_ylim(y.min() - 0.6, y.max() + 0.9)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0%", "50%", "100%"], fontsize=SMALL_PT)
    ax.set_xlabel("share of pathway members", fontsize=LABEL_PT, color=INK)
    ax.tick_params(length=0, pad=2)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.legend(handles=[
        Rectangle((0, 0), 1, 1, facecolor=LOST, label="lost on activation"),
        Rectangle((0, 0), 1, 1, facecolor=GAINED, label="gained on activation"),
        Rectangle((0, 0), 1, 1, facecolor=MIXED, label="mixed across regulators"),
        Rectangle((0, 0), 1, 1, facecolor=ZERO, edgecolor=GRID, label="no reversal")],
        loc="lower left", bbox_to_anchor=(0.0, 1.012), ncol=2, fontsize=SMALL_PT, frameon=False,
        handlelength=1.1, columnspacing=1.1, borderpad=0.2, labelspacing=0.3)


def panel_driver_matrix(ax, d: dict) -> None:
    drivers = d["cross_pathway_drivers"]
    targets = sorted(drivers, key=lambda t: (-drivers[t]["n_pathways"], t))
    rows = ordered_rows(d)
    terms = [r["term"] for r in rows]
    col_of = {t: i for i, t in enumerate(terms)}

    for yi, t in enumerate(targets):
        info = drivers[t]
        col = {"lost_on_activation": LOST, "gained_on_activation": GAINED,
               "mixed": MIXED}[info["orientation"]]
        hatch = "////" if info["orientation"] == "gained_on_activation" else None
        for term in info["pathways"]:
            xi = col_of[term]
            ax.add_patch(Rectangle((xi + 0.06, yi + 0.06), 0.88, 0.88, facecolor=col,
                                   edgecolor="white", lw=0.4, hatch=hatch))
    for xi in range(len(terms) + 1):
        ax.axvline(xi, color=GRID, lw=0.4)
    n_kegg = sum(1 for r in rows if r["library"] == "KEGG_Legacy")
    ax.axvline(n_kegg, color=INK, lw=1.0)
    ax.set_xlim(0, len(terms))
    ax.set_ylim(0, len(targets))
    ax.set_xticks([i + 0.5 for i in range(len(terms))])
    ax.set_xticklabels([str(i + 1) for i in range(len(terms))], fontsize=LABEL_PT)
    ax.set_yticks([i + 0.5 for i in range(len(targets))])
    # yi=0 is placed at the top by invert_yaxis(), so labels must remain in
    # the same order as the tile rows. Reversing them here mislabels every row.
    ax.set_yticklabels(targets, fontsize=SMALL_PT)
    ax.invert_yaxis()
    ax.tick_params(length=0, pad=2)
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(handles=[
        Rectangle((0, 0), 1, 1, facecolor=LOST, label="lost"),
        Rectangle((0, 0), 1, 1, facecolor=GAINED, hatch="////", label="gained"),
        Rectangle((0, 0), 1, 1, facecolor=MIXED, label="mixed")],
        loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3, fontsize=SMALL_PT, frameon=False,
        handlelength=1.1, columnspacing=1.1, borderpad=0.2)


def audit(fig) -> dict:
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    def in_view(ax, label, axis):
        lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
        try:
            v = float(label.get_position()[0 if axis == "x" else 1])
        except (TypeError, ValueError):
            return True
        return min(lo, hi) - 1e-9 <= v <= max(lo, hi) + 1e-9
    items = []
    for ax in fig.axes:
        items += list(ax.texts)
        items += [t for t in ax.get_xticklabels() if in_view(ax, t, "x")]
        items += [t for t in ax.get_yticklabels() if in_view(ax, t, "y")]
        items += [ax.xaxis.label, ax.yaxis.label]
        if ax.get_legend() is not None:
            items += list(ax.get_legend().get_texts())
    items += list(fig.texts)
    boxes, sizes, labels = [], [], []
    for t in items:
        if not t.get_text() or not t.get_visible():
            continue
        b = t.get_window_extent(renderer=r)
        if b.width <= 0 or b.height <= 0:
            continue
        boxes.append(b); sizes.append(t.get_fontsize()); labels.append(t.get_text()[:44])
    W, H = fig.get_window_extent().width, fig.get_window_extent().height
    off = [f"{labels[i]} x0={b.x0:.0f} x1={b.x1:.0f} y0={b.y0:.0f} y1={b.y1:.0f} (W={W:.0f} H={H:.0f})"
           for i, b in enumerate(boxes)
           if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5]
    pad, hits = 0.5, []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            if (a.x0 < b.x1 - pad and b.x0 < a.x1 - pad
                    and a.y0 < b.y1 - pad and b.y0 < a.y1 - pad):
                hits.append(f"{labels[i]} | {labels[j]}")
    return {"n_text": len(boxes), "min_font_pt": round(min(sizes), 2),
            "max_font_pt": round(max(sizes), 2),
            "n_off_page": len(off), "off_page": off,
            "n_collisions": len(hits), "collisions": hits[:10]}


def main() -> None:
    d = load()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")
    ax_a = fig.add_axes([0.510, 0.680, 0.220, 0.230])
    ax_b = fig.add_axes([0.135, 0.075, 0.585, 0.480])
    panel_pathways(ax_a, d)
    panel_driver_matrix(ax_b, d)

    fig.text(0.020, 0.978, "a", fontsize=HEAD_PT, fontweight="bold", color=INK)
    fig.text(0.048, 0.978,
             "Reversal pattern inside the corrected excess-burden\n"
             "enrichment result (9 KEGG statement-eligible + 10 GO\n"
             "globally maxT-stable; Supplementary Results R7;\n"
             "* = gene-set overlap, see caption)",
             fontsize=HEAD_PT, fontweight="bold", color=INK, va="top")
    fig.text(0.020, 0.578, "b", fontsize=HEAD_PT, fontweight="bold", color=INK)
    fig.text(0.058, 0.578,
             f"{len(d['cross_pathway_drivers'])} targets reversal-positive in "
             f"{d['shared_min_pathways']}+ of the 19 pathways (pathway numbers match panel a)",
             fontsize=HEAD_PT, fontweight="bold", color=INK, va="top")

    fig.text(0.020, 0.018,
             "Descriptive overlay on the corrected, independently audited enrichment result\n"
             "(results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2); no new statistic "
             "is computed here.\nThe exploratory v1 result (14 KEGG at within-library q < 0.10, "
             "weaker null and score) is superseded and is\nnot used here. "
             "* KEGG_PRION_DISEASES shares ubiquitin-proteasome genes with KEGG_PROTEASOME (R7).\n"
             "A target's orientation is its own dataset-wide direction across all its regulators, "
             "not restricted to\npathway-member regulators. R7 itself does not test direction.",
             fontsize=SMALL_PT, color=MUTED, va="bottom")

    a = audit(fig)
    assert a["min_font_pt"] >= 5.0 and a["max_font_pt"] <= 7.0, a
    assert a["n_off_page"] == 0, a["off_page"]
    assert a["n_collisions"] == 0, a["collisions"]
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{STEM}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    (OUT / "manifest.json").write_text(json.dumps({
        "_what": "reversal-pattern overlay on the corrected, promoted excess-burden enrichment (v2)",
        "_script": "scripts/fig_gsea_pathway_reversal_overlay.py",
        "_status": "MAIN-FIGURE CANDIDATE, not promoted",
        "inputs": {"gsea_pathway_reversal_overlay.json": sha256(SRC)},
        "audit": a,
    }, indent=2) + "\n")
    print(f"wrote {OUT}/{STEM}.pdf\naudit {a}")


if __name__ == "__main__":
    main()
