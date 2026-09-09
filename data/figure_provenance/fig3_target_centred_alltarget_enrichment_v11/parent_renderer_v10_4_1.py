#!/usr/bin/env python3
"""Build diagnostic Figure 3 v10.4.1 from the audited v10.3 layout.

This readability-only candidate preserves every audited v10.3 panel, datum,
order, scale, legend and short canvas while increasing panel d's x- and y-axis
tick labels from 5.0 pt to exactly 5.5 pt. Panel d's cells increase to 5.5-pt
true squares and panel e narrows so the inherited sans-serif labels remain
collision-free on the same page. The arrangement remains:

* top: smaller square target panels a and b, with the complete R49 x T279
  overview as panel c beside them;
* middle: the literal first-80-column enlargement d beside the frozen GO-BP
  slim descriptive summary, relabelled e; panel e displays the 23 nonzero
  fixed terms plus the outside-slim complement while retaining all 25 fixed
  terms in its bundled source;
* bottom: the former e, f and g panels relabelled f, g and h.

No source values, selected terms, denominators, encodings, tests, claim
boundaries or page geometry change. Run only through
scripts/run_memory_guarded.py. Canonical and previously retained artifacts are
protected by paired before/after hashes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
V9_RENDERER = ROOT / "scripts/fig3_target_centred_biological_summary_v9.py"
V8_RENDERER = ROOT / "scripts/fig3_target_centred_effect_zoom_v8.py"
V7_RENDERER = ROOT / "scripts/fig3_target_centred_effect_zoom_v7.py"
BIO_RENDERER = ROOT / "scripts/fig3c_go_slim_biological_summary_v2.py"
V10_3_RENDERER = ROOT / "scripts/fig3_target_centred_top_overview_layout_v10.py"
V10_3_DIR = (
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_3"
)
V10_4_PREVIEW_DIR = (
    ROOT
    / "figures_final/diagnostics/_archive/"
    "fig3h_target_centred_top_overview_layout_v10_4_failed_header_clip_20260904_1248"
)
DEFAULT_OUT = (
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_4_1"
)
SCHEMA_VERSION = "fig3-target-centred-top-overview-layout-v10.4.1"

FIGURE_WIDTH_IN = 210.0 / 25.4
FIGURE_HEIGHT_IN = 9.60
TOP_SQUARE_IN = 1.34
OVERVIEW_CELL_PT = 0.92
PANEL_D_CELL_PT = 5.5
DETAIL_VERTICAL_SHIFT_IN = 0.25
PANEL_D_HEAT_OFFSET_IN = 0.340
PANEL_D_ROW_TRACK_X_IN = 0.025
PANEL_D_PANEL_E_GAP_IN = 0.040
PANEL_D_INHERITED_AXIS_LABEL_PT = 5.0
PANEL_D_AXIS_LABEL_PT = 5.5
PDF_METADATA_DATE = datetime(2026, 9, 4, tzinfo=timezone.utc)
v9_module: Any = None

PROTECTED_PATHS = (
    ROOT / "figures_final/Fig3.pdf",
    ROOT / "figures_final/Fig3.png",
    ROOT / "MANUSCRIPT.md",
    ROOT / "FIGURE_LEGENDS.md",
    ROOT / "PAPER_FACTS.md",
    ROOT / "HumanCD4CoDEGNet_signflip.docx",
    ROOT / "HumanCD4CoDEGNet_signflip.pdf",
    ROOT / "scripts/active_claim_guard.py",
    V9_RENDERER,
    V8_RENDERER,
    V7_RENDERER,
    BIO_RENDERER,
    V10_3_RENDERER,
    ROOT / "figures_final/diagnostics/fig3_target_centred_effect_zoom_v8_2/manifest.json",
    ROOT / "figures_final/diagnostics/fig3_target_centred_effect_zoom_v8_2/Fig3_target_centred_effect_zoom_v8_2.pdf",
    ROOT / "figures_final/diagnostics/fig3_target_centred_effect_zoom_v8_2/Fig3_target_centred_effect_zoom_v8_2.png",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_biological_summary_v9_1/manifest.json",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_biological_summary_v9_1/Fig3_target_centred_biological_summary_v9_1.pdf",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_biological_summary_v9_1/Fig3_target_centred_biological_summary_v9_1.png",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10/manifest.json",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10/Fig3_target_centred_top_overview_layout_v10.pdf",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10/Fig3_target_centred_top_overview_layout_v10.png",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_1/manifest.json",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_1/Fig3_target_centred_top_overview_layout_v10_1.pdf",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_1/Fig3_target_centred_top_overview_layout_v10_1.png",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_2/manifest.json",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_2/Fig3_target_centred_top_overview_layout_v10_2.pdf",
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_2/Fig3_target_centred_top_overview_layout_v10_2.png",
    V10_3_DIR / "manifest.json",
    V10_3_DIR / "Fig3_target_centred_top_overview_layout_v10_3.pdf",
    V10_3_DIR / "Fig3_target_centred_top_overview_layout_v10_3.png",
    V10_4_PREVIEW_DIR / "manifest.json",
    V10_4_PREVIEW_DIR / "Fig3_target_centred_top_overview_layout_v10_4.pdf",
    V10_4_PREVIEW_DIR / "Fig3_target_centred_top_overview_layout_v10_4.png",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def protected_hashes() -> dict[str, str]:
    missing = [relative(path) for path in PROTECTED_PATHS if not path.is_file()]
    if missing:
        raise FileNotFoundError("protected artifact missing: " + ", ".join(missing))
    return {relative(path): sha256(path) for path in PROTECTED_PATHS}


def draw_overview_panel_c(
    parent: mpl.axes.Axes,
    fig: mpl.figure.Figure,
    v7: Any,
    bridge: Any,
    dense: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, mpl.axes.Axes]]:
    state = dense["state"]
    score = dense["score"]
    rows = dense["rows"]
    targets = dense["targets"]
    if state.shape != (49, 279) or score.shape != state.shape:
        raise RuntimeError("panel c complete matrix shape drift")

    parent.set_xticks([])
    parent.set_yticks([])
    parent.set_axis_off()
    bridge.panel_tag(parent, "c", x=0.000, y=1.01)
    parent.text(
        0.5, 0.975, "R49 × T279 effect-strength overview",
        transform=parent.transAxes, ha="center", va="top",
        fontsize=7.0, fontweight="bold", color=v7.INK,
    )
    parent.text(
        0.5, 0.120,
        "Frame = columns 1–80, enlarged literally in d; same values, order and scales.",
        transform=parent.transAxes, ha="center", va="bottom",
        fontsize=5.0, color=v7.MUTED,
    )

    slot = parent.get_position()
    slot_x0 = slot.x0 * fig.get_figwidth()
    slot_y0 = slot.y0 * fig.get_figheight()
    slot_right = slot.x1 * fig.get_figwidth()
    heat_width = 279 * OVERVIEW_CELL_PT / 72.0
    heat_height = 49 * OVERVIEW_CELL_PT / 72.0
    row_track_x = slot_x0 + 0.008
    heat_x = row_track_x + v7.ROW_TRACK_PT / 72.0 + v7.ROW_TRACK_GAP_PT / 72.0
    heat_y = slot_y0 + 0.575
    if heat_x + heat_width > slot_right + 1e-9:
        raise RuntimeError(
            f"panel c does not fit top slot: need {heat_x + heat_width:.4f}, "
            f"have {slot_right:.4f} inches"
        )

    heat = v7.add_absolute_axis(fig, heat_x, heat_y, heat_width, heat_height)
    state_cmap = ListedColormap([v7.OUTSIDE, v7.STABLE, v7.OUTSIDE, v7.OUTSIDE, v7.SELF])
    heat.pcolormesh(
        np.arange(280), np.arange(50), state,
        cmap=state_cmap,
        norm=BoundaryNorm(np.arange(-0.5, 5.5), 5),
        shading="flat", antialiased=False, edgecolors=v7.GRID,
        linewidth=0.035, rasterized=False, zorder=1,
    )
    heat.pcolormesh(
        np.arange(280), np.arange(50), np.ma.masked_invalid(score),
        cmap=v7.effect_colormap(), norm=Normalize(-1, 1),
        shading="flat", antialiased=False, rasterized=False, zorder=2,
    )
    heat.set_xlim(0, 279)
    heat.set_ylim(49, 0)
    heat.set_aspect("equal", adjustable="box")
    heat.set_anchor("SW")
    heat.set_xticks([])
    heat.set_yticks([])
    for spine in heat.spines.values():
        spine.set_color("#a7adb2")
        spine.set_linewidth(0.35)
    heat.add_patch(
        Rectangle(
            (0, 0), 80, 49, fill=False, edgecolor=v7.DENSE_FRAME,
            linewidth=0.85, zorder=5, clip_on=True,
        )
    )

    row_track = v7.add_absolute_axis(
        fig, row_track_x, heat_y, v7.ROW_TRACK_PT / 72.0, heat_height
    )
    row_track.pcolormesh(
        [0, 1], np.arange(50), rows["all_total"].to_numpy(dtype=float)[:, None],
        cmap="Blues", norm=Normalize(0, 113), shading="flat", rasterized=False,
    )
    row_track.set_xlim(0, 1)
    row_track.set_ylim(49, 0)
    row_track.set_axis_off()

    gap = v7.TRACK_GAP_PT / 72.0
    track_height = v7.TRACK_PT / 72.0
    orientation_y = heat_y + heat_height + gap
    orientation = v7.add_absolute_axis(fig, heat_x, orientation_y, heat_width, track_height)
    orientation_values = np.where(
        targets["orientation"].eq("lost_on_activation"), 0, 1
    )[None, :]
    orientation.pcolormesh(
        np.arange(280), [0, 1], orientation_values,
        cmap=ListedColormap([v7.LOST, v7.GAINED]),
        norm=BoundaryNorm([-0.5, 0.5, 1.5], 2),
        shading="flat", rasterized=False,
    )
    orientation.set_xlim(0, 279)
    orientation.set_ylim(0, 1)
    orientation.set_axis_off()

    burden_y = orientation_y + track_height + gap
    burden = v7.add_absolute_axis(fig, heat_x, burden_y, heat_width, track_height)
    burden.pcolormesh(
        np.arange(280), [0, 1], targets["n_reversals"].to_numpy(dtype=float)[None, :],
        cmap="Purples", norm=Normalize(3, 18), shading="flat", rasterized=False,
    )
    burden.set_xlim(0, 279)
    burden.set_ylim(0, 1)
    burden.set_axis_off()

    fig.canvas.draw()
    box = heat.get_window_extent()
    cell_width = box.width / fig.dpi * 72.0 / 279.0
    cell_height = box.height / fig.dpi * 72.0 / 49.0
    if abs(cell_width - cell_height) > 0.02 or abs(cell_width - OVERVIEW_CELL_PT) > 0.02:
        raise RuntimeError(f"panel c cell geometry drift: {cell_width} x {cell_height} pt")
    summary = {
        "shape": [49, 279],
        "cell_side_points": float(cell_width),
        "reversal_cells": int(np.isin(state, [2, 3]).sum()),
        "lost_cells": int((state == 2).sum()),
        "gained_cells": int((state == 3).sum()),
        "stable_cells": int((state == 1).sum()),
        "outside_cells": int((state == 0).sum()),
        "self_cells": int((state == 4).sum()),
        "regulators_with_reversal": int(np.isin(state, [2, 3]).any(axis=1).sum()),
        "targets_with_reversal": int(np.isin(state, [2, 3]).any(axis=0).sum()),
        "row_labels_shown": False,
        "target_labels_shown": False,
        "dendrogram": False,
        "ordering": "frozen dense-tail-v2 display_rank; no clustering",
    }
    return summary, {
        "heat": heat,
        "row_track": row_track,
        "orientation_track": orientation,
        "burden_track": burden,
    }


def measure_heat_axes(
    fig: mpl.figure.Figure,
    axes: dict[str, dict[str, mpl.axes.Axes]],
) -> dict[str, Any]:
    fig.canvas.draw()
    expected = {
        "c": (49, 279, OVERVIEW_CELL_PT),
        "d": (49, 80, PANEL_D_CELL_PT),
    }
    result: dict[str, Any] = {}
    for letter, (nrows, ncols, wanted) in expected.items():
        box = axes[letter]["heat"].get_window_extent()
        width_pt = box.width / fig.dpi * 72.0
        height_pt = box.height / fig.dpi * 72.0
        cell_width = width_pt / ncols
        cell_height = height_pt / nrows
        if abs(cell_width - cell_height) > 0.02 or abs(cell_width - wanted) > 0.02:
            raise RuntimeError(
                f"panel {letter} saved geometry drift: {cell_width:.5f} x {cell_height:.5f} pt"
            )
        result[letter] = {
            "matrix_width_points": width_pt,
            "matrix_height_points": height_pt,
            "cell_width_points": cell_width,
            "cell_height_points": cell_height,
        }
    return result


def set_panel_d_axis_label_typography(
    axes_d: dict[str, mpl.axes.Axes],
) -> dict[str, Any]:
    """Change only panel d's target/regulator tick-label size to 5.5 pt."""
    heat = axes_d["heat"]
    x_labels = heat.get_xticklabels()
    y_labels = heat.get_yticklabels()
    if len(x_labels) != 80 or len(y_labels) != 49:
        raise RuntimeError(
            "panel d axis-label count drift: "
            f"x={len(x_labels)}, y={len(y_labels)}"
        )
    if any(not label.get_text() for label in x_labels + y_labels):
        raise RuntimeError("panel d contains an empty x/y tick label")
    before_sizes = [float(label.get_fontsize()) for label in x_labels + y_labels]
    before_families = sorted(
        {tuple(label.get_fontfamily()) for label in x_labels + y_labels}
    )
    if any(
        abs(size - PANEL_D_INHERITED_AXIS_LABEL_PT) > 1e-8
        for size in before_sizes
    ):
        raise RuntimeError(
            "panel d inherited axis-label typography drift: "
            f"{sorted(set(before_sizes))}"
        )
    x_text_before = [label.get_text() for label in x_labels]
    y_text_before = [label.get_text() for label in y_labels]
    for label in x_labels + y_labels:
        label.set_fontsize(PANEL_D_AXIS_LABEL_PT)
    if x_text_before != [label.get_text() for label in x_labels] or y_text_before != [
        label.get_text() for label in y_labels
    ]:
        raise RuntimeError("panel d axis-label text changed while setting typography")
    after_sizes = [float(label.get_fontsize()) for label in x_labels + y_labels]
    if any(abs(size - PANEL_D_AXIS_LABEL_PT) > 1e-8 for size in after_sizes):
        raise RuntimeError("panel d axis-label typography did not reach exactly 5.5 pt")
    return {
        "x_tick_label_count": len(x_labels),
        "y_tick_label_count": len(y_labels),
        "inherited_font_size_points": PANEL_D_INHERITED_AXIS_LABEL_PT,
        "delivered_font_size_points": PANEL_D_AXIS_LABEL_PT,
        "inherited_and_delivered_font_families": [list(value) for value in before_families],
        "x_label_text_sha256": hash_bytes("\n".join(x_text_before).encode("utf-8")),
        "y_label_text_sha256": hash_bytes("\n".join(y_text_before).encode("utf-8")),
        "text_unchanged": True,
    }


def audit_panel_d_axis_label_typography(
    axes_d: dict[str, mpl.axes.Axes],
    expected: dict[str, Any],
) -> dict[str, Any]:
    """Verify the 5.5-pt setting survives clamping, drawing and saving."""
    heat = axes_d["heat"]
    x_labels = heat.get_xticklabels()
    y_labels = heat.get_yticklabels()
    observed = {
        "x_tick_label_count": len(x_labels),
        "y_tick_label_count": len(y_labels),
        "font_sizes_points": sorted(
            {float(label.get_fontsize()) for label in x_labels + y_labels}
        ),
        "font_families": sorted(
            {tuple(label.get_fontfamily()) for label in x_labels + y_labels}
        ),
        "x_label_text_sha256": hash_bytes(
            "\n".join(label.get_text() for label in x_labels).encode("utf-8")
        ),
        "y_label_text_sha256": hash_bytes(
            "\n".join(label.get_text() for label in y_labels).encode("utf-8")
        ),
    }
    if observed["x_tick_label_count"] != expected["x_tick_label_count"] or observed[
        "y_tick_label_count"
    ] != expected["y_tick_label_count"]:
        raise RuntimeError(f"panel d axis-label count changed: {observed}")
    if observed["font_sizes_points"] != [PANEL_D_AXIS_LABEL_PT]:
        raise RuntimeError(f"panel d axis labels are not exactly 5.5 pt: {observed}")
    if observed["font_families"] != [
        tuple(value) for value in expected["inherited_and_delivered_font_families"]
    ]:
        raise RuntimeError(f"panel d axis-label font-family drift: {observed}")
    if observed["x_label_text_sha256"] != expected["x_label_text_sha256"] or observed[
        "y_label_text_sha256"
    ] != expected["y_label_text_sha256"]:
        raise RuntimeError("panel d axis-label text drifted after typography change")
    observed["exact_5_5_points"] = True
    observed["text_unchanged"] = True
    return observed


def draw_panel_e_nonzero_go_slim(
    fig: mpl.figure.Figure,
    v7: Any,
    bridge: Any,
    source_rows: list[dict[str, Any]],
    *,
    x_in: float,
    y_in: float,
    width_in: float,
    height_in: float,
) -> tuple[dict[str, Any], dict[str, mpl.axes.Axes]]:
    """Draw panel e while retaining zero rows only in bundled provenance."""
    if len(source_rows) != 26:
        raise RuntimeError(f"panel e expects 26 frozen source rows, observed {len(source_rows)}")
    omitted = [
        row for row in source_rows
        if row["row_type"] == "official_bp_slim"
        and int(row["lost_targets"]) == 0
        and int(row["gained_targets"]) == 0
    ]
    omitted_identity = [(str(row["go_id"]), str(row["go_term"])) for row in omitted]
    expected_omitted = [
        ("GO:0071554", "cell wall organization or biogenesis"),
        ("GO:0019748", "secondary metabolic process"),
    ]
    if omitted_identity != expected_omitted:
        raise RuntimeError(f"panel e zero-row identity drift: {omitted_identity}")
    rows = [row for row in source_rows if row not in omitted]
    if len(rows) != 24:
        raise RuntimeError(f"panel e display row count drift: {len(rows)}")
    official_nonzero = [row for row in rows if row["row_type"] == "official_bp_slim"]
    complements = [row for row in rows if row["row_type"] == "explicit_complement"]
    if len(official_nonzero) != 23 or len(complements) != 1:
        raise RuntimeError(
            f"panel e expected 23 nonzero fixed terms plus one complement: "
            f"{len(official_nonzero)}, {len(complements)}"
        )
    if (
        int(complements[0]["lost_targets"]),
        int(complements[0]["gained_targets"]),
        int(complements[0]["total_targets"]),
    ) != (22, 17, 39):
        raise RuntimeError("panel e outside-slim complement counts drifted")

    panel = v7.add_absolute_axis(fig, x_in, y_in, width_in, height_in)
    panel.set_xticks([])
    panel.set_yticks([])
    panel.set_axis_off()
    bridge.panel_tag(panel, "e", x=-0.035, y=1.01)
    panel.text(
        0.5, 0.988, "GO-BP-slim summary",
        transform=panel.transAxes, ha="center", va="top",
        fontsize=7.0, fontweight="bold", color=v7.INK,
    )
    panel.text(
        0.5, 0.948, "240/279 mapped; 39 outside this BP slim",
        transform=panel.transAxes, ha="center", va="top",
        fontsize=5.0, color=v7.MUTED,
    )
    panel.text(
        0.5, 0.910, "23 nonzero fixed GO-BP-slim terms",
        transform=panel.transAxes, ha="center", va="top",
        fontsize=5.0, color=v7.MUTED,
    )
    panel.text(
        0.5, 0.878, "+ outside-slim complement;\nmemberships overlap",
        transform=panel.transAxes, ha="center", va="top",
        fontsize=5.0, color=v7.MUTED,
    )
    panel.text(
        0.5, 0.835, "shade = % within orientation; number = n",
        transform=panel.transAxes, ha="center", va="top",
        fontsize=5.0, color=v7.INK,
    )

    table_height = len(rows) * v9_module.PANEL_H_ROW_PT / 72.0
    # The two 5-pt, two-line orientation headers need at least 0.60 inches to
    # remain measurably disjoint in the delivered PDF.
    table_width = min(0.62, width_in * 0.320)
    table_x = x_in + width_in - table_width - 0.025
    table_top = y_in + height_in * 0.770
    table_y = table_top - table_height
    if table_y <= y_in + 1.25:
        raise RuntimeError("panel e table leaves insufficient room for the c/d key")
    heat = v7.add_absolute_axis(fig, table_x, table_y, table_width, table_height)

    percentages = np.array(
        [[row["lost_percent"], row["gained_percent"]] for row in rows], dtype=float
    )
    counts = np.array(
        [[row["lost_targets"], row["gained_targets"]] for row in rows], dtype=int
    )
    rgba = np.empty((len(rows), 2, 4), dtype=float)
    for row_index in range(len(rows)):
        for column_index, color in enumerate((v7.LOST, v7.GAINED)):
            value = percentages[row_index, column_index]
            intensity = 0.0 if value == 0 else 0.12 + 0.88 * value / 40.0
            rgba[row_index, column_index] = v9_module.mix_with_white(color, intensity)
    heat.imshow(rgba, interpolation="none", aspect="auto", origin="upper")
    heat.set_xlim(-0.5, 1.5)
    heat.set_ylim(len(rows) - 0.5, -0.5)
    heat.set_xticks([0, 1])
    heat.set_xticklabels(["lost\nn=161", "gained\nn=118"], fontsize=5.0)
    heat.xaxis.tick_top()
    heat.tick_params(axis="x", length=0, pad=1.5)
    heat.get_xticklabels()[0].set_color(v7.LOST)
    heat.get_xticklabels()[1].set_color(v7.GAINED)

    labels = [str(row["go_term"]) for row in rows]
    display_contractions = {
        "generation of precursor metabolites and energy": "precursor metabolites / energy",
    }
    labels = [display_contractions.get(label, label) for label in labels]
    labels[-1] = "outside this BP slim"
    heat.set_yticks(np.arange(len(rows)))
    heat.set_yticklabels(
        labels, fontsize=5.0, fontfamily="DejaVu Sans", fontstretch="condensed"
    )
    heat.tick_params(axis="y", length=0, pad=1.6)
    heat.get_yticklabels()[-1].set_fontstyle("italic")
    heat.get_yticklabels()[-1].set_color(v7.MUTED)
    for row_index in range(len(rows)):
        for column_index in range(2):
            heat.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5), 1.0, 1.0,
                    fill=False, edgecolor="#d5dadd", linewidth=0.35, clip_on=True,
                )
            )
            text_color = "#ffffff" if percentages[row_index, column_index] >= 22 else v7.INK
            heat.text(
                column_index, row_index, str(counts[row_index, column_index]),
                ha="center", va="center", fontsize=5.0,
                fontweight="bold" if percentages[row_index, column_index] >= 22 else "normal",
                color=text_color,
            )
    heat.axhline(len(rows) - 1.5, color="#7f878d", linewidth=0.60)
    for spine in heat.spines.values():
        spine.set_color("#9da4aa")
        spine.set_linewidth(0.40)

    key_top = table_y - 0.08
    key_height = 1.30
    key_bottom = key_top - key_height
    if key_bottom < y_in + 0.04:
        raise RuntimeError(
            f"compact c/d key extends below panel e: {key_bottom - y_in:.3f} inches"
        )
    if key_height < 1.05:
        raise RuntimeError(f"c/d key height too small: {key_height:.3f} inches")
    key = v9_module.draw_panel_d_key(
        fig, v7,
        x_in=x_in + 0.035,
        y_in=key_bottom,
        width_in=width_in - 0.070,
        height_in=key_height,
    )
    orientation_labels = [
        artist for artist in key.texts if artist.get_text() in {"lost", "gained"}
    ]
    if len(orientation_labels) != 2:
        raise RuntimeError("compact c/d key orientation-label identity drift")
    for artist in orientation_labels:
        artist.set_y(0.525)
    orientation_swatches = [
        patch for patch in key.patches
        if isinstance(patch, Rectangle) and abs(float(patch.get_y()) - 0.535) < 1e-9
    ]
    if len(orientation_swatches) != 2:
        raise RuntimeError("compact c/d key orientation-swatch identity drift")
    for patch in orientation_swatches:
        patch.set_y(0.495)
    return (
        {
            "source_rows": 26,
            "displayed_rows": 24,
            "source_official_bp_slim_terms": 25,
            "displayed_nonzero_official_bp_slim_terms": 23,
            "explicit_complement_rows": 1,
            "omitted_zero_rows": [
                {"go_id": go_id, "go_term": go_term}
                for go_id, go_term in omitted_identity
            ],
            "targets": 279,
            "mapped_targets": 240,
            "outside_this_slim_targets": 39,
            "lost_denominator": 161,
            "gained_denominator": 118,
            "memberships_overlap": True,
            "shade": "within-orientation percentage on a fixed 0--40% scale",
            "printed_cell_value": "unique-target count",
            "teal": v7.GAINED,
            "row_pitch_points": v9_module.PANEL_H_ROW_PT,
            "display_label_contractions": display_contractions,
            "claim_boundary": (
                "descriptive frozen GO-BP-slim membership; zero-total fixed terms "
                "omitted from display only; no enrichment, lost-versus-gained test, "
                "pathway activity, or mechanism"
            ),
        },
        {"parent": panel, "heat": heat, "panel_d_key": key},
    )


def overlap_dimensions(first: Any, second: Any) -> tuple[float, float]:
    return (
        max(0.0, min(first.x1, second.x1) - max(first.x0, second.x0)),
        max(0.0, min(first.y1, second.y1) - max(first.y0, second.y0)),
    )


def is_effectively_rendered_text(artist: mpl.text.Text) -> bool:
    if not artist.get_visible() or not artist.get_text().strip():
        return False
    axis = artist.axes
    if axis is None or axis.axison:
        return True
    # Axis-off parents still render explicit annotations and titles, but not
    # their latent default tick labels.
    return (
        artist in axis.texts
        or artist is axis.title
        or artist is getattr(axis, "_left_title", None)
        or artist is getattr(axis, "_right_title", None)
    )


def audit_all_text_bboxes(fig: mpl.figure.Figure) -> dict[str, Any]:
    """Fail on any measurable overlap among text artists actually rendered."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    artists = [
        artist for artist in fig.findobj(mpl.text.Text)
        if is_effectively_rendered_text(artist)
    ]
    boxes = [artist.get_window_extent(renderer=renderer) for artist in artists]
    collisions: list[dict[str, Any]] = []
    for left_index, left in enumerate(artists):
        for right_index in range(left_index + 1, len(artists)):
            width, height = overlap_dimensions(boxes[left_index], boxes[right_index])
            if width <= 0.50 or height <= 0.50:
                continue
            collisions.append(
                {
                    "first": left.get_text(),
                    "second": artists[right_index].get_text(),
                    "overlap_width_pixels": width,
                    "overlap_height_pixels": height,
                }
            )
    if collisions:
        raise RuntimeError(f"rendered text bbox collisions: {collisions[:12]}")
    return {
        "rendered_text_artists": len(artists),
        "pairwise_overlap_threshold_pixels_each_axis": 0.50,
        "collision_count": 0,
        "collisions": [],
    }


def audit_reported_text_relationships(
    fig: mpl.figure.Figure,
    axes_e: dict[str, mpl.axes.Axes],
    ax_f: mpl.axes.Axes,
    ax_g: mpl.axes.Axes,
) -> dict[str, Any]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    explanatory = [
        artist for artist in axes_e["parent"].texts
        if artist.get_text() == "shade = % within orientation; number = n"
    ]
    if len(explanatory) != 1:
        raise RuntimeError(f"panel e explanatory-line identity drift: {len(explanatory)}")
    orientation_headers = list(axes_e["heat"].get_xticklabels())
    if [header.get_text() for header in orientation_headers] != [
        "lost\nn=161", "gained\nn=118"
    ]:
        raise RuntimeError("panel e lost/gained header identity drift")

    relationships: list[dict[str, Any]] = []

    def record(label: str, first: mpl.text.Text, second: mpl.text.Text) -> None:
        width, height = overlap_dimensions(
            first.get_window_extent(renderer=renderer),
            second.get_window_extent(renderer=renderer),
        )
        relationships.append(
            {
                "relationship": label,
                "overlap_width_pixels": width,
                "overlap_height_pixels": height,
            }
        )
        if width > 0.50 and height > 0.50:
            raise RuntimeError(
                f"targeted text collision {label}: {width:.3f} x {height:.3f} pixels"
            )

    for header in orientation_headers:
        record(f"panel e explanatory vs {header.get_text()!r}", explanatory[0], header)

    f_headers = [artist for artist in ax_f.texts if artist.get_text() in {"CCNC", "TADA2B", "SUPT7L"}]
    g_headers = [artist for artist in ax_g.texts if artist.get_text() in {"AHR", "ARNT"}]
    if len(f_headers) != 3 or len(g_headers) != 2:
        raise RuntimeError(
            f"bottom regulator-header identity drift: f={len(f_headers)}, g={len(g_headers)}"
        )
    for header in f_headers:
        record(f"panel f title vs {header.get_text()}", ax_f.title, header)
    for header in g_headers:
        record(f"panel g title vs {header.get_text()}", ax_g.title, header)

    g_tags = [artist for artist in ax_g.texts if artist.get_text() == "g"]
    if len(g_tags) != 1:
        raise RuntimeError(f"panel g tag identity drift: {len(g_tags)}")
    record("panel g tag vs title", g_tags[0], ax_g.title)
    return {"threshold_pixels_each_axis": 0.50, "relationships": relationships}


def audit_panel_d_header_clearance(
    fig: mpl.figure.Figure,
    ax_c_parent: mpl.axes.Axes,
    ax_de_parent: mpl.axes.Axes,
    axes_d: dict[str, mpl.axes.Axes],
) -> dict[str, Any]:
    """Keep the raised panel-d header between c's footnote and d's top track."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def exactly_one(artists: list[mpl.text.Text], label: str) -> mpl.text.Text:
        if len(artists) != 1:
            raise RuntimeError(f"panel d header audit expected one {label}, found {len(artists)}")
        return artists[0]

    footnote = exactly_one(
        [
            artist for artist in ax_c_parent.texts
            if artist.get_text().startswith("Frame = columns 1–80")
        ],
        "panel-c footnote",
    )
    tag = exactly_one(
        [artist for artist in ax_de_parent.texts if artist.get_text() == "d"],
        "panel-d tag",
    )
    title = exactly_one(
        [
            artist for artist in ax_de_parent.texts
            if artist.get_text().startswith("Enlargement of framed c")
        ],
        "panel-d title",
    )
    subtitle = exactly_one(
        [
            artist for artist in ax_de_parent.texts
            if artist.get_text().startswith("R49 × first-80 intersection")
        ],
        "panel-d subtitle",
    )
    burden_box = axes_d["burden_track"].get_window_extent(renderer=renderer)
    footnote_box = footnote.get_window_extent(renderer=renderer)
    header_boxes = {
        "tag": tag.get_window_extent(renderer=renderer),
        "title": title.get_window_extent(renderer=renderer),
        "subtitle": subtitle.get_window_extent(renderer=renderer),
    }
    lower_clearances = {
        name: float(box.y0 - burden_box.y1) for name, box in header_boxes.items()
    }
    upper_clearances = {
        name: float(footnote_box.y0 - box.y1) for name, box in header_boxes.items()
    }
    if min(lower_clearances.values()) <= 0.50:
        raise RuntimeError(f"panel d header intersects top burden track: {lower_clearances}")
    if min(upper_clearances.values()) <= 0.50:
        raise RuntimeError(f"panel d header intersects panel c footnote: {upper_clearances}")
    return {
        "threshold_pixels": 0.50,
        "lower_clearance_from_burden_track_pixels": lower_clearances,
        "upper_clearance_to_panel_c_footnote_pixels": upper_clearances,
        "all_clear": True,
    }


def measure_content_bbox(fig: mpl.figure.Figure) -> dict[str, Any]:
    """Record the delivered-page content bounds and occupancy in inches."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = fig.get_tightbbox(renderer)
    page_width, page_height = fig.get_size_inches()
    values = {
        "x0_inches": float(box.x0),
        "y0_inches": float(box.y0),
        "x1_inches": float(box.x1),
        "y1_inches": float(box.y1),
        "width_inches": float(box.width),
        "height_inches": float(box.height),
        "page_width_inches": float(page_width),
        "page_height_inches": float(page_height),
        "width_occupancy_fraction": float(box.width / page_width),
        "height_occupancy_fraction": float(box.height / page_height),
        "area_occupancy_fraction": float(box.width * box.height / (page_width * page_height)),
        "left_margin_inches": float(box.x0),
        "right_margin_inches": float(page_width - box.x1),
        "bottom_margin_inches": float(box.y0),
        "top_margin_inches": float(page_height - box.y1),
    }
    if min(values[key] for key in (
        "left_margin_inches", "right_margin_inches",
        "bottom_margin_inches", "top_margin_inches",
    )) < -1e-6:
        raise RuntimeError(f"content bbox extends outside short canvas: {values}")
    return values


def build_figure(
    v9: Any,
    v7: Any,
    v6: Any,
    dense: dict[str, Any],
    biological_rows: list[dict[str, Any]],
) -> tuple[mpl.figure.Figure, dict[str, Any], dict[str, pd.DataFrame], dict[str, Any]]:
    global v9_module
    v9_module = v9
    # Substitute only the panel-e renderer.  The inherited middle-row routine
    # still draws panel d and enforces its exact c[:, :80] geometry.
    v9.draw_panel_h = draw_panel_e_nonzero_go_slim
    base = v6.load_module(v6.BASE_RENDERER, "fig3_v10_base")
    bridge = base.bridge if hasattr(base, "bridge") else v6.load_module(
        v6.BRIDGE_RENDERER, "fig3_v10_bridge"
    )
    bridge.require_inputs()
    arce = bridge.read_json(bridge.INPUTS["arce"])
    concentration = bridge.read_json(bridge.INPUTS["concentration"])
    unanimity = bridge.read_json(bridge.INPUTS["unanimity"])
    cofactor = pd.read_csv(bridge.INPUTS["cofactor"])
    sensor = pd.read_csv(bridge.INPUTS["sensor"])
    memberships = pd.read_csv(bridge.INPUTS["programme_memberships"])

    bridge.configure_style()
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
        }
    )
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white")
    grid = fig.add_gridspec(
        3, 1,
        height_ratios=[2.00, 4.80, 1.85],
        left=0.100, right=0.985, top=0.975, bottom=0.055, hspace=0.040,
    )
    top = grid[0, 0].subgridspec(1, 2, width_ratios=[3.45, 3.72], wspace=0.04)
    bottom = grid[2, 0].subgridspec(
        1, 3, width_ratios=base.BOTTOM_WIDTH_RATIOS, wspace=base.BOTTOM_WSPACE
    )
    ax_ab = fig.add_subplot(top[0, 0])
    ax_c_parent = fig.add_subplot(top[0, 1])
    ax_de_parent = fig.add_subplot(grid[1, 0])
    ax_f = fig.add_subplot(bottom[0, 0])
    ax_g = fig.add_subplot(bottom[0, 1])
    ax_h = fig.add_subplot(bottom[0, 2])
    for parent in (ax_ab, ax_c_parent, ax_de_parent):
        parent.set_xticks([])
        parent.set_yticks([])

    v7.TOP_SQUARE_IN = TOP_SQUARE_IN
    source_a, source_b, top_geometry, top_axes = v7.clear_and_square_target_panels(
        ax_ab, concentration, unanimity, bridge, fig
    )
    panel_b = top_axes[1]
    if panel_b.get_title() != "Direction unanimity":
        raise RuntimeError(f"panel b inherited title drift: {panel_b.get_title()!r}")
    panel_b.set_title("Single-orientation targets", fontsize=7.0,
                      fontweight="bold", pad=4)
    summary_c, axes_c = draw_overview_panel_c(
        ax_c_parent, fig, v7, bridge, dense
    )

    middle_box = ax_de_parent.get_position()
    ax_de_parent.set_axes_locator(None)
    ax_de_parent.set_position([0.015, middle_box.y0, 0.980, middle_box.height], which="both")
    real_middle_box = ax_de_parent.get_position()
    page_width_in = fig.get_figwidth()
    real_left_in = real_middle_box.x0 * page_width_in
    real_right_in = real_middle_box.x1 * page_width_in
    desired_heat_x_in = real_left_in + PANEL_D_HEAT_OFFSET_IN
    # The inherited middle renderer has a hard minimum width for its temporary
    # biological-summary panel. Give that temporary call a wider off-page
    # parent, then remove its temporary summary axes and redraw panel e in the
    # audited on-page slot below. Panel d itself is placed on-page throughout.
    fake_left_in = desired_heat_x_in - 0.500
    fake_right_in = desired_heat_x_in + 80 * PANEL_D_CELL_PT / 72.0 + 0.075 + 1.825
    ax_de_parent.set_position(
        [
            fake_left_in / page_width_in,
            real_middle_box.y0,
            (fake_right_in - fake_left_in) / page_width_in,
            real_middle_box.height,
        ],
        which="both",
    )
    v9.PANEL_D_CELL_PT = PANEL_D_CELL_PT
    summary_d, axes_d, summary_e, axes_e = v9.draw_panel_d_and_h(
        ax_de_parent, fig, v7, bridge, dense, biological_rows
    )
    ax_de_parent.set_position(real_middle_box, which="both")
    for temporary_e_axis in dict.fromkeys(axes_e.values()):
        temporary_e_axis.remove()
    detail_heat_right_in = axes_d["heat"].get_position().x1 * page_width_in
    panel_e_x_in = detail_heat_right_in + PANEL_D_PANEL_E_GAP_IN
    panel_e_width_in = real_right_in - panel_e_x_in
    if panel_e_width_in < 1.58:
        raise RuntimeError(f"panel e on-page width too small: {panel_e_width_in:.3f} inches")
    summary_e, axes_e = draw_panel_e_nonzero_go_slim(
        fig,
        v7,
        bridge,
        biological_rows,
        x_in=panel_e_x_in,
        y_in=real_middle_box.y0 * fig.get_figheight(),
        width_in=panel_e_width_in,
        height_in=real_middle_box.height * fig.get_figheight(),
    )
    row_track_box = axes_d["row_track"].get_position()
    axes_d["row_track"].set_position(
        [
            PANEL_D_ROW_TRACK_X_IN / page_width_in,
            row_track_box.y0,
            row_track_box.width,
            row_track_box.height,
        ],
        which="both",
    )
    detail_center_fraction = (
        (
            axes_d["heat"].get_position().x0
            + axes_d["heat"].get_position().x1
        )
        / 2.0
        - real_middle_box.x0
    ) / real_middle_box.width
    for artist in ax_de_parent.texts:
        if artist.get_text() == "d":
            artist.set_y(1.055)
        elif artist.get_text().startswith("Enlargement of framed c"):
            artist.set_x(detail_center_fraction)
            artist.set_y(1.058)
        elif artist.get_text().startswith("R49 × first-80 intersection"):
            artist.set_x(detail_center_fraction)
            artist.set_y(1.030)
    panel_d_axis_typography = set_panel_d_axis_label_typography(axes_d)
    # Move the heatmap and all three marginal tracks as one rigid block,
    # preserving their exact square-cell geometry, values and mutual alignment.
    for detail_axis in axes_d.values():
        detail_box = detail_axis.get_position()
        detail_axis.set_position(
            [
                detail_box.x0,
                detail_box.y0 + DETAIL_VERTICAL_SHIFT_IN / fig.get_figheight(),
                detail_box.width,
                detail_box.height,
            ],
            which="both",
        )
    e_tags = [artist for artist in axes_e["parent"].texts if artist.get_text() == "e"]
    if len(e_tags) != 1:
        raise RuntimeError(f"expected one biological-summary tag e, found {len(e_tags)}")
    key_titles = [artist for artist in axes_e["panel_d_key"].texts if artist.get_text() == "d key"]
    if len(key_titles) != 1:
        raise RuntimeError(f"expected one inherited d-key title, found {len(key_titles)}")
    key_titles[0].set_text("c/d key")
    if not np.array_equal(dense["dense_state"], dense["state"][:, :80]):
        raise RuntimeError("panel d states are not panel c[:, :80]")
    if not np.array_equal(
        np.isnan(dense["dense_score"]), np.isnan(dense["score"][:, :80])
    ) or not np.allclose(
        np.nan_to_num(dense["dense_score"]), np.nan_to_num(dense["score"][:, :80]),
        rtol=0, atol=0,
    ):
        raise RuntimeError("panel d effects are not panel c[:, :80]")
    if dense["targets"].iloc[:80]["target_ensg"].tolist() != list(
        dense["states_frame"].columns[:80]
    ):
        raise RuntimeError("panel d target order is not panel c first-80 order")

    cofactor_matrix, _ = bridge.shared_wide_to_matrix(cofactor, ["CCNC", "TADA2B", "SUPT7L"])
    sensor_matrix, _ = bridge.shared_wide_to_matrix(sensor, ["AHR", "ARNT"])
    maximum = max(float(np.abs(cofactor_matrix).max()), float(np.abs(sensor_matrix).max()))
    shared_norm = TwoSlopeNorm(vmin=-maximum, vcenter=0.0, vmax=maximum)
    with v6.force_tag(bridge, "f"):
        image = bridge.panel_d_cofactor(ax_f, cofactor, shared_norm, memberships)
    with v6.force_tag(bridge, "g"):
        bridge.panel_e_sensor(ax_g, sensor, shared_norm)
    source_h = base.panel_f_external(ax_h, arce)
    inherited_tags = [artist for artist in ax_h.texts if artist.get_text() == "f"]
    if len(inherited_tags) != 1:
        raise RuntimeError(f"expected one inherited external-panel tag f, found {len(inherited_tags)}")
    inherited_tags[0].set_text("h")

    ax_f.set_title("Cofactor trio\n16 direction-concordant targets",
                   fontsize=7.0, fontweight="bold", y=1.19, pad=0)
    ax_g.set_title("Sensor pair\n17 direction-concordant targets",
                   fontsize=7.0, fontweight="bold", y=1.19, pad=0)
    ax_h.set_title("External edge recurrence", fontsize=7.0, fontweight="bold", pad=7)
    panel_g_tags = [artist for artist in ax_g.texts if artist.get_text() == "g"]
    if len(panel_g_tags) != 1:
        raise RuntimeError(f"expected one panel-g tag, found {len(panel_g_tags)}")
    panel_g_tags[0].set_x(-0.42)
    ax_f.xaxis.labelpad = 0.0
    ax_g.xaxis.labelpad = 0.0
    colorbar_rect = (
        base.COLORBAR_RECT[0], 0.018,
        base.COLORBAR_RECT[2], base.COLORBAR_RECT[3],
    )
    colorbar_ax = fig.add_axes(colorbar_rect)
    colorbar = fig.colorbar(image, cax=colorbar_ax, orientation="horizontal")
    colorbar.set_ticks([-1.5, 0.0, 1.5])
    colorbar.set_label("net knockdown coefficient (log2 fold change)", labelpad=2)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.xaxis.set_ticks_position("bottom")
    colorbar.ax.tick_params(labelsize=5.0)
    colorbar.ax.xaxis.label.set_fontsize(6.0)
    colorbar.outline.set_linewidth(0.5)

    v7.clamp_visible_fonts(fig)
    fig.canvas.draw()
    panel_d_axis_typography["post_clamp_and_draw"] = (
        audit_panel_d_axis_label_typography(axes_d, panel_d_axis_typography)
    )
    panel_d_renderer = fig.canvas.get_renderer()
    target_label_audit = v9.label_overlap_audit(
        list(axes_d["heat"].get_xticklabels()), panel_d_renderer, "x"
    )
    regulator_label_audit = v9.label_overlap_audit(
        list(axes_d["heat"].get_yticklabels()), panel_d_renderer, "y"
    )
    if target_label_audit["overlap_pairs_above_0.2px_each_axis"] or (
        regulator_label_audit["overlap_pairs_above_0.2px_each_axis"]
    ):
        raise RuntimeError(
            "5.5-pt panel d label collision: "
            f"targets={target_label_audit}; regulators={regulator_label_audit}"
        )
    summary_d["inherited_label_points"] = PANEL_D_INHERITED_AXIS_LABEL_PT
    summary_d["label_points"] = PANEL_D_AXIS_LABEL_PT
    summary_d["label_font_families"] = panel_d_axis_typography[
        "inherited_and_delivered_font_families"
    ]
    summary_d["target_label_collision_audit"] = target_label_audit
    summary_d["regulator_label_collision_audit"] = regulator_label_audit
    top_after = v7.measure_top_squares(fig, top_axes)
    if top_after != top_geometry:
        raise RuntimeError("top square geometry changed during final draw")
    heat_geometry = measure_heat_axes(fig, {"c": axes_c, "d": axes_d})
    text_audit = v7.audit_text(fig)
    targeted_text_geometry = audit_reported_text_relationships(fig, axes_e, ax_f, ax_g)
    panel_d_header_clearance = audit_panel_d_header_clearance(
        fig, ax_c_parent, ax_de_parent, axes_d
    )
    all_text_bbox_audit = audit_all_text_bboxes(fig)
    content_bbox = measure_content_bbox(fig)
    rendered_strings = [
        artist.get_text() for artist in fig.findobj(mpl.text.Text)
        if artist.get_visible() and artist.get_text().strip()
    ]
    forbidden = [
        value for value in rendered_strings
        if "direction unanimity" in value.lower() or "unanimous" in value.lower()
    ]
    if forbidden:
        raise RuntimeError(f"retired terminology remains rendered: {forbidden}")

    panel_tag_contract = {
        "a": [artist.get_text() for artist in top_axes[0].texts if artist.get_text() == "a"],
        "b": [artist.get_text() for artist in top_axes[1].texts if artist.get_text() == "b"],
        "c": [artist.get_text() for artist in ax_c_parent.texts if artist.get_text() == "c"],
        "d": [artist.get_text() for artist in ax_de_parent.texts if artist.get_text() == "d"],
        "e": [artist.get_text() for artist in axes_e["parent"].texts if artist.get_text() == "e"],
        "f": [artist.get_text() for artist in ax_f.texts if artist.get_text() == "f"],
        "g": [artist.get_text() for artist in ax_g.texts if artist.get_text() == "g"],
        "h": [artist.get_text() for artist in ax_h.texts if artist.get_text() == "h"],
    }
    if any(values != [letter] for letter, values in panel_tag_contract.items()):
        raise RuntimeError(f"panel letter contract failed: {panel_tag_contract}")

    summaries = {
        "layout_inches": {
            "page": [FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN],
            "top_square_side": TOP_SQUARE_IN,
            "overview_cell_side_points": OVERVIEW_CELL_PT,
            "detail_cell_side_points": PANEL_D_CELL_PT,
            "detail_vertical_shift_inches": DETAIL_VERTICAL_SHIFT_IN,
            "detail_heat_x_inches": axes_d["heat"].get_position().x0 * fig.get_figwidth(),
            "detail_row_track_x_inches": (
                axes_d["row_track"].get_position().x0 * fig.get_figwidth()
            ),
            "panel_e_x_inches": axes_e["parent"].get_position().x0 * fig.get_figwidth(),
            "panel_e_width_inches": (
                axes_e["parent"].get_position().width * fig.get_figwidth()
            ),
            "detail_panel_e_gap_inches": PANEL_D_PANEL_E_GAP_IN,
            "top_row_height_inches": ax_ab.get_position().height * fig.get_figheight(),
            "middle_row_height_inches": ax_de_parent.get_position().height * fig.get_figheight(),
            "bottom_row_height_inches": ax_f.get_position().height * fig.get_figheight(),
            "top_ab_slot_width_inches": ax_ab.get_position().width * fig.get_figwidth(),
            "top_c_slot_width_inches": ax_c_parent.get_position().width * fig.get_figwidth(),
            "bottom_colorbar_rect": list(colorbar_rect),
        },
        "panel_ab_geometry": top_geometry,
        "panel_c": summary_c,
        "panel_d": summary_d,
        "panel_d_axis_label_typography": panel_d_axis_typography,
        "panel_e": summary_e,
        "heat_geometry": heat_geometry,
        "text_audit": text_audit,
        "targeted_text_geometry": targeted_text_geometry,
        "panel_d_header_clearance": panel_d_header_clearance,
        "all_text_bbox_audit": all_text_bbox_audit,
        "content_bbox": content_bbox,
        "panel_tag_contract": {letter: values[0] for letter, values in panel_tag_contract.items()},
        "literal_zoom_assertions": {
            "states_equal_c_first_80": True,
            "scores_equal_c_first_80": True,
            "targets_equal_c_first_80": True,
            "row_order_equal": True,
            "orientation_track_equal_c_first_80": True,
            "target_burden_track_equal_c_first_80": True,
            "row_burden_track_identical_full_panel_scale": True,
        },
        "terminology_contract": {
            "b": "Single-orientation targets",
            "f": "Cofactor trio; 16 direction-concordant targets",
            "g": "Sensor pair; 17 direction-concordant targets",
            "forbidden_matches": [],
        },
    }
    source_frames = {
        "source_panel_a_target_concentration.csv": source_a,
        "source_panel_b_target_unanimity.csv": source_b,
        "source_panel_e_go_bp_slim_summary.tsv": pd.DataFrame(biological_rows),
        "source_panel_f_cofactor_targets.csv": cofactor,
        "source_panel_g_sensor_targets.csv": sensor,
        "source_panel_h_arce_recurrence.csv": source_h,
    }
    handles = {
        "top_axes": top_axes,
        "matrix_axes": {"c": axes_c, "d": axes_d},
        "panel_e_axes": axes_e,
    }
    return fig, summaries, source_frames, handles


def render_to_bytes(
    v7: Any,
    fig: mpl.figure.Figure,
    summaries: dict[str, Any],
    handles: dict[str, Any],
) -> tuple[bytes, bytes, dict[str, Any]]:
    before_top = v7.measure_top_squares(fig, handles["top_axes"])
    before_heat = measure_heat_axes(fig, handles["matrix_axes"])
    before_axis_typography = audit_panel_d_axis_label_typography(
        handles["matrix_axes"]["d"], summaries["panel_d_axis_label_typography"]
    )
    pdf_buffer = io.BytesIO()
    png_buffer = io.BytesIO()
    fig.savefig(
        pdf_buffer,
        format="pdf",
        facecolor="white",
        metadata={
            "Title": "Figure 3 top-overview layout v10.4.1 diagnostic",
            "Creator": "scripts/fig3_target_centred_top_overview_layout_v10_4_1.py",
            "CreationDate": PDF_METADATA_DATE,
            "ModDate": PDF_METADATA_DATE,
        },
    )
    fig.savefig(
        png_buffer,
        format="png",
        dpi=300,
        facecolor="white",
        metadata={"Software": "CoDEG Figure 3 top-overview layout v10.4.1"},
    )
    after_top = v7.measure_top_squares(fig, handles["top_axes"])
    after_heat = measure_heat_axes(fig, handles["matrix_axes"])
    after_axis_typography = audit_panel_d_axis_label_typography(
        handles["matrix_axes"]["d"], summaries["panel_d_axis_label_typography"]
    )
    if (
        before_top != after_top
        or before_heat != after_heat
        or before_axis_typography != after_axis_typography
    ):
        raise RuntimeError("axes geometry or panel d axis typography changed during save")
    summaries["post_save_canvas_geometry"] = {
        "top_equal_pre_save": True,
        "heat_equal_pre_save": True,
        "panel_d_axis_labels_equal_pre_save": True,
    }
    return pdf_buffer.getvalue(), png_buffer.getvalue(), summaries


def pdf_page_geometry(pdf: Path) -> dict[str, Any]:
    result = subprocess.run(["pdfinfo", str(pdf)], check=True, capture_output=True, text=True)
    page_line = next(
        (line for line in result.stdout.splitlines() if line.startswith("Page size:")), None
    )
    if page_line is None:
        raise RuntimeError("pdfinfo did not report page size")
    fields = page_line.split()
    width_pt = float(fields[2])
    height_pt = float(fields[4])
    if abs(width_pt - FIGURE_WIDTH_IN * 72.0) > 0.03 or abs(
        height_pt - FIGURE_HEIGHT_IN * 72.0
    ) > 0.03:
        raise RuntimeError(
            f"saved PDF does not match custom canvas: {width_pt} x {height_pt} pt"
        )
    return {
        "page_width_points": width_pt,
        "page_height_points": height_pt,
        "page_width_inches": width_pt / 72.0,
        "page_height_inches": height_pt / 72.0,
        "page_standard": "custom short portrait",
        "tight_bbox": False,
    }


def write_frame(path: Path, frame: pd.DataFrame) -> None:
    if path.suffix == ".tsv":
        frame.to_csv(path, sep="\t", index=False)
    else:
        frame.to_csv(path, index=False)


def build(output_dir: Path) -> None:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite retained candidate: {output_dir}")
    if list(output_dir.parent.glob(f".{output_dir.name}.staging*")):
        raise FileExistsError("refusing to reuse retained staging material")

    protected_before = protected_hashes()
    predecessor_manifest = json.loads((V10_3_DIR / "manifest.json").read_text(encoding="utf-8"))
    if predecessor_manifest.get("schema_version") != (
        "fig3-target-centred-top-overview-layout-v10.3"
    ):
        raise RuntimeError("audited v10.3 predecessor schema drift")
    v9 = load_module(V9_RENDERER, "fig3_v10_v9")
    v7 = load_module(V7_RENDERER, "fig3_v10_v7")
    bio = load_module(BIO_RENDERER, "fig3_v10_bio")
    dense = v7.load_dense_sources()
    v6 = v7.load_module(v7.V6_RENDERER, "fig3_v10_v6")
    bio_manifest = bio.validate_sources()
    biological_rows, _ = bio.read_rows()

    first_fig, first_summaries, source_frames, first_handles = build_figure(
        v9, v7, v6, dense, biological_rows
    )
    first_pdf, first_png, first_summaries = render_to_bytes(
        v7, first_fig, first_summaries, first_handles
    )
    plt.close(first_fig)
    second_fig, second_summaries, _, second_handles = build_figure(
        v9, v7, v6, dense, biological_rows
    )
    second_pdf, second_png, second_summaries = render_to_bytes(
        v7, second_fig, second_summaries, second_handles
    )
    plt.close(second_fig)
    if first_pdf != second_pdf or first_png != second_png:
        raise RuntimeError("fresh-canvas deterministic rerender bytes differ")
    if first_summaries != second_summaries:
        raise RuntimeError("fresh-canvas semantic/geometry summaries differ")
    preserved_geometry_keys = (
        "panel_ab_geometry",
        "panel_c",
        "panel_e",
        "literal_zoom_assertions",
        "terminology_contract",
        "panel_tag_contract",
    )
    for key in preserved_geometry_keys:
        if first_summaries[key] != predecessor_manifest["geometry"][key]:
            raise RuntimeError(f"v10.3 preserved geometry/semantic field drift: {key}")
    preserved_layout_fields = (
        "page",
        "top_square_side",
        "overview_cell_side_points",
        "top_row_height_inches",
        "middle_row_height_inches",
        "bottom_row_height_inches",
        "top_ab_slot_width_inches",
        "top_c_slot_width_inches",
        "bottom_colorbar_rect",
    )
    for key in preserved_layout_fields:
        if first_summaries["layout_inches"][key] != predecessor_manifest["geometry"][
            "layout_inches"
        ][key]:
            raise RuntimeError(f"v10.3 preserved page/top/bottom layout field drift: {key}")
    if first_summaries["heat_geometry"]["c"] != predecessor_manifest["geometry"][
        "heat_geometry"
    ]["c"]:
        raise RuntimeError("v10.3 panel-c heat geometry drift")
    panel_d_semantic_fields = (
        "shape",
        "reversal_cells",
        "lost_cells",
        "gained_cells",
        "stable_cells",
        "outside_cells",
        "self_cells",
        "regulators_with_reversal",
        "targets_with_reversal",
        "row_labels_shown",
        "target_labels_shown",
        "dendrogram",
        "ordering",
    )
    for key in panel_d_semantic_fields:
        if first_summaries["panel_d"][key] != predecessor_manifest["geometry"]["panel_d"][key]:
            raise RuntimeError(f"v10.3 panel-d semantic field drift: {key}")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging_", dir=output_dir.parent))
    pdf = staging / "Fig3_target_centred_top_overview_layout_v10_4_1.pdf"
    png = staging / "Fig3_target_centred_top_overview_layout_v10_4_1.png"
    pdf.write_bytes(first_pdf)
    png.write_bytes(first_png)
    saved_page = pdf_page_geometry(pdf)

    copied: dict[str, dict[str, Any]] = {}

    def copy_closed(source: Path, destination_name: str) -> None:
        destination = staging / destination_name
        shutil.copyfile(source, destination)
        copied[destination_name] = {
            "source": relative(source),
            "sha256": sha256(destination),
            "bytes": destination.stat().st_size,
        }

    for name in ["manifest.json", *v7.DENSE_INPUT_HASHES.keys()]:
        copy_closed(v7.DENSE / name, f"input_dense_{name}")
    base = v6.load_module(v6.BASE_RENDERER, "fig3_v10_bundle_base")
    bridge = base.bridge if hasattr(base, "bridge") else v6.load_module(
        v6.BRIDGE_RENDERER, "fig3_v10_bundle_bridge"
    )
    for name, source in bridge.INPUTS.items():
        copy_closed(source, f"input_{name}{source.suffix}")
    for source in (bio.SOURCE_MEMBERSHIPS, bio.SOURCE_COVERAGE, bio.SOURCE_MANIFEST):
        copy_closed(source, f"input_go_slim_{source.name}")

    renderers = {
        "renderer_fig3_target_centred_top_overview_layout_v10_4_1.py": Path(__file__).resolve(),
        "predecessor_renderer_fig3_target_centred_top_overview_layout_v10.py": V10_3_RENDERER,
        "renderer_fig3_target_centred_biological_summary_v9.py": V9_RENDERER,
        "renderer_fig3_target_centred_effect_zoom_v8.py": V8_RENDERER,
        "renderer_fig3_target_centred_effect_zoom_v7.py": V7_RENDERER,
        "renderer_fig3_target_centred_v6.py": v7.V6_RENDERER,
        "renderer_base_r5.py": v6.BASE_RENDERER,
        "renderer_bridge.py": v6.BRIDGE_RENDERER,
        "renderer_fig3c_go_slim_biological_summary_v2.py": BIO_RENDERER,
    }
    for destination_name, source in renderers.items():
        copy_closed(source, destination_name)

    for name, frame in source_frames.items():
        write_frame(staging / name, frame)
    score_frame = pd.DataFrame(
        dense["score"],
        index=dense["rows"]["regulator_ensg"],
        columns=dense["targets"]["target_ensg"],
    )
    score_frame.index.name = "regulator_ensg"
    score_frame.to_csv(staging / "source_panel_c_effect_score_matrix.tsv", sep="\t", na_rep="")
    score_frame.iloc[:, :80].to_csv(
        staging / "source_panel_d_effect_score_matrix.tsv", sep="\t", na_rep=""
    )
    predecessor_source_identity: dict[str, dict[str, Any]] = {}
    for current in sorted(staging.glob("source_panel_*")):
        predecessor = V10_3_DIR / current.name
        if not predecessor.is_file():
            raise FileNotFoundError(f"v10.3 predecessor source missing: {predecessor}")
        current_hash = sha256(current)
        predecessor_hash = sha256(predecessor)
        if current_hash != predecessor_hash:
            raise RuntimeError(f"panel source bytes changed from v10.3: {current.name}")
        predecessor_source_identity[current.name] = {
            "v10_3_sha256": predecessor_hash,
            "v10_4_1_sha256": current_hash,
            "byte_identical": True,
        }

    protected_after = protected_hashes()
    if protected_after != protected_before:
        changed = sorted(
            path for path in protected_before if protected_before[path] != protected_after.get(path)
        )
        raise RuntimeError(f"protected artifacts changed during render: {changed}")

    output_records: dict[str, dict[str, Any]] = {}
    for path in sorted(staging.glob("source_panel_*")) + [pdf, png]:
        output_records[path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "renderer": {
            "path": "scripts/fig3_target_centred_top_overview_layout_v10_4_1.py",
            "sha256": sha256(Path(__file__).resolve()),
        },
        "layout_order": {
            "top": ["a", "b", "c"],
            "middle": ["d", "e"],
            "bottom": ["f", "g", "h"],
        },
        "panel_map": {
            "a": "unchanged target concentration, smaller true square",
            "b": "unchanged single-orientation target statistic, smaller true square",
            "c": "unchanged complete R49 x T279 fixed-order signed effect overview",
            "d": (
                "literal panel c[:, :80] enlargement; 5.5-pt true-square cells "
                "and 5.5-pt inherited-sans-serif axis tick labels"
            ),
            "e": (
                "former h: 23 nonzero fixed GO-BP-slim terms plus outside-slim "
                "complement; two zero rows retained in source only"
            ),
            "f": "former e: cofactor trio, 16 direction-concordant targets",
            "g": "former f: sensor pair, 17 direction-concordant targets",
            "h": "former g: external edge recurrence",
        },
        "semantic_counts": {
            "panel_c": dense["full_summary"],
            "panel_d": dense["zoom_summary"],
            "panel_d_zero_reversal_targets": dense["zero_zoom_targets"],
            "panel_e": first_summaries["panel_e"],
        },
        "revision_from_v10_3": {
            "predecessor_directory": relative(V10_3_DIR),
            "predecessor_manifest_sha256": sha256(V10_3_DIR / "manifest.json"),
            "predecessor_pdf_sha256": sha256(
                V10_3_DIR / "Fig3_target_centred_top_overview_layout_v10_3.pdf"
            ),
            "predecessor_png_sha256": sha256(
                V10_3_DIR / "Fig3_target_centred_top_overview_layout_v10_3.png"
            ),
            "readability_changes_only": (
                "panel d x-axis target and y-axis regulator tick-label font sizes "
                "increase from 5.0 to exactly 5.5 points; panel d cells increase "
                "from 5.0 to 5.5-point true squares and panel e narrows to retain "
                "the same short page"
            ),
            "panel_d_x_tick_labels": 80,
            "panel_d_y_tick_labels": 49,
            "page_canvas_top_bottom_layout_data_order_scales_legends_unchanged": True,
            "middle_panel_d_e_geometry_changed_only_for_readability": True,
            "panel_source_files_byte_identical": predecessor_source_identity,
            "preserved_geometry_and_semantic_fields": list(preserved_geometry_keys),
            "preserved_layout_fields": list(preserved_layout_fields),
            "preserved_panel_d_semantic_fields": list(panel_d_semantic_fields),
        },
        "superseded_preview_retained": {
            "directory": relative(V10_4_PREVIEW_DIR),
            "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
            "reason": "panel-d header overprinted by top tracks and panel-e line clipped",
            "manifest_sha256": sha256(V10_4_PREVIEW_DIR / "manifest.json"),
            "pdf_sha256": sha256(
                V10_4_PREVIEW_DIR / "Fig3_target_centred_top_overview_layout_v10_4.pdf"
            ),
            "png_sha256": sha256(
                V10_4_PREVIEW_DIR / "Fig3_target_centred_top_overview_layout_v10_4.png"
            ),
            "bytes_unchanged_by_v10_4_1_build": True,
        },
        "claim_boundary": (
            "Relative to audited v10.3, only panel d's 80 target and 49 regulator "
            "axis tick-label sizes change from 5.0 to exactly 5.5 points; its cells "
            "increase from 5.0 to 5.5-point true squares, and panel e narrows to keep "
            "the same short page. Panel e describes "
            "frozen broad GO "
            "Biological Process slim membership among 279 single-orientation targets; it "
            "is not enrichment, a lost-versus-gained test, pathway activity, community "
            "structure, or mechanism. Its 39-target complement means outside this fixed "
            "BP slim, not biologically unannotated. The two zero-total fixed terms, cell "
            "wall organization or biogenesis and secondary metabolic process, remain in "
            "the bundled 26-row source but are omitted from the 24-row display. Panel d "
            "remains the literal first-80 enlargement of panel c and is descriptive, not "
            "separately tested."
        ),
        "source_closure": {
            "dense_source_manifest_sha256": dense["manifest_sha256"],
            "go_slim_source_manifest_schema": bio_manifest["schema_version"],
            "go_slim_source_manifest_status": bio_manifest["status"],
            "panel_e_source_rows": 26,
            "panel_e_displayed_rows": 24,
            "panel_e_displayed_nonzero_fixed_terms": 23,
            "panel_e_omitted_zero_rows_retained_in_source": [
                "GO:0071554 cell wall organization or biogenesis",
                "GO:0019748 secondary metabolic process",
            ],
            "all_input_bytes_copied_into_bundle": True,
            "all_copied_input_and_renderer_hashes_recorded": True,
        },
        "geometry": first_summaries,
        "saved_pdf_geometry": saved_page,
        "export": {
            "page": "custom short portrait",
            "page_inches": [FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN],
            "pdf_vector_heatmaps": True,
            "pdf_fonttype": 42,
            "png_dpi": 300,
            "tight_bbox": False,
            "fresh_canvas_count": 2,
            "pdf_byte_identical": True,
            "png_byte_identical": True,
            "pdf_sha256": hash_bytes(first_pdf),
            "png_sha256": hash_bytes(first_png),
        },
        "no_write_guard": {
            "paired_before_after_hashes": True,
            "protected_before": protected_before,
            "protected_after": protected_after,
            "all_equal": True,
        },
        "bundled_inputs_and_renderers": copied,
        "outputs": output_records,
        "canonical_changes": False,
        "manuscript_changes": False,
        "legend_changes": False,
        "active_claim_guard_changes": False,
        "prior_candidate_changes": False,
        "trash_changes": False,
    }
    (staging / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(staging, output_dir)
    print(f"wrote {output_dir}")
    print(f"  PDF sha256 {sha256(output_dir / pdf.name)}")
    print(f"  PNG sha256 {sha256(output_dir / png.name)}")
    print(f"  manifest sha256 {sha256(output_dir / 'manifest.json')}")
    print("  v10.3 order a/b/c | d/e | f/g/h retained; all source panels byte-identical")
    print("  panel d labels/cells: 5.0 pt -> 5.5 pt; panel e narrowed on same short page")
    print("  superseded v10.4 preview retained unchanged")
    print("  canonical/manuscript/legend/claim-guard/prior candidates unchanged")


def audit_only() -> None:
    v9 = load_module(V9_RENDERER, "fig3_v10_audit_v9")
    v7 = load_module(V7_RENDERER, "fig3_v10_audit_v7")
    bio = load_module(BIO_RENDERER, "fig3_v10_audit_bio")
    dense = v7.load_dense_sources()
    v6 = v7.load_module(v7.V6_RENDERER, "fig3_v10_audit_v6")
    bio.validate_sources()
    biological_rows, _ = bio.read_rows()
    fig, summaries, _, handles = build_figure(v9, v7, v6, dense, biological_rows)
    _, _, summaries = render_to_bytes(v7, fig, summaries, handles)
    plt.close(fig)
    print(json.dumps({"status": "PASS", "geometry": summaries}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    if args.audit_only:
        audit_only()
    else:
        build(args.output_dir)


if __name__ == "__main__":
    main()
