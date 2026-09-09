#!/usr/bin/env python3
"""Build a diagnostic Figure 3 v11 with corrected all-target enrichment in panel e.

This candidate inherits the audited v10.4.1 layout and replaces only the
descriptive GO-BP-slim panel e.  The replacement is a compact rendering of the
corrected v2 all-target composition overlay: all nine repeatable-FDR KEGG terms
and all ten bank-stable global-maxT GO Biological Process terms.  Panels a-d
and f-h, the 5.5-point square-cell zoom, and the c/d key are preserved.

The script refuses to overwrite an existing candidate, renders twice from fresh
canvases, verifies unchanged panel sources and geometry against v10.4.1, and
hashes protected canonical and historical artifacts before and after.  It never
promotes the candidate or edits active manuscript files.

Run only through scripts/run_memory_guarded.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PARENT_RENDERER = ROOT / "scripts/fig3_target_centred_top_overview_layout_v10_4_1.py"
PARENT_DIR = (
    ROOT / "figures_final/diagnostics/fig3h_target_centred_top_overview_layout_v10_4_1"
)
OVERLAY_JSON = ROOT / "results_data/census/gsea_pathway_reversal_overlay.json"
AUTHORITATIVE_TERMS = (
    ROOT
    / "results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2"
    / "production_v2/01_terms.tsv"
)
FINAL_AUDIT = (
    ROOT
    / "results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2"
    / "supplement_package_final_audit.json"
)
CENSUS = ROOT / "results_data/census/census_master_edges.csv.gz"
OVERLAY_GENERATOR = ROOT / "scripts/gsea_pathway_reversal_overlay.py"
OVERLAY_RENDERER = ROOT / "scripts/fig_gsea_pathway_reversal_overlay.py"
DEFAULT_OUT = (
    ROOT / "figures_final/diagnostics/fig3_target_centred_alltarget_enrichment_v11"
)
SCHEMA_VERSION = "fig3-target-centred-alltarget-enrichment-v11"
STEM = "Fig3_target_centred_alltarget_enrichment_v11"
PDF_METADATA_DATE = datetime(2026, 9, 4, tzinfo=timezone.utc)

LOST = "#c1272d"
GAINED = "#2a9d8f"
MIXED = "#8a6fb0"
NO_REVERSAL = "#e3e6e9"
INK = "#22262B"
MUTED = "#666B70"
GRID = "#c7ccd1"

TITLE_PT = 7.0
TEXT_PT = 5.0
ROW_PITCH_PT = 7.0
KEY_HEIGHT_IN = 1.30
REFERENCE_PANEL_E_ROWS = 24
REFERENCE_PANEL_E_ROW_PITCH_PT = 6.05

DISPLAY_LABELS = {
    "Cytokine Cytokine Receptor Interaction": "Cytokine-receptor interaction",
    "Proteasome": "Proteasome",
    "Apoptosis": "Apoptosis",
    "Antigen Processing And Presentation": "Antigen presentation",
    "Jak Stat Signaling Pathway": "JAK-STAT signalling",
    "Prion Diseases": "Prion diseases",
    "Thyroid Cancer": "Thyroid cancer",
    "Chronic Myeloid Leukemia": "Chronic myeloid leukaemia",
    "P53 Signaling Pathway": "p53 signalling",
    "Apoptotic Process": "Apoptotic process",
    "Organelle Organization": "Organelle organisation",
    "Intrinsic Apoptotic Signaling Pathway In Response To Endoplasmic Reticulum Stress": (
        "ER-stress intrinsic apoptosis"
    ),
    "Cellular Response To Lipopolysaccharide": "Response to LPS",
    "Cellular Response To Cytokine Stimulus": "Response to cytokine",
    "Intrinsic Apoptotic Signaling Pathway": "Intrinsic apoptotic signalling",
    "Mitotic Cytokinesis": "Mitotic cytokinesis",
    "Regulation Of Type II Interferon Production": "Type II IFN regulation",
    "Positive Regulation Of Apoptotic Process": "Positive apoptosis regulation",
    "Proteasome-Mediated Ubiquitin-Dependent Protein Catabolic Process": (
        "Proteasomal protein turnover"
    ),
}

V10: Any = None
PANEL_ROWS: list[dict[str, Any]] = []
PANEL_SOURCE_FRAME: pd.DataFrame | None = None
PANEL_VALIDATION: dict[str, Any] = {}


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


def protected_paths() -> list[Path]:
    paths = [
        ROOT / "figures_final/Fig3.pdf",
        ROOT / "figures_final/Fig3.png",
        ROOT / "MANUSCRIPT.md",
        ROOT / "FIGURE_LEGENDS.md",
        ROOT / "PAPER_FACTS.md",
        ROOT / "HumanCD4CoDEGNet_signflip.docx",
        ROOT / "HumanCD4CoDEGNet_signflip.pdf",
        ROOT / "scripts/active_claim_guard.py",
        PARENT_RENDERER,
        OVERLAY_JSON,
        AUTHORITATIVE_TERMS,
        FINAL_AUDIT,
        CENSUS,
        OVERLAY_GENERATOR,
        OVERLAY_RENDERER,
    ]
    for directory in (
        PARENT_DIR,
        ROOT / "figures_final/diagnostics/gsea_pathway_reversal_overlay_v1",
        ROOT / "figures_final/diagnostics/gsea_pathway_reversal_overlay_v2",
    ):
        if directory.is_dir():
            paths.extend(path for path in directory.iterdir() if path.is_file())
    unique = sorted({path.resolve() for path in paths})
    missing = [relative(path) for path in unique if not path.is_file()]
    if missing:
        raise FileNotFoundError("protected artifact missing: " + ", ".join(missing))
    return unique


def file_hashes(paths: list[Path]) -> dict[str, str]:
    return {relative(path): sha256(path) for path in paths}


def read_authoritative_terms() -> list[dict[str, str]]:
    with AUTHORITATIVE_TERMS.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_panel_source() -> tuple[list[dict[str, Any]], pd.DataFrame, dict[str, Any]]:
    overlay = json.loads(OVERLAY_JSON.read_text(encoding="utf-8"))
    expected_inputs = {
        "01_terms.tsv": sha256(AUTHORITATIVE_TERMS),
        "supplement_package_final_audit.json": sha256(FINAL_AUDIT),
        "census_master_edges.csv.gz": sha256(CENSUS),
    }
    if overlay.get("inputs") != expected_inputs:
        raise RuntimeError(
            f"overlay input hashes do not close to authoritative files: {overlay.get('inputs')}"
        )

    audit = json.loads(FINAL_AUDIT.read_text(encoding="utf-8"))
    if audit.get("status") != "PASS":
        raise RuntimeError("broad-enrichment supplementary package audit is not PASS")
    authoritative = audit.get("authoritative_results", {})
    expected_authoritative = {
        "kegg_tested": 168,
        "kegg_statement_eligible": 9,
        "kegg_global_maxT_stable": 0,
        "go_bp_tested": 2425,
        "go_bp_statement_eligible": 161,
        "go_bp_pooled_global_maxT": 12,
        "go_bp_global_maxT_stable": 10,
    }
    if authoritative != expected_authoritative:
        raise RuntimeError(f"authoritative result-count drift: {authoritative}")

    terms = read_authoritative_terms()
    kegg_source = [
        row for row in terms
        if row["library"] == "KEGG_Legacy_v2026.1"
        and row["statement_eligible"] == "True"
    ]
    go_source = [
        row for row in terms
        if row["library"] == "GO_Biological_Process_2023"
        and row["global_familywise_stable"] == "True"
    ]
    if len(kegg_source) != 9 or len(go_source) != 10:
        raise RuntimeError("authoritative 9-KEGG/10-GO selection drift")
    if any(row["global_familywise_stable"] == "True" for row in kegg_source):
        raise RuntimeError("a displayed KEGG term unexpectedly passes global maxT")
    if any(row["global_familywise_stable"] != "True" for row in go_source):
        raise RuntimeError("a displayed GO term is not bank-stable under global maxT")

    source_by_name = {row["term_name"]: row for row in kegg_source + go_source}
    overlay_rows = list(overlay.get("pathways", []))
    if len(overlay_rows) != 19 or set(source_by_name) != {
        str(row["term"]) for row in overlay_rows
    }:
        raise RuntimeError("overlay term identities differ from authoritative v2 selection")

    kegg_rows = sorted(
        (row for row in overlay_rows if row["library"] == "KEGG_Legacy"),
        key=lambda row: float(row["q"]),
    )
    go_rows = sorted(
        (row for row in overlay_rows if row["library"].startswith("GO_Biological")),
        key=lambda row: float(row["q"]),
    )
    rows = kegg_rows + go_rows
    if len(kegg_rows) != 9 or len(go_rows) != 10:
        raise RuntimeError("overlay library partition is not 9 KEGG plus 10 GO-BP")

    records: list[dict[str, Any]] = []
    for display_order, row in enumerate(rows, start=1):
        term = str(row["term"])
        source = source_by_name[term]
        if term not in DISPLAY_LABELS:
            raise RuntimeError(f"missing compact display label for {term}")
        expected_library = (
            "KEGG_Legacy"
            if source["library"] == "KEGG_Legacy_v2026.1"
            else "GO_Biological_Process"
        )
        checks = {
            "library": row["library"] == expected_library,
            "nes": math.isclose(
                float(row["NES"]), round(float(source["observed_nes"]), 4),
                rel_tol=0.0, abs_tol=1e-12,
            ),
            "q": math.isclose(
                float(row["q"]), round(float(source["bh_q_within_library"]), 4),
                rel_tol=0.0, abs_tol=1e-12,
            ),
            "members": int(row["n_members_in_universe"]) == int(source["member_targets_n"]),
            "partition": (
                int(row["n_lost"])
                + int(row["n_gained"])
                + int(row["n_mixed"])
                == int(row["n_members_with_reversal"])
            ),
        }
        if not all(checks.values()):
            raise RuntimeError(f"overlay/source mismatch for {term}: {checks}")
        n_members = int(row["n_members_in_universe"])
        n_no_reversal = n_members - int(row["n_members_with_reversal"])
        if n_no_reversal < 0:
            raise RuntimeError(f"negative no-reversal remainder for {term}")
        evidence = (
            "repeatable within-library FDR q<0.10; not global maxT"
            if expected_library == "KEGG_Legacy"
            else "bank-stable global maxT"
        )
        records.append(
            {
                "display_order": display_order,
                "library": expected_library,
                "evidence_tier": evidence,
                "term": term,
                "display_label": DISPLAY_LABELS[term],
                "n_members_in_universe": n_members,
                "n_lost": int(row["n_lost"]),
                "n_gained": int(row["n_gained"]),
                "n_mixed": int(row["n_mixed"]),
                "n_no_reversal": n_no_reversal,
                "lost_share": int(row["n_lost"]) / n_members,
                "gained_share": int(row["n_gained"]) / n_members,
                "mixed_share": int(row["n_mixed"]) / n_members,
                "no_reversal_share": n_no_reversal / n_members,
                "network_null_normalised_es": float(row["NES"]),
                "bh_q_within_library": float(row["q"]),
                "global_familywise_stable": bool(row["global_familywise_stable"]),
            }
        )

    # The retired asterisk/footnote asserted an overlap that the authoritative
    # collection does not contain. Keep that correction machine-checked.
    source_all = {row["term_name"]: row for row in terms}
    prion = source_all["Prion Diseases"]
    proteasome = source_all["Proteasome"]
    overlap_counts: dict[str, int] = {}
    for field in (
        "member_targets",
        "leading_edge_targets",
        "leading_edge_targets_tie_expanded",
    ):
        left = set(prion[field].split("|")) - {""}
        right = set(proteasome[field].split("|")) - {""}
        overlap_counts[field] = len(left & right)
    if any(overlap_counts.values()):
        raise RuntimeError(f"Prion/Proteasome overlap unexpectedly nonzero: {overlap_counts}")
    if any("*" in record["display_label"] for record in records):
        raise RuntimeError("retired asterisk remains in a panel-e display label")

    validation = {
        "all_target_universe": 9554,
        "targets_with_zero_observed_reversals": 6977,
        "displayed_terms": 19,
        "kegg_repeatable_fdr_terms": 9,
        "kegg_global_maxT_stable_terms": 0,
        "go_bp_global_maxT_stable_terms": 10,
        "go_bp_statement_eligible_terms_total": 161,
        "term_selection_matches_authoritative_table": True,
        "composition_values_match_overlay_json": True,
        "prion_proteasome_overlap_counts": overlap_counts,
        "false_overlap_asterisk_or_footnote_rendered": False,
        "upstream_inputs": expected_inputs,
    }
    return rows, pd.DataFrame.from_records(records), validation


def draw_term_group(
    ax: mpl.axes.Axes,
    rows: list[dict[str, Any]],
    *,
    show_x: bool,
) -> None:
    y_positions = np.arange(len(rows))[::-1]
    for y_value, row in zip(y_positions, rows):
        n_members = int(row["n_members_in_universe"])
        values = (
            int(row["n_lost"]),
            int(row["n_gained"]),
            int(row["n_mixed"]),
            n_members - int(row["n_members_with_reversal"]),
        )
        left = 0.0
        for value, color in zip(values, (LOST, GAINED, MIXED, NO_REVERSAL)):
            fraction = value / n_members
            ax.barh(
                y_value,
                fraction,
                left=left,
                height=0.68,
                color=color,
                edgecolor="white",
                linewidth=0.20,
                clip_on=True,
            )
            left += fraction
        if not math.isclose(left, 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"composition does not sum to one for {row['term']}: {left}")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [DISPLAY_LABELS[str(row["term"])] for row in rows],
        fontsize=TEXT_PT,
        fontfamily="DejaVu Sans",
        fontstretch="condensed",
        color=INK,
    )
    ax.tick_params(axis="y", length=0, pad=1.2)
    if show_x:
        ax.set_xticks([0.0, 1.0])
        ax.set_xticklabels(["0%", "100%"], fontsize=TEXT_PT, color=MUTED)
        ax.get_xticklabels()[0].set_horizontalalignment("left")
        ax.get_xticklabels()[1].set_horizontalalignment("right")
        ax.tick_params(axis="x", length=0, pad=1.0)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_color(GRID)
        ax.spines["bottom"].set_linewidth(0.45)
    else:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def draw_panel_e_alltarget_terms(
    fig: mpl.figure.Figure,
    v7: Any,
    bridge: Any,
    _inherited_go_slim_rows: list[dict[str, Any]],
    *,
    x_in: float,
    y_in: float,
    width_in: float,
    height_in: float,
) -> tuple[dict[str, Any], dict[str, mpl.axes.Axes]]:
    if len(PANEL_ROWS) != 19 or PANEL_SOURCE_FRAME is None:
        raise RuntimeError("validated panel-e source was not initialized")
    kegg_rows = [row for row in PANEL_ROWS if row["library"] == "KEGG_Legacy"]
    go_rows = [
        row for row in PANEL_ROWS if row["library"].startswith("GO_Biological")
    ]
    if len(kegg_rows) != 9 or len(go_rows) != 10:
        raise RuntimeError("panel-e 9/10 term partition drift")

    panel = v7.add_absolute_axis(fig, x_in, y_in, width_in, height_in)
    panel.set_xticks([])
    panel.set_yticks([])
    panel.set_axis_off()
    bridge.panel_tag(panel, "e", x=-0.035, y=1.01)
    panel.text(
        0.5,
        0.988,
        "All-target enriched terms",
        transform=panel.transAxes,
        ha="center",
        va="top",
        fontsize=TITLE_PT,
        fontweight="bold",
        color=INK,
    )
    panel.text(
        0.5,
        0.948,
        "all 9,554 targets",
        transform=panel.transAxes,
        ha="center",
        va="top",
        fontsize=TEXT_PT,
        color=MUTED,
    )
    legend = panel.legend(
        handles=[
            Rectangle((0, 0), 1, 1, facecolor=LOST, label="lost"),
            Rectangle((0, 0), 1, 1, facecolor=GAINED, label="gained"),
            Rectangle((0, 0), 1, 1, facecolor=MIXED, label="mixed"),
            Rectangle((0, 0), 1, 1, facecolor=NO_REVERSAL, edgecolor=GRID, label="no reversal"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        fontsize=TEXT_PT,
        frameon=False,
        handlelength=1.0,
        handleheight=0.7,
        columnspacing=0.75,
        borderpad=0.0,
        labelspacing=0.20,
    )
    legend.set_zorder(10)

    # Preserve the c/d key at the exact v10.4.1 position. Its position was
    # defined from the former 24-row table, not from the replacement's 19 rows.
    reference_table_top = y_in + height_in * 0.770
    reference_table_height = (
        REFERENCE_PANEL_E_ROWS * REFERENCE_PANEL_E_ROW_PITCH_PT / 72.0
    )
    reference_table_y = reference_table_top - reference_table_height
    key_top = reference_table_y - 0.08
    key_bottom = key_top - KEY_HEIGHT_IN
    if key_bottom < y_in + 0.04:
        raise RuntimeError("preserved c/d key does not fit the panel-e slot")

    bar_width = min(0.48, width_in * 0.30)
    bar_x = x_in + width_in - bar_width - 0.025
    go_height = len(go_rows) * ROW_PITCH_PT / 72.0
    kegg_height = len(kegg_rows) * ROW_PITCH_PT / 72.0
    go_y = key_top + 0.220
    go_header_y = go_y + go_height + 0.052
    kegg_y = go_header_y + 0.165
    kegg_header_y = kegg_y + kegg_height + 0.052
    if kegg_header_y > y_in + height_in * 0.855:
        raise RuntimeError("panel-e term groups intrude into the legend/header")
    if go_y <= key_top + 0.10:
        raise RuntimeError("panel-e GO terms leave insufficient clearance above c/d key")

    panel.text(
        0.02,
        (kegg_header_y - y_in) / height_in,
        "KEGG · repeatable FDR (9; 0 maxT)",
        transform=panel.transAxes,
        ha="left",
        va="center",
        fontsize=TEXT_PT,
        fontweight="bold",
        color=INK,
    )
    panel.text(
        0.02,
        (go_header_y - y_in) / height_in,
        "GO-BP · global maxT-stable (10)",
        transform=panel.transAxes,
        ha="left",
        va="center",
        fontsize=TEXT_PT,
        fontweight="bold",
        color=INK,
    )

    kegg_ax = v7.add_absolute_axis(fig, bar_x, kegg_y, bar_width, kegg_height)
    go_ax = v7.add_absolute_axis(fig, bar_x, go_y, bar_width, go_height)
    draw_term_group(kegg_ax, kegg_rows, show_x=False)
    draw_term_group(go_ax, go_rows, show_x=True)

    key = V10.v9_module.draw_panel_d_key(
        fig,
        v7,
        x_in=x_in + 0.035,
        y_in=key_bottom,
        width_in=width_in - 0.070,
        height_in=KEY_HEIGHT_IN,
    )
    orientation_labels = [
        artist for artist in key.texts if artist.get_text() in {"lost", "gained"}
    ]
    if len(orientation_labels) != 2:
        raise RuntimeError("preserved c/d key orientation-label identity drift")
    for artist in orientation_labels:
        artist.set_y(0.525)
    orientation_swatches = [
        patch
        for patch in key.patches
        if isinstance(patch, Rectangle) and abs(float(patch.get_y()) - 0.535) < 1e-9
    ]
    if len(orientation_swatches) != 2:
        raise RuntimeError("preserved c/d key orientation-swatch identity drift")
    for patch in orientation_swatches:
        patch.set_y(0.495)

    summary = {
        **PANEL_VALIDATION,
        "display": "lost/gained/mixed/no-reversal share of census members",
        "displayed_per_term_n_nes_q": False,
        "statistics_retained_in_source_table": True,
        "row_pitch_points": ROW_PITCH_PT,
        "term_label_points": TEXT_PT,
        "heading_points": TITLE_PT,
        "bar_width_inches": bar_width,
        "evidence_tiers_visibly_separated": True,
        "labels_call_go_entries_terms_not_pathways": True,
        "c_d_key_preserved_at_v10_4_1_position": True,
        "claim_boundary": (
            "all-target network-adjusted pooled reversal-burden enrichment; KEGG and GO-BP "
            "use explicitly separated evidence tiers. Bar colours are a descriptive partition "
            "of each term's census members by dataset-wide target orientation and do not test "
            "lost-versus-gained enrichment, pathway activity, mechanism, individual edges, or "
            "enrichment among the 279 targets in panels c-d. Terms overlap."
        ),
    }
    axes = {
        "parent": panel,
        # The inherited temporary middle-row constructor audits the minimum
        # label x-position through a historical `heat` handle before removing
        # that temporary panel. KEGG contains the widest compact label.
        "heat": kegg_ax,
        "kegg": kegg_ax,
        "go_bp": go_ax,
        "panel_d_key": key,
    }
    return summary, axes


def audit_replacement_text_relationships(
    fig: mpl.figure.Figure,
    axes_e: dict[str, mpl.axes.Axes],
    ax_f: mpl.axes.Axes,
    ax_g: mpl.axes.Axes,
) -> dict[str, Any]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    parent_strings = [artist.get_text() for artist in axes_e["parent"].texts]
    required = {
        "e",
        "All-target enriched terms",
        "all 9,554 targets",
        "KEGG · repeatable FDR (9; 0 maxT)",
        "GO-BP · global maxT-stable (10)",
    }
    if not required.issubset(parent_strings):
        raise RuntimeError(f"panel-e required text missing: {required - set(parent_strings)}")
    term_labels = [
        label.get_text()
        for axis_name in ("kegg", "go_bp")
        for label in axes_e[axis_name].get_yticklabels()
    ]
    if term_labels != [DISPLAY_LABELS[str(row["term"])] for row in PANEL_ROWS]:
        raise RuntimeError("panel-e rendered term order or compact labels drifted")
    if any("*" in value for value in parent_strings + term_labels):
        raise RuntimeError("retired Prion/Proteasome asterisk remains rendered")

    relationships: list[dict[str, Any]] = []

    def record(label: str, first: mpl.text.Text, second: mpl.text.Text) -> None:
        width, height = V10.overlap_dimensions(
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
            raise RuntimeError(f"targeted text collision {label}: {width} x {height}")

    f_headers = [
        artist for artist in ax_f.texts
        if artist.get_text() in {"CCNC", "TADA2B", "SUPT7L"}
    ]
    g_headers = [
        artist for artist in ax_g.texts if artist.get_text() in {"AHR", "ARNT"}
    ]
    if len(f_headers) != 3 or len(g_headers) != 2:
        raise RuntimeError("unchanged bottom-panel regulator header identity drift")
    for header in f_headers:
        record(f"panel f title vs {header.get_text()}", ax_f.title, header)
    for header in g_headers:
        record(f"panel g title vs {header.get_text()}", ax_g.title, header)
    return {
        "threshold_pixels_each_axis": 0.50,
        "panel_e_term_labels": len(term_labels),
        "panel_e_required_text_present": True,
        "retired_overlap_asterisks_present": False,
        "relationships": relationships,
    }


def render_to_bytes(
    v7: Any,
    fig: mpl.figure.Figure,
    summaries: dict[str, Any],
    handles: dict[str, Any],
) -> tuple[bytes, bytes, dict[str, Any]]:
    before_top = v7.measure_top_squares(fig, handles["top_axes"])
    before_heat = V10.measure_heat_axes(fig, handles["matrix_axes"])
    before_typography = V10.audit_panel_d_axis_label_typography(
        handles["matrix_axes"]["d"], summaries["panel_d_axis_label_typography"]
    )
    pdf_buffer = io.BytesIO()
    png_buffer = io.BytesIO()
    fig.savefig(
        pdf_buffer,
        format="pdf",
        facecolor="white",
        metadata={
            "Title": "Figure 3 all-target enrichment panel v11 diagnostic",
            "Subject": "Panel e replacement candidate; not promoted",
            "Creator": "scripts/fig3_target_centred_alltarget_enrichment_v11.py",
            "CreationDate": PDF_METADATA_DATE,
            "ModDate": PDF_METADATA_DATE,
        },
    )
    fig.savefig(
        png_buffer,
        format="png",
        dpi=300,
        facecolor="white",
        metadata={"Software": "CoDEG Figure 3 all-target enrichment v11 diagnostic"},
    )
    after_top = v7.measure_top_squares(fig, handles["top_axes"])
    after_heat = V10.measure_heat_axes(fig, handles["matrix_axes"])
    after_typography = V10.audit_panel_d_axis_label_typography(
        handles["matrix_axes"]["d"], summaries["panel_d_axis_label_typography"]
    )
    if before_top != after_top or before_heat != after_heat or before_typography != after_typography:
        raise RuntimeError("unchanged axes geometry or panel-d typography changed during save")
    summaries["post_save_canvas_geometry"] = {
        "top_equal_pre_save": True,
        "heat_equal_pre_save": True,
        "panel_d_axis_labels_equal_pre_save": True,
    }
    return pdf_buffer.getvalue(), png_buffer.getvalue(), summaries


def assert_parent_geometry(summaries: dict[str, Any], parent: dict[str, Any]) -> dict[str, Any]:
    parent_geometry = parent["geometry"]

    def json_value(value: Any) -> Any:
        """Compare runtime tuples/NumPy scalars as their manifest JSON values."""
        return json.loads(json.dumps(value, sort_keys=True))

    exact_keys = (
        "panel_ab_geometry",
        "panel_c",
        "panel_d",
        "panel_d_axis_label_typography",
        "heat_geometry",
        "literal_zoom_assertions",
        "panel_d_header_clearance",
        "panel_tag_contract",
        "terminology_contract",
    )
    for key in exact_keys:
        if json_value(summaries[key]) != parent_geometry[key]:
            raise RuntimeError(f"v10.4.1 unchanged geometry/semantic field drift: {key}")
    layout_keys = (
        "page",
        "top_square_side",
        "overview_cell_side_points",
        "detail_cell_side_points",
        "detail_vertical_shift_inches",
        "detail_heat_x_inches",
        "detail_row_track_x_inches",
        "panel_e_x_inches",
        "panel_e_width_inches",
        "detail_panel_e_gap_inches",
        "top_row_height_inches",
        "middle_row_height_inches",
        "bottom_row_height_inches",
        "top_ab_slot_width_inches",
        "top_c_slot_width_inches",
        "bottom_colorbar_rect",
    )
    for key in layout_keys:
        if json_value(summaries["layout_inches"][key]) != parent_geometry["layout_inches"][key]:
            raise RuntimeError(f"v10.4.1 layout field drift: {key}")
    return {
        "geometry_fields_byte_value_equal": list(exact_keys),
        "layout_fields_byte_value_equal": list(layout_keys),
        "panel_e_slot_width_inches": summaries["layout_inches"]["panel_e_width_inches"],
        "panel_e_slot_x_inches": summaries["layout_inches"]["panel_e_x_inches"],
    }


def raster_preservation_audit(
    candidate_png: bytes,
    panel_axis: mpl.axes.Axes,
    key_axis: mpl.axes.Axes,
    figure_dpi: float,
) -> dict[str, Any]:
    with Image.open(io.BytesIO(candidate_png)) as candidate_image:
        candidate = np.asarray(candidate_image.convert("RGBA"))
    with Image.open(PARENT_DIR / "Fig3_target_centred_top_overview_layout_v10_4_1.png") as parent_image:
        parent = np.asarray(parent_image.convert("RGBA"))
    if candidate.shape != parent.shape:
        raise RuntimeError(f"candidate/parent PNG shape drift: {candidate.shape} vs {parent.shape}")
    height, width = candidate.shape[:2]
    scale = 300.0 / figure_dpi

    def top_down_box(axis: mpl.axes.Axes, pad: int = 0) -> tuple[int, int, int, int]:
        box = axis.get_window_extent()
        x0 = max(0, int(math.floor(box.x0 * scale)) - pad)
        x1 = min(width, int(math.ceil(box.x1 * scale)) + pad)
        y0 = max(0, height - int(math.ceil(box.y1 * scale)) - pad)
        y1 = min(height, height - int(math.floor(box.y0 * scale)) + pad)
        return x0, y0, x1, y1

    panel_box = top_down_box(panel_axis, pad=16)
    mask = np.ones((height, width), dtype=bool)
    x0, y0, x1, y1 = panel_box
    mask[y0:y1, x0:x1] = False
    different = np.any(candidate != parent, axis=2)
    outside_differences = int(np.count_nonzero(different & mask))
    if outside_differences:
        raise RuntimeError(
            f"candidate differs from v10.4.1 outside panel-e mask: {outside_differences} pixels"
        )

    key_box = top_down_box(key_axis, pad=0)
    kx0, ky0, kx1, ky1 = key_box
    key_differences = int(
        np.count_nonzero(
            np.any(
                candidate[ky0:ky1, kx0:kx1] != parent[ky0:ky1, kx0:kx1], axis=2
            )
        )
    )
    if key_differences:
        raise RuntimeError(f"c/d key raster changed from v10.4.1: {key_differences} pixels")
    return {
        "png_shape": list(candidate.shape),
        "panel_e_difference_mask_pixels_top_down": list(panel_box),
        "outside_panel_e_changed_pixels": outside_differences,
        "c_d_key_crop_pixels_top_down": list(key_box),
        "c_d_key_changed_pixels": key_differences,
        "all_other_rendered_pixels_identical_to_v10_4_1": True,
        "c_d_key_pixels_identical_to_v10_4_1": True,
    }


def write_frame(path: Path, frame: pd.DataFrame) -> None:
    if path.suffix == ".tsv":
        frame.to_csv(path, sep="\t", index=False)
    else:
        frame.to_csv(path, index=False)


def build(output_dir: Path) -> None:
    global V10, PANEL_ROWS, PANEL_SOURCE_FRAME, PANEL_VALIDATION
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite retained candidate: {output_dir}")
    if list(output_dir.parent.glob(f".{output_dir.name}.staging*")):
        raise FileExistsError("refusing to reuse retained staging material")

    protected = protected_paths()
    protected_before = file_hashes(protected)
    PANEL_ROWS, PANEL_SOURCE_FRAME, PANEL_VALIDATION = validate_panel_source()
    parent_manifest = json.loads((PARENT_DIR / "manifest.json").read_text(encoding="utf-8"))
    if parent_manifest.get("schema_version") != "fig3-target-centred-top-overview-layout-v10.4.1":
        raise RuntimeError("parent v10.4.1 candidate schema drift")

    V10 = load_module(PARENT_RENDERER, "fig3_v11_parent")
    V10.draw_panel_e_nonzero_go_slim = draw_panel_e_alltarget_terms
    V10.audit_reported_text_relationships = audit_replacement_text_relationships
    v9 = V10.load_module(V10.V9_RENDERER, "fig3_v11_v9")
    v7 = V10.load_module(V10.V7_RENDERER, "fig3_v11_v7")
    bio = V10.load_module(V10.BIO_RENDERER, "fig3_v11_bio")
    dense = v7.load_dense_sources()
    v6 = v7.load_module(v7.V6_RENDERER, "fig3_v11_v6")
    bio.validate_sources()
    inherited_rows, _ = bio.read_rows()

    first_fig, first_summaries, first_frames, first_handles = V10.build_figure(
        v9, v7, v6, dense, inherited_rows
    )
    geometry_identity = assert_parent_geometry(first_summaries, parent_manifest)
    first_pdf, first_png, first_summaries = render_to_bytes(
        v7, first_fig, first_summaries, first_handles
    )
    raster_identity = raster_preservation_audit(
        first_png,
        first_handles["panel_e_axes"]["parent"],
        first_handles["panel_e_axes"]["panel_d_key"],
        first_fig.dpi,
    )
    plt.close(first_fig)

    second_fig, second_summaries, second_frames, second_handles = V10.build_figure(
        v9, v7, v6, dense, inherited_rows
    )
    assert_parent_geometry(second_summaries, parent_manifest)
    second_pdf, second_png, second_summaries = render_to_bytes(
        v7, second_fig, second_summaries, second_handles
    )
    plt.close(second_fig)
    if first_pdf != second_pdf or first_png != second_png:
        raise RuntimeError("fresh-canvas deterministic rerender bytes differ")
    if first_summaries != second_summaries:
        raise RuntimeError("fresh-canvas semantic/geometry summaries differ")
    if set(first_frames) != set(second_frames):
        raise RuntimeError("fresh-canvas inherited source-frame names differ")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging_", dir=output_dir.parent))
    pdf = staging / f"{STEM}.pdf"
    png = staging / f"{STEM}.png"
    pdf.write_bytes(first_pdf)
    png.write_bytes(first_png)
    saved_page = V10.pdf_page_geometry(pdf)

    source_frames = dict(first_frames)
    source_frames.pop("source_panel_e_go_bp_slim_summary.tsv")
    source_frames["source_panel_e_alltarget_term_composition.tsv"] = PANEL_SOURCE_FRAME
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

    unchanged_sources: dict[str, dict[str, Any]] = {}
    changed_source_name = "source_panel_e_alltarget_term_composition.tsv"
    for current in sorted(staging.glob("source_panel_*")):
        if current.name == changed_source_name:
            continue
        parent_source = PARENT_DIR / current.name
        if not parent_source.is_file():
            raise FileNotFoundError(f"parent source missing: {parent_source}")
        current_hash = sha256(current)
        parent_hash = sha256(parent_source)
        if current_hash != parent_hash:
            raise RuntimeError(f"unchanged panel source drift: {current.name}")
        unchanged_sources[current.name] = {
            "v10_4_1_sha256": parent_hash,
            "v11_sha256": current_hash,
            "byte_identical": True,
        }

    copied: dict[str, dict[str, Any]] = {}

    def copy_closed(source: Path, destination_name: str) -> None:
        destination = staging / destination_name
        shutil.copyfile(source, destination)
        copied[destination_name] = {
            "source": relative(source),
            "sha256": sha256(destination),
            "bytes": destination.stat().st_size,
        }

    copy_closed(PARENT_DIR / "manifest.json", "parent_v10_4_1_manifest.json")
    copy_closed(PARENT_RENDERER, "parent_renderer_v10_4_1.py")
    copy_closed(OVERLAY_JSON, "input_panel_e_overlay_v2.json")
    copy_closed(AUTHORITATIVE_TERMS, "input_panel_e_authoritative_terms.tsv")
    copy_closed(FINAL_AUDIT, "input_panel_e_supplement_package_final_audit.json")
    copy_closed(OVERLAY_GENERATOR, "source_panel_e_overlay_generator.py")
    copy_closed(OVERLAY_RENDERER, "source_panel_e_standalone_renderer.py")
    copy_closed(Path(__file__).resolve(), "renderer_fig3_target_centred_alltarget_enrichment_v11.py")

    protected_after = file_hashes(protected)
    if protected_after != protected_before:
        changed = [
            path for path, before in protected_before.items()
            if protected_after.get(path) != before
        ]
        raise RuntimeError(f"protected artifacts changed during v11 build: {changed}")

    outputs: dict[str, dict[str, Any]] = {}
    for path in sorted(staging.glob("source_panel_*")) + [pdf, png]:
        outputs[path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "decision_scope": (
            "Replace only Figure 3 panel e with the corrected all-target v2 composition overlay; "
            "preserve panels a-d and f-h and the c/d key. No canonical or active-text write."
        ),
        "renderer": {
            "path": "scripts/fig3_target_centred_alltarget_enrichment_v11.py",
            "sha256": sha256(Path(__file__).resolve()),
        },
        "parent": {
            "directory": relative(PARENT_DIR),
            "manifest_sha256": sha256(PARENT_DIR / "manifest.json"),
            "pdf_sha256": sha256(PARENT_DIR / "Fig3_target_centred_top_overview_layout_v10_4_1.pdf"),
            "png_sha256": sha256(PARENT_DIR / "Fig3_target_centred_top_overview_layout_v10_4_1.png"),
            "renderer_sha256": sha256(PARENT_RENDERER),
        },
        "layout_order": {
            "top": ["a", "b", "c"],
            "middle": ["d", "e"],
            "bottom": ["f", "g", "h"],
        },
        "panel_map": {
            "a": "unchanged target concentration",
            "b": "unchanged single-orientation target statistic",
            "c": "unchanged complete R49 x T279 fixed-order signed effect overview",
            "d": "unchanged literal panel c[:, :80] enlargement",
            "e": "corrected v2 all-target enriched-term composition; 9 KEGG plus 10 GO-BP",
            "f": "unchanged cofactor-trio direction-concordant targets",
            "g": "unchanged sensor-pair direction-concordant targets",
            "h": "unchanged Arce edge recurrence",
        },
        "panel_e": first_summaries["panel_e"],
        "claim_boundary": first_summaries["panel_e"]["claim_boundary"],
        "scientific_caption_boundary": (
            "Panel e uses all 9,554 census targets, including 6,977 with no observed reversal, "
            "not the 279 targets in panels c-d. Nine of 168 KEGG Legacy terms met the repeatable "
            "within-library FDR criterion in the pooled permutations and both independent banks; "
            "none passed global maxT across all 2,593 terms. Ten of 2,425 GO Biological Process "
            "terms passed global maxT in the pooled analysis and both banks. Bars descriptively "
            "partition term members into lost, gained, mixed, or no observed reversal; the "
            "enrichment test pools direction. Terms overlap and do not establish pathway activity."
        ),
        "unchanged_content_proof": {
            "panel_sources": unchanged_sources,
            "geometry": geometry_identity,
            "saved_png_raster": raster_identity,
            "panels_a_d_and_f_h_preserved": True,
            "c_d_key_preserved": True,
        },
        "geometry": first_summaries,
        "saved_pdf_geometry": saved_page,
        "export": {
            "page_inches": [V10.FIGURE_WIDTH_IN, V10.FIGURE_HEIGHT_IN],
            "pdf_fonttype": 42,
            "png_dpi": 300,
            "tight_bbox": False,
            "fresh_canvas_count": 2,
            "pdf_byte_identical": True,
            "png_byte_identical": True,
            "pdf_sha256": hash_bytes(first_pdf),
            "png_sha256": hash_bytes(first_png),
        },
        "source_closure": {
            "overlay_json_sha256": sha256(OVERLAY_JSON),
            "authoritative_terms_sha256": sha256(AUTHORITATIVE_TERMS),
            "supplement_package_final_audit_sha256": sha256(FINAL_AUDIT),
            "census_master_edges_sha256": sha256(CENSUS),
            "overlay_generator_sha256": sha256(OVERLAY_GENERATOR),
            "standalone_overlay_renderer_sha256": sha256(OVERLAY_RENDERER),
            "panel_source_rows": len(PANEL_SOURCE_FRAME),
            "all_display_values_retained_in_source_panel": True,
            "copied_inputs_and_renderers": copied,
        },
        "outputs": outputs,
        "no_write_guard": {
            "paired_before_after_hashes": True,
            "protected_before": protected_before,
            "protected_after": protected_after,
            "all_equal": True,
        },
        "canonical_changes": False,
        "manuscript_changes": False,
        "legend_changes": False,
        "docx_changes": False,
        "active_claim_guard_changes": False,
        "prior_artifact_changes": False,
        "promotion_performed": False,
        "permanent_deletion_performed": False,
    }
    (staging / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(staging, output_dir)
    print(f"wrote {output_dir}")
    print(f"  PDF sha256 {sha256(output_dir / f'{STEM}.pdf')}")
    print(f"  PNG sha256 {sha256(output_dir / f'{STEM}.png')}")
    print(f"  manifest sha256 {sha256(output_dir / 'manifest.json')}")
    print("  panel e: 9 KEGG repeatable-FDR terms (0 maxT) + 10 GO-BP maxT-stable terms")
    print("  panels a-d/f-h and c/d key unchanged; canonical and active text unchanged")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    build(args.output_dir)


if __name__ == "__main__":
    main()
