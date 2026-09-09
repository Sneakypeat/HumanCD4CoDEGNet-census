#!/usr/bin/env python3
"""Render the Figure 4 v6 presentation with its Figure 3 cross-reference relabelled.

This is a presentation-only diagnostic successor.  It delegates the complete
render to the canonical Figure 4 renderer and changes exactly one visible text
artist: ``◆ also in Fig. 3d trio`` becomes ``◆ also in Fig. 3e trio``.  The
scientific data, panel geometry, typography, colours and all other text remain
unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt

import fig4_horizontal_transposed_diagnostic_v4 as PREDECESSOR
import fig4_production_v5_1 as P


DEFAULT_OUTPUT_DIR = (
    P.PAPER
    / "figures_final/diagnostics/fig4_horizontal_transposed_candidate_v7_fig3e_crossref"
)
OLD_LABEL = "◆ also in Fig. 3d trio"
NEW_LABEL = "◆ also in Fig. 3e trio"
PREDECESSOR_RENDERER = (
    P.PAPER / "scripts/response_shape/fig4_horizontal_transposed_diagnostic_v4.py"
)
PREDECESSOR_BUNDLE = (
    P.PAPER
    / "figures_final/diagnostics/fig4_horizontal_transposed_candidate_v6_detection_05bbaa"
)
CANONICAL_FILES = (
    P.PAPER / "figures_final/Fig4.pdf",
    P.PAPER / "figures_final/Fig4.png",
    P.PAPER / "figures_final/Fig4.presentation.json",
    P.PAPER / "figures_final/Fig4.panel_f_source.tsv",
)
_ORIGINAL_PANEL_TARGET_DOMINANCE = PREDECESSOR.panel_target_dominance_wide


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    """Write JSON atomically inside an already-created diagnostic bundle."""
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp",
        delete=False, encoding="utf-8",
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def panel_target_dominance_fig3e(
    ax: plt.Axes,
    bundle: P.ProductionBundle,
    n_targets: int = 28,
    min_edges: int = 3,
) -> None:
    """Delegate panel e, then relabel its sole Figure 3 cross-reference."""
    _ORIGINAL_PANEL_TARGET_DOMINANCE(
        ax, bundle, n_targets=n_targets, min_edges=min_edges
    )
    matches = [text for text in ax.texts if text.get_text() == OLD_LABEL]
    if len(matches) != 1:
        raise AssertionError(
            f"expected one predecessor cross-reference, found {len(matches)}"
        )
    matches[0].set_text(NEW_LABEL)
    if any(text.get_text() == OLD_LABEL for text in ax.texts):
        raise AssertionError("predecessor Figure 3d label remains visible")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty successor directory: {output_dir}"
        )

    canonical_before = {
        str(path.relative_to(P.PAPER)): sha256(path) for path in CANONICAL_FILES
    }
    predecessor_before = {
        str(path.relative_to(P.PAPER)): sha256(path)
        for path in (
            PREDECESSOR_RENDERER,
            PREDECESSOR_BUNDLE / "Fig4.pdf",
            PREDECESSOR_BUNDLE / "Fig4.png",
            PREDECESSOR_BUNDLE / "Fig4.presentation.json",
            PREDECESSOR_BUNDLE / "Fig4.panel_f_source.tsv",
        )
    }

    PREDECESSOR.panel_target_dominance_wide = panel_target_dominance_fig3e
    predecessor_argv = sys.argv
    try:
        sys.argv = [str(Path(__file__).resolve()), "--output-dir", str(output_dir)]
        PREDECESSOR.main()
    finally:
        sys.argv = predecessor_argv
        PREDECESSOR.panel_target_dominance_wide = _ORIGINAL_PANEL_TARGET_DOMINANCE

    presentation_path = output_dir / "Fig4.presentation.json"
    record = json.loads(presentation_path.read_text(encoding="utf-8"))
    record.update({
        "schema": "response-shape-figure4-fig3e-crossref-diagnostic-v7",
        "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
        "canonical_figure_modified": False,
        "scientific_panels_or_data_changed": False,
        "layout_change_only": False,
        "presentation_only_change": True,
        "cross_reference_only_successor": {
            "old_visible_label": OLD_LABEL,
            "new_visible_label": NEW_LABEL,
            "reason": (
                "Figure 3 merged its overview and enlargement as panel c and "
                "relettered the former cofactor-trio panel d as panel e"
            ),
            "all_other_visible_content_intended_unchanged": True,
        },
        "diagnostic_wrapper": {
            "path": str(Path(__file__).resolve().relative_to(P.PAPER)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    })
    connector = dict(record.get("connector", {}))
    connector.update({
        "visible_key": NEW_LABEL,
        "role": (
            "backward descriptive display cross-reference from Figure 4e to "
            "the Figure 3e cofactor-trio intersection"
        ),
    })
    record["connector"] = connector
    atomic_json(presentation_path, record)

    renderer_copy = output_dir / Path(__file__).name
    shutil.copyfile(Path(__file__).resolve(), renderer_copy)

    canonical_after = {
        str(path.relative_to(P.PAPER)): sha256(path) for path in CANONICAL_FILES
    }
    predecessor_after = {
        str(path.relative_to(P.PAPER)): sha256(path)
        for path in (
            PREDECESSOR_RENDERER,
            PREDECESSOR_BUNDLE / "Fig4.pdf",
            PREDECESSOR_BUNDLE / "Fig4.png",
            PREDECESSOR_BUNDLE / "Fig4.presentation.json",
            PREDECESSOR_BUNDLE / "Fig4.panel_f_source.tsv",
        )
    }
    if canonical_after != canonical_before:
        raise AssertionError("canonical Figure 4 changed during diagnostic render")
    if predecessor_after != predecessor_before:
        raise AssertionError("sealed predecessor changed during diagnostic render")

    output_files = (
        output_dir / "Fig4.pdf",
        output_dir / "Fig4.png",
        output_dir / "Fig4.presentation.json",
        output_dir / "Fig4.panel_f_source.tsv",
        renderer_copy,
    )
    manifest = {
        "schema": "response-shape-figure4-fig3e-crossref-successor-v1",
        "status": "COMPLETE_DIAGNOSTIC_ONLY_NOT_PROMOTED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "change_scope": "one visible cross-reference label only",
        "old_visible_label": OLD_LABEL,
        "new_visible_label": NEW_LABEL,
        "scientific_panels_or_data_changed": False,
        "layout_fonts_colours_changed": False,
        "canonical_before_and_after_identical": canonical_before,
        "predecessor_before_and_after_identical": predecessor_before,
        "outputs": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in output_files
        },
        "verification_pending": "independent deterministic replay and visual audit",
    }
    atomic_json(output_dir / "manifest.json", manifest)
    print(json.dumps({
        "status": manifest["status"],
        "output_dir": str(output_dir),
        "pdf_sha256": manifest["outputs"]["Fig4.pdf"]["sha256"],
        "png_sha256": manifest["outputs"]["Fig4.png"]["sha256"],
        "old_visible_label": OLD_LABEL,
        "new_visible_label": NEW_LABEL,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
