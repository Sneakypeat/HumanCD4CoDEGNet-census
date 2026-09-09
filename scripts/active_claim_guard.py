#!/usr/bin/env python3
"""Fail closed when active paper sources drift from the canonical estimands.

Historical/provenance files may retain superseded analyses when explicitly labelled. This guard
targets only sources that feed the manuscript, submission claims, or supplementary package.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import re
import subprocess
import unicodedata
import zipfile
from collections import Counter
from pathlib import Path

from fig3_unified_zoom_guard_v13 import (
    audit_candidate_v13,
    check_fig3_unified_zoom_v13_promotion,
)
from fig4_fig3e_crossref_guard_v7 import (
    audit_candidate_v7,
    check_fig4_fig3e_crossref_promotion,
)
from legend_concision_guard_v1 import check_legend_concision_successor

ROOT = Path(__file__).resolve().parents[1]

LEGEND_CONCISION_GUARD = (
    "scripts/legend_concision_guard_v1.py",
    "5d8509a10ece69d02682e7c77dc84aefb5642cffeff16dfa8fe0cdb2ec1ceeab",
)

ACTIVE_FIGURE_SUCCESSORS = {
    "active_figure_3": {
        "figure": "Figure 3 v13",
        "promotion_path": (
            "data/census/"
            "fig3_target_centred_unified_zoom_v13_promotion_20260906.json"
        ),
        "promotion_schema": "fig3-target-centred-unified-zoom-v13-promotion-v1",
        "candidate_directory": (
            "data/figure_provenance/fig3_target_centred_unified_zoom_v13"
        ),
        "candidate_manifest": (
            "data/figure_provenance/fig3_target_centred_unified_zoom_v13/"
            "manifest.json"
        ),
        "candidate_manifest_sha256": (
            "f93af0ff8c29a9a7c2861f0783e4f8a10e66e60b83fa12c6d9725bcfc270ca06"
        ),
        "builder": "scripts/fig3_target_centred_unified_zoom_v13.py",
        "builder_sha256": (
            "5839e53bec9d09290b98988a0d963c017594c3692cf7954620c570a17a619ca5"
        ),
        "guard": "scripts/fig3_unified_zoom_guard_v13.py",
        "promoter": "scripts/promote_fig3_unified_zoom_v13.py",
        "canonical": {
            "pdf": (
                "figures/Fig3.pdf",
                "3b1868d0948b066072c58213d8b7838b424d900955bdf9d42a67898f8144776d",
            ),
            "png": (
                "figures/Fig3.png",
                "7195775b051b965acceb1397f3cea932341b75582735ec561636b06fb78851ae",
            ),
        },
    },
    "active_figure_4": {
        "figure": "Figure 4 v7",
        "promotion_path": (
            "data/census/fig4_fig3e_crossref_v7_promotion_20260906.json"
        ),
        "promotion_schema": "fig4-fig3e-crossref-v7-promotion-v1",
        "candidate_directory": (
            "data/figure_provenance/"
            "fig4_horizontal_transposed_candidate_v7_fig3e_crossref"
        ),
        "candidate_manifest": (
            "data/figure_provenance/"
            "fig4_horizontal_transposed_candidate_v7_fig3e_crossref/manifest.json"
        ),
        "candidate_manifest_sha256": (
            "157e41a599a7f71ccc111c8827c01c2934b2d869cd01b7176a9f39a1bfca5b76"
        ),
        "builder": (
            "scripts/response_shape/"
            "fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7.py"
        ),
        "builder_sha256": (
            "ca75dd009b0fc3395020e62e5bd08e63cbbb7d1e886e5a1037b511cf3ae5c06a"
        ),
        "guard": "scripts/fig4_fig3e_crossref_guard_v7.py",
        "promoter": "scripts/promote_fig4_fig3e_crossref_v7.py",
        "canonical": {
            "pdf": (
                "figures/Fig4.pdf",
                "88001c154dda2213b344f3d37ca5886a24fb28b627672a91e889203238582c66",
            ),
            "png": (
                "figures/Fig4.png",
                "d621a797a9e0eba979868c0f66a5f103e569a0ab8c70bd15ca32860652cdbcce",
            ),
            "presentation": (
                "figures/Fig4.presentation.json",
                None,
            ),
            "panel_f_source": (
                "figures/Fig4.panel_f_source.tsv",
                "63f055b7584248aba9495df7fd546980973a7e4a650d1f398eed9f5a25c0b521",
            ),
        },
    },
}

HISTORICAL_FIGURE_PROMOTIONS = {
    "figure_3_v12": (
        "data/census/"
        "fig3_target_centred_compact_zoom_v12_promotion_20260906.json",
        "00d5415a80f66c6465851820ff1a69a53dd122018ec9bf70be612bfaf4541538",
    ),
    "figure_4_cyan_v2": (
        "data/census/fig4_horizontal_transposed_cyan_promotion_v2.json",
        "ed233e1abec6b13623223540620613bc87750e771db3900a8c11d9f5b4b59c4f",
    ),
}

FIG2_DIAGNOSTIC_RECONCILIATION = (
    "data/census/fig2_eight_panel_live_diagnostic_reconciliation_20260831.json"
)
FIG2_DIAGNOSTIC_RECONCILIATION_SHA256 = (
    "aa385a7374633b78f8bfbfcea9340e4abbab5f9e0c74795527d903ffdbe07cf9"
)

FIG2_SIGNED_KEGG_V11_PROMOTION = (
    "data/census/fig2_signed_kegg_v11_overlap_clearance_promotion_20260906.json"
)
FIG2_SIGNED_KEGG_V11_CANDIDATE = (
    "data/figure_provenance/"
    "fig2_signed_kegg_side_by_side_candidate_v11_final_review"
)
FIG2_SIGNED_KEGG_V11_EXPECTED = {
    "pdf": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/"
        "Fig2_signed_kegg_side_by_side_candidate_v11_final_review.pdf",
        "c0f9a50eb371dbff80ef11313bd9dbefd5d66fed87eddce8859fd0c863dc3bc5",
    ),
    "png": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/"
        "Fig2_signed_kegg_side_by_side_candidate_v11_final_review.png",
        "8a096a13fe2ab1b21057e214bf9b99c5a91f644a875f40fa6033116cdb4768c5",
    ),
    "manifest": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/manifest.json",
        "9992863bbb3cff21f21fc26c905f878908af11a835ae6e10d06fefcbe4c3baf2",
    ),
    "provenance": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/PROVENANCE.md",
        "0aacad796be168beba4d558057abe53ae5ac800c7880ab8cfaf452808c2e67b4",
    ),
    "render_guard": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/render_guard.log",
        "c393fd2b5a39c61275aef7162ea964f47c943f83b64eb4cb1eb929984f932cef",
    ),
    "determinism_guard": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/determinism_guard.log",
        "3293f3177359c9fcdc98b658a652d395b8183dd4c89e0fecadf130abc475c217",
    ),
    "determinism_receipt": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/determinism_receipt.json",
        "c371ec7357fd4f523cfe6e07afce3070eb6d8600e343bdb5aebe2b9b7b4a1a4b",
    ),
    "panel_b_source": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/source_panel_b_regulators.tsv",
        "4013602f4d0284a33fdc4e9ecb201302837b1c15020709efd129fdefeb02db4e",
    ),
    "panel_d_source": (
        f"{FIG2_SIGNED_KEGG_V11_CANDIDATE}/source_panel_d_signed_kegg.tsv",
        "240f4fa497685be13a433584d5fe046b0e592e964865ca5d3c8aaa16ba95afb0",
    ),
    "renderer": (
        "scripts/fig2_signed_kegg_side_by_side_candidate_v11.py",
        "c1dabf58102a0afa26ee04a81a3fe0b742c218b2f3bc8954afc3d5606b17bbc2",
    ),
    "audit_script": (
        "scripts/audit_fig2_signed_kegg_side_by_side_candidate_v11.py",
        "d238d3994db73961a643c7eb39281687bdf048f600be9fb8221e5b0b380d9ab9",
    ),
}

FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_SCHEMA = (
    "fig2-v11-active-guard-normalization-reseal-v1"
)
FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_SCRIPT = (
    "scripts/reseal_fig2_signed_kegg_v11_guard_normalization.py"
)
FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE = (
    "figures/repository/superseded_figure2_presentations_2026-09-06/"
    "fig2_signed_kegg_v11_pre_guard_normalization_reseal_20260906"
)
FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS = {
    "prior_promotion_manifest": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "promotion.sealed.pre_guard_normalization.json"
        ),
        "sha256": "7e70417c1c98b9f8be5e7e685f3f288d5f464bf72be64d54ce55813bd5abb72b",
        "bytes": 28864,
    },
    "prior_status_registry": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "status.sealed.pre_guard_normalization.json"
        ),
        "sha256": "81a3859c8a3a2f17ad50e91d07c40a08f18530ea4eb6f42a6af35d121ec291e5",
        "bytes": 25026,
    },
    "active_guard_before": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "active_claim_guard.pre_normalization.py"
        ),
        "sha256": "6e638e67606d04afaf90ab6bcf9156dfdbffd3e3a6b3d145955f454b02417d8b",
        "bytes": 254053,
    },
    "prior_docx": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "HumanCD4CoDEGNet_signflip.pre_guard_normalization.docx"
        ),
        "sha256": "eca5fb58184914ae445348328a1a57271d6df1ffccdb53bb886bdbd9c79971f8",
        "bytes": 8207999,
    },
    "prior_pdf": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "HumanCD4CoDEGNet_signflip.pre_guard_normalization.pdf"
        ),
        "sha256": "30484585314515e7f0031d743f4653c05f3b4fa8f66ebef76d8d8311990c1912",
        "bytes": 3151802,
    },
    "prior_supplementary_pdf": {
        "path": (
            f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
            "Supplementary_Information.pre_guard_normalization.pdf"
        ),
        "sha256": "5cfeb96ae38b0f6c18a4d86a210b1788fd0fb80d0d04bf0624080c12aa6dc509",
        "bytes": 1143610,
    },
}

KINETICS_SIGN_STATE_BINDINGS = {
    "scripts/kinetics_8hr_sign_states.py":
        "68639d57de496a8039e120b77f9ad9e1d75850174a9f11a454095244a669b83e",
    "data/kinetics/kinetics_8hr_sign_states.csv.gz":
        "9d34bdb23acdc1a1556854bea0e148a8f1945fb5dfe14c41ad5162a4cd11c222",
    "data/kinetics/kinetics_8hr_sign_states.paired_regulators.tsv":
        "a8bc41f6745335396fdd51726011533b558f3a7ced452b67fd878174c97ec3ba",
    "data/kinetics/kinetics_8hr_sign_states.summary.json":
        "af034de06b5b0de8a414e8f8bd3c692377de5faa8a3d7899b6db83d0f2a70a9a",
    "data/kinetics/kinetics_8hr_sign_states.manifest.json":
        "2b26eb8dd424eff79a2912cd99c045ed562628f3d12848e29360c178e46ac212",
}

ACTIVE_TEXT = [
    "README.md",
    "MANUSCRIPT.md",
    "FIGURE_LEGENDS.md",
    "PAPER_FACTS.md",
    "REPLICATION_SCOPE.md",
    "COVER_LETTER_MSB.md",
    "SUPPLEMENTARY_NOTE_negatives.md",
    "SUPPLEMENTARY_METHODS.md",
    "SUPPLEMENTARY_RESULTS.md",
    "submission/COVER_LETTER.md",
    "submission/SIGNIFICANCE_STATEMENT.md",
]

FORBIDDEN = {
    r"4,379\s+(?:of|/)\s+84,044": "inclusive count used as the reversal denominator",
    r"84,044\s+trans(?:-|\s)*(?:regulatory\s+)?(?:edges|effects|interactions)":
        "inclusive count labelled as trans",
    r"trans(?:-|\s)*(?:regulatory\s+)?(?:edges|effects|interactions)[^\n.]{0,80}84,044":
        "inclusive count labelled as trans",
    r"79,665\s+(?:keep|stable|concordant)": "inclusive stable count in an active claim",
    r"83,490": "symbol-only trans count retained after Ensembl identity correction",
    r"79,111": "symbol-only stable trans count retained after Ensembl identity correction",
    r"554\s+(?:on-target|self)": "symbol-only self-edge count retained after Ensembl identity correction",
    r"9\.6\s*(?:-?fold|×|x)": "withdrawn state-shifter enrichment",
    r"8\.8\s*(?:-?fold|×|x)[^\n]{0,80}(?:state|cell-count)":
        "withdrawn state-shifter sensitivity",
    r"Figure\s+6\s*[dD]\s*[–-]\s*f": "superseded discovery-panel reference",
    r"accrue(?:s|d)?\s+late": "unsupported late-firing claim",
    r"late-versus-early\s+split.{0,120}even\s+split":
        "obsolete pooled 8 h timing test",
    r"regulator-level\s+test\s+over\s+50\s+regulators.{0,120}p\s*=\s*0\.004":
        "superseded pooled 8 h timing result",
    r"3,366\s+reversing\s+edges": "superseded Figure 3 composition coverage",
    r"0\.0117\s+log2": "superseded Figure 3 composition magnitude",
    r"0\.0045\s+log2": "superseded Figure 3 random-state floor",
    r"(?:roughly|about|~)\s*19(?:-fold|\s*[x×])":
        "superseded Figure 3 endpoint-gap estimate",
    r"decomposition\s+shows\s+(?:that\s+)?(?:the\s+)?reversals\s+are\s+(?:a\s+)?within-cell":
        "composition falsification overclaimed as a demonstrated within-cell mechanism",
    r"true\s+edge-level[^\n.]{0,100}(?:holdout|leave-one)[^\n.]{0,100}impossible":
        "stale claim that donor/guide edge holdout is impossible",
}

TARGET_TRAJECTORY_TEXT = {
    "README.md": (
        r"of 910 targets[^\n]{0,160}903 mapped",
        r"post-primary specificity analysis",
    ),
    "MANUSCRIPT.md": (
        r"631/910",
        r"of 910 eligible targets, 903 mapped",
        r"of 633 targets[^\n]{0,160}629 mapped",
        r"22 donors contribute[^\n]{0,100}except 12 h, where 19 contribute",
        r"prewritten post-primary specificity analysis",
    ),
    "FIGURE_LEGENDS.md": (
        r"all 903 mapped targets",
        r"629[^\n]{0,100}633",
        r"22 donors at every time except 12 h, where 19 contribute",
        r"prewritten but post-primary specificity analysis",
    ),
    "PAPER_FACTS.md": (
        r"631/910",
        r"910 targets[\s\S]{0,250}903[^\n]{0,40}map",
        r"633[\s\S]{0,200}629/633 map",
        r"prewritten post-primary specificity",
    ),
    "REFRAME.md": (
        r"903/910",
        r"of 633 unanimous[^\n]{0,160}629 mapped",
        r"22 donors paired at 0 and 48 h",
        r"post-primary",
    ),
}

EXPECTED_SUPPLEMENTARY_FIGURES = {
    f"Figure_S{i}.{ext}" for i in range(1, 5) for ext in ("pdf", "png")
}
EXPECTED_SUPPLEMENTARY_DATA = {
    "Supplementary_Data_1_sign_reversing_edges.csv",
    "Supplementary_Data_2_per_regulator_rates.csv",
    "Supplementary_Data_3_kinetics_8h_sign_states.csv.gz",
    "Supplementary_Data_4_arce_recurrence_edges.csv.gz",
    "Supplementary_Data_5_donor_holdout_edges.csv.gz",
    "Supplementary_Data_6_edge_state_composition.csv.gz",
    "Supplementary_Data_7_target_trajectory_targets.tsv",
    "Supplementary_Data_8_target_trajectory_donors.tsv",
    "Supplementary_Data_9_target_trajectory_timecourse.tsv",
    "Supplementary_Data_10_target_trajectory_nulls.json",
    "Supplementary_Data_11_figure3_all49_programme_counts.csv",
    "Supplementary_Data_12_figure3_all49_programme_memberships.csv",
    "Supplementary_Data_13_figure3_cofactor_shared_targets.csv",
    "Supplementary_Data_14_figure3_AHR_ARNT_shared_targets.csv",
    "Supplementary_Data_15_figureS2_sign_calibration.csv",
    "Supplementary_Data_16_figureS2_composition.csv",
    "Supplementary_Data_17_figure2h_activation_movement_regulators.tsv",
    "Supplementary_Data_18_figure2h_activation_movement_summary.json",
    "Supplementary_Data_19_broad_kegg_go_rank_enrichment_terms.tsv",
    "Supplementary_Data_20_broad_kegg_go_target_scores.tsv.gz",
    "Supplementary_Data_21_broad_kegg_go_leading_edges.tsv.gz",
}
EXPECTED_SUPPLEMENTARY_PACKAGE_FILES = {
    "README.md",
    "Supplementary_Information.pdf",
    *{f"figures/{name}" for name in EXPECTED_SUPPLEMENTARY_FIGURES},
    *{f"data/{name}" for name in EXPECTED_SUPPLEMENTARY_DATA},
    "tables/Supplementary_Table_S1.csv",
    "tables/Supplementary_Table_S1.md",
    "tables/Supplementary_Table_S2_druggable_census.csv",
    "notes/Supplementary_Methods_M1-M20.md",
    "notes/Supplementary_Results_R1-R7.md",
    "notes/Supplementary_Note_S1_negative_results.md",
}


def _normalise_text(text: str) -> str:
    """Collapse renderer-specific whitespace for resilient derived-artifact checks."""
    normalized = unicodedata.normalize("NFKC", text)
    normalized = re.sub(r"\s+", " ", normalized).strip().lower()
    normalized = re.sub(r"\s*=\s*", "=", normalized)
    return re.sub(r"\s+([;,:])", r"\1", normalized)


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _historical_binding_record(
        rel: object, expected_hash: object, label: str,
        failures: list[str]) -> None:
    """Validate an immutable manifest record without binding it to mutable live bytes."""
    if not isinstance(rel, str) or not rel or Path(rel).is_absolute() or \
            ".." in Path(rel).parts:
        failures.append(f"historical {label} has unsafe path {rel!r}")
    if not isinstance(expected_hash, str) or not re.fullmatch(
            r"[0-9a-f]{64}", expected_hash):
        failures.append(f"historical {label} has invalid SHA-256 {expected_hash!r}")


def _active_successor_registry_state(
        biology: dict, registry_key: str, *, check_derived: bool,
        failures: list[str]) -> str:
    """Return ``pending``, ``bound`` or ``invalid`` for an active figure successor.

    Pending is intentionally accepted only by the source guard so the document
    builds can refresh bytes that older Figure 2 provenance records historical-
    ly describe.  Post-build mode requires one atomic registry rebind to the two
    sealed successor manifests.
    """
    config = ACTIVE_FIGURE_SUCCESSORS[registry_key]
    label = config["figure"]
    active = biology.get(registry_key)
    if not isinstance(active, dict):
        failures.append(f"biology-first registry lacks {registry_key}")
        return "invalid"
    expected_fields = {
        "promotion_manifest_path": config["promotion_path"],
        "promotion_schema_version": config["promotion_schema"],
        "candidate_directory": config["candidate_directory"],
        "candidate_manifest_path": config["candidate_manifest"],
        "candidate_manifest_sha256": config["candidate_manifest_sha256"],
        "builder_path": config["builder"],
        "builder_sha256": config["builder_sha256"],
        "guard_helper_path": config["guard"],
        "promotion_script_path": config["promoter"],
    }
    for key, expected in expected_fields.items():
        if active.get(key) != expected:
            failures.append(
                f"{label} registry {key}={active.get(key)!r}; expected {expected!r}")

    for rel, expected_hash, item_label in (
        (config["candidate_manifest"], config["candidate_manifest_sha256"],
         "candidate manifest"),
        (config["builder"], config["builder_sha256"], "builder"),
    ):
        observed = _sha256(ROOT / rel)
        if observed != expected_hash:
            failures.append(
                f"{label} {item_label} hash {observed} != {expected_hash}")
    for key in ("guard", "promoter"):
        if not (ROOT / config[key]).is_file():
            failures.append(f"{label} {key} is missing: {config[key]}")

    historical = biology.get("historical_figure_promotions", {})
    for historical_key, (expected_path, expected_hash) in \
            HISTORICAL_FIGURE_PROMOTIONS.items():
        record = historical.get(historical_key, {}) \
            if isinstance(historical, dict) else {}
        if record.get("path") != expected_path or \
                record.get("sha256") != expected_hash or \
                record.get("status_at_supersession") != \
                "PROMOTED_AND_POSTBUILD_GUARDED" or \
                record.get("historical_without_rewrite") is not True:
            failures.append(
                f"historical figure promotion {historical_key} binding changed")
        if _sha256(ROOT / expected_path) != expected_hash:
            failures.append(
                f"historical figure promotion {historical_key} bytes changed")

    state = active.get("binding_status")
    if state == "PENDING_PROMOTION_BINDING":
        for key in (
            "promotion_manifest_sha256", "guard_helper_sha256",
            "promotion_script_sha256",
        ):
            if active.get(key) is not None:
                failures.append(f"{label} pending registry {key} must be null")
        if check_derived:
            failures.append(
                f"{label} registry is still pending; run the sealed two-figure rebind")
        return "pending"

    if state != "BOUND_POSTBUILD":
        failures.append(f"{label} has invalid binding_status {state!r}")
        return "invalid"
    if active.get("promotion_status") != "PROMOTED_AND_POSTBUILD_GUARDED":
        failures.append(f"{label} bound registry does not name a sealed promotion")
    promotion_path = ROOT / config["promotion_path"]
    promotion_hash = _sha256(promotion_path)
    if promotion_hash != active.get("promotion_manifest_sha256"):
        failures.append(
            f"{label} promotion hash {promotion_hash} != registry "
            f"{active.get('promotion_manifest_sha256')}")
    if promotion_path.is_file():
        try:
            promotion = json.loads(promotion_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"{label} promotion manifest is unreadable: {exc}")
        else:
            if promotion.get("schema_version") != config["promotion_schema"] or \
                    promotion.get("status") != "PROMOTED_AND_POSTBUILD_GUARDED":
                failures.append(f"{label} promotion manifest is not the sealed successor")

    for record_key, expected_path in (
        ("guard_helper", config["guard"]),
        ("promotion_script", config["promoter"]),
        ("canonical_renderer", config["builder"]),
    ):
        record = active.get(record_key, {})
        if not isinstance(record, dict) or record.get("path") != expected_path or \
                _sha256(ROOT / expected_path) != record.get("sha256"):
            failures.append(f"{label} bound {record_key} record changed")
    canonical = active.get("canonical_outputs", {})
    for key, (expected_path, expected_hash) in config["canonical"].items():
        record = canonical.get(key, {}) if isinstance(canonical, dict) else {}
        recorded_hash = record.get("sha256")
        if record.get("path") != expected_path or \
                not isinstance(recorded_hash, str) or \
                not re.fullmatch(r"[0-9a-f]{64}", recorded_hash) or \
                (expected_hash is not None and recorded_hash != expected_hash) or \
                _sha256(ROOT / expected_path) != recorded_hash:
            failures.append(f"{label} bound canonical {key} record changed")
    return "bound"


def _pdf_text(path: Path, failures: list[str]) -> str:
    if not path.is_file():
        failures.append(f"{path.relative_to(ROOT)}: missing derived PDF")
        return ""
    try:
        proc = subprocess.run(
            ["pdftotext", str(path), "-"], capture_output=True, text=True, check=False)
    except FileNotFoundError:
        failures.append("post-build inspection requires pdftotext on PATH")
        return ""
    if proc.returncode != 0:
        failures.append(
            f"{path.relative_to(ROOT)}: pdftotext failed with exit {proc.returncode}")
        return ""
    return _normalise_text(proc.stdout)


def _docx_text(path: Path, failures: list[str]) -> str:
    if not path.is_file():
        failures.append(f"{path.relative_to(ROOT)}: missing derived DOCX")
        return ""
    try:
        with zipfile.ZipFile(path) as archive:
            xml = archive.read("word/document.xml").decode("utf-8")
    except (KeyError, OSError, zipfile.BadZipFile, UnicodeDecodeError) as exc:
        failures.append(f"{path.relative_to(ROOT)}: cannot inspect DOCX: {exc}")
        return ""
    # Word stores ordinary text in ``w:t`` and equation text in OMML ``m:t`` runs.
    # Extract both in document order so a styled mathematical p-value is not silently
    # dropped before claim checking.
    text_runs = re.findall(
        r"<(?P<namespace>w|m):t(?:\s[^>]*)?>(.*?)</(?P=namespace):t>", xml)
    return _normalise_text(" ".join(html.unescape(value) for _, value in text_runs))


def _docx_media_hashes(path: Path, failures: list[str]) -> set[str]:
    """Return hashes for every embedded Word media object, failing closed on bad ZIPs."""
    if not path.is_file():
        failures.append(f"{path.relative_to(ROOT)}: missing derived DOCX")
        return set()
    try:
        with zipfile.ZipFile(path) as archive:
            return {
                hashlib.sha256(archive.read(name)).hexdigest()
                for name in archive.namelist()
                if name.startswith("word/media/") and not name.endswith("/")
            }
    except (OSError, zipfile.BadZipFile) as exc:
        failures.append(f"{path.relative_to(ROOT)}: cannot inspect DOCX media: {exc}")
        return set()


def _docx_picture_extents(
        path: Path, failures: list[str]) -> dict[str, tuple[int, int]]:
    """Return OOXML drawing extents keyed by the embedded source basename."""
    if not path.is_file():
        failures.append(f"{path.relative_to(ROOT)}: missing derived DOCX")
        return {}
    try:
        with zipfile.ZipFile(path) as archive:
            xml = archive.read("word/document.xml").decode("utf-8")
    except (KeyError, OSError, zipfile.BadZipFile, UnicodeDecodeError) as exc:
        failures.append(f"{path.relative_to(ROOT)}: cannot inspect DOCX extents: {exc}")
        return {}
    extents = {}
    for block in re.findall(r"<wp:inline>.*?</wp:inline>", xml):
        source = re.search(r'<pic:cNvPr[^>]+descr="([^"]+)"', block)
        extent = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"\s*/>', block)
        if source is not None and extent is not None:
            extents[Path(source.group(1)).name] = tuple(map(int, extent.groups()))
    return extents


def _check_required_phrases(
        rel: str, text: str, phrases: tuple[str, ...], failures: list[str]) -> None:
    for phrase in phrases:
        if _normalise_text(phrase) not in text:
            failures.append(f"{rel}: derived artifact lacks {phrase!r}")


def _check_derived_artifacts(failures: list[str]) -> None:
    """Inspect generated deliverables only after a build has completed.

    This is deliberately excluded from the default source guard so stale generated files cannot
    prevent the build that is needed to refresh them.
    """
    docx_path = ROOT / "HumanCD4CoDEGNet_signflip.docx"
    main_phrases = (
        "532 of the 4,379 census reversals",
        "170 regulators",
        "382 targets",
        "all 49 enriched regulators are shown against these 279 targets",
        "across all 9,554 census targets",
        "nine of 168 kegg legacy pathways",
        "ten of 2,425 go biological process terms",
        "none passed global maxt correction across all 2,593 terms",
        "the enrichment test pools reversal direction",
        "all 17 reversal targets shared by ahr and arnt",
        "48.1% of 241 catalogue reversals",
        "22.0% of 3,256 stable edges",
        "mantel–haenszel odds ratio of 3.23",
        "eight of the 16 exhaustive ccnc",
        "at the cutoff, 29 targets tie for eight display positions",
        "prewritten post-primary specificity analysis",
        "extended the comparison to 903 mapped targets carrying at least two reversals and at least one stable edge",
        "631/910",
        "of 910 eligible targets, 903 mapped",
        "of 633 targets whose reversals all shared one orientation, 629 mapped",
        "22 donors contribute at every plotted time except 12 h, where 19 contribute",
        "529/704 (75.1%)",
        "175/704 (24.9%)",
        "295/575 (51.3%)",
        "280/575 (48.7%)",
        "3,098 remained unresolved",
        "rest and 48 h defined the reversal set and 8 h classified it afterwards",
    )
    for rel in ("HumanCD4CoDEGNet_signflip.pdf", "HumanCD4CoDEGNet_signflip.docx"):
        path = ROOT / rel
        text = _pdf_text(path, failures) if path.suffix == ".pdf" else _docx_text(path, failures)
        if text:
            _check_required_phrases(rel, text, main_phrases, failures)
            if re.search(r"supplementary figure(?:s)? s[6-9]", text):
                failures.append(f"{rel}: derived manuscript cites a retired supplementary figure")

    picture_extents = _docx_picture_extents(docx_path, failures)
    expected_picture_names = {
        "Fig1.png", "Fig2.png", "Fig3.png", "Fig4.png",
        "FigS1.png", "FigS2.png", "FigS3.png", "FigS4.png",
    }
    if set(picture_extents) != expected_picture_names:
        failures.append(
            "HumanCD4CoDEGNet_signflip.docx: expected exactly the eight canonical "
            f"figure drawings; observed {sorted(picture_extents)!r}")
    for name in expected_picture_names:
        extent = picture_extents.get(name)
        if name == "Fig4.png":
            expected_width = round(6.5 * 914400)
        elif name == "FigS4.png":
            expected_width = round(7.2 * 914400)
        else:
            expected_width = 5334000
        if extent is not None and extent[0] != expected_width:
            failures.append(
                f"HumanCD4CoDEGNet_signflip.docx: {name} width is {extent[0]} EMU; "
                f"expected {expected_width} EMU")
        if name == "FigS4.png" and extent is not None and \
                extent[1] > 9 * 914400:
            failures.append(
                f"HumanCD4CoDEGNet_signflip.docx: {name} height exceeds the "
                "9-inch physical content height")

    promotion_path = (ROOT / "data/census/fig3_biological_bridge/"
                       "promotion_manifest_20260901_scope_link_and_functional_groups.json")
    v11_promotion_path = (
        ROOT / "data/census/"
        "fig3_target_centred_alltarget_enrichment_v11_promotion_20260904.json")
    v12_promotion_path = (
        ROOT / "data/census/"
        "fig3_target_centred_compact_zoom_v12_promotion_20260906.json")
    v13_promotion_path = (
        ROOT / "data/census/"
        "fig3_target_centred_unified_zoom_v13_promotion_20260906.json")
    if promotion_path.is_file():
        promotion = json.loads(promotion_path.read_text())
        media_hashes = _docx_media_hashes(docx_path, failures)
        if v13_promotion_path.is_file():
            v13_promotion = json.loads(v13_promotion_path.read_text())
            expected_fig3_png = v13_promotion.get("canonical_outputs", {}).get(
                "png", {}).get("sha256")
        elif v12_promotion_path.is_file():
            v12_promotion = json.loads(v12_promotion_path.read_text())
            expected_fig3_png = v12_promotion.get("canonical_outputs", {}).get(
                "png", {}).get("sha256")
        elif v11_promotion_path.is_file():
            v11_promotion = json.loads(v11_promotion_path.read_text())
            expected_fig3_png = v11_promotion.get("canonical", {}).get(
                "png", {}).get("sha256")
        else:
            expected_fig3_png = promotion.get("figure3", {}).get(
                "canonical_png", {}).get("sha256")
        fig4_cyan_promotion_path = (
            ROOT / "data/census/fig4_horizontal_transposed_cyan_promotion_v2.json")
        fig4_v7_promotion_path = (
            ROOT / "data/census/fig4_fig3e_crossref_v7_promotion_20260906.json")
        fig4f_promotion_path = (
            ROOT / "data/census/fig4f_landmark_census_promotion_v1.json")
        if fig4_v7_promotion_path.is_file():
            fig4f_promotion = json.loads(fig4_v7_promotion_path.read_text())
            expected_fig4_png = fig4f_promotion.get(
                "canonical_outputs", {}).get("png", {}).get("sha256")
        elif fig4_cyan_promotion_path.is_file():
            fig4f_promotion = json.loads(fig4_cyan_promotion_path.read_text())
            expected_fig4_png = fig4f_promotion.get(
                "canonical_outputs", {}).get("png", {}).get("sha256")
        elif fig4f_promotion_path.is_file():
            fig4f_promotion = json.loads(fig4f_promotion_path.read_text())
            expected_fig4_png = fig4f_promotion.get(
                "canonical_outputs", {}).get("png", {}).get("sha256")
        else:
            expected_fig4_png = promotion.get("figure4", {}).get(
                "canonical_png", {}).get("sha256")
        if expected_fig3_png not in media_hashes:
            failures.append(
                "HumanCD4CoDEGNet_signflip.docx: promoted canonical Figure 3 PNG "
                f"{expected_fig3_png!r} is not embedded")
        if expected_fig4_png not in media_hashes:
            failures.append(
                "HumanCD4CoDEGNet_signflip.docx: promoted canonical Figure 4 PNG "
                f"{expected_fig4_png!r} is not embedded")
        build_integration = promotion.get("build_integration", {})
        if build_integration.get("status") != "COMPLETE":
            failures.append("joint Figure 3/Figure 4 promotion is not post-build complete")
        # These whole-manuscript hashes seal the historical Figure 3/4 promotion
        # build.  Later, independently guarded figure promotions legitimately
        # change the manuscript bytes, so current integration is checked through
        # the embedded canonical Figure 3/4 media and the required claim text
        # above rather than by comparing against the historical whole-file hash.
        for key, expected_path in (
                ("manuscript_pdf", "HumanCD4CoDEGNet_signflip.pdf"),
                ("manuscript_docx", "HumanCD4CoDEGNet_signflip.docx")):
            record = build_integration.get(key, {})
            if record.get("path") != expected_path:
                failures.append(
                    f"build integration {key} path {record.get('path')!r}; "
                    f"expected {expected_path!r}")
                continue
            historical_hash = record.get("sha256")
            if not isinstance(historical_hash, str) or not re.fullmatch(
                    r"[0-9a-f]{64}", historical_hash):
                failures.append(
                    f"build integration {key} has invalid historical SHA-256 "
                    f"{historical_hash!r}")

    fig2_registry_path = ROOT / "data/reframe_biology_status_20260814.json"
    fig2_promotion_path = None
    if fig2_registry_path.is_file():
        try:
            active_fig2 = json.loads(fig2_registry_path.read_text()).get(
                "active_figure_2", {})
            promotion_rel = active_fig2.get("promotion_manifest_path")
            if isinstance(promotion_rel, str) and promotion_rel:
                fig2_promotion_path = ROOT / promotion_rel
        except json.JSONDecodeError:
            failures.append("active Figure 2 registry is invalid JSON")
    if fig2_promotion_path is not None and fig2_promotion_path.is_file():
        fig2_promotion = json.loads(fig2_promotion_path.read_text())
        if fig2_promotion.get("schema_version") == \
                "fig2-signed-kegg-side-by-side-v11-promotion-v1":
            signed_kegg_phrases = (
                "a conventional binary analysis of the 2,577 targets receiving at least one reversal",
                "14 of 113 filtered kegg legacy pathways",
                "global p = 0.184; 0/113 maxt-supported",
                "descriptive sensitivity rather than evidence of pathway enrichment or activity",
                "12 lost and 43 gained reversal edges",
                "complementary margins and do not constitute a joint test",
            )
            for rel in (
                    "HumanCD4CoDEGNet_signflip.pdf",
                    "HumanCD4CoDEGNet_signflip.docx"):
                path = ROOT / rel
                text = (_pdf_text(path, failures) if path.suffix == ".pdf"
                        else _docx_text(path, failures))
                if text:
                    _check_required_phrases(
                        rel, text, signed_kegg_phrases, failures)
        expected_fig2_png = fig2_promotion.get("canonical_outputs", {}).get(
            "png", {}).get("sha256")
        if expected_fig2_png not in _docx_media_hashes(docx_path, failures):
            failures.append(
                "HumanCD4CoDEGNet_signflip.docx: promoted canonical Figure 2 PNG "
                f"{expected_fig2_png!r} is not embedded")

    si_rel = "submission/supplementary/Supplementary_Information.pdf"
    si_text = _pdf_text(ROOT / si_rel, failures)
    if si_text:
        _check_required_phrases(
            si_rel, si_text,
            ("supplementary figure s4", "supplementary_data_3_kinetics_8h_sign_states.csv.gz",
             "supplementary_data_16_figures2_composition.csv",
             "supplementary_data_17_figure2h_activation_movement_regulators.tsv",
             "supplementary_data_18_figure2h_activation_movement_summary.json",
             "supplementary_data_19_broad_kegg_go_rank_enrichment_terms.tsv",
             "supplementary_data_20_broad_kegg_go_target_scores.tsv.gz",
             "supplementary_data_21_broad_kegg_go_leading_edges.tsv.gz"),
            failures)
        if re.search(r"supplementary figure(?:s)? s[6-9]", si_text):
            failures.append(f"{si_rel}: SI contains a retired supplementary figure reference")

    readme_rel = "submission/supplementary/README.md"
    readme_path = ROOT / readme_rel
    if not readme_path.is_file():
        failures.append(f"{readme_rel}: missing derived package manifest")
    else:
        readme_text = _normalise_text(readme_path.read_text())
        _check_required_phrases(
            readme_rel, readme_text,
            ("supplementary figures s1–s4", "supplementary_data_3_kinetics_8h_sign_states.csv.gz",
             "supplementary_data_16_figures2_composition.csv",
             "supplementary_data_17_figure2h_activation_movement_regulators.tsv",
             "supplementary_data_18_figure2h_activation_movement_summary.json",
             "supplementary_data_19_broad_kegg_go_rank_enrichment_terms.tsv",
             "supplementary_data_20_broad_kegg_go_target_scores.tsv.gz",
             "supplementary_data_21_broad_kegg_go_leading_edges.tsv.gz"),
            failures)

    figures_dir = ROOT / "submission/supplementary/figures"
    actual_figures = {
        path.name for path in figures_dir.iterdir()
        if path.is_file() and re.fullmatch(r"Figure_S\d+\.(?:pdf|png)", path.name)
    } if figures_dir.is_dir() else set()
    if actual_figures != EXPECTED_SUPPLEMENTARY_FIGURES:
        failures.append(
            "submission/supplementary/figures: S1-S4 contract mismatch; "
            f"missing={sorted(EXPECTED_SUPPLEMENTARY_FIGURES - actual_figures)!r}, "
            f"unexpected={sorted(actual_figures - EXPECTED_SUPPLEMENTARY_FIGURES)!r}")

    data_dir = ROOT / "submission/supplementary/data"
    actual_data = {
        path.name for path in data_dir.iterdir()
        if path.is_file() and re.fullmatch(r"Supplementary_Data_\d+.*", path.name)
    } if data_dir.is_dir() else set()
    if actual_data != EXPECTED_SUPPLEMENTARY_DATA:
        failures.append(
            "submission/supplementary/data: Data1-21 contract mismatch; "
            f"missing={sorted(EXPECTED_SUPPLEMENTARY_DATA - actual_data)!r}, "
            f"unexpected={sorted(actual_data - EXPECTED_SUPPLEMENTARY_DATA)!r}")

    package_root = ROOT / "submission/supplementary"
    actual_package_files = {
        str(path.relative_to(package_root))
        for path in package_root.rglob("*") if path.is_file()
    } if package_root.is_dir() else set()
    if actual_package_files != EXPECTED_SUPPLEMENTARY_PACKAGE_FILES:
        failures.append(
            "submission/supplementary: exact 37-file reader-facing contract mismatch; "
            f"missing={sorted(EXPECTED_SUPPLEMENTARY_PACKAGE_FILES - actual_package_files)!r}, "
            f"unexpected={sorted(actual_package_files - EXPECTED_SUPPLEMENTARY_PACKAGE_FILES)!r}")


def check_active_sources(*, check_derived: bool = False) -> None:
    failures: list[str] = []

    legend_guard_path, legend_guard_hash = LEGEND_CONCISION_GUARD
    observed_legend_guard_hash = _sha256(ROOT / legend_guard_path)
    if observed_legend_guard_hash != legend_guard_hash:
        failures.append(
            "legend-concision guard hash "
            f"{observed_legend_guard_hash} != {legend_guard_hash}"
        )
        legend_concision_state = "invalid"
    else:
        legend_concision_state = check_legend_concision_successor(
            ROOT, check_derived=check_derived, failures=failures
        )
    legend_concision_active = legend_concision_state in {"pending", "sealed"}

    # Formal Supplementary Figure S4 is the audited broad all-target enrichment
    # display formerly numbered S5.  Bind its immutable artwork directly here,
    # require the renumbering provenance, and keep the retired S5 alias out of
    # the active namespace.  The transaction may be source-pending only until
    # the document/package rebuild has completed.
    broad_promotion_path = (
        ROOT / "data/census/broad_enrichment_s5_to_s4_promotion_v1.json"
    )
    broad_integration_path = (
        ROOT / "data/census/"
        "broad_enrichment_s5_to_s4_source_integration_v1.json"
    )
    broad_expected = {
        "figures/FigS4.pdf":
            "32cc2210083f94467c284ec83b9bfbdb2fe9e83e06c025ec2608d6411581f291",
        "figures/FigS4.png":
            "33b0d9a64be5cee4f6a8f45a6ee5a1e9055bc54c63a6f6e09555b884bf4696d5",
    }
    for rel, expected_hash in broad_expected.items():
        path = ROOT / rel
        observed_hash = (
            hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        )
        if observed_hash != expected_hash:
            failures.append(
                f"formal broad-enrichment S4 {rel!r} hash "
                f"{observed_hash!r} != {expected_hash!r}"
            )
    for suffix in ("pdf", "png", "provenance.json"):
        rel = f"figures/FigS{5}.{suffix}"
        if (ROOT / rel).exists():
            failures.append(f"retired formal S5 remains active: {rel}")

    broad_provenance_path = ROOT / "figures/FigS4.provenance.json"
    if not broad_provenance_path.is_file():
        failures.append("formal broad-enrichment S4 provenance is missing")
    else:
        try:
            broad_provenance = json.loads(broad_provenance_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"formal broad-enrichment S4 provenance is unreadable: {exc}")
        else:
            if broad_provenance.get("schema_version") != \
                    "formal-supplementary-figure-s4-broad-gsea-v2-renumber-v1":
                failures.append("formal broad-enrichment S4 provenance schema changed")
            if broad_provenance.get("status") != \
                    "PROMOTED_BY_BYTE_IDENTICAL_S5_TO_S4_RENUMBER":
                failures.append("formal broad-enrichment S4 provenance status changed")
            if broad_provenance.get("canonical_outputs") != {
                    "FigS4.pdf": broad_expected["figures/FigS4.pdf"],
                    "FigS4.png": broad_expected["figures/FigS4.png"],
            }:
                failures.append("formal broad-enrichment S4 provenance output hashes changed")
            if broad_provenance.get("permanent_deletion_performed") is not False:
                failures.append("formal broad-enrichment S4 lost its no-delete contract")

    if not broad_integration_path.is_file():
        failures.append("broad-enrichment S5-to-S4 source integration is missing")
    else:
        try:
            broad_integration = json.loads(broad_integration_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"broad-enrichment source integration is unreadable: {exc}")
        else:
            if broad_integration.get("schema_version") != \
                    "broad-enrichment-s5-to-s4-source-integration-v1" or \
                    broad_integration.get("status") != \
                    "READY_FOR_BYTE_IDENTICAL_RENUMBER_AND_REBUILD":
                failures.append("broad-enrichment source integration contract changed")

    if not broad_promotion_path.is_file():
        failures.append("broad-enrichment S5-to-S4 promotion manifest is missing")
    else:
        try:
            broad_promotion = json.loads(broad_promotion_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"broad-enrichment promotion manifest is unreadable: {exc}")
        else:
            broad_pending = "PROMOTED_SOURCE_POSTBUILD_PENDING"
            broad_final = "PROMOTED_AND_POSTBUILD_GUARDED"
            if broad_promotion.get("schema_version") != \
                    "broad-enrichment-s5-to-s4-promotion-v1":
                failures.append("broad-enrichment promotion schema changed")
            if broad_promotion.get("status") not in {broad_pending, broad_final}:
                failures.append("broad-enrichment promotion status is invalid")
            if check_derived and broad_promotion.get("status") != broad_final:
                failures.append(
                    "broad-enrichment S4 promotion is not final in --post-build mode"
                )
            canonical = broad_promotion.get("canonical", {})
            for key, rel in (("pdf", "figures/FigS4.pdf"),
                             ("png", "figures/FigS4.png")):
                record = canonical.get(key, {})
                if record.get("path") != rel or \
                        record.get("sha256") != broad_expected[rel]:
                    failures.append(
                        f"broad-enrichment promotion canonical {key} record changed"
                    )
            if canonical.get("promoted_by_byte_copy") is not True or \
                    canonical.get("former_s5_and_new_s4_bytes_identical") is not True:
                failures.append("broad-enrichment promotion lost byte-copy provenance")
            if check_derived:
                integration = broad_promotion.get("build_integration", {})
                package_record = integration.get("outputs", {}).get(
                    "supplementary_package", {}
                )
                if integration.get("status") != "COMPLETE" or \
                        package_record.get("files") != 37:
                    failures.append(
                        "broad-enrichment S4 final build/package seal is incomplete"
                    )

    for rel in ACTIVE_TEXT:
        path = ROOT / rel
        if not path.exists():
            failures.append(f"{rel}: missing active source")
            continue
        text = path.read_text()
        for pattern, reason in FORBIDDEN.items():
            for match in re.finditer(pattern, text, flags=re.I):
                line = text.count("\n", 0, match.start()) + 1
                failures.append(f"{rel}:{line}: {reason}: {match.group(0)!r}")

    # Only S1-S4 are submission figures. Historical captions live after an explicit repository
    # boundary in FIGURE_LEGENDS.md and are intentionally excluded from this check.
    active_supplement_sources = [
        "README.md", "MANUSCRIPT.md", "PAPER_FACTS.md", "REPLICATION_SCOPE.md",
        "COVER_LETTER_MSB.md", "SUPPLEMENTARY_NOTE_negatives.md",
        "submission/COVER_LETTER.md", "submission/SIGNIFICANCE_STATEMENT.md",
    ]
    retired_ref = re.compile(r"Supplementary Figure(?:s)?\s+S[6-9]", re.I)
    for rel in active_supplement_sources:
        text = (ROOT / rel).read_text()
        for match in retired_ref.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            failures.append(
                f"{rel}:{line}: active source cites retired supplementary figure: "
                f"{match.group(0)!r}")
    formal_legend_text = (ROOT / "FIGURE_LEGENDS.md").read_text().split(
        "# Repository-only figure legends", 1)[0]
    for match in retired_ref.finditer(formal_legend_text):
        line = formal_legend_text.count("\n", 0, match.start()) + 1
        failures.append(
            f"FIGURE_LEGENDS.md:{line}: formal legend cites retired supplementary figure: "
            f"{match.group(0)!r}")

    # REFRAME.md is the governance record and intentionally documents withdrawn analyses.
    # It must use the corrected census, but phrases such as "former 9.6-fold" are evidence of
    # retirement, not active claims.
    reframe = (ROOT / "REFRAME.md").read_text()
    for stale in ("83,490", "79,111", "554 self", "554 on-target"):
        if stale in reframe:
            line = reframe.count("\n", 0, reframe.index(stale)) + 1
            failures.append(f"REFRAME.md:{line}: stale canonical count {stale!r}")

    # Target-trajectory provenance is easy to corrupt because three related universes are used:
    # 910 eligible targets, 903 externally mapped specificity targets, and a 633-target display
    # universe of which 629 map. The 631 count belongs only to the eligible target-conditioned
    # census; four of those targets are unmapped, leaving 627 unanimous mapped specificity rows.
    target_theory = json.loads(
        (ROOT / "data/census/target_theory_decisive.json").read_text())
    target_min2 = target_theory.get(
        "D1_conditional_on_target_sign_composition", {}).get("min_rev_2", {})
    for key, expected in {"n_targets": 910, "n_unanimous_observed": 631}.items():
        if target_min2.get(key) != expected:
            failures.append(
                f"target_theory_decisive.json min_rev_2.{key}={target_min2.get(key)!r}; "
                f"expected {expected}")

    specificity = json.loads((ROOT / "data/target_gating/"
                              "gse140244_reversal_specificity_results.json").read_text())
    if specificity.get("status") != "POST_PRIMARY_REVERSAL_SPECIFICITY_SENSITIVITY_COMPLETE":
        failures.append("target specificity status does not preserve its post-primary designation")
    provenance = specificity.get("design_provenance", "").lower()
    if "execution was required after independent review of the unanimous-label auc" not in provenance:
        failures.append("target specificity provenance no longer records post-primary execution")
    specificity_universe = specificity.get("universe", {})
    for key, expected in {
            "identity_safe_targets_at_least_2_reversals_and_1_stable": 910,
            "mapped_targets": 903,
            "unanimous_mapped_targets": 627,
            "non_unanimous_mapped_targets": 276}.items():
        if specificity_universe.get(key) != expected:
            failures.append(
                f"target specificity universe.{key}={specificity_universe.get(key)!r}; "
                f"expected {expected}")
    primary_specificity = specificity.get(
        "primary_target_conditioned", {}).get("result", {})
    if primary_specificity.get("n_targets") != 903:
        failures.append(
            "target specificity primary result does not contain all 903 mapped targets")

    gating = json.loads((ROOT / "data/target_gating/"
                        "gse140244_target_gating_results.json").read_text())
    if gating.get("status") != "GSE140244_FROZEN_TARGET_GATING_COMPLETE":
        failures.append("frozen target-gating result is not complete")
    if gating.get("processing", {}).get("donor_time_pseudobulks") != 173:
        failures.append("target-gating result does not lock the 173 available donor-time pairs")
    gating_mapping = gating.get("mapping", {})
    for key, expected in {"frozen_targets": 633, "mapped_targets": 629}.items():
        if gating_mapping.get(key) != expected:
            failures.append(
                f"target-gating mapping.{key}={gating_mapping.get(key)!r}; expected {expected}")
    gating_min2 = gating.get("primary_and_sensitivities", {}).get("minimum_2", {})
    if gating_min2.get("n_mapped_targets") != 629:
        failures.append("target-gating minimum_2 does not contain all 629 mapped display targets")
    for estimand in ("donor_auc", "donor_stratified_auc", "median_delta_difference"):
        if gating_min2.get(estimand, {}).get("n") != 22:
            failures.append(f"target-gating minimum_2.{estimand}.n is not 22 donors")

    if 631 - (633 - 629) != specificity_universe.get("unanimous_mapped_targets"):
        failures.append("631/633/629 provenance does not reconcile to 627 mapped unanimous targets")
    if specificity_universe.get("mapped_targets") != (
            specificity_universe.get("unanimous_mapped_targets", -1) +
            specificity_universe.get("non_unanimous_mapped_targets", -1)):
        failures.append("903-target specificity universe does not partition into 627 + 276")

    row_contract = {
        "data/target_gating/gse140244_target_gating_universe.tsv": 633,
        "data/target_gating/gse140244_target_activation_results.tsv": 629,
        "data/target_gating/gse140244_reversal_specificity_targets.tsv": 903,
    }
    for rel, expected in row_contract.items():
        path = ROOT / rel
        if not path.is_file():
            failures.append(f"{rel}: missing canonical target-trajectory table")
            continue
        with path.open(newline="") as handle:
            observed = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
        if observed != expected:
            failures.append(f"{rel}: {observed} data rows; expected {expected}")

    timecourse_path = ROOT / "data/target_gating/gse140244_timecourse_donor_auc.tsv"
    if not timecourse_path.is_file():
        failures.append("target-trajectory donor time-course table is missing")
    else:
        with timecourse_path.open(newline="") as handle:
            time_rows = list(csv.DictReader(handle, delimiter="\t"))
        time_counts = Counter(int(row["time"]) for row in time_rows)
        time_donor_pairs = {(int(row["time"]), row["donor"]) for row in time_rows}
        if len(time_donor_pairs) != len(time_rows):
            failures.append("target-trajectory time course contains duplicate time-donor rows")
        expected_time_counts = {0: 22, 2: 22, 4: 22, 8: 22, 12: 19,
                                24: 22, 48: 22, 72: 22}
        if dict(sorted(time_counts.items())) != expected_time_counts:
            failures.append(
                f"target-trajectory donor counts={dict(sorted(time_counts.items()))!r}; "
                f"expected {expected_time_counts!r}")

    for rel, phrases in TARGET_TRAJECTORY_TEXT.items():
        path = ROOT / rel
        if not path.is_file():
            failures.append(f"{rel}: missing target-trajectory claim source")
            continue
        source_text = path.read_text().lower()
        for phrase in phrases:
            if not re.search(phrase, source_text):
                failures.append(
                    f"{rel}: target specificity must remain explicitly post-primary "
                    f"(missing /{phrase}/)")

    formal_legends = (ROOT / "FIGURE_LEGENDS.md").read_text().split(
        "# Repository-only figure legends", 1)[0]
    main_numbers = [int(x) for x in re.findall(r"^## Figure (\d+)\.", formal_legends, re.M)]
    supp_numbers = [int(x) for x in re.findall(
        r"^## Supplementary Figure S(\d+)\.", formal_legends, re.M)]
    if main_numbers != [1, 2, 3, 4]:
        failures.append(
            f"formal main-figure contract is {main_numbers!r}; expected [1, 2, 3, 4]")
    if supp_numbers != [1, 2, 3, 4]:
        failures.append(
            f"formal supplementary-figure contract is {supp_numbers!r}; expected S1-S4")
    formal_legend_normalized = _normalise_text(formal_legends)
    _check_required_phrases(
        "FIGURE_LEGENDS.md formal block",
        formal_legend_normalized,
        (
            "supplementary figure s4. network-adjusted all-target kegg and go biological process rank enrichment.",
            "9/168 kegg legacy pathways and 161/2,425 go biological process terms",
            "the first 12 of 92 go representatives",
            "does not test reversal direction, pathway activity, mechanism, individual-edge validity or enrichment among the 279 targets",
        ),
        failures,
    )

    # Figure 1a uses the shared complete-census renderer with the standalone figure's exact
    # top-14 candidate pool and deterministic collision filter.  Lock the promoted artwork,
    # its identity-safe authority, and the dedicated legend paragraph without hashing the
    # entire multi-figure legends file: unrelated figure edits must not invalidate this seal.
    # 2026-09-03: the independently audited v2 arrangement enlarged panel a, narrowed panels
    # b-e and set their internal/title typography to 5/7 pt.  The scientific panel-a contract,
    # all inputs, axis limits and panels f-h remained unchanged.
    fig1_promotion_path = (
        ROOT / "data/census/fig1_compact_arrangement_v2_promotion_20260903.json")
    if not fig1_promotion_path.is_file():
        failures.append("Figure 1 compact-arrangement v2 promotion manifest is missing")
    else:
        fig1_promotion = json.loads(fig1_promotion_path.read_text())
        if fig1_promotion.get("schema_version") != \
                "fig1-compact-arrangement-v2-promotion-v1":
            failures.append("Figure 1 compact-arrangement v2 promotion schema changed")
        if fig1_promotion.get("status") != "PROMOTED_AFTER_INDEPENDENT_AUDIT_GO":
            failures.append("Figure 1 compact-arrangement v2 promotion status changed")

        def check_fig1_record(record: dict, label: str) -> None:
            rel = record.get("path")
            expected_hash = record.get("sha256")
            path = ROOT / rel if rel else None
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest()
                if path is not None and path.is_file() else None)
            if observed_hash != expected_hash:
                failures.append(
                    f"{label} {rel!r} hash {observed_hash} != {expected_hash}")
            expected_size = record.get("size_bytes")
            if (expected_size is not None and path is not None and path.is_file()
                    and path.stat().st_size != expected_size):
                failures.append(
                    f"{label} {rel!r} size {path.stat().st_size} != {expected_size}")

        check_fig1_record(
            fig1_promotion.get("candidate_manifest", {}),
            "Figure 1 candidate manifest")
        for name, record in fig1_promotion.get("canonical_outputs", {}).items():
            if isinstance(record, dict):
                check_fig1_record(record, f"canonical Figure 1 {name}")
        for name, record in fig1_promotion.get("authority", {}).items():
            check_fig1_record(record, f"Figure 1 authority {name}")
        fig1_audit = fig1_promotion.get("independent_audit", {})
        check_fig1_record(fig1_audit, "Figure 1 independent audit")
        if fig1_audit.get("verdict") != "GO":
            failures.append("Figure 1 compact-arrangement v2 lacks independent GO")
        archive = fig1_promotion.get("preceding_canonical_archive", {})
        for name in ("pdf", "png", "readme"):
            check_fig1_record(archive.get(name, {}), f"archived preceding Figure 1 {name}")
        if archive.get("recoverable") is not True or \
                archive.get("permanent_deletion_performed") is not False:
            failures.append("preceding Figure 1 archive lost its recoverable no-delete contract")

        fig1_manifest_record = fig1_promotion.get("candidate_manifest", {})
        fig1_manifest_path = ROOT / fig1_manifest_record.get("path", "")
        if not fig1_manifest_path.is_file():
            failures.append("Figure 1 compact-arrangement v2 candidate manifest is missing")
        else:
            fig1_manifest = json.loads(fig1_manifest_path.read_text())
            if fig1_manifest.get("schema_version") != \
                    "fig1-compact-arrangement-candidate-v2":
                failures.append("active Figure 1 candidate is not compact arrangement v2")
            for closure_name in ("source_closure", "input_closure"):
                for rel, expected_hash in fig1_promotion.get(closure_name, {}).items():
                    path = ROOT / rel
                    observed_hash = (
                        hashlib.sha256(path.read_bytes()).hexdigest()
                        if path.is_file() else None)
                    if observed_hash != expected_hash:
                        failures.append(
                            f"Figure 1 {closure_name} {rel!r} hash "
                            f"{observed_hash} != {expected_hash}")

        expected_fig1_labels = [
            "CCNC", "TADA2B", "ARNT", "SGF29", "AHR",
            "NCKAP1L", "ATAD5", "ARHGAP30", "ELOF1", "SENP5"]
        expected_fig1_contract = {
            "heading": "4,379 sign-reversing edges connect 324 regulators to 2,577 targets",
            "displayed_regulator_labels": expected_fig1_labels,
            "n_reg_labels_candidate_pool": 14,
            "edges": 4379,
            "regulators": 324,
            "targets": 2577,
            "lost_edges": 2298,
            "gained_edges": 2081,
            "fixed_effector_display_set_size": 25,
            "fixed_effector_targets_present": 16,
            "new_inference": False,
        }
        if fig1_promotion.get("panel_a_contract") != expected_fig1_contract:
            failures.append("Figure 1 complete-census panel-a contract changed")
        fig1_layout = fig1_promotion.get("layout_contract", {})
        if fig1_layout.get("page_size_in") != [8.27, 8.75]:
            failures.append("Figure 1 compact-arrangement v2 page size changed")
        if fig1_layout.get("upper_width_ratios") != [1.28, 0.92]:
            failures.append("Figure 1 compact-arrangement v2 width ratios changed")
        if fig1_layout.get("panel_a_size_in") != [4.18705, 3.68603]:
            failures.append("Figure 1 compact-arrangement v2 panel-a geometry changed")
        if fig1_layout.get("right_panel_internal_font_pt") != 5.0 or \
                fig1_layout.get("right_panel_title_and_tag_font_pt") != 7.0:
            failures.append("Figure 1 compact-arrangement v2 typography changed")
        if fig1_layout.get("panel_a_axis_limits_changed") is not False or \
                fig1_layout.get("panels_f_to_h_geometry_or_typography_changed") is not False:
            failures.append("Figure 1 compact-arrangement v2 exceeded its presentation scope")
        legend_contract = fig1_promotion.get("legend_contract", {})
        legend_paragraph = legend_contract.get("dedicated_panel_a_paragraph", "")
        # The sealed Figure 1 manifest owns the former paragraph verbatim.  The
        # legend-concision successor instead binds an exact new whole-file hash
        # and checks the same panel-a claims semantically in its dedicated guard.
        if not legend_concision_active and (
                not legend_paragraph or legend_paragraph not in formal_legends):
            failures.append("formal Figure 1 legend lost the dedicated panel-a disclosure")
        if legend_contract.get("whole_file_hash_is_not_a_live_lock") is not True:
            failures.append("Figure 1 legend contract reverted to a brittle whole-file lock")

    figure3_legend_match = re.search(
        r"^## Figure 3\..*?(?=^## Figure 4\.)", formal_legends, flags=re.M | re.S)
    figure3_legend = figure3_legend_match.group(0) if figure3_legend_match else ""
    if not figure3_legend:
        failures.append("formal Figure 3 legend is missing")
    else:
        for required in (
                "(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"):
            if required not in figure3_legend:
                failures.append(
                    f"formal Figure 3 legend lacks required seven-panel token {required!r}")
        panel_markers = list(re.finditer(r"\*\*\(([a-g])\)\*\*", figure3_legend))
        panel_letters = [match.group(1) for match in panel_markers]
        figure3_panel_text: dict[str, str] = {}
        if panel_letters != list("abcdefg"):
            failures.append(
                f"formal Figure 3 panel order is {panel_letters!r}; expected a--g exactly once")
        else:
            for index, match in enumerate(panel_markers):
                stop = panel_markers[index + 1].start() if index + 1 < len(panel_markers) \
                    else len(figure3_legend)
                figure3_panel_text[match.group(1)] = figure3_legend[match.start():stop]
            panel_requirements = {
                "b": (
                    "Single-orientation means that the different regulators with a reversing "
                    "effect on one target all agree on its reversal orientation",
                    "at least three reversals",
                ),
                "c": (
                    "Red and teal reversal tiles",
                    "square root of the weaker endpoint",
                    "No z-score or clipping is used",
                    "fixed 3–18 scale",
                    "fixed 0–113 scale",
                    "outside the both-endpoint census",
                    "no clustering or dendrogram",
                    "displays the 37 regulators and 77 targets",
                    "All 208 reversal cells are retained",
                    "12 zero-reversal rows, three zero-reversal columns and 16 comparable non-reversals are omitted",
                ),
                "d": (
                    "all 9,554 eligible census targets",
                    "6,977 with no observed reversal",
                    "Nine of 168 KEGG Legacy pathways",
                    "none passed global maxT correction",
                    "Ten of 2,425 GO Biological Process terms",
                    "enrichment test pools reversal direction",
                    "Terms overlap", "pathway activity",
                ),
                "e": ("CCNC", "TADA2B", "SUPT7L", "direction-concordant"),
                "f": ("AHR", "ARNT", "direction-concordant"),
                "g": ("Arce", "Mantel–Haenszel", "241", "3,256", "3.23"),
            }
            if not legend_concision_active:
                panel_requirements["c"] += (
                    "complete 49-by-80 source matrix is retained",
                    "5 pt",
                )
            for letter, required_tokens in panel_requirements.items():
                panel_text = figure3_panel_text.get(letter, "")
                for required in required_tokens:
                    if required not in panel_text:
                        failures.append(
                            f"formal Figure 3{letter} legend lacks owned token {required!r}")
        stale_main_patterns = {
            r"regulators\s+with\s+(?:at\s+least|≥)\s*50\s+comparable\s+(?:trans\s+)?edges":
                "retired at-least-50-comparable-edge panel-a definition",
            r"≥\s*50[-\s]+comparable[-\s]+edge":
                "retired >=50-comparable-edge panel-a shorthand",
            r"only\s+31\s+of\s+49\s+meet\s+panel\s+a":
                "retired 31-of-49 panel-a overlap",
        }
        manuscript_text = (ROOT / "MANUSCRIPT.md").read_text()
        manuscript_figure1_match = re.search(
            r"^### Figure 1\..*?(?=^### Figure 2\.)",
            manuscript_text, flags=re.M | re.S)
        manuscript_figure1 = (
            manuscript_figure1_match.group(0) if manuscript_figure1_match else "")
        manuscript_figure3_match = re.search(
            r"^### Figure 3\..*?(?=^### Figure 4\.)",
            manuscript_text, flags=re.M | re.S)
        manuscript_figure3 = (
            manuscript_figure3_match.group(0) if manuscript_figure3_match else "")
        if not manuscript_figure3:
            failures.append("MANUSCRIPT.md Figure 3 Results section is missing")
        high_magnitude_figure1_sentence = (
            "The high-magnitude subnetwork contained 532 of the 4,379 census reversals, "
            "each with an absolute log2 fold change above 0.5 at both endpoints, linking "
            "170 regulators to 382 targets; 318 effects were lost and 214 gained after "
            "stimulation (Figure 1a).")
        if high_magnitude_figure1_sentence not in manuscript_figure1:
            failures.append(
                "MANUSCRIPT.md Figure 1 Results lost the 532-edge high-magnitude sentence")
        if "The high-magnitude subnetwork contained 532" in manuscript_figure3:
            failures.append(
                "MANUSCRIPT.md Figure 3 Results still owns the Figure 1 high-magnitude sentence")
        active_figure3_claim_text = {
            "MANUSCRIPT.md Figure 3 Results": manuscript_figure3,
            "FIGURE_LEGENDS.md formal Figure 3": figure3_legend,
        }
        for rel, text in active_figure3_claim_text.items():
            for pattern, reason in stale_main_patterns.items():
                match = re.search(pattern, text, flags=re.I)
                if match:
                    failures.append(f"{rel}: {reason}: {match.group(0)!r}")
            retired_term = re.search(
                r"direction\s+unanimity|\bunanimous\b", text, flags=re.I)
            if retired_term:
                failures.append(
                    f"{rel}: retired Figure 3 reader-facing term: "
                    f"{retired_term.group(0)!r}")
        if not re.search(r"Figure\s+3d\b", manuscript_figure3, flags=re.I):
            failures.append("MANUSCRIPT.md does not cite Figure 3d all-target enrichment")
        if not re.search(r"Figure\s+3e\b", manuscript_figure3, flags=re.I):
            failures.append("MANUSCRIPT.md does not cite the Figure 3e cofactor-trio panel")
        if not re.search(r"Figure\s+3f\b", manuscript_figure3, flags=re.I):
            failures.append("MANUSCRIPT.md does not cite the Figure 3f AHR-ARNT panel")
        if not re.search(r"Figure\s+3g\b", manuscript_figure3, flags=re.I):
            failures.append("MANUSCRIPT.md does not cite Figure 3g targeted Arce recurrence")
        for rel, text in active_figure3_claim_text.items():
            if re.search(r"Figure\s+4e\b|filled\s+(?:circle|diamond)|navy\s+diamond", text,
                         flags=re.I):
                failures.append(
                    f"{rel}: Figure 3 still owns the downstream Figure 4e connector")

    manuscript_text = (ROOT / "MANUSCRIPT.md").read_text()
    manuscript_figure4_match = re.search(
        r"^### Figure 4\..*?(?=^## Discussion)", manuscript_text, flags=re.M | re.S)
    manuscript_figure4 = (
        manuscript_figure4_match.group(0) if manuscript_figure4_match else "")
    figure4_legend_match = re.search(
        r"^## Figure 4\..*?(?=^## Fig_complete)", formal_legends, flags=re.M | re.S)
    figure4_legend = figure4_legend_match.group(0) if figure4_legend_match else ""
    connector_requirements = {
        "MANUSCRIPT.md Figure 4 Results": (
            manuscript_figure4,
            ("Figure 3e", "Eight of the 16", "Membership did not affect"),
        ),
        "FIGURE_LEGENDS.md formal Figure 4": (
            figure4_legend,
            ("Figure 3e", "Figure 3f", "BCL2L1", "DDIT4", "GPR108", "IL2RA",
             "PSMB9", "PSME1", "PSME2", "SDCBP", "29 targets tie", "11/16",
             "1/17"),
        ),
    }
    for rel, (text, required_phrases) in connector_requirements.items():
        for required in required_phrases:
            if required not in text:
                failures.append(f"{rel} lacks downstream connector token {required!r}")
        if not re.search(r"(?:not an (?:additional )?test|rather than an overlap test)",
                         text, flags=re.I):
            failures.append(f"{rel} does not preserve the connector claim boundary")
    figure4_facts = (ROOT / "PAPER_FACTS.md").read_text()
    corrected_figure4_ci_text = {
        "MANUSCRIPT.md Figure 4 Results": (
            manuscript_figure4,
            ("confidence interval 85.5–92.2",
             "95% confidence interval 14.0–19.7"),
        ),
        "FIGURE_LEGENDS.md formal Figure 4": (
            figure4_legend,
            ("2,000 regulator-cluster bootstrap draws",
             "retain the multiplicity of every sampled regulator"),
        ),
        "PAPER_FACTS.md": (
            figure4_facts,
            ("regulator-cluster 95% CI 85.5–92.2%",
             "656/3,905 (16.8%; 95% CI 14.0–19.7%)"),
        ),
    }
    for rel, (text, required_phrases) in corrected_figure4_ci_text.items():
        for phrase in required_phrases:
            if phrase not in text:
                failures.append(f"{rel} lacks corrected Figure 4 interval text {phrase!r}")
        for obsolete in ("86.4–91.8", "14.4–18.8"):
            if obsolete in text:
                failures.append(f"{rel} retains obsolete Figure 4 interval {obsolete!r}")
    if "Labels name the three estimated edges furthest from the origin" in figure4_legend:
        failures.append("FIGURE_LEGENDS.md falsely claims that Figure 4b labels three edges")

    fig = json.loads((ROOT / "data/census/figure_data.json").read_text())["census"]
    wanted = {"n_eligibility_gate_passed": 618, "n_census_regulators": 617,
              "n_gate_passed_without_shared_edge": 1,
              "n_candidate_regulators": 617, "edges_sig_both": 83489,
              "same_sign": 79110, "inverting": 4379, "self_effects_excluded": 555,
              "n_inverting_regulators": 324}
    for key, expected in wanted.items():
        if fig.get(key) != expected:
            failures.append(f"figure_data.json census.{key}={fig.get(key)!r}; expected {expected}")
    if fig.get("trans_edge_definition") != \
            "regulator_ensg != target_ensg (master is_self == false)":
        failures.append("figure_data.json does not lock the Ensembl-ID trans-edge definition")

    manifest = json.loads(
        (ROOT / "data/census/census_master_manifest.json").read_text())
    trans = manifest.get("trans_totals", {})
    for key, expected in {"edges": 83489, "same_sign": 79110, "reversing": 4379,
                          "regulators": 617, "self_edges_excluded": 555}.items():
        if trans.get(key) != expected:
            failures.append(
                f"census_master_manifest.json trans_totals.{key}={trans.get(key)!r}; "
                f"expected {expected}")
    identity = manifest.get("identity_rule", {})
    if identity.get("self_definition") != "regulator_ensg == target_ensg":
        failures.append("census master does not define self effects by stable Ensembl identity")
    if identity.get("identity_authority") != \
            "raw atlas obs/target_contrast and var/gene_ids":
        failures.append("census master Ensembl identity is not sourced from the raw atlas")
    if identity.get("symbol_alias_self_edges") != ["QNG1->C9orf64"]:
        failures.append("census master does not record the QNG1->C9orf64 self-alias edge")
    if identity.get("frozen_self_mask_agrees") is not True:
        failures.append("raw-atlas and frozen identity maps do not agree on self membership")
    if identity.get("frozen_historical_id_mismatch_symbols") != \
            {"regulator": 15, "target": 309}:
        failures.append("census master does not record the frozen historical-ID drift")
    for key, expected in {"trans_edges_without_8h": 74,
                          "trans_reversing_edges_without_8h": 2}.items():
        if manifest.get(key) != expected:
            failures.append(
                f"census_master_manifest.json {key}={manifest.get(key)!r}; expected {expected}")
    scope = manifest.get("regulator_scope", {})
    for key, expected in {"state_specific_gate_passed": 618,
                          "represented_in_trans_census": 617,
                          "gate_passed_without_a_both_significant_edge": 1}.items():
        if scope.get(key) != expected:
            failures.append(
                f"census_master_manifest.json regulator_scope.{key}={scope.get(key)!r}; "
                f"expected {expected}")
    if scope.get("gate_passed_without_a_both_significant_edge_ids") != ["CARD6"]:
        failures.append("census master regulator scope does not identify CARD6 as the "
                        "gate-passing perturbation without a shared edge")

    # The active 8 h export uses sign states conditional on endpoint-defined reversals.
    # The older t-based MID/LATE/EARLY file remains historical provenance but must not feed the
    # submission package or manuscript claim.
    for rel, expected_hash in KINETICS_SIGN_STATE_BINDINGS.items():
        path = ROOT / rel
        if not path.is_file():
            failures.append(f"active 8 h sign-state binding is missing: {rel}")
            continue
        observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed_hash != expected_hash:
            failures.append(
                f"active 8 h sign-state binding {rel} hash {observed_hash} != {expected_hash}")

    sign_summary_path = ROOT / "data/kinetics/kinetics_8hr_sign_states.summary.json"
    if sign_summary_path.is_file():
        sign_summary = json.loads(sign_summary_path.read_text())
        if sign_summary.get("analysis_status") != "POST_HOC_EXPLORATORY":
            failures.append("8 h sign-state summary lost its post-hoc status")
        coverage = sign_summary.get("coverage", {})
        expected_coverage = {
            "n_endpoint_defined_reversals_in_8h_input": 4377,
            "n_resolved": 1279,
            "n_Rest_sign_retaining": 824,
            "n_Stim48_sign_matching": 455,
            "n_unresolved": 3098,
        }
        for key, expected in expected_coverage.items():
            if coverage.get(key) != expected:
                failures.append(f"8 h sign-state coverage.{key} changed")
        direction = sign_summary.get("direction_description", {})
        expected_direction = {
            "lost_on_activation": (529, 704, 75.1),
            "gained_on_activation": (295, 575, 51.3),
        }
        for key, (n_rest, n_resolved, percentage) in expected_direction.items():
            record = direction.get(key, {})
            if (record.get("n_Rest_sign_retaining"), record.get("n_resolved"),
                    record.get("percentage_Rest_sign_retaining")) != \
                    (n_rest, n_resolved, percentage):
                failures.append(f"8 h sign-state direction result changed for {key}")
        paired = sign_summary.get("paired_regulator_result", {})
        if (paired.get("n_regulators_with_both_directions"), paired.get("n_positive"),
                paired.get("n_negative"), paired.get("n_tied")) != (75, 40, 14, 21):
            failures.append("8 h paired-regulator counts changed")
        if abs(float(paired.get("median_difference", float("nan"))) -
               0.08928571428571429) > 1e-15 or abs(float(
                   paired.get("p_value", float("nan"))) - 0.0049970904) > 1e-12:
            failures.append("8 h paired-regulator effect or p-value changed")
        if paired.get("statistic") != 417.0 or \
                "lost_minus_gained_exact_fraction" not in paired.get(
                    "authoritative_persisted_test_input", ""):
            failures.append("8 h paired-regulator exact-tie input changed")
        calibration = sign_summary.get("stable_calibration", {})
        if (calibration.get("n_8h_sign_matching_shared_endpoint_sign"),
                calibration.get("n_resolved_stable_edges")) != (47248, 47650):
            failures.append("8 h stable-edge sign calibration changed")
        if "does not define edges" not in sign_summary.get("claim_boundary", "") or \
                "exact switch times" not in sign_summary.get("claim_boundary", ""):
            failures.append("8 h sign-state claim boundary changed")

    sign_table_path = ROOT / "data/kinetics/kinetics_8hr_sign_states.csv.gz"
    if sign_table_path.is_file():
        with gzip.open(sign_table_path, "rt", encoding="utf-8", newline="") as handle:
            header = next(csv.reader(handle), [])
        required_columns = {
            "direction", "eight_hour_resolved", "eight_hour_sign_state",
            "lfc_Rest", "lfc_Stim8hr", "lfc_Stim48hr", "padj_Stim8hr",
        }
        if not required_columns.issubset(header):
            failures.append("active 8 h sign-state table is missing required columns")
        if {"t", "kinetic_class"}.intersection(header):
            failures.append("active 8 h sign-state table retains obsolete continuous classes")

    discovery = json.loads(
        (ROOT / "data/atlas/atlas_discovery_corrected.json").read_text())
    if discovery.get("main_text_authorized") is not False:
        failures.append("atlas discovery overlap is not explicitly fail-closed")

    # Figure 3 composition is rebuilt from hash-pinned all-guide derived arrays. The canonical
    # table must be the exact file recorded in its manifest; a stale July table previously sat at
    # the same path and cannot be distinguished safely by filename alone.
    comp_manifest_path = ROOT / "data/atlas/atlas_edge_composition_manifest.json"
    comp_table_path = ROOT / "data/atlas/atlas_edge_composition.csv.gz"
    comp_final_path = ROOT / "data/atlas/atlas_composition_FINAL.json"
    if not all(path.exists() for path in (comp_manifest_path, comp_table_path, comp_final_path)):
        failures.append("Figure 3 composition canonical table, manifest, or summary is missing")
    else:
        comp = json.loads(comp_manifest_path.read_text())
        if comp.get("status") != "PASS" or comp.get("analysis_version") != \
                "fig3-composition-v2.0.0":
            failures.append("Figure 3 composition manifest is not the passed v2 analysis")
        ident = comp.get("identity", {})
        expected_identity = {
            "authority": "raw atlas obs/target_contrast and var/gene_ids",
            "raw_atlas_ids_used": True,
            "target_array_id_set_exact": True,
            "historical_ids_used": False,
        }
        for key, expected in expected_identity.items():
            if ident.get(key) != expected:
                failures.append(
                    f"Figure 3 composition identity.{key}={ident.get(key)!r}; "
                    f"expected {expected!r}")
        comp_scope = comp.get("scope", {})
        expected_scope = {
            "n_edges": 83415, "n_reversals": 4377,
            "n_reversing_regulators": 323, "n_regulators": 615,
            "n_targets": 9553, "parent_trans_edges": 83489,
            "parent_reversals": 4379, "missing_8h_edges": 74,
            "missing_8h_reversals": 2,
            "regulators_absent_from_8h": ["KAT7", "STARD10"],
        }
        for key, expected in expected_scope.items():
            if comp_scope.get(key) != expected:
                failures.append(
                    f"Figure 3 composition scope.{key}={comp_scope.get(key)!r}; "
                    f"expected {expected!r}")
        coverage = comp.get("coverage", {})
        for key, expected in {
                "n_edges_with_paired_composition": 79002,
                "n_reversals_with_paired_composition": 4102,
                "n_reversals_missing_paired_composition": 275}.items():
            if coverage.get(key) != expected:
                failures.append(
                    f"Figure 3 composition coverage.{key}={coverage.get(key)!r}; "
                    f"expected {expected}")
        mcat = comp.get("regulator_checks", {}).get("MCAT", {})
        for key, expected in {"scope_edges": 18, "scope_reversals": 10,
                              "paired_edges": 18, "paired_reversals": 10,
                              "fully_resolved": True}.items():
            if mcat.get(key) != expected:
                failures.append(
                    f"Figure 3 composition MCAT.{key}={mcat.get(key)!r}; expected {expected}")
        expected_hash = comp.get("outputs", {}).get("complete_table", {}).get("sha256")
        observed_hash = hashlib.sha256(comp_table_path.read_bytes()).hexdigest()
        if expected_hash != observed_hash:
            failures.append(
                f"Figure 3 canonical table hash {observed_hash} != manifest {expected_hash}")
        primary = comp.get("primary_reversal_magnitude_check", {})
        final = json.loads(comp_final_path.read_text())
        mapping = {
            "n_rev_edges_tested": "n_rev_edges_tested",
            "real_median_abs_composition": "real_median_abs_composition",
            "zeroinfo_median_abs_composition": "zeroinfo_median_abs_composition",
            "real_over_zeroinfo": "real_over_zeroinfo",
            "median_rev_minmag": "median_rev_minmag",
            "composition_x_below_endpoint": "composition_x_below_endpoint",
            "n_that_could_flip": "n_that_could_flip",
        }
        for final_key, primary_key in mapping.items():
            got, expected = final.get(final_key), primary.get(primary_key)
            if got != expected:
                failures.append(
                    f"atlas_composition_FINAL.json {final_key}={got!r}; "
                    f"manifest has {expected!r}")

    # Figure 3 and Figure 4 are promoted as one cross-figure presentation contract.
    # Figure 3 owns the biological source pattern; Figure 4e owns the backward
    # connector when eight of those targets recur in its membership-blind top 28.
    fig3_source_paths = [
        ROOT / "scripts/fig3_biological_bridge_candidate.py",
        ROOT / "scripts/fig3_biological_bridge_r5_ahr_arce_candidate.py",
        ROOT / "scripts/fig3_biological_bridge_r5_ahr_arce_candidate_v5.py",
        ROOT / "scripts/fig3_target_centred_v6.py",
        ROOT / "scripts/fig3_target_centred_effect_zoom_v7.py",
        ROOT / "scripts/fig3_target_centred_effect_zoom_v8.py",
        ROOT / "scripts/fig3_target_centred_biological_summary_v9.py",
        ROOT / "scripts/fig3c_go_slim_biological_summary_v2.py",
        ROOT / "scripts/fig3_target_centred_top_overview_layout_v10.py",
    ]
    for fig3_source_path in fig3_source_paths:
        if not fig3_source_path.is_file():
            failures.append(f"Figure 3 source is missing: {fig3_source_path}")
            continue
        fig3_source = fig3_source_path.read_text()
        if "fig_atlas_discovery" in fig3_source or "panel_shifters" in fig3_source:
            failures.append(
                f"Figure 3 source still imports or draws the withdrawn discovery arm: "
                f"{fig3_source_path}")
    active_fig3_renderer_path = ROOT / "scripts/fig3_target_centred_top_overview_layout_v10.py"
    if active_fig3_renderer_path.is_file():
        active_fig3_renderer = active_fig3_renderer_path.read_text()
        for forbidden_token in (
                "show_shifters", "state_shifter", "STATE_SHIFTER", "occupancy_map"):
            if forbidden_token in active_fig3_renderer:
                failures.append(
                    "active Figure 3 renderer contains a forbidden occupancy-overlay hook: "
                    f"{forbidden_token!r}")

    promotion_path = (
        ROOT / "data/census/fig3_biological_bridge/"
        "promotion_manifest_20260901_scope_link_and_functional_groups.json")
    promotion: dict = {}
    if not promotion_path.is_file():
        failures.append("joint Figure 3/Figure 4 connector promotion manifest is missing")
    else:
        promotion = json.loads(promotion_path.read_text())
        if promotion.get("schema_version") !=                 "fig3-scope-link-and-functional-groups-promotion-v1":
            failures.append("joint Figure 3/Figure 4 promotion schema changed")
        joint_promotion_status = promotion.get("status")
        if joint_promotion_status not in {
                "PROMOTED_SOURCE_POSTBUILD_PENDING", "PROMOTED_AND_POSTBUILD_GUARDED"}:
            failures.append("joint Figure 3/Figure 4 promotion status is invalid")
        if check_derived and joint_promotion_status != "PROMOTED_AND_POSTBUILD_GUARDED":
            failures.append(
                "joint Figure 3/Figure 4 promotion is not final in --post-build mode")
        joint_build_integration = promotion.get("build_integration", {})
        expected_joint_build_status = (
            "FIGURE3_V8_2_POSTBUILD_PENDING"
            if joint_promotion_status == "PROMOTED_SOURCE_POSTBUILD_PENDING" else "COMPLETE")
        expected_joint_effect_status = (
            "POSTBUILD_PENDING"
            if joint_promotion_status == "PROMOTED_SOURCE_POSTBUILD_PENDING"
            else "POSTBUILD_GUARDED")
        if joint_build_integration.get("status") != expected_joint_build_status:
            failures.append(
                "joint Figure 3/Figure 4 build-integration stage is inconsistent with status")
        if joint_build_integration.get("active_effect_zoom_v8_2") != {
                "promotion_manifest": (
                    "data/census/fig3_effect_zoom_v8_2_promotion_20260903.json"),
                "status": expected_joint_effect_status,
        }:
            failures.append("joint promotion has an inconsistent Figure 3 v8.2 build stage")
        if promotion.get("supersedes") != {
                "path": (
                    "data/census/fig3_biological_bridge/"
                    "promotion_manifest_20260828_connector_in_fig4e.json"),
                "sha256": (
                    "e9af65e39d457bc943c893bbac8616cdee0e93936b8d9c535591641f14ecda91"),
        }:
            failures.append("joint promotion does not preserve the exact preceding Figure 3 record")

        def check_bound_record(record: dict, label: str) -> None:
            rel = record.get("path")
            expected_hash = record.get("sha256")
            if not isinstance(rel, str) or not rel or Path(rel).is_absolute() or \
                    ".." in Path(rel).parts:
                failures.append(f"{label} has invalid relative path {rel!r}")
                return
            if not isinstance(expected_hash, str) or not re.fullmatch(
                    r"[0-9a-f]{64}", expected_hash):
                failures.append(f"{label} has invalid SHA-256 {expected_hash!r}")
                return
            path = ROOT / rel
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest()
                if path.is_file() else None)
            if observed_hash != expected_hash:
                failures.append(
                    f"{label} {rel!r} hash {observed_hash} != {expected_hash}")

        fig3_promotion = promotion.get("figure3", {})
        fig4_promotion = promotion.get("figure4", {})
        # The joint record is historical for Figure 3 after the v10.3 promotion.
        # Its former canonical and active-text paths are preserved in the new
        # immediate-predecessor archive and are validated by the live v10.3
        # promotion block below.  Continue to bind its immutable candidate,
        # renderer, provenance and Figure 4 records here.
        for label, record in (
                ("promoted Figure 3 renderer", fig3_promotion.get("renderer", {})),
                ("promoted Figure 3 manifest", fig3_promotion.get("candidate_manifest", {})),
                ("promoted Figure 3 bridge table",
                 fig3_promotion.get("bridge_provenance_table", {})),
                ("active Figure 3 effect-zoom promotion",
                 fig3_promotion.get("active_effect_zoom_promotion", {})),
                ("promoted Figure 4 renderer", fig4_promotion.get("renderer", {})),
                ("promoted Figure 4 manifest", fig4_promotion.get("candidate_manifest", {})),
                ("promoted Figure 4 panel-e source", fig4_promotion.get("panel_e_source", {})),
        ):
            check_bound_record(record, label)
        for figure_key, figure_record in (
                ("Figure 3", fig3_promotion), ("Figure 4", fig4_promotion)):
            for index, record in enumerate(figure_record.get("independent_audits", []), start=1):
                check_bound_record(record, f"{figure_key} independent audit {index}")
        if len(fig3_promotion.get("independent_audits", [])) != 1 or \
                fig3_promotion.get("independent_audits", [{}])[0].get("decision") != "GO":
            failures.append("joint promotion does not bind one GO audit for Figure 3 v8.2")
        expected_joint_fig3_terminology = {
            "panel_b": "single-orientation",
            "panels_e_f": "direction-concordant",
            "retired_reader_facing_terms": ["direction unanimity", "unanimous"],
        }
        if fig3_promotion.get("terminology_contract") != expected_joint_fig3_terminology:
            failures.append("joint promotion Figure 3 terminology contract changed")

        expected_fig3_dir = (
            "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2")
        if fig3_promotion.get("candidate_directory") != expected_fig3_dir:
            failures.append("joint promotion points to the wrong Figure 3 candidate")
        fig3_manifest_path = ROOT / expected_fig3_dir / "manifest.json"
        if fig3_manifest_path.is_file():
            fig3_manifest = json.loads(fig3_manifest_path.read_text())
            if fig3_manifest.get("schema_version") != \
                    "fig3-target-centred-effect-zoom-v8.2":
                failures.append("active Figure 3 candidate is not effect-zoom v8.2")
            if fig3_manifest.get("status") != "DIAGNOSTIC_ONLY_NOT_PROMOTED":
                failures.append("Figure 3 v8.2 immutable candidate status changed")
            if sorted(fig3_manifest.get("panel_map", {})) != list("abcdefg"):
                failures.append("Figure 3 v8.2 panel map is not exactly a--g")

            def check_fig3_candidate_member(name: str, record: dict, label: str) -> None:
                path = ROOT / expected_fig3_dir / name
                observed_hash = (
                    hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
                if observed_hash != record.get("sha256"):
                    failures.append(
                        f"Figure 3 v8.2 {label} {name!r} hash "
                        f"{observed_hash} != {record.get('sha256')}")
                if path.is_file() and path.stat().st_size != record.get("bytes"):
                    failures.append(
                        f"Figure 3 v8.2 {label} {name!r} size "
                        f"{path.stat().st_size} != {record.get('bytes')}")

            bundled = fig3_manifest.get("bundled_inputs_and_renderers", {})
            expected_bundled_names = {
                "input_all49.csv", "input_arce.json", "input_bridge_manifest.json",
                "input_cofactor.csv", "input_concentration.json",
                "input_dense_cell_state_matrix.tsv",
                "input_dense_displayed_edge_values.tsv", "input_dense_manifest.json",
                "input_dense_regulator_order.tsv", "input_dense_target_order.tsv",
                "input_programme_manifest.json", "input_programme_memberships.csv",
                "input_programmes.csv", "input_sensor.csv", "input_unanimity.json",
                "renderer_base_r5.py", "renderer_base_v6.py", "renderer_bridge.py",
                "renderer_fig3_target_centred_effect_zoom_v7.py",
                "renderer_fig3_target_centred_effect_zoom_v8.py",
            }
            if set(bundled) != expected_bundled_names:
                failures.append(
                    "Figure 3 v8.2 bundled input/renderer inventory changed: "
                    f"{sorted(bundled)!r}")
            for name, record in bundled.items():
                check_fig3_candidate_member(name, record, "bundled artifact")
            outputs = fig3_manifest.get("outputs", {})
            expected_output_names = {
                "Fig3_target_centred_effect_zoom_v8_2.pdf",
                "Fig3_target_centred_effect_zoom_v8_2.png",
                "source_panel_a_target_concentration.csv",
                "source_panel_b_target_unanimity.csv",
                "source_panel_c_effect_score_matrix.tsv",
                "source_panel_d_effect_score_matrix.tsv",
                "source_panel_e_cofactor_targets.csv",
                "source_panel_f_sensor_targets.csv",
                "source_panel_g_arce_recurrence.csv",
            }
            if set(outputs) != expected_output_names:
                failures.append(
                    "Figure 3 v8.2 output inventory changed: "
                    f"{sorted(outputs)!r}")
            for name, record in outputs.items():
                check_fig3_candidate_member(name, record, "output")

            for label, record in (
                    ("renderer", fig3_manifest.get("renderer", {})),
                    ("v6 base renderer", fig3_manifest.get("base_renderer", {})),
                    ("v7 implementation renderer",
                     fig3_manifest.get("implementation_base_renderer_v7", {}))):
                rel = record.get("path")
                path = ROOT / rel if rel else None
                observed = (
                    hashlib.sha256(path.read_bytes()).hexdigest()
                    if path is not None and path.is_file() else None)
                if observed != record.get("sha256"):
                    failures.append(
                        f"Figure 3 v8.2 {label} {rel!r} hash "
                        f"{observed} != {record.get('sha256')}")
            if fig3_manifest.get("base_renderer", {}).get("edited") is not False or \
                    fig3_manifest.get("implementation_base_renderer_v7", {}).get(
                        "edited") is not False:
                failures.append("Figure 3 v8.2 base-renderer custody changed")

            dense_authority = fig3_manifest.get("authoritative_dense_source", {})
            dense_rel = dense_authority.get("path")
            dense_path = ROOT / dense_rel if dense_rel else None
            dense_manifest_path = dense_path / "manifest.json" if dense_path else None
            dense_manifest_hash = (
                hashlib.sha256(dense_manifest_path.read_bytes()).hexdigest()
                if dense_manifest_path is not None and dense_manifest_path.is_file() else None)
            if dense_manifest_hash != dense_authority.get("manifest_sha256"):
                failures.append("Figure 3 v8.2 dense-source manifest hash changed")
            verified_dense_outputs = dense_authority.get("verified_output_hashes", {})
            expected_dense_output_names = {
                "cell_state_matrix.tsv", "displayed_edge_values.tsv",
                "regulator_order.tsv", "target_order.tsv",
            }
            if set(verified_dense_outputs) != expected_dense_output_names:
                failures.append(
                    "Figure 3 v8.2 dense-source output inventory changed: "
                    f"{sorted(verified_dense_outputs)!r}")
            for name, expected_hash in verified_dense_outputs.items():
                path = dense_path / name if dense_path is not None else None
                observed = (
                    hashlib.sha256(path.read_bytes()).hexdigest()
                    if path is not None and path.is_file() else None)
                if observed != expected_hash:
                    failures.append(
                        f"Figure 3 v8.2 dense-source {name!r} hash "
                        f"{observed} != {expected_hash}")

            expected_panel_c = {
                "cells": 13671, "reversal": 537, "lost": 330, "gained": 207,
                "stable": 430, "outside": 12703, "self": 1,
                "positive_rows": 43, "positive_targets": 257,
            }
            expected_panel_d = {
                "cells": 3920, "reversal": 208, "lost": 137, "gained": 71,
                "stable": 128, "outside": 3584, "self": 0,
                "positive_rows": 37, "positive_targets": 77,
            }
            semantic = fig3_manifest.get("semantic_counts", {})
            if semantic.get("panel_c") != expected_panel_c:
                failures.append("Figure 3 v8.2 panel-c semantic counts changed")
            if semantic.get("panel_d") != expected_panel_d:
                failures.append("Figure 3 v8.2 panel-d semantic counts changed")
            if semantic.get("panel_d_zero_reversal_targets") != [
                    "CD47", "NSMCE1", "SPIN1"]:
                failures.append("Figure 3 v8.2 zero-reversal zoom targets changed")
            if semantic.get("t279_orientation") != {"gained": 118, "lost": 161} or \
                    semantic.get("first80_orientation") != {"gained": 28, "lost": 52}:
                failures.append("Figure 3 v8.2 orientation composition changed")

            geometry = fig3_manifest.get("geometry", {})
            panel_c = geometry.get("panel_c", {})
            panel_d = geometry.get("panel_d", {})
            expected_panel_geometry = {
                "c": ([49, 279], 537, 330, 207, 430, 12703, 1, 43, 257),
                "d": ([49, 80], 208, 137, 71, 128, 3584, 0, 37, 77),
            }
            for letter, panel, expected in (
                    ("c", panel_c, expected_panel_geometry["c"]),
                    ("d", panel_d, expected_panel_geometry["d"])):
                observed = (
                    panel.get("shape"), panel.get("reversal_cells"),
                    panel.get("lost_cells"), panel.get("gained_cells"),
                    panel.get("stable_cells"), panel.get("outside_cells"),
                    panel.get("self_cells"), panel.get("regulators_with_reversal"),
                    panel.get("targets_with_reversal"),
                )
                if observed != expected:
                    failures.append(
                        f"Figure 3 v8.2 panel-{letter} geometry/count contract changed")
                if panel.get("dendrogram") is not False or \
                        "no clustering" not in panel.get("ordering", ""):
                    failures.append(
                        f"Figure 3 v8.2 panel-{letter} acquired clustering/dendrogram")
            if panel_c.get("row_labels_shown") is not False or \
                    panel_c.get("target_labels_shown") is not False:
                failures.append("Figure 3 v8.2 overview label policy changed")
            if panel_d.get("row_labels_shown") is not True or \
                    panel_d.get("target_labels_shown") is not True:
                failures.append("Figure 3 v8.2 zoom label policy changed")
            if panel_d.get("regulator_label_points") != 5.5 or \
                    panel_d.get("target_label_points") != 5.5:
                failures.append("Figure 3 v8.2 panel-d labels are not exactly 5.5 pt")

            literal_zoom = geometry.get("literal_zoom_assertions", {})
            expected_literal_zoom = {
                "states_equal_c_first_80": True,
                "scores_equal_c_first_80": True,
                "targets_equal_c_first_80": True,
                "row_order_equal": True,
                "orientation_track_equal_c_first_80": True,
                "target_burden_track_equal_c_first_80": True,
                "row_burden_track_identical_full_panel_scale": True,
            }
            if literal_zoom != expected_literal_zoom:
                failures.append("Figure 3 v8.2 literal c[:, :80] assertions changed")
            c_matrix_path = ROOT / expected_fig3_dir / \
                "source_panel_c_effect_score_matrix.tsv"
            d_matrix_path = ROOT / expected_fig3_dir / \
                "source_panel_d_effect_score_matrix.tsv"
            if c_matrix_path.is_file() and d_matrix_path.is_file():
                with c_matrix_path.open(newline="") as handle:
                    c_matrix = list(csv.reader(handle, delimiter="\t"))
                with d_matrix_path.open(newline="") as handle:
                    d_matrix = list(csv.reader(handle, delimiter="\t"))
                literal_prefix = (
                    len(c_matrix) == len(d_matrix) == 50
                    and all(d_row == c_row[:81]
                            for c_row, d_row in zip(c_matrix, d_matrix)))
                if not literal_prefix:
                    failures.append("Figure 3 v8.2 panel-d score table is not literal c[:, :80]")

            heat = geometry.get("heat_geometry", {})
            if heat.get("c", {}).get("cell_width_points") != 1.59 or \
                    abs(float(heat.get("c", {}).get("cell_height_points", 0)) - 1.59) > 1e-12:
                failures.append("Figure 3 v8.2 panel-c square-cell geometry changed")
            if heat.get("d", {}).get("cell_width_points") != 5.5 or \
                    heat.get("d", {}).get("cell_height_points") != 5.5:
                failures.append("Figure 3 v8.2 panel-d cells are not 5.5-pt squares")
            label_audit = geometry.get("panel_d_label_audit", {})
            if any(label_audit.get(key) != 5.5 for key in (
                    "requested_points", "cell_pitch_points",
                    "observed_min_points", "observed_max_points")):
                failures.append("Figure 3 v8.2 delivered panel-d typography changed")
            for axis_name, expected_labels in (
                    ("target_labels", 80), ("regulator_labels", 49)):
                axis_audit = label_audit.get(axis_name, {})
                if axis_audit.get("labels") != expected_labels or \
                        axis_audit.get("overlap_pairs_above_0.2px_each_axis") != 0 or \
                        axis_audit.get("maximum_overlap_area_pixels2") != 0.0:
                    failures.append(
                        f"Figure 3 v8.2 {axis_name} collision audit changed")

            expected_graphical_key = {
                "prose_legend_removed": True,
                "shared_for_panels": ["c", "d"],
                "location": "beside panel d; no duplicate overview key",
                "orientation": {"lost": "#c1272d", "gained": "#2a9d8f"},
                "target_burden": {"cmap": "Purples", "minimum": 3, "maximum": 18},
                "row_burden": {"cmap": "Blues", "minimum": 0, "maximum": 113},
                "outside": {"fill": "#ffffff", "bordered": True},
                "stable_self": {
                    "stable": "#e4e7ea", "self": "#7f8a94", "split_swatch": True},
                "effect_ticks_raw_signed_log2fc": [-2.0, -0.5, 0.5, 2.0],
            }
            if geometry.get("graphical_heatmap_key") != expected_graphical_key:
                failures.append("Figure 3 v8.2 graphical c/d key changed")
            encoding = fig3_manifest.get("encoding", {})
            if encoding.get("effect") != (
                    "sign * sqrt(min(abs(lfc_Rest), abs(lfc_Stim48)) / "
                    "2.21400360540582)") or encoding.get("z_score") is not None or \
                    encoding.get("clipping") is not False:
                failures.append("Figure 3 v8.2 effect transform/z-score/clipping changed")
            graphical_encoding = encoding.get("graphical_legend", {})
            if graphical_encoding.get("shared_for_panels") != ["c", "d"] or \
                    graphical_encoding.get("location") != \
                    "beside panel d; no duplicate overview key":
                failures.append("Figure 3 v8.2 shared graphical-key placement changed")

            terminology = fig3_manifest.get("rendered_terminology_contract", {})
            expected_terminology = {
                "panel_b": "Single-orientation targets",
                "panel_e": "Cofactor trio\n16 direction-concordant targets",
                "panel_f": "Sensor pair\n17 direction-concordant targets",
                "forbidden_rendered_terms": ["direction unanimity", "unanimous"],
                "forbidden_matches": [],
            }
            if terminology != expected_terminology:
                failures.append("Figure 3 v8.2 rendered terminology contract changed")
            title_audit = geometry.get("panel_title_terminology_guard", {})
            if title_audit.get("rendered") != {
                    "b": "Single-orientation targets",
                    "e": "Cofactor trio\n16 direction-concordant targets",
                    "f": "Sensor pair\n17 direction-concordant targets",
            } or title_audit.get("forbidden_matches") != [] or \
                    title_audit.get("e_f_title_overlap") is not False or \
                    title_audit.get("title_points") != 7.0:
                failures.append("Figure 3 v8.2 title terminology/fit audit changed")

            text_audit = geometry.get("text_audit", {})
            if (text_audit.get("min_visible_text_points") or 0) < 5.0 or \
                    (text_audit.get("max_visible_text_points") or 99) > 7.0 or \
                    text_audit.get("off_page_text") != []:
                failures.append("Figure 3 v8.2 visible text violates the 5--7 pt/page contract")
            deterministic = fig3_manifest.get("export", {}).get(
                "deterministic_rerender", {})
            if deterministic.get("fresh_canvas_count") != 2 or \
                    deterministic.get("pdf_byte_identical") is not True or \
                    deterministic.get("png_byte_identical") is not True or \
                    deterministic.get("pdf_sha256") != outputs.get(
                        "Fig3_target_centred_effect_zoom_v8_2.pdf", {}).get("sha256") or \
                    deterministic.get("png_sha256") != outputs.get(
                        "Fig3_target_centred_effect_zoom_v8_2.png", {}).get("sha256"):
                failures.append("Figure 3 v8.2 deterministic-render closure changed")
            export = fig3_manifest.get("export", {})
            if export.get("page") != "A4" or \
                    export.get("pdf_vector_heatmaps") is not True or \
                    export.get("pdf_fonttype") != 42 or export.get("png_dpi") != 300 or \
                    export.get("tight_bbox") is not False:
                failures.append("Figure 3 v8.2 export contract changed")
            if "not a cluster" not in fig3_manifest.get("claim_boundary", "") or \
                    "were not separately tested" not in fig3_manifest.get(
                        "claim_boundary", ""):
                failures.append("Figure 3 v8.2 descriptive/no-cluster claim boundary changed")
            if fig3_manifest.get("canonical_changes") is not False or \
                    fig3_manifest.get("manuscript_changes") is not False or \
                    fig3_manifest.get("legend_changes") is not False or \
                    fig3_manifest.get("prior_candidate_changes") is not False or \
                    fig3_manifest.get("trash_changes") is not False:
                failures.append("Figure 3 v8.2 immutable pre-promotion custody changed")
        else:
            failures.append("Figure 3 v8.2 candidate manifest is missing")

        expected_fig4_dir = (
            "data/figure_provenance/fig4_v5_1_fig3d_connector_candidate_v1")
        if fig4_promotion.get("candidate_directory") != expected_fig4_dir:
            failures.append("joint promotion points to the wrong Figure 4 candidate")
        fig4_manifest_path = ROOT / expected_fig4_dir / "manifest.json"
        expected_marked = [
            "BCL2L1", "DDIT4", "GPR108", "IL2RA",
            "PSMB9", "PSME1", "PSME2", "SDCBP"]
        if fig4_manifest_path.is_file():
            fig4_manifest = json.loads(fig4_manifest_path.read_text())
            if fig4_manifest.get("schema") !=                     "response-shape-figure4-fig3d-connector-candidate-v1":
                failures.append("active Figure 4 connector candidate schema changed")
            if fig4_manifest.get("status") != "COMPLETE":
                failures.append("active Figure 4 connector candidate is not complete")
            if fig4_manifest.get("marked_target_count") != 8 or                     fig4_manifest.get("marked_targets") != expected_marked:
                failures.append("Figure 4 connector candidate changed the exact eight targets")
            if fig4_manifest.get("ahr_arnt_top28_overlap") != 0:
                failures.append("Figure 4 connector candidate changed AHR-ARNT overlap")
            if fig4_manifest.get("new_inference") is not False or                     fig4_manifest.get("resampling_performed") is not False:
                failures.append("Figure 4 connector candidate reran or introduced inference")
            if fig4_manifest.get("top28_rule") != (
                    "estimated reversal edges; target size >=3; nlargest(28, size); "
                    "sort_values(size)"):
                failures.append("Figure 4 connector candidate changed the top-28 rule")
            expected_outputs = {
                "pdf": "b514e8ee88dac9321bfc4d54f9e97c0c8df6e94a05b7094ae764089cb0a7a37e",
                "png": "67eebe559facdf8309744212ac0e384f6cf588a8b414586923e63a13b116e08e",
            }
            for key, expected_hash in expected_outputs.items():
                if fig4_manifest.get("outputs", {}).get(key, {}).get("sha256") != expected_hash:
                    failures.append(f"Figure 4 connector manifest changed {key} output")
            bridge_records = fig4_manifest.get("bridge_table_records", [])
            if len(bridge_records) != 28:
                failures.append(
                    f"Figure 4 connector table has {len(bridge_records)} rows; expected 28")
            marked_records = [
                row["target"] for row in bridge_records if row.get("fig3d_trio_shared") is True]
            if sorted(marked_records) != expected_marked:
                failures.append(
                    f"Figure 4 connector table marks {sorted(marked_records)!r}; "
                    f"expected {expected_marked!r}")
            for row in bridge_records:
                if row.get("supported_reversals") != (
                        row.get("detection_dominant", 0) + row.get("intensity_dominant", 0)):
                    failures.append(
                        f"Figure 4 connector row {row.get('target')!r} is not additive")
        else:
            failures.append("Figure 4 connector candidate manifest is missing")

        expected_connector_contract = {
            "encoding_owner": "Figure 4e",
            "encoding": "navy diamond plus bold target label and direct in-panel key",
            "cofactor_display_overlap_n": 8,
            "cofactor_display_overlap_targets": expected_marked,
            "sensor_display_overlap_n": 0,
            "display_cutoff_supported_reversals": 4,
            "targets_above_cutoff": 20,
            "targets_tied_at_cutoff": 29,
            "display_slots_at_cutoff": 8,
            "cofactor_inclusive_at_cutoff_n": 11,
            "sensor_inclusive_at_cutoff_n": 1,
            "all_marked_targets_above_cutoff": True,
            "top28_selection_uses_connector_membership": False,
            "top28_membership_order_and_bar_values_unchanged": True,
            "new_statistics_rerun": False,
            "new_inference": False,
        }
        if promotion.get("connector_contract") != expected_connector_contract:
            failures.append("joint promotion connector contract changed")

        fig4_presentation_path = ROOT / "figures/Fig4.presentation.json"
        if fig4_presentation_path.is_file():
            fig4_presentation = json.loads(fig4_presentation_path.read_text())
            presentation_connector = fig4_presentation.get("connector", {})
            presentation_expected = {
                "role": (
                    "backward descriptive display cross-reference from Figure 4e to the "
                    "exhaustive Figure 3d cofactor-trio intersection"),
                "encoding": "navy diamond plus bold target label and direct in-panel key",
                "marked_target_count": 8,
                "marked_targets": expected_marked,
                "ahr_arnt_displayed_overlap": 0,
                "top28_selection_uses_connector_membership": False,
                "top28_membership_order_and_bar_values_unchanged": True,
                "display_cutoff_supported_reversals": 4,
                "targets_above_cutoff": 20,
                "targets_tied_at_cutoff": 29,
                "display_slots_at_cutoff": 8,
                "cofactor_inclusive_at_cutoff": 11,
                "sensor_inclusive_at_cutoff": 1,
                "interpretation": (
                    "descriptive rendered-display cross-reference; not an overlap test, "
                    "enrichment test, target test or mechanism claim"),
            }
            if fig4_presentation.get("schema") == \
                    "response-shape-figure4-horizontal-transposed-canonical-v3":
                presentation_expected["role"] = (
                    "backward descriptive display cross-reference from Figure 4e to the "
                    "Figure 3e cofactor-trio intersection")
                presentation_expected["visible_key"] = "◆ also in Fig. 3e trio"
            if presentation_connector != presentation_expected:
                failures.append("canonical Figure 4 presentation connector record changed")
            if fig4_presentation.get("schema") not in {
                    "response-shape-figure4-presentation-v5.1.1",
                    "response-shape-figure4-presentation-v5.2-landmark-census",
                    "response-shape-figure4-horizontal-transposed-canonical-v2",
                    "response-shape-figure4-horizontal-transposed-canonical-v3",
            }:
                failures.append("canonical Figure 4 presentation schema is unsupported")
            expected_bootstrap = {
                "seed": 20260824,
                "draws": 2000,
                "cluster_unit": "regulator_ensg",
                "sampling": "regulator clusters sampled with replacement",
                "pooling": (
                    "all edges in each sampled cluster, retaining sample multiplicity"),
                "interval": "2.5th and 97.5th percentiles",
            }
            if fig4_presentation.get("alignment_bootstrap") != expected_bootstrap:
                failures.append("canonical Figure 4 bootstrap method or frozen seed changed")
            expected_alignment = {
                "matched stable": {
                    "ci_low": 0.14013058032670647,
                    "ci_high": 0.1967170500869265,
                    "pooled": 656 / 3905,
                    "edges": 3905,
                    "regulators": 104,
                },
                "reversals": {
                    "ci_low": 0.8547043806483036,
                    "ci_high": 0.9220900905896241,
                    "pooled": 1694 / 1899,
                    "edges": 1899,
                    "regulators": 101,
                },
            }
            observed_alignment = fig4_presentation.get("alignment_null", {})
            for group, expected in expected_alignment.items():
                observed = observed_alignment.get(group, {})
                for key, value in expected.items():
                    if isinstance(value, float):
                        if abs(float(observed.get(key, float("nan"))) - value) > 1e-15:
                            failures.append(
                                f"canonical Figure 4 {group}.{key} changed from {value!r}")
                    elif observed.get(key) != value:
                        failures.append(
                            f"canonical Figure 4 {group}.{key}={observed.get(key)!r}; "
                            f"expected {value!r}")
        else:
            failures.append("canonical Figure 4 presentation record is missing")

        correction = fig4_promotion.get("cluster_bootstrap_correction_2026_09_01", {})
        if correction.get("status") != "CORRECTED_RERENDERED_AND_LOCALLY_GUARDED":
            failures.append("Figure 4 cluster-bootstrap correction is not sealed")
        if (correction.get("seed"), correction.get("draws")) != (20260824, 2000):
            failures.append("Figure 4 cluster-bootstrap correction changed seed or draw count")
        if correction.get("reversal_exact_ci") != [
                0.8547043806483036, 0.9220900905896241] or \
                correction.get("matched_stable_exact_ci") != [
                    0.14013058032670647, 0.1967170500869265]:
            failures.append("Figure 4 correction manifest does not preserve the exact intervals")
        for label, record in (
                ("Figure 4 bootstrap regression test", correction.get("regression_test", {})),
                ("Figure 4 bootstrap focused-test log", correction.get("focused_test_log", {})),
                ("Figure 4 bootstrap render log", correction.get("canonical_render_log", {})),
        ):
            check_bound_record(record, label)
        correction_archive = correction.get("preceding_canonical_archive", {})
        correction_archive_path = ROOT / correction_archive.get("path", "")
        for name, key in (
                ("Fig4.pdf", "pdf_sha256"),
                ("Fig4.png", "png_sha256"),
                ("Fig4.presentation.json", "presentation_record_sha256"),
                ("fig4_presentation_v5_1.py", "renderer_sha256"),
        ):
            path = correction_archive_path / name
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != correction_archive.get(key):
                failures.append(
                    f"preceding Figure 4 correction archive {name} hash "
                    f"{observed_hash} != {correction_archive.get(key)}")
        if correction_archive.get("recoverable") is not True or \
                correction_archive.get("permanent_deletion_performed") is not False:
            failures.append("preceding Figure 4 correction archive is not recoverable")

        # The Figure 3 bridge table independently locks both the rendered and inclusive counts.
        # It belongs to the retained v6 predecessor; v8.2 changes c/d to an effect map but
        # preserves the e-trio -> Figure 4e connector without rewriting this historical table.
        historical_fig3_bridge_dir = (
            ROOT / "data/figure_provenance/fig3_target_centred_v6")
        bridge_table_path = (
            historical_fig3_bridge_dir / "source_fig3d_e_to_fig4e_target_bridge.csv")
        if bridge_table_path.is_file():
            with bridge_table_path.open(newline="") as handle:
                bridge_rows = list(csv.DictReader(handle))
            d_rows = [row for row in bridge_rows if row["panel"] == "d"]
            e_rows = [row for row in bridge_rows if row["panel"] == "e"]
            d_displayed = sorted(
                row["target"] for row in d_rows if row["in_fig4e_top28"] == "True")
            e_displayed = sorted(
                row["target"] for row in e_rows if row["in_fig4e_top28"] == "True")
            d_inclusive = sum(
                bool(row["supported_reversals"])
                and float(row["supported_reversals"]) >= 4 for row in d_rows)
            e_inclusive = sum(
                bool(row["supported_reversals"])
                and float(row["supported_reversals"]) >= 4 for row in e_rows)
            if d_displayed != expected_marked or e_displayed != []:
                failures.append(
                    f"bridge table displayed overlaps d={d_displayed!r}, e={e_displayed!r}")
            if (d_inclusive, e_inclusive) != (11, 1):
                failures.append(
                    f"bridge table inclusive overlaps {(d_inclusive, e_inclusive)!r}; "
                    "expected (11, 1)")

        # These hashes preserve the exact text snapshot at the time of connector promotion.
        # They are not live whole-file locks: MANUSCRIPT.md, FIGURE_LEGENDS.md and PAPER_FACTS.md
        # also contain unrelated figures and must be allowed to evolve.  The extracted active
        # Figure 3/Figure 4 sections above enforce the connector's current semantic contract.
        for label, record in promotion.get("active_text", {}).items():
            rel = record.get("path")
            if not rel or not (ROOT / rel).is_file():
                failures.append(
                    f"active connector text {label} source is missing: {rel!r}")

        # Retained assets and archives are part of the no-delete contract.
        preserved = promotion.get("preserved_related_assets", {})
        for rel, key in (
                ("figures/repository/uncited_retained/FigS7.pdf", "repository_figure_r5_pdf_sha256"),
                ("figures/repository/uncited_retained/FigS7.png", "repository_figure_r5_png_sha256"),
                ("figures/repository/uncited_retained/Fig3_target_axis_candidate.pdf",
                 "target_axis_candidate_pdf_sha256"),
                ("figures/repository/uncited_retained/Fig3_target_axis_candidate.png",
                 "target_axis_candidate_png_sha256"),
        ):
            path = ROOT / rel
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != preserved.get(key):
                failures.append(
                    f"preserved related asset {rel!r} hash {observed_hash} "
                    f"!= {preserved.get(key)}")
        for archive_key in ("previous_canonical_archive",):
            for figure_name, figure_record in (
                    ("Figure 3", fig3_promotion), ("Figure 4", fig4_promotion)):
                archive = figure_record.get(archive_key, {})
                archive_path = ROOT / archive.get("path", "")
                if not archive_path.is_dir():
                    failures.append(
                        f"{figure_name} recoverable archive is missing: {archive_path}")
        old_fig3_archive = (
            ROOT / fig3_promotion.get("previous_canonical_archive", {}).get("path", ""))
        for name, key in (
                ("fig3_biological_bridge_r5_ahr_arce_candidate_v5_downstream_connector/manifest.json", "manifest_sha256"),
                ("Fig3.pdf", "pdf_sha256"),
                ("Fig3.png", "png_sha256"),
                ("promotion_manifest_20260828_connector_in_fig4e.json", "promotion_manifest_sha256"),
                ("README.md", "readme_sha256"),
        ):
            expected_hash = fig3_promotion.get(
                "previous_canonical_archive", {}).get(key)
            path = old_fig3_archive / name
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != expected_hash:
                failures.append(
                    f"archived displaced Figure 3 {name} hash {observed_hash} != {expected_hash}")
        old_fig4_archive = (
            ROOT / fig4_promotion.get("previous_canonical_archive", {}).get("path", ""))
        for name, key in (
                ("Fig4.pdf", "pdf_sha256"),
                ("Fig4.png", "png_sha256"),
                ("Fig4.presentation.json", "presentation_record_sha256"),
                ("fig4_presentation_v5_1.py", "renderer_sha256"),
        ):
            expected_hash = fig4_promotion.get(
                "previous_canonical_archive", {}).get(key)
            path = old_fig4_archive / name
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != expected_hash:
                failures.append(
                    f"archived Figure 4 v5.1 {name} hash {observed_hash} != {expected_hash}")

        diagnostics_root = ROOT / "data/figure_provenance"
        observed_fig3_directories = sorted(
            path.name for path in diagnostics_root.iterdir()
            if path.is_dir()
            and path.name.lstrip(".").startswith(("fig3_", "fig_ge50_")))
        # Retain the promoted v8.2 candidate, its v6 predecessor and the source-closed
        # v11 candidate. Superseded drafts live in recoverable quarantine and must not
        # silently re-enter the active diagnostic surface.
        required_fig3_directories = {
            "fig3_target_centred_alltarget_enrichment_v11",
            "fig3_target_centred_compact_zoom_v12",
            "fig3_target_centred_compact_zoom_v12_promotion",
            "fig3_target_centred_effect_zoom_v8_2",
            "fig3_target_centred_unified_zoom_v13",
            "fig3_target_centred_v6",
        }
        allowed_fig3_directories = required_fig3_directories | {
            "fig3_target_centred_unified_zoom_v13_promotion",
        }
        observed_fig3_set = set(observed_fig3_directories)
        if not required_fig3_directories.issubset(observed_fig3_set) or \
                not observed_fig3_set.issubset(allowed_fig3_directories):
            failures.append(
                "Figure 3 diagnostics cleanup mismatch; "
                f"observed={observed_fig3_directories!r}, "
                f"required={sorted(required_fig3_directories)!r}, "
                f"allowed={sorted(allowed_fig3_directories)!r}")

        bridge_manifest_path = (
            ROOT / "data/census/fig3_biological_bridge/manifest.json")
        if bridge_manifest_path.is_file():
            bridge_manifest = json.loads(bridge_manifest_path.read_text())
            if "nineteen KEGG terms" not in bridge_manifest.get(
                    "programme_annotation_policy", ""):
                failures.append(
                    "Figure 3 biological bridge manifest lost the nineteen-term dictionary")
            table_expectations = {
                "shared_reversal_targets_CCNC_TADA2B_SUPT7L.csv": (16, 6),
                "shared_reversal_targets_AHR_ARNT.csv": (17, 4),
            }
            for name, (expected_rows, expected_classified) in table_expectations.items():
                path = bridge_manifest_path.parent / name
                with path.open(newline="") as handle:
                    rows = list(csv.DictReader(handle))
                classified = sum(
                    row.get("program_annotation") != "unclassified" for row in rows)
                if len(rows) != expected_rows or classified != expected_classified:
                    failures.append(
                        f"Figure 3 {name} has {len(rows)} rows/{classified} classified; "
                        f"expected {expected_rows}/{expected_classified}")

    # Figure 3 v8.2 has its own promotion record because it superseded only Figure 3
    # inside the older joint Figure 3/Figure 4 presentation contract.  Bind both the
    # immutable diagnostic bundle and the recoverable promotion/cleanup custody here.
    effect_promotion_rel = (
        "data/census/fig3_effect_zoom_v8_2_promotion_20260903.json")
    effect_promotion_path = ROOT / effect_promotion_rel
    if not effect_promotion_path.is_file():
        failures.append("Figure 3 v8.2 standalone promotion record is missing")
    else:
        effect_promotion = json.loads(effect_promotion_path.read_text())
        expected_effect_schema = "fig3-effect-zoom-v8.2-promotion-v1"
        if effect_promotion.get("schema_version") != expected_effect_schema:
            failures.append("Figure 3 v8.2 standalone promotion schema changed")
        effect_status = effect_promotion.get("status")
        effect_pending_status = "PROMOTED_SOURCE_POSTBUILD_PENDING"
        effect_final_status = (
            "PROMOTED_AFTER_INDEPENDENT_AUDIT_GO_AND_POSTBUILD_GUARDED")
        if effect_status not in {effect_pending_status, effect_final_status}:
            failures.append("Figure 3 v8.2 standalone promotion status is invalid")
        if check_derived and effect_status != effect_final_status:
            failures.append(
                "Figure 3 v8.2 standalone promotion is not final in --post-build mode")

        def check_effect_record(
                record: dict, label: str, expected_path: str | None = None) -> str | None:
            if not isinstance(record, dict):
                failures.append(f"Figure 3 v8.2 {label} record is not an object")
                return None
            rel = record.get("path")
            if expected_path is not None and rel != expected_path:
                failures.append(
                    f"Figure 3 v8.2 {label} path {rel!r} != {expected_path!r}")
            path = ROOT / rel if isinstance(rel, str) else None
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest()
                if path is not None and path.is_file() else None)
            if observed_hash != record.get("sha256"):
                failures.append(
                    f"Figure 3 v8.2 {label} {rel!r} hash "
                    f"{observed_hash} != {record.get('sha256')}")
            return observed_hash

        expected_supersedes = {
            "path": "data/census/fig3_target_centred_promotion_20260903.json",
            "sha256": "36fe72e1aaf03d1a78a8c8a2cc41f54d09e60696d543064b0193ee0728a2be11",
        }
        if effect_promotion.get("supersedes") != expected_supersedes:
            failures.append("Figure 3 v8.2 promotion lost its exact v6 predecessor record")
        check_effect_record(
            effect_promotion.get("supersedes", {}), "superseded v6 promotion",
            expected_supersedes["path"])

        expected_effect_candidate = {
            "directory": "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2",
            "manifest": {
                "path": (
                    "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2/"
                    "manifest.json"),
                "sha256": (
                    "3f602f7193f113527027fb4576047a9c94fad343eace0b9bff165426194a2bc6"),
            },
            "pdf": {
                "path": (
                    "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2/"
                    "Fig3_target_centred_effect_zoom_v8_2.pdf"),
                "sha256": (
                    "48cbf1202c51604d94d6fa7ffb5420d6d6699edd5f9d3cdd2a07a708ff2afa91"),
            },
            "png": {
                "path": (
                    "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2/"
                    "Fig3_target_centred_effect_zoom_v8_2.png"),
                "sha256": (
                    "df9622531b3c74dec29bda73425096ee071bb3a2ef2db27893511a2279321618"),
            },
            "renderer": {
                "path": "scripts/fig3_target_centred_effect_zoom_v8.py",
                "sha256": (
                    "d709c21e15ae0d3e450cca1c8ca4e1f46e50651b8815021a72b20c1c61486937"),
            },
            "independent_audit": {
                "path": (
                    "data/figure_provenance/fig3_target_centred_effect_zoom_v8_2/"
                    "independent_audit.log"),
                "sha256": (
                    "ffafdec2bf77e268adc158f80317f1e1145f10229165055f5933e2971eb4c254"),
                "decision": "GO",
            },
            "independent_rerender_guard": {
                "path": (
                    "data/figure_provenance/"
                    "fig3_target_centred_effect_zoom_v8_2.independent_rerender.guard.log"),
                "sha256": (
                    "7420a91add23659e1f274938bd9d6d126f9c94d52ef40f6e16356a6d82f5b557"),
                "pdf_byte_identical": True,
                "png_byte_identical": True,
            },
        }
        effect_candidate = effect_promotion.get("candidate", {})
        if effect_candidate != expected_effect_candidate:
            failures.append("Figure 3 v8.2 standalone candidate custody changed")
        for key, expected_path in (
                ("manifest", expected_effect_candidate["manifest"]["path"]),
                ("pdf", expected_effect_candidate["pdf"]["path"]),
                ("png", expected_effect_candidate["png"]["path"]),
                ("renderer", expected_effect_candidate["renderer"]["path"]),
                ("independent_audit",
                 expected_effect_candidate["independent_audit"]["path"]),
                ("independent_rerender_guard",
                 expected_effect_candidate["independent_rerender_guard"]["path"]),
        ):
            check_effect_record(
                effect_candidate.get(key, {}), f"candidate {key}", expected_path)

        expected_effect_canonical = {
            "pdf": {
                "path": "figures/Fig3.pdf",
                "sha256": (
                    "48cbf1202c51604d94d6fa7ffb5420d6d6699edd5f9d3cdd2a07a708ff2afa91"),
            },
            "png": {
                "path": "figures/Fig3.png",
                "sha256": (
                    "df9622531b3c74dec29bda73425096ee071bb3a2ef2db27893511a2279321618"),
            },
            "promoted_by_byte_copy": True,
            "candidate_and_canonical_bytes_identical": True,
        }
        effect_canonical = effect_promotion.get("canonical", {})
        if effect_canonical != expected_effect_canonical:
            failures.append("Figure 3 v8.2 canonical byte-copy contract changed")
        for key in ("pdf", "png"):
            if effect_canonical.get(key, {}).get("sha256") != \
                    effect_candidate.get(key, {}).get("sha256"):
                failures.append(
                    f"Figure 3 v8.2 canonical {key} is not byte-identical to candidate")

        expected_immediate_archive = {
            "path": (
                "figures/repository/superseded_figure3_candidates_2026-09-03/"
                "pre_effect_zoom_v8_20260903"),
            "readme_sha256": (
                "1dfcd528dcf1584b855d01e3df89940b552191cd623df136f129166a2c0a466a"),
            "sha256sums_sha256": (
                "d0903374c0d5a1f0775765b505654510b72f664b56ffcfb7dfe6ea6054afce5c"),
            "old_canonical_pdf_sha256": (
                "421b95137fed0788b562fcd5f738f4d15a1ea6ecc57508eca961ea03e6ac844a"),
            "old_canonical_png_sha256": (
                "1f321049a45d417de021ec5fafc6117d3d455e2bb03be929d3c498b26b697427"),
            "recoverable": True,
        }
        immediate_archive = effect_promotion.get("immediate_predecessor_archive", {})
        if immediate_archive != expected_immediate_archive:
            failures.append("Figure 3 v8.2 immediate-predecessor archive contract changed")
        immediate_archive_path = ROOT / expected_immediate_archive["path"]
        for name, hash_key in (
                ("README.md", "readme_sha256"),
                ("SHA256SUMS.txt", "sha256sums_sha256"),
                ("Fig3.pdf", "old_canonical_pdf_sha256"),
                ("Fig3.png", "old_canonical_png_sha256"),
        ):
            path = immediate_archive_path / name
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != expected_immediate_archive[hash_key]:
                failures.append(
                    f"Figure 3 v8.2 immediate archive {name} hash "
                    f"{observed_hash} != {expected_immediate_archive[hash_key]}")
        archive_sums_path = immediate_archive_path / "SHA256SUMS.txt"
        if archive_sums_path.is_file():
            for line_number, line in enumerate(
                    archive_sums_path.read_text().splitlines(), start=1):
                match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
                if match is None:
                    failures.append(
                        "Figure 3 v8.2 archive SHA256SUMS.txt has malformed line "
                        f"{line_number}")
                    continue
                expected_hash, member_rel = match.groups()
                member_rel_path = Path(member_rel)
                if member_rel_path.is_absolute() or ".." in member_rel_path.parts:
                    failures.append(
                        f"Figure 3 v8.2 archive checksum has unsafe path {member_rel!r}")
                    continue
                member_path = immediate_archive_path / member_rel_path
                observed_hash = (
                    hashlib.sha256(member_path.read_bytes()).hexdigest()
                    if member_path.is_file() else None)
                if observed_hash != expected_hash:
                    failures.append(
                        f"Figure 3 v8.2 archived member {member_rel!r} hash "
                        f"{observed_hash} != {expected_hash}")

        expected_moved_items = [
            "fig3_target_centred_effect_zoom_v7",
            "fig3_target_centred_effect_zoom_v7_1",
            "fig3_target_centred_effect_zoom_v7_2",
            "fig3_target_centred_effect_zoom_v7_3",
            "fig3_target_centred_effect_zoom_v7_4",
            "fig3_target_centred_effect_zoom_v7_5",
            "fig3_target_centred_effect_zoom_v8",
            "fig3_target_centred_effect_zoom_v8_1",
            "fig3_target_centred_effect_zoom_v7.guard.log",
            "fig3_target_centred_effect_zoom_v7_1.guard.log",
            "fig3_target_centred_effect_zoom_v7_2.guard.log",
            "fig3_target_centred_effect_zoom_v7_3.guard.log",
            "fig3_target_centred_effect_zoom_v7_4.guard.log",
            "fig3_target_centred_effect_zoom_v7_5.guard.log",
            "fig3_target_centred_effect_zoom_v8.guard.log",
            "fig3_target_centred_effect_zoom_v8_1.guard.log",
        ]
        expected_cleanup = {
            "recoverable_quarantine": "_legacy/quarantine_20260903_fig3_effect_zoom_drafts",
            "why_quarantined_sha256": (
                "b66e15d53d8d732ecde830cc26f967646a388d767728c461bbdeb3158db7da8a"),
            "inventory_sha256": (
                "db526a5e22a33deeb48e5cf78d19fbe81a37ea267d8838dcb7a00a439e3d960f"),
            "moved_items": expected_moved_items,
            "permanent_deletion_performed": False,
        }
        cleanup = effect_promotion.get("cleanup", {})
        if cleanup != expected_cleanup:
            failures.append("Figure 3 v8.2 recoverable cleanup contract changed")
        quarantine_path = ROOT / expected_cleanup["recoverable_quarantine"]
        for name, hash_key in (
                ("WHY_QUARANTINED.md", "why_quarantined_sha256"),
                ("INVENTORY.tsv", "inventory_sha256"),
        ):
            path = quarantine_path / name
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != expected_cleanup[hash_key]:
                failures.append(
                    f"Figure 3 v8.2 quarantine {name} hash "
                    f"{observed_hash} != {expected_cleanup[hash_key]}")
        observed_quarantine_items = (
            {path.name for path in quarantine_path.iterdir()}
            if quarantine_path.is_dir() else set())
        expected_quarantine_items = set(expected_moved_items) | {
            "WHY_QUARANTINED.md", "INVENTORY.tsv"}
        if observed_quarantine_items != expected_quarantine_items:
            failures.append(
                "Figure 3 v8.2 quarantine inventory changed; "
                f"observed={sorted(observed_quarantine_items)!r}")
        quarantine_inventory_path = quarantine_path / "INVENTORY.tsv"
        if quarantine_inventory_path.is_file():
            with quarantine_inventory_path.open(newline="") as handle:
                quarantine_reader = csv.DictReader(handle, delimiter="\t")
                quarantine_rows = list(quarantine_reader)
            expected_inventory_fields = [
                "kind", "sha256", "path", "files", "allocated_bytes"]
            if quarantine_reader.fieldnames != expected_inventory_fields:
                failures.append("Figure 3 v8.2 quarantine inventory columns changed")
            inventory_paths = [row.get("path") for row in quarantine_rows]
            if len(quarantine_rows) != 16 or inventory_paths != expected_moved_items or \
                    len(set(inventory_paths)) != 16:
                failures.append("Figure 3 v8.2 quarantine inventory rows changed")
            for row in quarantine_rows:
                member_name = row.get("path")
                member_path = quarantine_path / member_name \
                    if isinstance(member_name, str) else None
                expected_hash = row.get("sha256")
                if row.get("kind") == "file":
                    observed_hash = (
                        hashlib.sha256(member_path.read_bytes()).hexdigest()
                        if member_path is not None and member_path.is_file() else None)
                elif row.get("kind") == "tree" and member_path is not None and \
                        member_path.is_dir():
                    tree_digest_lines = []
                    for tree_file in sorted(
                            (path for path in member_path.rglob("*")
                             if path.is_file() and not path.is_symlink()),
                            key=lambda path: path.relative_to(member_path).as_posix()):
                        tree_rel = tree_file.relative_to(member_path).as_posix()
                        tree_hash = hashlib.sha256(tree_file.read_bytes()).hexdigest()
                        tree_digest_lines.append(f"{tree_hash}  ./{tree_rel}\n")
                    observed_hash = hashlib.sha256(
                        "".join(tree_digest_lines).encode()).hexdigest()
                else:
                    observed_hash = None
                if observed_hash != expected_hash:
                    failures.append(
                        f"Figure 3 v8.2 quarantined {member_name!r} hash "
                        f"{observed_hash} != {expected_hash}")

        expected_encoding_contract = {
            "orientation": {"lost": "red #c1272d", "gained": "teal #2a9d8f"},
            "target_burden": "Purples, fixed 3-18",
            "row_burden": "Blues, fixed 0-113",
            "outside_census": "bordered white",
            "stable_and_self": "split light grey #e4e7ea and dark grey #7f8a94",
            "effect_score": (
                "sign * sqrt(min(abs(lfc_Rest), abs(lfc_Stim48)) / "
                "2.21400360540582)"),
            "z_scored": False,
            "clipped": False,
            "shared_graphical_key_for_panels_c_d": True,
            "d_is_literal_c_first_80_columns": True,
            "ordering": "frozen dense-tail-v2 display order; no clustering or dendrogram",
            "panel_d_axis_font_points": 5.5,
            "all_visible_text_points_min": 5.0,
            "all_visible_text_points_max": 7.0,
            "rejected_5_25_point_trial": {
                "target_adjacent_overlap_pairs": 79,
                "regulator_adjacent_overlap_pairs": 48,
                "decision": "retain 5.5-point square cells",
            },
        }
        if effect_promotion.get("encoding_contract") != expected_encoding_contract:
            failures.append("Figure 3 v8.2 promotion encoding contract changed")
        expected_terminology_contract = {
            "panel_b": {
                "term": "single-orientation",
                "definition": (
                    "For a target at the downstream at-least-three threshold, every distinct "
                    "regulator reversing on it agrees on direction."),
            },
            "panels_e_f": {
                "term": "direction-concordant",
                "definition": (
                    "For the fixed CCNC/TADA2B/SUPT7L trio or AHR/ARNT pair, each shared "
                    "target moves the same way under those specific perturbations."),
            },
            "forbidden_reader_facing_terms": ["direction unanimity", "unanimous"],
        }
        if effect_promotion.get("terminology_contract") != expected_terminology_contract:
            failures.append("Figure 3 v8.2 promotion terminology contract changed")
        expected_scientific_scope = {
            "new_statistics_rerun": False,
            "new_clustering": False,
            "panel_c_shape": [49, 279],
            "panel_c_reversals": 537,
            "panel_d_shape": [49, 80],
            "panel_d_reversals": 208,
            "panels_a_b_e_f_g_scientific_artists_unchanged": True,
        }
        if effect_promotion.get("scientific_scope") != expected_scientific_scope:
            failures.append("Figure 3 v8.2 promotion scientific scope changed")

        expected_active_text_paths = {
            "manuscript": "MANUSCRIPT.md",
            "figure_legends": "FIGURE_LEGENDS.md",
            "paper_facts": "PAPER_FACTS.md",
        }
        active_effect_text = effect_promotion.get("active_text", {})
        if set(active_effect_text) != set(expected_active_text_paths):
            failures.append("Figure 3 v8.2 active-text record inventory changed")
        for key, expected_path in expected_active_text_paths.items():
            if active_effect_text.get(key, {}).get("path") != expected_path:
                failures.append(
                    f"Figure 3 v8.2 historical active text {key} path changed")

        expected_retention_policy = {
            "permanent_deletion_performed": False,
            "superseded_artifacts_archived_recoverably": True,
            "routine_cleanup_may_not_delete_archived_bundles": True,
        }
        if effect_promotion.get("retention_policy") != expected_retention_policy:
            failures.append("Figure 3 v8.2 retention policy changed")
        effect_guard_record = effect_promotion.get("integration_controls", {}).get(
            "active_claim_guard", {})
        if effect_guard_record.get("path") != "scripts/active_claim_guard.py":
            failures.append("Figure 3 v8.2 active-guard path changed")
        if effect_status == effect_pending_status:
            if effect_guard_record.get("sha256") is not None:
                failures.append("pending Figure 3 v8.2 promotion prematurely seals the guard")

        derived = effect_promotion.get("derived_artifacts", {})
        expected_derived_paths = {
            "pdf": "HumanCD4CoDEGNet_signflip.pdf",
            "docx": "HumanCD4CoDEGNet_signflip.docx",
        }
        if effect_status == effect_pending_status:
            if derived.get("status") != "PENDING_SOURCE_AND_POSTBUILD_GUARDS":
                failures.append("pending Figure 3 v8.2 promotion has wrong derived status")
            for key, expected_path in expected_derived_paths.items():
                if derived.get(key) != {"path": expected_path, "sha256": None}:
                    failures.append(
                        f"pending Figure 3 v8.2 promotion has invalid derived {key} record")
        elif effect_status == effect_final_status:
            if derived.get("status") != "POSTBUILD_GUARDED":
                failures.append("final Figure 3 v8.2 promotion is not post-build guarded")
            for key, expected_path in expected_derived_paths.items():
                if derived.get(key, {}).get("path") != expected_path or not re.fullmatch(
                        r"[0-9a-f]{64}", str(derived.get(key, {}).get("sha256", ""))):
                    failures.append(
                        f"historical Figure 3 v8.2 derived {key} record changed")

        # The joint record must name exactly this promoted Figure 3 while retaining its
        # untouched Figure 4 half and the older historical archive pointer checked above.
        if promotion:
            joint_fig3 = promotion.get("figure3", {})
            observed_effect_promotion_hash = hashlib.sha256(
                effect_promotion_path.read_bytes()).hexdigest()
            if joint_fig3.get("active_effect_zoom_promotion") != {
                    "path": effect_promotion_rel,
                    "sha256": observed_effect_promotion_hash,
            }:
                failures.append("joint promotion does not bind the active v8.2 promotion record")
            if joint_fig3.get("candidate_directory") != effect_candidate.get("directory") or \
                    joint_fig3.get("candidate_manifest") != effect_candidate.get("manifest") or \
                    joint_fig3.get("renderer") != effect_candidate.get("renderer") or \
                    joint_fig3.get("canonical_pdf") != effect_canonical.get("pdf") or \
                    joint_fig3.get("canonical_png") != effect_canonical.get("png") or \
                    joint_fig3.get("independent_audits") != [
                        effect_candidate.get("independent_audit")]:
                failures.append("joint promotion Figure 3 does not point exactly to v8.2")
            if joint_fig3.get("immediate_predecessor_archive") != immediate_archive:
                failures.append(
                    "joint promotion lost the v8.2 immediate-predecessor archive link")
            if promotion.get("active_text") != active_effect_text:
                failures.append("joint and standalone Figure 3 active-text records differ")
            joint_guard_record = promotion.get("integration_controls", {}).get(
                "active_claim_guard", {})
            if joint_guard_record.get("path") != "scripts/active_claim_guard.py":
                failures.append("joint promotion active-guard path changed")

        if check_derived and promotion.get("build_integration", {}).get(
                "status") != "COMPLETE":
            failures.append(
                "joint Figure 3/Figure 4 build integration is not complete in --post-build mode")

    # Figure 3 v10.4.1 is historical after the v11 promotion.  Continue to bind its
    # immutable candidate, promotion record and rollback archive, but do not compare its
    # former canonical or active-text paths with the current v11 files.  It was superseded
    # before its own post-build stage, so its pending derived status is historical truth.
    top_promotion_rel = (
        "data/census/fig3_top_overview_v10_4_1_promotion_20260904.json")
    top_promotion_path = ROOT / top_promotion_rel
    if not top_promotion_path.is_file():
        failures.append("Figure 3 v10.4.1 promotion record is missing")
    else:
        expected_top_promotion_hash = (
            "d53f5b7af556c3909adf3b5f1bdfd5d31c3989e8c0f4f40b7899c264bc6ca68f")
        observed_top_promotion_hash = hashlib.sha256(
            top_promotion_path.read_bytes()).hexdigest()
        if observed_top_promotion_hash != expected_top_promotion_hash:
            failures.append(
                "historical Figure 3 v10.4.1 promotion record hash changed")
        top_promotion = json.loads(top_promotion_path.read_text())
        expected_top_schema = "fig3-target-centred-top-overview-v10.4.1-promotion-v1"
        if top_promotion.get("schema_version") != expected_top_schema:
            failures.append("Figure 3 v10.3 promotion schema changed")
        top_status = top_promotion.get("status")
        top_pending_status = "PROMOTED_SOURCE_POSTBUILD_PENDING"
        top_final_status = "PROMOTED_AND_POSTBUILD_GUARDED"
        if top_status not in {top_pending_status, top_final_status}:
            failures.append("Figure 3 v10.3 promotion status is invalid")

        def check_top_record(
                record: dict, label: str, expected_path: str | None = None,
                expected_hash: str | None = None) -> str | None:
            if not isinstance(record, dict):
                failures.append(f"Figure 3 v10.3 {label} record is not an object")
                return None
            rel = record.get("path")
            if expected_path is not None and rel != expected_path:
                failures.append(
                    f"Figure 3 v10.3 {label} path {rel!r} != {expected_path!r}")
            if not isinstance(rel, str) or not rel or Path(rel).is_absolute() or \
                    ".." in Path(rel).parts:
                failures.append(
                    f"Figure 3 v10.3 {label} has unsafe relative path {rel!r}")
                return None
            recorded_hash = record.get("sha256")
            if not isinstance(recorded_hash, str) or not re.fullmatch(
                    r"[0-9a-f]{64}", recorded_hash):
                failures.append(
                    f"Figure 3 v10.3 {label} has invalid SHA-256 {recorded_hash!r}")
                return None
            if expected_hash is not None and recorded_hash != expected_hash:
                failures.append(
                    f"Figure 3 v10.3 {label} recorded hash {recorded_hash} "
                    f"!= {expected_hash}")
            path = ROOT / rel
            observed_hash = (
                hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
            if observed_hash != recorded_hash:
                failures.append(
                    f"Figure 3 v10.3 {label} {rel!r} hash "
                    f"{observed_hash} != {recorded_hash}")
            return observed_hash

        expected_top_supersedes = {
            "path": "data/census/fig3_effect_zoom_v8_2_promotion_20260903.json",
            "sha256": "a332c154b68032a4a75e207d0c7ee37bf2ad21a18d791d19d4da94f2b4368f94",
        }
        if top_promotion.get("supersedes") != expected_top_supersedes:
            failures.append("Figure 3 v10.3 lost its exact v8.2 predecessor record")
        check_top_record(
            top_promotion.get("supersedes", {}), "superseded v8.2 promotion",
            expected_top_supersedes["path"], expected_top_supersedes["sha256"])

        top_candidate_dir_rel = (
            "data/figure_provenance/fig3h_target_centred_top_overview_layout_v10_4_1")
        top_candidate = top_promotion.get("candidate", {})
        if top_candidate.get("directory") != top_candidate_dir_rel:
            failures.append("Figure 3 v10.3 promotion points to the wrong candidate directory")
        expected_top_candidate_records = {
            "manifest": (
                f"{top_candidate_dir_rel}/manifest.json",
                "e90b58662bc7dd3eed709323bfa55ca44daed84a800325bc77e348b47120020b"),
            "pdf": (
                f"{top_candidate_dir_rel}/Fig3_target_centred_top_overview_layout_v10_4_1.pdf",
                "8652a5130cbe57c00a2c34337eabc2c177312d8b91c803195862baca6d2e4627"),
            "png": (
                f"{top_candidate_dir_rel}/Fig3_target_centred_top_overview_layout_v10_4_1.png",
                "8e3e4ab040e55cda511637c024d739970b06df686c4d48ff776c52149bff57de"),
            "renderer": (
                "scripts/fig3_target_centred_top_overview_layout_v10_4_1.py",
                "79c9747af744bd94cc2c12d6c4f00279ef2651c117048e2b377c6f7b88867858"),
            "render_guard": (
                "data/figure_provenance/"
                "fig3h_target_centred_top_overview_layout_v10_4_1.guard.log",
                "83dba762310d6565a0812ed8a7e25d38474402189bc1c0dbe6076baa59688391"),
            "geometry_audit": (
                "data/figure_provenance/"
                "fig3h_target_centred_top_overview_layout_v10_4_1_geometry_final.json",
                "723c23b52264a200c2bd69867a98c4af5c67db0802c3dd53feeac00498f31047"),
            "geometry_guard": (
                "data/figure_provenance/"
                "fig3h_target_centred_top_overview_layout_v10_4_1_geometry_final.guard.log",
                "524ad5afcacbb63b71d0721b939755c232c6d57e7b37eafe8c8f3f3d9bf708d5"),
            "source_pdf_audit": (
                "data/figure_provenance/"
                "fig3h_target_centred_top_overview_layout_v10_4_1_source_pdf_final.json",
                "48363529ceb1866b246ce89fa60f420fd2e28fbaab25b3864efcee0f6d959794"),
            "source_pdf_guard": (
                "data/figure_provenance/"
                "fig3h_target_centred_top_overview_layout_v10_4_1_source_pdf_final.guard.log",
                "234145ad0ce9b34a4d89ce7053a0e19f2464173ccf4ed160022f606adab12f33"),
            "audit_renderer": (
                "scripts/audit_fig3_top_overview_v10_4_1.py",
                "6a5978e377cbb54a517fc1fc50460797fec91e521913cc5f9f0faf5a657a070d"),
        }
        if set(top_candidate) != {"directory", *expected_top_candidate_records}:
            failures.append("Figure 3 v10.3 candidate record inventory changed")
        for key, (expected_path, expected_hash) in expected_top_candidate_records.items():
            check_top_record(
                top_candidate.get(key, {}), f"candidate {key}",
                expected_path, expected_hash)

        top_canonical = top_promotion.get("canonical", {})
        expected_top_canonical_paths = {
            "pdf": "figures/Fig3.pdf",
            "png": "figures/Fig3.png",
        }
        for key, expected_path in expected_top_canonical_paths.items():
            candidate_record = top_candidate.get(key, {})
            canonical_record = top_canonical.get(key, {})
            if canonical_record.get("path") != expected_path or \
                    canonical_record.get("sha256") != candidate_record.get("sha256"):
                failures.append(
                    f"historical Figure 3 v10.4.1 canonical {key} record changed")
        if top_canonical.get("promoted_by_byte_copy") is not True or \
                top_canonical.get("candidate_and_canonical_bytes_identical") is not True:
            failures.append("Figure 3 v10.3 canonical byte-copy contract changed")

        expected_top_archive_rel = (
            "figures/repository/superseded_figure3_candidates_2026-09-04/"
            "pre_top_overview_v10_3_20260904")
        expected_top_archive_hashes = {
            "FIGURE_LEGENDS.md": (
                "524e5fb5d662298da013b363713126d53014a7812dc51fbf021b06622274905b"),
            "Fig3.pdf": (
                "48cbf1202c51604d94d6fa7ffb5420d6d6699edd5f9d3cdd2a07a708ff2afa91"),
            "Fig3.png": (
                "df9622531b3c74dec29bda73425096ee071bb3a2ef2db27893511a2279321618"),
            "MANUSCRIPT.md": (
                "314ec095387413a52e6f9bfb536d1582f1f4b2173a676dbefc0b54a4c798a5cf"),
            "PAPER_FACTS.md": (
                "683e0fc8d3d982e7677a1c54f977c8226c1fd1eeb6a42d3a78e12d63dafa7440"),
            "README.md": (
                "74fb971514814acd842d2808d68826b5bc1c0c5884a1829f34d59f7d85065bf3"),
            "SHA256SUMS.txt": (
                "e7d9f3468b13872f463c4728ac0cc81b6f20d47e6638209109d445d7bdd929ae"),
            "active_claim_guard.py": (
                "97e46133c78bca8ca1f17bad636894c8706c6116c6a433fad60ae69bb6c42632"),
            "fig3_effect_zoom_v8_2_promotion_20260903.json": (
                "a332c154b68032a4a75e207d0c7ee37bf2ad21a18d791d19d4da94f2b4368f94"),
            "promotion_manifest_20260901_scope_link_and_functional_groups.json": (
                "132d0a2d9a10e3aa9d417e89dd1d421a26ba1e31a1565fc79b42680759756b78"),
        }
        top_archive = top_promotion.get("immediate_predecessor_archive", {})
        if top_archive.get("path") != expected_top_archive_rel or \
                top_archive.get("recoverable") is not True or \
                top_archive.get("permanent_deletion_performed") is not False:
            failures.append("Figure 3 v10.3 immediate-predecessor archive contract changed")
        archive_records = top_archive.get("records", {})
        if set(archive_records) != set(expected_top_archive_hashes):
            failures.append("Figure 3 v10.3 archive record inventory changed")
        for name, expected_hash in expected_top_archive_hashes.items():
            check_top_record(
                archive_records.get(name, {}), f"archived {name}",
                f"{expected_top_archive_rel}/{name}", expected_hash)
        top_archive_path = ROOT / expected_top_archive_rel
        top_archive_sums = top_archive_path / "SHA256SUMS.txt"
        if top_archive_sums.is_file():
            for line_number, line in enumerate(
                    top_archive_sums.read_text().splitlines(), start=1):
                match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
                if match is None:
                    failures.append(
                        "Figure 3 v10.3 archive SHA256SUMS.txt has malformed line "
                        f"{line_number}")
                    continue
                expected_hash, member_rel = match.groups()
                member_path_fragment = Path(member_rel)
                if member_path_fragment.is_absolute() or \
                        ".." in member_path_fragment.parts:
                    failures.append(
                        f"Figure 3 v10.3 archive checksum has unsafe path {member_rel!r}")
                    continue
                member_path = top_archive_path / member_path_fragment
                observed_hash = (
                    hashlib.sha256(member_path.read_bytes()).hexdigest()
                    if member_path.is_file() else None)
                if observed_hash != expected_hash:
                    failures.append(
                        f"Figure 3 v10.3 archived member {member_rel!r} hash "
                        f"{observed_hash} != {expected_hash}")

        top_manifest_path = ROOT / top_candidate_dir_rel / "manifest.json"
        if top_manifest_path.is_file():
            top_manifest = json.loads(top_manifest_path.read_text())
            if top_manifest.get("schema_version") != \
                    "fig3-target-centred-top-overview-layout-v10.4.1" or \
                    top_manifest.get("status") != "DIAGNOSTIC_ONLY_NOT_PROMOTED":
                failures.append("Figure 3 v10.3 immutable candidate identity changed")
            if top_manifest.get("layout_order") != {
                    "top": ["a", "b", "c"],
                    "middle": ["d", "e"],
                    "bottom": ["f", "g", "h"]}:
                failures.append("Figure 3 v10.3 layout order changed")
            if sorted(top_manifest.get("panel_map", {})) != list("abcdefgh"):
                failures.append("Figure 3 v10.3 panel map is not exactly a--h")
            expected_top_panel_counts = {
                "panel_c": {
                    "cells": 13671, "reversal": 537, "lost": 330, "gained": 207,
                    "stable": 430, "outside": 12703, "self": 1,
                    "positive_rows": 43, "positive_targets": 257,
                },
                "panel_d": {
                    "cells": 3920, "reversal": 208, "lost": 137, "gained": 71,
                    "stable": 128, "outside": 3584, "self": 0,
                    "positive_rows": 37, "positive_targets": 77,
                },
            }
            top_semantic = top_manifest.get("semantic_counts", {})
            for key, expected in expected_top_panel_counts.items():
                if top_semantic.get(key) != expected:
                    failures.append(f"Figure 3 v10.3 {key} semantic counts changed")
            top_panel_e = top_semantic.get("panel_e", {})
            expected_zero_rows = [
                {"go_id": "GO:0071554",
                 "go_term": "cell wall organization or biogenesis"},
                {"go_id": "GO:0019748",
                 "go_term": "secondary metabolic process"},
            ]
            expected_panel_e_values = {
                "targets": 279,
                "mapped_targets": 240,
                "outside_this_slim_targets": 39,
                "source_official_bp_slim_terms": 25,
                "source_rows": 26,
                "displayed_nonzero_official_bp_slim_terms": 23,
                "displayed_rows": 24,
                "explicit_complement_rows": 1,
                "lost_denominator": 161,
                "gained_denominator": 118,
                "memberships_overlap": True,
                "printed_cell_value": "unique-target count",
                "teal": "#2a9d8f",
            }
            for key, expected in expected_panel_e_values.items():
                if top_panel_e.get(key) != expected:
                    failures.append(
                        f"Figure 3 v10.3 panel-e {key}={top_panel_e.get(key)!r}; "
                        f"expected {expected!r}")
            if top_panel_e.get("omitted_zero_rows") != expected_zero_rows:
                failures.append("Figure 3 v10.3 panel-e omitted-zero rows changed")
            for required_boundary in (
                    "descriptive frozen GO-BP-slim membership", "no enrichment",
                    "lost-versus-gained test", "pathway activity", "mechanism"):
                if required_boundary not in top_panel_e.get("claim_boundary", ""):
                    failures.append(
                        f"Figure 3 v10.3 panel-e boundary lacks {required_boundary!r}")

            top_geometry = top_manifest.get("geometry", {})
            if top_geometry.get("literal_zoom_assertions") != {
                    "states_equal_c_first_80": True,
                    "scores_equal_c_first_80": True,
                    "targets_equal_c_first_80": True,
                    "row_order_equal": True,
                    "orientation_track_equal_c_first_80": True,
                    "target_burden_track_equal_c_first_80": True,
                    "row_burden_track_identical_full_panel_scale": True}:
                failures.append("Figure 3 v10.3 literal c[:, :80] contract changed")
            heat_geometry = top_geometry.get("heat_geometry", {})
            for letter, expected_side in (("c", 0.92), ("d", 5.5)):
                observed = heat_geometry.get(letter, {})
                if abs(float(observed.get("cell_width_points", 0)) - expected_side) > 1e-12 or \
                        abs(float(observed.get("cell_height_points", 0)) - expected_side) > 1e-12:
                    failures.append(
                        f"Figure 3 v10.3 panel-{letter} square-cell geometry changed")
            panel_d_geometry = top_geometry.get("panel_d", {})
            if panel_d_geometry.get("label_points") != 5.5 or \
                    panel_d_geometry.get("target_label_collision_audit", {}).get(
                        "overlap_pairs_above_0.2px_each_axis") != 0 or \
                    panel_d_geometry.get("regulator_label_collision_audit", {}).get(
                        "overlap_pairs_above_0.2px_each_axis") != 0:
                failures.append(
                    "Figure 3 v10.4.1 panel-d 5.5-pt label/collision contract changed")
            if top_geometry.get("panel_d_header_clearance", {}).get(
                    "all_clear") is not True or \
                    float(top_geometry.get("content_bbox", {}).get(
                        "right_margin_inches", -1)) <= 0:
                failures.append(
                    "Figure 3 v10.4.1 header clearance or page-fit contract changed")
            text_audit = top_geometry.get("text_audit", {})
            if text_audit.get("min_visible_text_points") != 5.0 or \
                    text_audit.get("max_visible_text_points") != 7.0 or \
                    text_audit.get("off_page_text") != [] or \
                    top_geometry.get("all_text_bbox_audit", {}).get(
                        "collision_count") != 0:
                failures.append("Figure 3 v10.3 violates the 5--7 pt/text-fit contract")
            top_export = top_manifest.get("export", {})
            if top_export.get("page") != "custom short portrait" or \
                    top_export.get("page_inches") != [8.267716535433072, 9.6] or \
                    top_export.get("fresh_canvas_count") != 2 or \
                    top_export.get("pdf_byte_identical") is not True or \
                    top_export.get("png_byte_identical") is not True or \
                    top_export.get("pdf_sha256") != expected_top_candidate_records[
                        "pdf"][1] or \
                    top_export.get("png_sha256") != expected_top_candidate_records[
                        "png"][1]:
                failures.append("Figure 3 v10.3 compact deterministic export contract changed")
            if top_manifest.get("no_write_guard", {}).get("all_equal") is not True or \
                    top_manifest.get("canonical_changes") is not False or \
                    top_manifest.get("manuscript_changes") is not False or \
                    top_manifest.get("legend_changes") is not False or \
                    top_manifest.get("trash_changes") is not False:
                failures.append("Figure 3 v10.3 diagnostic no-write custody changed")

            top_outputs = top_manifest.get("outputs", {})
            expected_top_outputs = {
                "Fig3_target_centred_top_overview_layout_v10_4_1.pdf",
                "Fig3_target_centred_top_overview_layout_v10_4_1.png",
                "source_panel_a_target_concentration.csv",
                "source_panel_b_target_unanimity.csv",
                "source_panel_c_effect_score_matrix.tsv",
                "source_panel_d_effect_score_matrix.tsv",
                "source_panel_e_go_bp_slim_summary.tsv",
                "source_panel_f_cofactor_targets.csv",
                "source_panel_g_sensor_targets.csv",
                "source_panel_h_arce_recurrence.csv",
            }
            if set(top_outputs) != expected_top_outputs:
                failures.append("Figure 3 v10.3 output inventory changed")
            for name, record in top_outputs.items():
                path = ROOT / top_candidate_dir_rel / name
                observed_hash = (
                    hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
                observed_size = path.stat().st_size if path.is_file() else None
                if observed_hash != record.get("sha256") or \
                        observed_size != record.get("bytes"):
                    failures.append(
                        f"Figure 3 v10.3 output {name!r} changed")
            top_bundled = top_manifest.get("bundled_inputs_and_renderers", {})
            if len(top_bundled) != 27:
                failures.append("Figure 3 v10.3 bundled source inventory changed")
            for name, record in top_bundled.items():
                path = ROOT / top_candidate_dir_rel / name
                observed_hash = (
                    hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
                observed_size = path.stat().st_size if path.is_file() else None
                if observed_hash != record.get("sha256") or \
                        observed_size != record.get("bytes"):
                    failures.append(
                        f"Figure 3 v10.3 bundled source {name!r} changed")

            go_source_path = (
                ROOT / top_candidate_dir_rel / "source_panel_e_go_bp_slim_summary.tsv")
            if go_source_path.is_file():
                with go_source_path.open(newline="") as handle:
                    go_rows = list(csv.DictReader(handle, delimiter="\t"))
                zero_rows = [
                    {"go_id": row.get("go_id"), "go_term": row.get("go_term")}
                    for row in go_rows
                    if row.get("row_type") == "official_bp_slim"
                    and int(row.get("lost_targets", "-1")) == 0
                    and int(row.get("gained_targets", "-1")) == 0]
                official_rows = [
                    row for row in go_rows if row.get("row_type") == "official_bp_slim"]
                complement_rows = [
                    row for row in go_rows if row.get("row_type") == "explicit_complement"]
                if len(go_rows) != 26 or len(official_rows) != 25 or \
                        len(complement_rows) != 1 or zero_rows != expected_zero_rows:
                    failures.append("Figure 3 v10.3 panel-e source-table closure changed")

        expected_top_panel_order = {
            "a": "target concentration",
            "b": "single-orientation",
            "c": "complete 49 x 279 effect map",
            "d": "literal first-80 enlargement",
            "e": "descriptive fixed GO-BP-slim membership",
            "f": "CCNC/TADA2B/SUPT7L direction-concordant targets",
            "g": "AHR/ARNT direction-concordant targets",
            "h": "Arce edge recurrence",
        }
        if top_promotion.get("panel_order") != expected_top_panel_order:
            failures.append("Figure 3 v10.3 promotion panel-order contract changed")
        expected_top_claim_boundary = {
            "panel_b_term": "single-orientation",
            "panels_f_g_term": "direction-concordant",
            "panel_e_is_enrichment": False,
            "panel_e_is_directional_test": False,
            "panel_e_is_pathway_activity": False,
            "panel_e_is_mechanism": False,
            "panel_e_memberships_overlap": True,
            "panel_e_outside_means_biologically_unannotated": False,
        }
        if top_promotion.get("claim_boundary") != expected_top_claim_boundary:
            failures.append("Figure 3 v10.3 promotion claim boundary changed")

        expected_top_active_paths = {
            "manuscript": "MANUSCRIPT.md",
            "figure_legends": "FIGURE_LEGENDS.md",
            "paper_facts": "PAPER_FACTS.md",
        }
        active_top_text = top_promotion.get("active_text", {})
        if set(active_top_text) != set(expected_top_active_paths):
            failures.append("Figure 3 v10.3 active-text inventory changed")
        for key, expected_path in expected_top_active_paths.items():
            record_value = active_top_text.get(key, {})
            if record_value.get("path") != expected_path or not re.fullmatch(
                    r"[0-9a-f]{64}", str(record_value.get("sha256", ""))):
                failures.append(
                    f"historical Figure 3 v10.4.1 active text {key} record changed")
        historical_top_guard = top_promotion.get("integration_controls", {}).get(
            "active_claim_guard", {})
        if historical_top_guard.get("path") != "scripts/active_claim_guard.py" or \
                not re.fullmatch(
                    r"[0-9a-f]{64}", str(historical_top_guard.get("sha256", ""))):
            failures.append("historical Figure 3 v10.4.1 guard record changed")

        expected_top_retention = {
            "permanent_deletion_performed": False,
            "superseded_artifacts_archived_recoverably": True,
            "routine_cleanup_may_not_delete_archived_bundles": True,
            "candidate_bundle_retained": True,
        }
        if top_promotion.get("retention_policy") != expected_top_retention:
            failures.append("Figure 3 v10.3 retention policy changed")
        top_derived = top_promotion.get("derived_artifacts", {})
        expected_derived_paths = {
            "pdf": "HumanCD4CoDEGNet_signflip.pdf",
            "docx": "HumanCD4CoDEGNet_signflip.docx",
        }
        if top_status == top_pending_status:
            if top_derived.get("status") != "NOT_REBUILT_FOR_THIS_PROMOTION":
                failures.append("pending Figure 3 v10.3 derived status changed")
            for key, expected_path in expected_derived_paths.items():
                if top_derived.get(key) != {"path": expected_path, "sha256": None}:
                    failures.append(
                        f"pending Figure 3 v10.3 derived {key} record changed")
        elif top_status == top_final_status:
            if top_derived.get("status") != "POSTBUILD_GUARDED":
                failures.append("final Figure 3 v10.3 derived status changed")
            for key, expected_path in expected_derived_paths.items():
                record = top_derived.get(key, {})
                if record.get("path") != expected_path or not re.fullmatch(
                        r"[0-9a-f]{64}", str(record.get("sha256", ""))):
                    failures.append(
                        f"final Figure 3 v10.3 derived {key} record changed")

    # Figure 3 v11 is the live Figure 3 contract.  It replaces only panel e with
    # network-adjusted all-target functional rank enrichment, keeps a--d/f--h
    # byte-identical at the rendered-pixel/source level, and preserves v10.4.1
    # without rewriting its historical promotion record.
    # Figure 4f is a later, independently guarded promotion. Once its manifest exists,
    # the Figure 3 v11 records below become historical for the three active text files,
    # the active guard, build_docx and the main manuscript artifacts. The immutable v11
    # candidate copies remain fully verified; the new Figure 4f manifest directly binds
    # every superseding live byte.
    fig4f_supersedes_v11_live_state = (
        ROOT / "data/census/fig4f_landmark_census_promotion_v1.json"
    ).is_file()
    fig2_signed_kegg_supersedes_v11_live_state = (
        ROOT / FIG2_SIGNED_KEGG_V11_PROMOTION
    ).is_file()
    later_promotion_supersedes_v11_live_state = (
        fig4f_supersedes_v11_live_state or
        fig2_signed_kegg_supersedes_v11_live_state
    )
    v12_supersedes_v11_canonical = (
        ROOT / "data/census/"
        "fig3_target_centred_compact_zoom_v12_promotion_20260906.json"
    ).is_file()
    v11_promotion_rel = (
        "data/census/"
        "fig3_target_centred_alltarget_enrichment_v11_promotion_20260904.json")
    v11_promotion_path = ROOT / v11_promotion_rel
    if not v11_promotion_path.is_file():
        failures.append("Figure 3 v11 promotion record is missing")
    else:
        v11 = json.loads(v11_promotion_path.read_text())
        v11_pending = "PROMOTED_SOURCE_POSTBUILD_PENDING"
        v11_final = "PROMOTED_AND_POSTBUILD_GUARDED"
        if v11.get("schema_version") != \
                "fig3-target-centred-alltarget-enrichment-v11-promotion-v1":
            failures.append("Figure 3 v11 promotion schema changed")
        if v11.get("status") not in {v11_pending, v11_final}:
            failures.append("Figure 3 v11 promotion status is invalid")
        if check_derived and v11.get("status") != v11_final:
            failures.append("Figure 3 v11 is not final in --post-build mode")

        def check_v11_record(
                record_value: dict, label: str, expected_path: str | None = None,
                expected_hash: str | None = None,
                live: bool = True) -> str | None:
            if not isinstance(record_value, dict):
                failures.append(f"Figure 3 v11 {label} record is not an object")
                return None
            rel = record_value.get("path")
            if expected_path is not None and rel != expected_path:
                failures.append(
                    f"Figure 3 v11 {label} path {rel!r} != {expected_path!r}")
            if not isinstance(rel, str) or not rel or Path(rel).is_absolute() or \
                    ".." in Path(rel).parts:
                failures.append(f"Figure 3 v11 {label} has unsafe path {rel!r}")
                return None
            recorded_hash = record_value.get("sha256")
            if not isinstance(recorded_hash, str) or not re.fullmatch(
                    r"[0-9a-f]{64}", recorded_hash):
                failures.append(
                    f"Figure 3 v11 {label} has invalid SHA-256 {recorded_hash!r}")
                return None
            if expected_hash is not None and recorded_hash != expected_hash:
                failures.append(
                    f"Figure 3 v11 {label} recorded hash {recorded_hash} "
                    f"!= {expected_hash}")
            path = ROOT / rel
            observed_hash = recorded_hash
            if live:
                observed_hash = (
                    hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None)
                if observed_hash != recorded_hash:
                    failures.append(
                        f"Figure 3 v11 {label} {rel!r} hash "
                        f"{observed_hash} != {recorded_hash}")
                recorded_bytes = record_value.get("bytes")
                if recorded_bytes is not None:
                    observed_bytes = path.stat().st_size if path.is_file() else None
                    if observed_bytes != recorded_bytes:
                        failures.append(
                            f"Figure 3 v11 {label} {rel!r} size "
                            f"{observed_bytes} != {recorded_bytes}")
            return observed_hash

        expected_v10_link = {
            "path": top_promotion_rel,
            "sha256": (
                "d53f5b7af556c3909adf3b5f1bdfd5d31c3989e8c0f4f40b7899c264bc6ca68f"),
            "status_at_supersession": "PROMOTED_SOURCE_POSTBUILD_PENDING",
            "historical_without_rewrite": True,
        }
        if v11.get("supersedes") != expected_v10_link:
            failures.append("Figure 3 v11 lost its exact historical v10.4.1 link")
        check_v11_record(
            v11.get("supersedes", {}), "historical v10.4.1 promotion",
            top_promotion_rel, expected_v10_link["sha256"])

        v11_dir_rel = (
            "data/figure_provenance/fig3_target_centred_alltarget_enrichment_v11")
        v11_candidate = v11.get("candidate", {})
        expected_v11_records = {
            "manifest": (
                f"{v11_dir_rel}/manifest.json",
                "54e5cbe31003f81d06918682c3ffce90fc69c01c9995b97de8f917c38a45607c"),
            "pdf": (
                f"{v11_dir_rel}/Fig3_target_centred_alltarget_enrichment_v11.pdf",
                "3582a9dc256e9d070701667a1d67a24699fd5bf43375a92abaeada119216b21a"),
            "png": (
                f"{v11_dir_rel}/Fig3_target_centred_alltarget_enrichment_v11.png",
                "67ae02e5d0d9f0ad7e20382b2bd45b56e10b367fe84c403a7f594142c7fded46"),
            "renderer": (
                "scripts/fig3_target_centred_alltarget_enrichment_v11.py",
                "2b0936bc9855e2eba5689e236982473974f3f0cce71f491f7463d78a29c2e845"),
        }
        if v11_candidate.get("directory") != v11_dir_rel:
            failures.append("Figure 3 v11 points to the wrong candidate directory")
        for key, (expected_path, expected_hash) in expected_v11_records.items():
            check_v11_record(
                v11_candidate.get(key, {}), f"candidate {key}",
                expected_path, expected_hash)

        expected_v11_files = {
            "manifest.json",
            "Fig3_target_centred_alltarget_enrichment_v11.pdf",
            "Fig3_target_centred_alltarget_enrichment_v11.png",
            "input_panel_e_authoritative_terms.tsv",
            "input_panel_e_overlay_v2.json",
            "input_panel_e_supplement_package_final_audit.json",
            "parent_renderer_v10_4_1.py",
            "parent_v10_4_1_manifest.json",
            "renderer_fig3_target_centred_alltarget_enrichment_v11.py",
            "source_panel_a_target_concentration.csv",
            "source_panel_b_target_unanimity.csv",
            "source_panel_c_effect_score_matrix.tsv",
            "source_panel_d_effect_score_matrix.tsv",
            "source_panel_e_alltarget_term_composition.tsv",
            "source_panel_e_overlay_generator.py",
            "source_panel_e_standalone_renderer.py",
            "source_panel_f_cofactor_targets.csv",
            "source_panel_g_sensor_targets.csv",
            "source_panel_h_arce_recurrence.csv",
        }
        if set(v11_candidate.get("exact_file_inventory", [])) != expected_v11_files or \
                v11_candidate.get("unique_file_count") != 19 or \
                v11_candidate.get("manifest_outputs") != 12 or \
                v11_candidate.get("manifest_copied_inputs_and_renderers") != 8:
            failures.append("Figure 3 v11 promotion candidate inventory changed")
        v11_dir = ROOT / v11_dir_rel
        observed_v11_files = {
            str(path.relative_to(v11_dir)) for path in v11_dir.rglob("*")
            if path.is_file()} if v11_dir.is_dir() else set()
        if observed_v11_files != expected_v11_files:
            failures.append(
                "Figure 3 v11 on-disk candidate inventory changed; "
                f"extra={sorted(observed_v11_files - expected_v11_files)!r}, "
                f"missing={sorted(expected_v11_files - observed_v11_files)!r}")

        v11_manifest_path = v11_dir / "manifest.json"
        if v11_manifest_path.is_file():
            v11_manifest = json.loads(v11_manifest_path.read_text())
            if v11_manifest.get("schema_version") != \
                    "fig3-target-centred-alltarget-enrichment-v11" or \
                    v11_manifest.get("status") != "DIAGNOSTIC_ONLY_NOT_PROMOTED":
                failures.append("Figure 3 v11 immutable candidate identity changed")
            if v11_manifest.get("layout_order") != {
                    "top": ["a", "b", "c"],
                    "middle": ["d", "e"],
                    "bottom": ["f", "g", "h"]}:
                failures.append("Figure 3 v11 a--h layout changed")
            if sorted(v11_manifest.get("panel_map", {})) != list("abcdefgh"):
                failures.append("Figure 3 v11 panel map is not exactly a--h")

            outputs = v11_manifest.get("outputs", {})
            copied = v11_manifest.get("source_closure", {}).get(
                "copied_inputs_and_renderers", {})
            if len(outputs) != 12 or len(copied) != 8:
                failures.append("Figure 3 v11 manifest source inventory count changed")
            for label, records in (("output", outputs), ("copied source", copied)):
                for name, record_value in records.items():
                    path = v11_dir / name
                    observed_hash = (
                        hashlib.sha256(path.read_bytes()).hexdigest()
                        if path.is_file() else None)
                    observed_size = path.stat().st_size if path.is_file() else None
                    if observed_hash != record_value.get("sha256") or \
                            observed_size != record_value.get("bytes"):
                        failures.append(
                            f"Figure 3 v11 {label} {name!r} changed")
                    if label == "copied source" and not \
                            later_promotion_supersedes_v11_live_state:
                        source_rel = record_value.get("source")
                        source_path = ROOT / source_rel if isinstance(source_rel, str) else None
                        source_hash = (
                            hashlib.sha256(source_path.read_bytes()).hexdigest()
                            if source_path is not None and source_path.is_file() else None)
                        if source_hash != record_value.get("sha256"):
                            failures.append(
                                f"Figure 3 v11 upstream source {source_rel!r} changed")

            panel_e = v11_manifest.get("panel_e", {})
            expected_panel_e = {
                "all_target_universe": 9554,
                "targets_with_zero_observed_reversals": 6977,
                "displayed_terms": 19,
                "kegg_repeatable_fdr_terms": 9,
                "kegg_global_maxT_stable_terms": 0,
                "go_bp_statement_eligible_terms_total": 161,
                "go_bp_global_maxT_stable_terms": 10,
                "term_label_points": 5.0,
                "heading_points": 7.0,
                "false_overlap_asterisk_or_footnote_rendered": False,
                "evidence_tiers_visibly_separated": True,
                "c_d_key_preserved_at_v10_4_1_position": True,
            }
            for key, expected in expected_panel_e.items():
                if panel_e.get(key) != expected:
                    failures.append(
                        f"Figure 3 v11 panel-e {key}={panel_e.get(key)!r}; "
                        f"expected {expected!r}")
            for required_boundary in (
                    "all-target network-adjusted pooled reversal-burden enrichment",
                    "do not test lost-versus-gained enrichment", "pathway activity",
                    "mechanism", "individual edges", "279 targets in panels c-d",
                    "Terms overlap"):
                if required_boundary not in v11_manifest.get("claim_boundary", ""):
                    failures.append(
                        f"Figure 3 v11 boundary lacks {required_boundary!r}")

            geometry = v11_manifest.get("geometry", {})
            text_audit = geometry.get("text_audit", {})
            if text_audit.get("min_visible_text_points") != 5.0 or \
                    text_audit.get("max_visible_text_points") != 7.0 or \
                    text_audit.get("off_page_text") != [] or \
                    geometry.get("all_text_bbox_audit", {}).get(
                        "collision_count") != 0:
                failures.append("Figure 3 v11 violates the 5--7 pt/text-fit contract")
            panel_d = geometry.get("panel_d", {})
            if panel_d.get("label_points") != 5.5 or \
                    panel_d.get("shape") != [49, 80]:
                failures.append("Figure 3 v11 panel-d zoom/typography changed")
            raster = v11_manifest.get("unchanged_content_proof", {}).get(
                "saved_png_raster", {})
            if raster.get("outside_panel_e_changed_pixels") != 0 or \
                    raster.get("c_d_key_changed_pixels") != 0 or \
                    v11_manifest.get("unchanged_content_proof", {}).get(
                        "panels_a_d_and_f_h_preserved") is not True:
                failures.append("Figure 3 v11 unchanged-panel proof changed")
            no_write = v11_manifest.get("no_write_guard", {})
            if no_write.get("all_equal") is not True or \
                    no_write.get("protected_before") != no_write.get("protected_after") or \
                    v11_manifest.get("canonical_changes") is not False or \
                    v11_manifest.get("manuscript_changes") is not False or \
                    v11_manifest.get("legend_changes") is not False:
                failures.append("Figure 3 v11 diagnostic no-write custody changed")

            source_e = v11_dir / "source_panel_e_alltarget_term_composition.tsv"
            if source_e.is_file():
                with source_e.open(newline="") as handle:
                    rows_e = list(csv.DictReader(handle, delimiter="\t"))
                if len(rows_e) != 19 or \
                        [int(row["display_order"]) for row in rows_e] != list(range(1, 20)) or \
                        [row["library"] for row in rows_e[:9]] != ["KEGG_Legacy"] * 9 or \
                        [row["library"] for row in rows_e[9:]] != \
                        ["GO_Biological_Process"] * 10:
                    failures.append("Figure 3 v11 panel-e row/group order changed")
                for row in rows_e:
                    orientation_total = sum(int(row[key]) for key in (
                        "n_lost", "n_gained", "n_mixed", "n_no_reversal"))
                    if orientation_total != int(row["n_members_in_universe"]):
                        failures.append(
                            f"Figure 3 v11 panel-e counts do not close for {row['term']!r}")

        v11_canonical = v11.get("canonical", {})
        for key, expected_path in (
                ("pdf", "figures/Fig3.pdf"),
                ("png", "figures/Fig3.png")):
            check_v11_record(
                v11_canonical.get(key, {}), f"canonical {key}", expected_path,
                expected_v11_records[key][1], live=not v12_supersedes_v11_canonical)
        if v11_canonical.get("promoted_by_byte_copy") is not True or \
                v11_canonical.get("candidate_and_canonical_bytes_identical") is not True:
            failures.append("Figure 3 v11 canonical byte-copy contract changed")

        def check_v11_archive(archive: dict, label: str, expected_path: str) -> None:
            if archive.get("path") != expected_path or \
                    archive.get("recoverable") is not True or \
                    archive.get("permanent_deletion_performed") is not False:
                failures.append(f"Figure 3 v11 {label} archive contract changed")
                return
            check_v11_record(
                archive.get("inventory", {}), f"{label} inventory",
                f"{expected_path}/inventory.json")
            check_v11_record(
                archive.get("checksums", {}), f"{label} checksums",
                f"{expected_path}/SHA256SUMS.txt")
            archive_path = ROOT / expected_path
            sums_path = archive_path / "SHA256SUMS.txt"
            ledger_rows = []
            if sums_path.is_file():
                for line_number, line in enumerate(
                        sums_path.read_text().splitlines(), start=1):
                    match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
                    if match is None:
                        failures.append(
                            f"Figure 3 v11 {label} checksum line {line_number} malformed")
                        continue
                    expected_hash, member_rel = match.groups()
                    member_fragment = Path(member_rel)
                    if member_fragment.is_absolute() or ".." in member_fragment.parts:
                        failures.append(
                            f"Figure 3 v11 {label} checksum has unsafe path {member_rel!r}")
                        continue
                    member_path = archive_path / member_fragment
                    observed_hash = (
                        hashlib.sha256(member_path.read_bytes()).hexdigest()
                        if member_path.is_file() else None)
                    if observed_hash != expected_hash:
                        failures.append(
                            f"Figure 3 v11 {label} member {member_rel!r} changed")
                    ledger_rows.append(member_rel)
            if archive.get("record_count") != len(ledger_rows):
                failures.append(f"Figure 3 v11 {label} archive count changed")

        check_v11_archive(
            v11.get("immediate_predecessor_archive", {}), "predecessor",
            "figures/repository/superseded_figure3_candidates_2026-09-04/"
            "pre_alltarget_enrichment_v11_20260904")
        check_v11_archive(
            v11.get("audit_byproduct_archive", {}), "audit-byproduct",
            "figures/repository/fig3_v11_audit_byproducts_2026-09-04")
        check_v11_archive(
            v11.get("prebuild_script_archive", {}), "prebuild-script",
            "figures/repository/fig3_v11_prebuild_scripts_2026-09-04")
        check_v11_archive(
            v11.get("supplementary_cleanup_archive", {}),
            "supplementary-cleanup",
            "figures/repository/"
            "fig3_v11_supplementary_nonreader_files_2026-09-04")

        expected_v11_panel_order = {
            "a": "target concentration",
            "b": "single-orientation",
            "c": "complete 49 x 279 effect map",
            "d": "literal first-80 enlargement",
            "e": "network-adjusted all-target functional rank enrichment",
            "f": "CCNC/TADA2B/SUPT7L direction-concordant targets",
            "g": "AHR/ARNT direction-concordant targets",
            "h": "Arce edge recurrence",
        }
        if v11.get("panel_order") != expected_v11_panel_order:
            failures.append("Figure 3 v11 panel-order contract changed")
        scope_e = v11.get("panel_e_scope", {})
        expected_scope_e = {
            "eligible_census_targets": 9554,
            "targets_with_no_observed_reversal": 6977,
            "displayed_rows": 19,
            "displayed_kegg_rows": 9,
            "displayed_go_bp_rows": 10,
            "kegg_tested_terms": 168,
            "kegg_repeatable_within_library_fdr_terms": 9,
            "kegg_global_maxT_terms": 0,
            "go_bp_tested_terms": 2425,
            "go_bp_repeatable_within_library_fdr_terms": 161,
            "go_bp_pooled_global_maxT_terms": 12,
            "go_bp_both_bank_global_maxT_terms": 10,
            "global_maxT_family_terms": 2593,
            "permutation_banks": 2,
            "permutations_per_bank": 10000,
        }
        for key, expected in expected_scope_e.items():
            if scope_e.get(key) != expected:
                failures.append(
                    f"Figure 3 v11 promotion panel-e {key} changed")
        if len(scope_e.get("display_order", [])) != 19:
            failures.append("Figure 3 v11 promotion panel-e display order changed")
        for required in (
                "all-target network-adjusted pooled reversal-burden enrichment",
                "do not test lost-versus-gained enrichment", "pathway activity",
                "mechanism", "individual edges", "279 targets in panels c-d",
                "Terms overlap"):
            if required not in v11.get("claim_boundary", ""):
                failures.append(
                    f"Figure 3 v11 promotion claim boundary lacks {required!r}")

        expected_v11_active_paths = {
            "manuscript": "MANUSCRIPT.md",
            "figure_legends": "FIGURE_LEGENDS.md",
            "paper_facts": "PAPER_FACTS.md",
            "supplementary_methods": "SUPPLEMENTARY_METHODS.md",
            "supplementary_results": "SUPPLEMENTARY_RESULTS.md",
        }
        active_v11_text = v11.get("active_text", {})
        if set(active_v11_text) != set(expected_v11_active_paths):
            failures.append("Figure 3 v11 active-text inventory changed")
        superseded_v11_text = (
            set(expected_v11_active_paths)
            if fig2_signed_kegg_supersedes_v11_live_state
            else {"manuscript", "figure_legends", "paper_facts"}
        )
        for key, expected_path in expected_v11_active_paths.items():
            if not (later_promotion_supersedes_v11_live_state and
                    key in superseded_v11_text):
                check_v11_record(
                    active_v11_text.get(key, {}), f"active text {key}", expected_path)
        if not later_promotion_supersedes_v11_live_state:
            check_v11_record(
                v11.get("integration_controls", {}).get("active_claim_guard", {}),
                "active claim guard", "scripts/active_claim_guard.py")
        expected_v11_build_scripts = {
            "docx": "scripts/build_docx.py",
            "pdf": "scripts/build_pdf.py",
            "supplementary": "scripts/build_supplementary.py",
        }
        v11_build_scripts = v11.get("build_scripts", {})
        if set(v11_build_scripts) != set(expected_v11_build_scripts):
            failures.append("Figure 3 v11 build-script inventory changed")
        for key, expected_path in expected_v11_build_scripts.items():
            superseded_build = (
                fig2_signed_kegg_supersedes_v11_live_state or
                (fig4f_supersedes_v11_live_state and key == "docx")
            )
            if not superseded_build:
                check_v11_record(
                    v11_build_scripts.get(key, {}), f"build script {key}", expected_path)

        expected_v11_retention = {
            "permanent_deletion_performed": False,
            "superseded_artifacts_archived_recoverably": True,
            "audit_byproducts_archived_recoverably": True,
            "routine_cleanup_may_not_delete_archived_bundles": True,
            "v10_candidate_bundle_retained": True,
            "v11_candidate_bundle_retained": True,
        }
        if v11.get("retention_policy") != expected_v11_retention:
            failures.append("Figure 3 v11 retention policy changed")

        derived_v11 = v11.get("derived_artifacts", {})
        expected_v11_derived_paths = {
            "pdf": "HumanCD4CoDEGNet_signflip.pdf",
            "docx": "HumanCD4CoDEGNet_signflip.docx",
            "supplementary_pdf": (
                "submission/supplementary/Supplementary_Information.pdf"),
        }
        if v11.get("status") == v11_pending:
            if derived_v11.get("status") != "NOT_REBUILT_FOR_THIS_PROMOTION":
                failures.append("pending Figure 3 v11 derived status changed")
            for key, expected_path in expected_v11_derived_paths.items():
                if derived_v11.get(key) != {"path": expected_path, "sha256": None}:
                    failures.append(
                        f"pending Figure 3 v11 derived {key} record changed")
        elif v11.get("status") == v11_final:
            if derived_v11.get("status") != "POSTBUILD_GUARDED":
                failures.append("final Figure 3 v11 derived status changed")
            for key, expected_path in expected_v11_derived_paths.items():
                superseded_derived = (
                    (fig2_signed_kegg_supersedes_v11_live_state and
                     key in {"docx", "pdf", "supplementary_pdf"}) or
                    (fig4f_supersedes_v11_live_state and key in {"docx", "pdf"})
                )
                if not superseded_derived:
                    check_v11_record(
                        derived_v11.get(key, {}), f"derived {key}", expected_path)
            package_record = derived_v11.get("supplementary_package", {})
            package_root = ROOT / "submission/supplementary"
            package_rows = []
            if package_root.is_dir():
                for path in sorted(
                        item for item in package_root.rglob("*") if item.is_file()):
                    package_rows.append(
                        f"{path.relative_to(package_root)}\0"
                        f"{hashlib.sha256(path.read_bytes()).hexdigest()}\0"
                        f"{path.stat().st_size}\n")
            package_hash = hashlib.sha256("".join(package_rows).encode()).hexdigest()
            if not fig2_signed_kegg_supersedes_v11_live_state and (
                    package_record.get("path") != "submission/supplementary" or
                    package_record.get("tree_sha256") != package_hash or
                    package_record.get("files") != len(package_rows) or
                    package_record.get("files") != 37):
                failures.append("Figure 3 v11 supplementary package seal changed")
            expected_build_logs = {
                "source_guard": (
                    "data/figure_provenance/"
                    "fig3_target_centred_alltarget_enrichment_v11_source_guard.log"),
                "docx_build": (
                    "data/figure_provenance/"
                    "fig3_target_centred_alltarget_enrichment_v11_docx_build.guard.log"),
                "pdf_build": (
                    "data/figure_provenance/"
                    "fig3_target_centred_alltarget_enrichment_v11_pdf_build.guard.log"),
                "supplementary_build": (
                    "data/figure_provenance/"
                    "fig3_target_centred_alltarget_enrichment_v11_supplementary_build.guard.log"),
            }
            build_logs = v11.get("build_logs", {})
            if set(build_logs) != set(expected_build_logs):
                failures.append("Figure 3 v11 build-log inventory changed")
            for key, expected_path in expected_build_logs.items():
                check_v11_record(
                    build_logs.get(key, {}), f"build log {key}", expected_path)
    # The three high-ceiling reframe analyses finished as outcome-blind NO-GO decisions.
    # Component-level engineering or marginal-calibration PASS files are provenance only and
    # must never be able to authorise a scientific run or manuscript claim.
    reframe_status_path = ROOT / "data/reframe_method_status_20260814.json"
    if not reframe_status_path.exists():
        failures.append("binding three-analysis reframe status is missing")
    else:
        reframe_status = json.loads(reframe_status_path.read_text())
        if reframe_status.get("decision") != "ALL_THREE_CURRENT_DESIGNS_NO_GO":
            failures.append("binding reframe status does not stop all three current designs")
        for key in ("scientific_run_authorized", "manuscript_consumption_authorized"):
            if reframe_status.get(key) is not False:
                failures.append(f"binding reframe status has {key} != false")
        expected_decisions = {
            "analysis1": "NO_GO_REPLACEMENT_V1_FAILS_VALIDITY_AND_POWER",
            "analysis2": "NO_GO_TERMINAL_CURRENT_DESIGN",
            "analysis3": "STOP_ANALYSIS3_CURRENT_DESIGN",
        }
        for analysis, expected_decision in expected_decisions.items():
            record = reframe_status.get("analyses", {}).get(analysis, {})
            rel = record.get("canonical_record")
            if record.get("decision") != expected_decision:
                failures.append(
                    f"binding reframe {analysis} decision={record.get('decision')!r}; "
                    f"expected {expected_decision!r}")
                continue
            path = ROOT / rel if rel else None
            if path is None or not path.exists():
                failures.append(f"binding reframe {analysis} canonical record is missing")
                continue
            observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            if observed_hash != record.get("sha256"):
                failures.append(
                    f"binding reframe {analysis} hash {observed_hash} != "
                    f"status {record.get('sha256')}")
            canonical = json.loads(path.read_text())
            if canonical.get("decision") != expected_decision:
                failures.append(
                    f"binding reframe {analysis} canonical decision "
                    f"{canonical.get('decision')!r} != {expected_decision!r}")

        crossed = json.loads((ROOT / "data/analysis3_preflight/"
                              "crossed_auprc_calibration_consolidated_20260814.json").read_text())
        if crossed.get("decision") != \
                "CROSSED_AUPRC_NOT_CERTIFIED_FAMILYWISE_AND_POWER_NO_GO":
            failures.append("Analysis 3 crossed-AUPRC status is not fail-closed")
        a3 = reframe_status.get("analyses", {}).get("analysis3", {})
        hpc_source = ROOT / a3.get("exact_hpc_calibration_source", "")
        if not hpc_source.is_file():
            failures.append("exact Analysis 3 HPC calibration source is not preserved locally")
        else:
            observed_hash = hashlib.sha256(hpc_source.read_bytes()).hexdigest()
            if observed_hash != a3.get("exact_hpc_calibration_source_sha256"):
                failures.append(
                    f"exact Analysis 3 HPC source hash {observed_hash} != "
                    f"status {a3.get('exact_hpc_calibration_source_sha256')}")

    # The biology-first branch is separately locked: the target-trajectory association passed,
    # whereas the orthogonal promoter-accessibility explanation failed. Neither result may be
    # silently promoted to mechanism, and the active Figure 2 must be the audited integration.
    biology_status_path = ROOT / "data/reframe_biology_status_20260814.json"
    figure3_successor_state = "invalid"
    figure4_successor_state = "invalid"
    fig2_shared_binding_not_live = False
    if not biology_status_path.is_file():
        failures.append("binding biology-first reframe status is missing")
    else:
        biology = json.loads(biology_status_path.read_text())
        if biology.get("status") != "BIOLOGY_FIRST_ASSOCIATION_PASS_MECHANISM_UNRESOLVED":
            failures.append("biology-first status no longer preserves the unresolved mechanism")
        figure3_successor_state = _active_successor_registry_state(
            biology, "active_figure_3", check_derived=check_derived,
            failures=failures)
        figure4_successor_state = _active_successor_registry_state(
            biology, "active_figure_4", check_derived=check_derived,
            failures=failures)
        # Figure 2's manifest faithfully records the source/build bytes that were
        # active when it was sealed.  Those shared records become historical once
        # Figure 3 v13 owns the active registry.  Source mode also permits the
        # explicit pending state so the successor builds can create their new
        # bytes before the final two-figure registry rebind.
        fig2_shared_binding_is_historical = figure3_successor_state == "bound"
        fig2_shared_binding_is_transitioning = (
            not check_derived and figure3_successor_state == "pending"
        )
        fig2_shared_binding_not_live = (
            fig2_shared_binding_is_historical or
            fig2_shared_binding_is_transitioning
        )
        target = biology.get("target_trajectory", {})
        movement = biology.get("activation_movement", {})
        timing = biology.get("eight_hour_sign_state", {})
        accessibility = biology.get("promoter_accessibility", {})
        if target.get("decision") != "PASS_BOUNDED_ASSOCIATION":
            failures.append("biology-first target-trajectory decision is not bounded PASS")
        if target.get("provenance") != "prewritten post-primary specificity analysis":
            failures.append("biology-first target-trajectory provenance is not post-primary")
        if movement.get("decision") != "NO_DETECTED_DIFFERENCE_NOT_EQUIVALENCE":
            failures.append("activation-movement result no longer preserves the negative boundary")
        if movement.get("design_provenance") != "post-hoc paired-regulator comparison":
            failures.append("activation-movement result is not explicitly post hoc")
        if timing.get("decision") != "POST_HOC_ORIENTATION_DEPENDENT_SIGN_RETENTION":
            failures.append("8 h sign-state decision is not the bounded post-hoc result")
        if (timing.get("resolved"), timing.get("unresolved"),
                timing.get("paired_regulators")) != (1279, 3098, 75):
            failures.append("8 h sign-state status lost its resolved/unresolved regulator counts")
        if abs(float(timing.get("wilcoxon_two_sided_p", float("nan"))) -
               0.0049970904) > 5e-11:
            failures.append("8 h sign-state status lost the exact paired result")
        if accessibility.get("decision") != "FAIL_PRIMARY_ACCESSIBILITY_COMPATIBILITY":
            failures.append("biology-first promoter-accessibility decision is not failed closed")
        for section, record, pairs in (
                ("target_trajectory", target, (
                    ("result_path", "result_sha256"),
                    ("production_seal_path", "production_seal_sha256"),
                    ("independent_audit_seal_path", "independent_audit_seal_sha256"))),
                ("promoter_accessibility", accessibility, (
                    ("result_path", "result_sha256"),
                    ("independent_audit_path", "independent_audit_sha256"),
                    ("final_seal_path", "final_seal_sha256"))),
                ("activation_movement", movement, (
                    ("result_path", "result_sha256"),
                    ("source_table_path", "source_table_sha256"))),
                ("eight_hour_sign_state", timing, (
                    ("summary_path", "summary_sha256"),
                    ("source_table_path", "source_table_sha256"))),
                ("active_figure_2", biology.get("active_figure_2", {}), (
                    ("pdf_path", "pdf_sha256"),
                    ("png_path", "png_sha256")))):
            for path_key, hash_key in pairs:
                rel = record.get(path_key)
                path = ROOT / rel if rel else None
                if path is None or not path.is_file():
                    failures.append(f"biology-first {section}.{path_key} is missing")
                    continue
                observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
                if observed_hash != record.get(hash_key):
                    failures.append(
                        f"biology-first {section} hash {observed_hash} != "
                        f"status {record.get(hash_key)}")

        fig2 = biology.get("active_figure_2", {})
        for path_key, hash_key in (
                ("builder_path", "builder_sha256"),
                ("promotion_manifest_path", "promotion_manifest_sha256")):
            rel = fig2.get(path_key)
            path = ROOT / rel if rel else None
            if path is None or not path.is_file():
                failures.append(f"biology-first active_figure_2.{path_key} is missing")
                continue
            observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            if observed_hash != fig2.get(hash_key):
                failures.append(
                    f"biology-first active_figure_2 {path_key} hash {observed_hash} != "
                    f"status {fig2.get(hash_key)}")

        promotion_rel = fig2.get("promotion_manifest_path")
        promotion_file = ROOT / promotion_rel if promotion_rel else None
        if promotion_file is not None and promotion_file.is_file():
            promotion = json.loads(promotion_file.read_text())
            fig2_schema = promotion.get("schema_version")
            fig2_signed_kegg_schema = (
                "fig2-signed-kegg-side-by-side-v11-promotion-v1"
            )
            if fig2_schema not in {
                    "fig2-panel-b-main-scale-visibility-promotion-v1",
                    "fig2-panel-b-no-zoom-docx-typography-promotion-v5",
                    fig2_signed_kegg_schema}:
                failures.append("Figure 2 promotion manifest schema is not an accepted audited successor")
            expected_fig2_statuses = {
                "fig2-panel-b-main-scale-visibility-promotion-v1": {
                    "PROMOTED_AFTER_INDEPENDENT_SCIENTIFIC_AND_VISUAL_AUDIT"},
                "fig2-panel-b-no-zoom-docx-typography-promotion-v5": {
                    "PROMOTED_AFTER_INDEPENDENT_NO_ZOOM_DOCX_TYPOGRAPHY_AUDIT"},
                fig2_signed_kegg_schema: {
                    "PROMOTED_SOURCE_POSTBUILD_PENDING",
                    "PROMOTED_AND_POSTBUILD_GUARDED"},
            }
            if promotion.get("status") not in expected_fig2_statuses.get(
                    fig2_schema, set()):
                failures.append("Figure 2 promotion manifest is not in its audited promoted state")
            if check_derived and fig2_schema == fig2_signed_kegg_schema and \
                    promotion.get("status") != "PROMOTED_AND_POSTBUILD_GUARDED":
                failures.append("Figure 2 signed-KEGG promotion is not post-build complete")
            if fig2_schema == fig2_signed_kegg_schema:
                expected_panel_order = {
                    "a": "observed versus out-degree-matched expected reversal burden",
                    "b": "all 49 enriched regulators with rev/test, rate, fold and L/G",
                    "c": "reversal rate versus out-degree",
                    "d": "full-census conventional reversal-target KEGG ORA with null disclosure",
                    "e": "all 39 reversal targets of six one-direction panel-b regulators",
                    "f": "Rest-sign retention at Rest, 8 h and 48 h with 8 h classification",
                    "g": "0-72 h donor AUC trajectory for frozen unanimous target groups",
                    "h": "absolute 0-48 h responses showing that both target groups rise",
                    "i": "903-target reversal-orientation association against two nulls",
                }
                if promotion.get("panel_order") != expected_panel_order:
                    failures.append(
                        "Figure 2 promotion manifest a--i panel semantics changed")
                if promotion_rel != FIG2_SIGNED_KEGG_V11_PROMOTION:
                    failures.append(
                        "active Figure 2 registry points to the wrong signed-KEGG promotion")
                expected_registry = {
                    "builder_path": FIG2_SIGNED_KEGG_V11_EXPECTED["renderer"][0],
                    "builder_sha256": FIG2_SIGNED_KEGG_V11_EXPECTED["renderer"][1],
                    "pdf_path": "figures/Fig2.pdf",
                    "pdf_sha256": FIG2_SIGNED_KEGG_V11_EXPECTED["pdf"][1],
                    "png_path": "figures/Fig2.png",
                    "png_sha256": FIG2_SIGNED_KEGG_V11_EXPECTED["png"][1],
                    "page_format": "A4_PORTRAIT_FULL_PAGE",
                }
                for key, expected in expected_registry.items():
                    if fig2.get(key) != expected:
                        failures.append(
                            f"active Figure 2 registry {key}={fig2.get(key)!r}; "
                            f"expected {expected!r}")
                lower_layout = fig2.get("lower_layout", "")
                main_contract = fig2.get("main_figure_contract", "")
                regulator_contract = fig2.get("regulator_panel_b", "")
                kegg_contract = fig2.get("panel_d_kegg", "")
                if not lower_layout.startswith("Nine panels.") or \
                        "d the full-census conventional KEGG sensitivity" not in lower_layout or \
                        "i the 903-target orientation analysis" not in lower_layout:
                    failures.append("active Figure 2 registry lower-layout contract is stale")
                if "Figure 2 contains panels a-i" not in main_contract or \
                        "global p=0.184; 0/113 maxT-supported" not in main_contract or \
                        "Panels e-i retain their prior analyses" not in main_contract:
                    failures.append("active Figure 2 registry main contract is stale")
                if "All 49 reversal-enriched regulators" not in regulator_contract or \
                        "rev/test, rate, fold and L/G" not in regulator_contract or \
                        "duplicate inset and direction donut are absent" not in regulator_contract:
                    failures.append("active Figure 2 registry panel-b contract is stale")
                if "2,577 reversal targets against all 9,554 eligible targets" not in \
                        kegg_contract or "14/113 KEGG Legacy rows" not in kegg_contract or \
                        "0/113 terms passed maxT" not in kegg_contract:
                    failures.append("active Figure 2 registry panel-d contract is stale")
                if fig2.get("typography") != {
                        "delivered_font_range_pt": [5.0, 7.0],
                        "panel_letter_and_title_pt": 7.0,
                        "axis_label_pt": 6.0,
                        "other_text_pt": 5.0,
                        "panel_b_metadata_pt": 5.0,
                }:
                    failures.append("active Figure 2 registry typography contract changed")
            elif list(promotion.get("panel_order", {})) != list("abcdefgh"):
                failures.append(
                    "Figure 2 promotion manifest panel order is not the active contract")

            def _check_fig2_bound_file(rel: str | None, expected: str | None, label: str) -> None:
                path = ROOT / rel if rel else None
                if path is None or not path.is_file():
                    failures.append(f"Figure 2 promotion binding missing {label}: {rel!r}")
                    return
                observed = hashlib.sha256(path.read_bytes()).hexdigest()
                if observed != expected:
                    failures.append(
                        f"Figure 2 promotion binding {label} hash {observed} != {expected}")

            def _check_fig2_shared_record(
                    rel: str | None, expected: str | None, label: str) -> None:
                if fig2_shared_binding_not_live:
                    _historical_binding_record(rel, expected, label, failures)
                else:
                    _check_fig2_bound_file(rel, expected, label)

            canonical_outputs = promotion.get("canonical_outputs", {})
            for label in ("pdf", "png"):
                record = canonical_outputs.get(label, {})
                _check_fig2_bound_file(record.get("path"), record.get("sha256"), label)
            if fig2_schema == fig2_signed_kegg_schema:
                for label, expected_path, expected_hash in (
                        ("pdf", "figures/Fig2.pdf",
                         FIG2_SIGNED_KEGG_V11_EXPECTED["pdf"][1]),
                        ("png", "figures/Fig2.png",
                         FIG2_SIGNED_KEGG_V11_EXPECTED["png"][1])):
                    item = canonical_outputs.get(label, {})
                    if item.get("path") != expected_path or \
                            item.get("sha256") != expected_hash:
                        failures.append(
                            f"Figure 2 canonical {label} is not the exact selected v11 byte copy")
                if canonical_outputs.get("promoted_by_byte_copy") is not True or \
                        canonical_outputs.get(
                            "candidate_and_canonical_bytes_identical") is not True:
                    failures.append("Figure 2 v11 lost its canonical byte-copy contract")
            supersedes = promotion.get("supersedes", {})
            _check_fig2_bound_file(
                supersedes.get("path"), supersedes.get("sha256"),
                "superseded panel-b inset visibility promotion manifest")
            candidate_manifest = promotion.get("candidate_manifest", {})
            _check_fig2_bound_file(
                candidate_manifest.get("path"), candidate_manifest.get("sha256"),
                "Figure 2 selected candidate manifest")
            if fig2_schema == fig2_signed_kegg_schema:
                candidate_records = promotion.get("candidate", {})
                candidate_payload = {}
                expected_candidate_keys = {
                    "directory", "manifest", "pdf", "png", "provenance",
                    "render_guard", "determinism_guard", "determinism_receipt",
                    "panel_b_source", "panel_d_source", "renderer", "audit_script",
                }
                if set(candidate_records) != expected_candidate_keys or \
                        candidate_records.get("directory") != \
                        FIG2_SIGNED_KEGG_V11_CANDIDATE:
                    failures.append("Figure 2 v11 candidate record inventory changed")
                for key, (expected_path, expected_hash) in \
                        FIG2_SIGNED_KEGG_V11_EXPECTED.items():
                    item = candidate_records.get(key, {})
                    if item.get("path") != expected_path or \
                            item.get("sha256") != expected_hash:
                        failures.append(
                            f"Figure 2 v11 candidate {key} binding changed")
                    _check_fig2_bound_file(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 v11 candidate {key}")
                if candidate_manifest != candidate_records.get("manifest"):
                    failures.append("Figure 2 v11 duplicated candidate-manifest bindings differ")

                candidate_manifest_path = ROOT / \
                    FIG2_SIGNED_KEGG_V11_EXPECTED["manifest"][0]
                if candidate_manifest_path.is_file():
                    candidate_payload = json.loads(candidate_manifest_path.read_text())
                    if candidate_payload.get("schema_version") != \
                            "fig2-signed-kegg-side-by-side-candidate-v11-final-review" or \
                            candidate_payload.get("status") != \
                            "NONCANONICAL_DIAGNOSTIC_CANDIDATE_FOR_REVIEW":
                        failures.append("Figure 2 v11 candidate identity changed")
                    if candidate_payload.get("canonical_figure_2", {}).get(
                            "written") is not False or \
                            candidate_payload.get("canonical_figure_2", {}).get("before") != \
                            candidate_payload.get("canonical_figure_2", {}).get("after"):
                        failures.append("Figure 2 v11 candidate no-write custody changed")
                    tags = []
                    stack = [candidate_payload]
                    while stack:
                        value = stack.pop()
                        if isinstance(value, dict):
                            if isinstance(value.get("panel_tags"), dict):
                                tags.append(value["panel_tags"])
                            stack.extend(value.values())
                        elif isinstance(value, list):
                            stack.extend(value)
                    if {letter: 1 for letter in "abcdefghi"} not in tags:
                        failures.append("Figure 2 v11 candidate lacks exactly one a--i tag set")
                    candidate_audit = candidate_payload.get("audit", {})
                    inherited = candidate_audit.get("predecessor_v9_audit", {})
                    metadata_audit = inherited.get("panel_b_metadata", {})
                    label_audit = inherited.get("panel_b_regulator_labels", {})
                    panel_d_audit = inherited.get("panel_d", {})
                    expected_candidate_typography = {
                        "axis_title_artists": 10,
                        "delivered_values_pt": [5.0, 6.0, 7.0],
                        "panel_title_artists": 10,
                        "visible_text_artists": 660,
                    }
                    shared_legend = candidate_audit.get(
                        "shared_functional_legend", {})
                    panel_c_annotation = candidate_audit.get(
                        "panel_c_annotation", {})
                    if candidate_audit.get("delta_only") is not True or \
                            candidate_audit.get(
                                "scientific_values_and_panels_preserved") is not True or \
                            candidate_audit.get("typography") != \
                            expected_candidate_typography or \
                            metadata_audit.get("headers") != [
                                "rev/test", "rate", "fold", "L/G"] or \
                            metadata_audit.get("delivered_font_pt") != 5.0 or \
                            float(metadata_audit.get(
                                "minimum_horizontal_gap_pt", -1)) <= 0 or \
                            label_audit.get("labels") != 49 or \
                            label_audit.get("fontweight") != "normal" or \
                            label_audit.get("lost_gained_colour_collisions") != [] or \
                            len(label_audit.get("functional_colours", {})) != 49 or \
                            float(panel_d_audit.get("width_gain_fraction", -1)) <= 0 or \
                            candidate_audit.get("inherited_geometry", {}).get(
                                "functional_key_in_artwork") is not True or \
                            shared_legend != {
                                "present": True,
                                "applies_to_panels": ["b", "c"],
                                "categories": 10,
                                "inside_panel_b": True,
                                "bar_overlap_count": 0,
                                "delivered_font_pt": 5.0,
                            } or \
                            panel_c_annotation.get("arrow_visible") is not False or \
                            panel_c_annotation.get(
                                "highlighted_label_overlap_count") != 0 or \
                            float(candidate_audit.get(
                                "panel_b_metadata_minimum_internal_gap_pt", -1)) < 0.35 or \
                            float(candidate_audit.get(
                                "panel_d_width_gain_from_v10_fraction", -1)) <= 0 or \
                            min(float(value) for value in candidate_audit.get(
                                "lower_panel_tag_to_title_gaps_pt", {}).values()) < 2.0 or \
                            any(int(value) for value in candidate_audit.get(
                                "lower_panel_tag_tick_overlap_counts", {}).values()) or \
                            min(float(candidate_audit.get(key, -1)) for key in (
                                "panel_a_xlabel_to_c_title_gap_pt",
                                "panel_h_to_i_ylabel_gap_pt",
                                "panel_b_title_to_subtitle_gap_pt",
                                "kegg_header_minimum_gap_pt",
                                "kegg_count_to_q_minimum_gap_pt",
                                "kegg_q_to_marker_minimum_gap_pt",
                                "panel_b_to_d_pathway_minimum_gap_pt",
                                "panel_c_xlabel_to_e_title_gap_pt",
                                "panel_g_title_to_h_ylabel_gap_pt")) < 2.0:
                        failures.append("Figure 2 v11 candidate typography/geometry audit changed")
                    receipt_path = ROOT / \
                        FIG2_SIGNED_KEGG_V11_EXPECTED["determinism_receipt"][0]
                    receipt = json.loads(receipt_path.read_text()) \
                        if receipt_path.is_file() else {}
                    expected_receipt = {
                        "all_named_cross_panel_clearances_at_least_2_pt": True,
                        "bundle_source_hashes": {
                            "panel_b": FIG2_SIGNED_KEGG_V11_EXPECTED[
                                "panel_b_source"][1],
                            "panel_d": FIG2_SIGNED_KEGG_V11_EXPECTED[
                                "panel_d_source"][1],
                        },
                        "candidate_hashes": {
                            "pdf": FIG2_SIGNED_KEGG_V11_EXPECTED["pdf"][1],
                            "png": FIG2_SIGNED_KEGG_V11_EXPECTED["png"][1],
                        },
                        "canonical_and_predecessors_unchanged": True,
                        "delivered_typography_pt": [5.0, 6.0, 7.0],
                        "manifest_schema_and_status_valid": True,
                        "page_text_overflow_count": 0,
                        "panel_b_metadata_delivered_pt": 5.0,
                        "panel_b_metadata_minimum_internal_gap_pt": 0.3540983999999753,
                        "panel_b_rows": 49,
                        "panel_c_annotation_leader_removed": True,
                        "panel_d_rows": 14,
                        "panel_d_width_gain_from_v10_fraction": 0.06870229007633588,
                        "rerenders_match_saved_candidate": True,
                        "saved_pdf_required_text_present": True,
                        "schema_version": "fig2-v11-determinism-receipt-v1",
                        "shared_b_c_functional_legend_categories": 10,
                        "shared_b_c_functional_legend_present": True,
                        "lower_panel_tag_tick_overlap_counts": {
                            "f": 0, "g": 0, "h": 0, "i": 0,
                        },
                        "two_independent_rerenders_identical": True,
                    }
                    if receipt != expected_receipt:
                        failures.append("Figure 2 v11 determinism receipt changed")
            for label, expected in promotion.get("renderer_closure", {}).items():
                _check_fig2_bound_file(label, expected, label)
            for label, expected in promotion.get("source_closure", {}).items():
                _check_fig2_bound_file(label, expected, label)
            if fig2_schema == fig2_signed_kegg_schema:
                expected_renderer_closure = candidate_payload.get(
                    "renderer_closure", {})
                if promotion.get("renderer_closure") != expected_renderer_closure or \
                        set(expected_renderer_closure) == set():
                    failures.append("Figure 2 v11 renderer closure changed")
                expected_source_closure = {
                    FIG2_SIGNED_KEGG_V11_EXPECTED["panel_b_source"][0]:
                        FIG2_SIGNED_KEGG_V11_EXPECTED["panel_b_source"][1],
                    FIG2_SIGNED_KEGG_V11_EXPECTED["panel_d_source"][0]:
                        FIG2_SIGNED_KEGG_V11_EXPECTED["panel_d_source"][1],
                }
                expected_source_closure.update({
                    path: item.get("sha256")
                    for path, item in candidate_payload.get(
                        "scientific_source_closure", {}).items()
                })
                if promotion.get("source_closure") != expected_source_closure:
                    failures.append("Figure 2 v11 scientific/source-table closure changed")
                expected_active_text = {
                    "MANUSCRIPT.md", "FIGURE_LEGENDS.md", "PAPER_FACTS.md",
                    "SUPPLEMENTARY_METHODS.md", "SUPPLEMENTARY_RESULTS.md",
                }
                active_text = promotion.get("active_text", {})
                if set(active_text) != expected_active_text:
                    failures.append(
                        "Figure 2 signed-KEGG active-text inventory changed")
                for label, item in active_text.items():
                    if item.get("path") != label:
                        failures.append(
                            f"Figure 2 signed-KEGG active text {label} path changed")
                    _check_fig2_shared_record(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 signed-KEGG active text {label}")
                expected_guard_paths = {
                    "active_guard": "scripts/active_claim_guard.py",
                    "fig3_guard_helper": "scripts/fig3_compact_zoom_guard_v12.py",
                    "fig4_guard_helper":
                        "scripts/fig4_horizontal_transposed_cyan_guard_v2.py",
                    "promotion_script": "scripts/promote_fig2_signed_kegg_v11.py",
                }
                for label, expected_path in expected_guard_paths.items():
                    item = promotion.get(label, {})
                    if item.get("path") != expected_path:
                        failures.append(
                            f"Figure 2 signed-KEGG {label} path changed")
                    checker = (
                        _check_fig2_bound_file
                        if label == "promotion_script"
                        else _check_fig2_shared_record
                    )
                    checker(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 signed-KEGG {label}")

                amendment = promotion.get("guard_normalization_reseal", {})
                if promotion.get("presentation_revision") == \
                        "overlap-clearance-and-shared-bc-key-v1":
                    if amendment != {}:
                        failures.append(
                            "Figure 2 overlap-clearance revision carries a stale reseal amendment")
                elif not isinstance(amendment, dict) or amendment.get(
                        "schema_version") != \
                        FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_SCHEMA:
                    failures.append(
                        "Figure 2 v11 lacks the audited active-guard normalization reseal")
                else:
                    archive_binding = amendment.get("pre_reseal_archive", {})
                    expected_archive_paths = {
                        **FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS,
                        "inventory": {
                            "path": (
                                f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
                                "inventory.json"
                            ),
                        },
                        "readme": {
                            "path": (
                                f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
                                "README.md"
                            ),
                        },
                        "sha256sums": {
                            "path": (
                                f"{FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE}/"
                                "SHA256SUMS.txt"
                            ),
                        },
                    }
                    if not isinstance(archive_binding, dict) or \
                            archive_binding.get("directory") != \
                            FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_ARCHIVE or \
                            archive_binding.get("recoverable") is not True or \
                            archive_binding.get(
                                "permanent_deletion_performed") is not False:
                        failures.append(
                            "Figure 2 v11 guard-reseal archive contract changed")
                    else:
                        for label, expected in expected_archive_paths.items():
                            item = archive_binding.get(label, {})
                            if item.get("path") != expected["path"] or any(
                                    item.get(key) != value
                                    for key, value in expected.items()
                                    if key != "path"):
                                failures.append(
                                    "Figure 2 v11 guard-reseal archive binding changed: "
                                    f"{label}")
                            _check_fig2_bound_file(
                                item.get("path"), item.get("sha256"),
                                f"Figure 2 v11 guard-reseal archive {label}")

                    if amendment.get("active_guard_before") != \
                            FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                                "active_guard_before"] or \
                            amendment.get("active_guard_after") != \
                            promotion.get("active_guard"):
                        failures.append(
                            "Figure 2 v11 guard-reseal before/after guard binding changed")
                    amendment_script = amendment.get("amendment_script", {})
                    if amendment_script.get("path") != \
                            FIG2_SIGNED_KEGG_V11_GUARD_RESEAL_SCRIPT:
                        failures.append(
                            "Figure 2 v11 guard-reseal transaction path changed")
                    _check_fig2_bound_file(
                        amendment_script.get("path"), amendment_script.get("sha256"),
                        "Figure 2 v11 guard-reseal transaction")

                    expected_normalization = {
                        "unicode_compatibility_normalization": "NFKC",
                        "collapse_whitespace": True,
                        "normalize_spaces_around_equals": True,
                        "remove_whitespace_before_punctuation": [";", ",", ":"],
                        "docx_text_namespaces": ["w:t", "m:t"],
                        "required_numbers_and_conditions_unchanged": True,
                    }
                    if amendment.get("normalization_scope") != expected_normalization or \
                            _normalise_text(
                                "global 𝑝 = 0.184 ; 0/113 maxT-supported") != \
                            "global p=0.184; 0/113 maxt-supported":
                        failures.append(
                            "Figure 2 v11 guard extraction normalization contract changed")

                    reseal_complete = promotion.get("status") == \
                        "PROMOTED_AND_POSTBUILD_GUARDED"
                    expected_change_boundary = {
                        "manuscript_or_legend_changed": False,
                        "canonical_figure_2_changed": False,
                        "scientific_numbers_or_conditions_changed": False,
                        "build_scripts_changed": False,
                        "derived_rebuild_required": True,
                        "derived_rebuild_performed": reseal_complete,
                        "original_build_integration_preserved_in_amendment": True,
                    }
                    if amendment.get("change_boundary") != expected_change_boundary:
                        failures.append(
                            "Figure 2 v11 guard-reseal change boundary changed")
                    expected_fresh_logs = {
                        "source_guard": (
                            "data/figure_provenance/"
                            "fig2_signed_kegg_v11_guard_normalization_reseal/"
                            "source_guard.guard.log"
                        ),
                        "docx_build": (
                            "data/figure_provenance/"
                            "fig2_signed_kegg_v11_guard_normalization_reseal/"
                            "docx_build.guard.log"
                        ),
                        "pdf_build": (
                            "data/figure_provenance/"
                            "fig2_signed_kegg_v11_guard_normalization_reseal/"
                            "pdf_build.guard.log"
                        ),
                        "supplementary_build": (
                            "data/figure_provenance/"
                            "fig2_signed_kegg_v11_guard_normalization_reseal/"
                            "supplementary_build.guard.log"
                        ),
                        "postbuild_guard": (
                            "data/figure_provenance/"
                            "fig2_signed_kegg_v11_guard_normalization_reseal/"
                            "postbuild_guard.guard.log"
                        ),
                    }
                    if amendment.get("fresh_guard_log_paths") != expected_fresh_logs:
                        failures.append(
                            "Figure 2 v11 guard-reseal fresh-log paths changed")
                    if amendment.get("build_sequence") != [
                            "source_guard", "docx_build", "pdf_build",
                            "supplementary_build"] or \
                            amendment.get("guard_rebound_utc") != \
                            promotion.get("guard_rebound_utc"):
                        failures.append(
                            "Figure 2 v11 guard-reseal build order/timestamp changed")

                    prior_manifest_path = ROOT / \
                        FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                            "prior_promotion_manifest"]["path"]
                    if prior_manifest_path.is_file():
                        prior_promotion = json.loads(prior_manifest_path.read_text())
                        allowed_added_keys = {
                            "guard_normalization_reseal", "guard_rebound_utc"}
                        if reseal_complete:
                            allowed_added_keys.add("guard_resealed_utc")
                        if set(promotion) != set(prior_promotion) | allowed_added_keys:
                            failures.append(
                                "Figure 2 v11 guard reseal changed the promotion key set")
                        for key, value in prior_promotion.items():
                            if key not in {
                                    "active_guard", "status", "sealed_utc",
                                    "build_integration"} and \
                                    promotion.get(key) != value:
                                failures.append(
                                    "Figure 2 v11 guard reseal changed protected promotion "
                                    f"content: {key}")
                        if promotion.get("active_guard", {}).get("path") != \
                                "scripts/active_claim_guard.py" or \
                                not isinstance(promotion.get("guard_rebound_utc"), str):
                            failures.append(
                                "Figure 2 v11 guard reseal lacks its live guard/timestamp")
                        prior_outputs = prior_promotion.get(
                            "build_integration", {}).get("outputs", {})
                        archived_outputs = amendment.get(
                            "derived_outputs_before_rebuild_archive", {})
                        expected_archived_outputs = {
                            "docx": FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                                "prior_docx"],
                            "pdf": FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                                "prior_pdf"],
                            "supplementary_pdf":
                                FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                                    "prior_supplementary_pdf"],
                        }
                        same_prior_output_bytes = all(
                            archived_outputs.get(key, {}).get("sha256") ==
                            prior_outputs.get(key, {}).get("sha256") and
                            archived_outputs.get(key, {}).get("bytes") ==
                            prior_outputs.get(key, {}).get("bytes")
                            for key in expected_archived_outputs
                        )
                        if amendment.get("prior_build_integration") != \
                                prior_promotion.get("build_integration") or \
                                archived_outputs != expected_archived_outputs or \
                                not same_prior_output_bytes:
                            failures.append(
                                "Figure 2 v11 guard reseal lost the prior build closure")
                        current_build = promotion.get("build_integration", {})
                        if reseal_complete:
                            if promotion.get("guard_resealed_utc") != \
                                    promotion.get("sealed_utc") or \
                                    amendment.get("guard_resealed_utc") != \
                                    promotion.get("sealed_utc") or \
                                    not isinstance(promotion.get("sealed_utc"), str) or \
                                    amendment.get("derived_outputs_after_rebuild") != \
                                    current_build.get("outputs"):
                                failures.append(
                                    "Figure 2 v11 guard reseal lacks its fresh sealed outputs")
                            observed_build_log_paths = {
                                key: current_build.get("logs", {}).get(key, {}).get("path")
                                for key in (
                                    "source_guard", "docx_build", "pdf_build",
                                    "supplementary_build")
                            }
                            if observed_build_log_paths != {
                                    key: expected_fresh_logs[key]
                                    for key in observed_build_log_paths}:
                                failures.append(
                                    "Figure 2 v11 guard reseal is not bound to the fresh logs")
                        elif promotion.get("sealed_utc") is not None or \
                                current_build != {
                                    "status": "PENDING",
                                    "sequence": [
                                        "source_guard", "docx_build", "pdf_build",
                                        "supplementary_build"],
                                    "logs": {},
                                    "outputs": {},
                                } or "guard_resealed_utc" in promotion or \
                                "derived_outputs_after_rebuild" in amendment:
                            failures.append(
                                "Figure 2 v11 guard rebind is not cleanly post-build pending")

                    prior_status_path = ROOT / \
                        FIG2_SIGNED_KEGG_V11_PRE_RESEAL_RECORDS[
                            "prior_status_registry"]["path"]
                    if prior_status_path.is_file():
                        prior_status = json.loads(prior_status_path.read_text())
                        expected_active_fig2 = prior_status.get("active_figure_2", {})
                        expected_active_fig2 = json.loads(json.dumps(
                            expected_active_fig2))
                        expected_active_fig2["promotion_manifest_sha256"] = \
                            hashlib.sha256(promotion_file.read_bytes()).hexdigest()
                        if fig2 != expected_active_fig2:
                            failures.append(
                                "Figure 2 v11 active registry changed beyond the resealed "
                                "promotion pointer")
                expected_build_scripts = {
                    "docx": "scripts/build_docx.py",
                    "pdf": "scripts/build_pdf.py",
                    "supplementary": "scripts/build_supplementary.py",
                }
                build_scripts = promotion.get("build_scripts", {})
                if set(build_scripts) != set(expected_build_scripts):
                    failures.append(
                        "Figure 2 signed-KEGG build-script inventory changed")
                for label, expected_path in expected_build_scripts.items():
                    item = build_scripts.get(label, {})
                    if item.get("path") != expected_path:
                        failures.append(
                            f"Figure 2 signed-KEGG {label} build-script path changed")
                    _check_fig2_shared_record(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 signed-KEGG {label} build script")

            layout = promotion.get("layout_preservation", {})
            if fig2_schema == fig2_signed_kegg_schema:
                new_panel_range = layout.get("new_panel_delivered_font_range_pt", [])
                typography_evidence = layout.get("typography_evidence", {})
                bad_layout = (
                    layout.get("panels_retained") != list("abcdefghi") or
                    "former panels d-h relettered e-i" not in
                    layout.get("outside_panels_b_and_d", "") or
                    layout.get("panel_letter_and_title_delivered_pt") != 7.0 or
                    layout.get("axis_label_delivered_pt") != 6.0 or
                    layout.get("other_text_delivered_pt") != 5.0 or
                    layout.get("delivered_pdf_text_sizes_pt") != [5.0, 6.0, 7.0] or
                    new_panel_range != [5.0, 7.0] or
                    layout.get("panel_b_metadata_delivered_pt") != 5.0 or
                    layout.get("panel_b_metadata_compacted") is not True or
                    layout.get("panel_d_widened_for_readability") is not True or
                    layout.get("panel_d_width_gain_from_v10_fraction") !=
                    0.06870229007633588 or
                    layout.get("shared_panel_b_c_functional_colour_legend") is not True or
                    layout.get("panel_ac_axes_narrowed_for_clearance") is not True or
                    layout.get("lower_hi_ylabels_cleared") is not True or
                    layout.get("harmful_cross_panel_collision_count") != 0 or
                    layout.get("independent_deterministic_rerender_passed") is not True or
                    typography_evidence.get("candidate_manifest_typography") != {
                        "axis_title_artists": 10,
                        "delivered_values_pt": [5.0, 6.0, 7.0],
                        "panel_title_artists": 10,
                        "visible_text_artists": 660,
                    } or
                    typography_evidence.get(
                        "saved_candidate_byte_matched_by_two_independent_rerenders") is not True or
                    typography_evidence.get("determinism_receipt") !=
                    promotion.get("candidate", {}).get("determinism_receipt") or
                    layout.get("off_page_text") != [])
            elif fig2_schema == "fig2-panel-b-main-scale-visibility-promotion-v1":
                bad_layout = (
                    layout.get("panels_a_and_c_to_h") !=
                    "pixel-identical to the superseded canonical Figure 2" or
                    "0–300 linear" not in layout.get("panel_b_main", "") or
                    "22 rows" not in layout.get("panel_b_inset", "") or
                    "0–15 linear" not in layout.get("panel_b_inset", "") or
                    "visibility only" not in layout.get("panel_b_outline", "") or
                    layout.get("pixel_differences") !=
                    "confined to panel b main-scale origin/mark strip" or
                    "ten one-edge" not in layout.get("main_scale_change", "") or
                    "0–300" not in layout.get("main_scale_axis", "") or
                    "49 bold category" not in layout.get("label_region", "") or
                    layout.get("independent_audit_passed") is not True)
            else:
                bad_layout = (
                    layout.get("panels_a_and_c_to_g") !=
                    "data artists and axes geometry unchanged; typography sizes only" or
                    layout.get("panel_h_geometry") !=
                    "left boundary moved right 0.012 page fraction; right boundary and data unchanged" or
                    "0–300" not in layout.get("panel_b_main", "") or
                    layout.get("panel_b_inset_removed") is not True or
                    "duplicate inset" not in layout.get("removed_element", "") or
                    layout.get("pixel_differences") !=
                    "duplicate inset removal, bounded text-glyph size changes and panel h label-clearance shift" or
                    "ten one-edge" not in layout.get("main_scale_change", "") or
                    "0–300" not in layout.get("main_scale_axis", "") or
                    "49 bold category" not in layout.get("label_region", "") or
                    layout.get("independent_audit_passed") is not True or
                    layout.get("typography_audit_passed") is not True or
                    layout.get("visible_text_artists") != 324 or
                    layout.get("source_font_range_pt") != [7.09, 9.92] or
                    layout.get("delivered_width_in") != 5334000 / 914400 or
                    not (5.0 <= layout.get("delivered_font_min_pt", 0.0) <=
                         layout.get("delivered_font_max_pt", 99.0) <= 7.0) or
                    layout.get("off_page_text") != [] or
                    layout.get("data_and_layout_invariants_passed") is not True or
                    layout.get("panel_h_ylabel_visibility_audit_passed") is not True)
            if bad_layout:
                failures.append("Figure 2 successor lost its audited panel-b visibility contract")

            panel_b = promotion.get("panel_b_display", {})
            expected_panel_b = {
                "regulators": 49,
                "both_orientations": 43,
                "lost_only_regulators": 3,
                "gained_only_regulators": 3,
                "total_edges": 2075,
                "lost_edges": 1154,
                "gained_edges": 921,
                "metadata_fields": "rev/test; rate; fold; L/G",
            }
            for key, expected in expected_panel_b.items():
                if panel_b.get(key) != expected:
                    failures.append(
                        f"Figure 2b display {key}={panel_b.get(key)!r}; expected {expected!r}")
            if fig2_schema == fig2_signed_kegg_schema:
                bad_panel_b = (
                    panel_b.get("tested_edges") != 14434 or
                    panel_b.get("metadata_visible_for_all_rows") is not True or
                    panel_b.get("metadata_delivered_pt") != 5.0 or
                    panel_b.get("regulator_label_weight") != "regular" or
                    panel_b.get("functional_label_colours_retained") is not True or
                    panel_b.get("functional_colour_key_owner") !=
                    "shared in-artwork legend for panels b and c" or
                    panel_b.get("functional_colour_key_categories") != 10 or
                    panel_b.get("duplicate_low_count_inset_displayed") is not False or
                    panel_b.get("direction_summary_donut_displayed") is not False or
                    panel_b.get("one_edge_visibility_ticks_retained") is not True or
                    panel_b.get("bar_endpoints_unchanged") is not True or
                    "linear 0–300" not in panel_b.get("main_scale", ""))
                b_source_path = ROOT / \
                    FIG2_SIGNED_KEGG_V11_EXPECTED["panel_b_source"][0]
                b_rows = []
                b_fields = []
                if b_source_path.is_file():
                    with b_source_path.open(newline="") as handle:
                        reader = csv.DictReader(handle, delimiter="\t")
                        b_fields = list(reader.fieldnames or [])
                        b_rows = list(reader)
                expected_b_fields = [
                    "display_rank", "gene", "n_inv", "n_both", "rate",
                    "enrichment_vs_decile", "p_two_sided", "q_two_sided",
                    "lost_on_activation", "gained_on_activation",
                ]
                if b_fields != expected_b_fields or len(b_rows) != 49 or \
                        [int(row["display_rank"]) for row in b_rows] != \
                        list(range(1, 50)) or \
                        panel_b.get("row_order") != [row["gene"] for row in b_rows] or \
                        sum(int(row["n_inv"]) for row in b_rows) != 2075 or \
                        sum(int(row["n_both"]) for row in b_rows) != 14434 or \
                        sum(int(row["lost_on_activation"]) for row in b_rows) != 1154 or \
                        sum(int(row["gained_on_activation"]) for row in b_rows) != 921 or \
                        any(int(row["n_inv"]) !=
                            int(row["lost_on_activation"]) +
                            int(row["gained_on_activation"]) for row in b_rows):
                    failures.append("Figure 2b exact 49-row source/order/count closure changed")
                inherited_label_audit = candidate_payload.get("audit", {}).get(
                    "predecessor_v9_audit", {}).get(
                        "panel_b_regulator_labels", {})
                if panel_b.get("functional_label_colours") != \
                        inherited_label_audit.get("functional_colours") or \
                        panel_b.get("lost_gained_colour_collisions") != [] or \
                        panel_b.get("fold_definition") != \
                        "observed reversal rate divided by out-degree-decile expectation":
                    failures.append("Figure 2b label-colour/fold-definition contract changed")
            elif fig2_schema == "fig2-panel-b-main-scale-visibility-promotion-v1":
                expected_inset = {
                    "inset_rows": 22, "inset_total_range": [3, 10],
                    "segment_outline_colour": "#3f454b",
                    "segment_outline_width_pt": 0.55,
                    "outline_role": "visibility aid only"}
                bad_panel_b = any(panel_b.get(key) != expected
                                  for key, expected in expected_inset.items()) or (
                    "linear 0–300" not in panel_b.get("main_scale", "") or
                    "linear 0–15" not in panel_b.get("inset_scale", "") or
                    "at most 10" not in panel_b.get("inset_rule", "") or
                    "one lost edge" not in panel_b.get("rasal3_fill", ""))
            else:
                bad_panel_b = (
                    "linear 0–300" not in panel_b.get("main_scale", "") or
                    panel_b.get("duplicate_low_count_inset_removed") is not True or
                    panel_b.get("remaining_child_axis") != "direction-summary donut" or
                    panel_b.get("one_edge_visibility_ticks_retained") is not True or
                    panel_b.get("bar_endpoints_unchanged") is not True or
                    panel_b.get("all_other_panel_data_artists_and_layout_unchanged") is not True or
                    panel_b.get("small_count_readout") !=
                    "exact adjacent lost/gained counts; no duplicate inset")
            if bad_panel_b:
                failures.append("Figure 2b scale, no-zoom or RASAL3 visibility contract changed")
            one_edge = panel_b.get("main_scale_one_edge_ticks", {})
            if one_edge.get("main_scale") != "linear 0–300 count scale, unchanged" or \
                    one_edge.get("bar_endpoints_unchanged") is not True or \
                    one_edge.get("left_spine_position_unchanged") is not True or \
                    one_edge.get("one_edge_segment_count") != 10 or \
                    one_edge.get("one_edge_rule") != \
                    "slim colour tick at the exact segment midpoint" or \
                    one_edge.get("tick_role") != \
                    "nonzero visibility only; bar length remains the count encoding" or \
                    one_edge.get("rasal3") != {
                        "lost": 1, "gained": 3, "lost_tick_x": 0.5}:
                failures.append("Figure 2b main-scale one-edge visibility contract changed")
            expected_one_edge = {
                ("ACTR8", "lost", 0.5),
                ("ATP5ME", "gained", 11.5),
                ("GSS", "lost", 0.5),
                ("ATP5PO", "gained", 9.5),
                ("ATP5F1A", "gained", 8.5),
                ("CCNF", "lost", 0.5),
                ("NDUFS8", "lost", 0.5),
                ("ZAP70", "lost", 0.5),
                ("RASAL3", "lost", 0.5),
                ("WAPL", "gained", 2.5),
            }
            observed_one_edge = {
                (record.get("gene"), record.get("direction"), record.get("x"))
                for record in one_edge.get("one_edge_segments", [])}
            if observed_one_edge != expected_one_edge:
                failures.append("Figure 2b one-edge tick membership or midpoint changed")

            if fig2_schema == fig2_signed_kegg_schema:
                panel_d_kegg = promotion.get("panel_d_kegg", {})
                expected_kegg = {
                    "analysis":
                        "conventional binary unique-reversal-target KEGG Legacy ORA",
                    "reversal_targets": 2577,
                    "eligible_census_targets": 9554,
                    "filtered_family_terms": 113,
                    "fisher_bh_q_lt_0_10_terms": 14,
                    "regulator_opportunity_maxT_supported_terms": 0,
                    "library": "MSigDB c2.cp.kegg_legacy.v2026.1.Hs",
                }
                for key, expected in expected_kegg.items():
                    if panel_d_kegg.get(key) != expected:
                        failures.append(
                            f"Figure 2d KEGG {key}={panel_d_kegg.get(key)!r}; "
                            f"expected {expected!r}")
                kegg_global_p = float(panel_d_kegg.get(
                    "regulator_opportunity_global_p", float("nan")))
                if not math.isfinite(kegg_global_p) or abs(
                        kegg_global_p - 0.18439078046097696) > 1e-15:
                    failures.append("Figure 2d KEGG global opportunity p-value changed")
                if "does not establish pathway enrichment" not in \
                        panel_d_kegg.get("claim_boundary", "") or \
                        "pathway activity" not in panel_d_kegg.get("claim_boundary", ""):
                    failures.append("Figure 2d KEGG claim boundary changed")
                antigen = panel_d_kegg.get("antigen_processing_orientation", {})
                if antigen != {
                        "lost_edges": 12,
                        "gained_edges": 43,
                        "target_conditioned_maxT_p": 0.0034498275086245686,
                        "regulator_margin_maxT_p": 4.999750012499375e-05,
                        "complementary_single_margin_tests_not_joint": True,
                }:
                    failures.append(
                        "Figure 2d antigen-processing orientation sensitivity changed")
                if panel_d_kegg.get("reversal_dyads") != 4379 or \
                        panel_d_kegg.get("fixed_family_covered_targets") != 796 or \
                        panel_d_kegg.get("fixed_family_covered_dyads") != 1443 or \
                        panel_d_kegg.get("pathway_memberships_overlap") is not True or \
                        panel_d_kegg.get("sole_dual_supported_orientation_term") != \
                        "Antigen processing and presentation":
                    failures.append("Figure 2d KEGG coverage/overlap disclosure changed")
                d_source_path = ROOT / \
                    FIG2_SIGNED_KEGG_V11_EXPECTED["panel_d_source"][0]
                d_rows = []
                d_fields = []
                if d_source_path.is_file():
                    with d_source_path.open(newline="") as handle:
                        reader = csv.DictReader(handle, delimiter="\t")
                        d_fields = list(reader.fieldnames or [])
                        d_rows = list(reader)
                required_d_fields = {
                    "display_order", "term", "kegg_display_name", "query_overlap",
                    "fisher_bh_q", "n_reversal_dyads", "n_lost_dyads",
                    "n_gained_dyads", "target_observed_T", "target_maxT_p",
                    "regulator_observed_T", "regulator_maxT_p",
                    "regulator_sign_relation", "both_single_margin_maxT_supported",
                    "support_marker",
                }
                supported_rows = [
                    row for row in d_rows
                    if row.get("both_single_margin_maxT_supported") == "True"
                ]
                if not required_d_fields.issubset(d_fields) or len(d_rows) != 14 or \
                        [int(row["display_order"]) for row in d_rows] != \
                        list(range(1, 15)) or \
                        panel_d_kegg.get("displayed_term_order") != \
                        [row["term"] for row in d_rows] or \
                        any(row["conventional_BH_descriptive"] != "True"
                            for row in d_rows) or \
                        any(int(row["n_reversal_dyads"]) !=
                            int(row["n_lost_dyads"]) + int(row["n_gained_dyads"])
                            for row in d_rows) or \
                        len(supported_rows) != 1 or \
                        supported_rows[0].get("term") != \
                        "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION":
                    failures.append("Figure 2d exact 14-row KEGG source/order changed")

            timing_letter = "f" if fig2_schema == fig2_signed_kegg_schema else "e"
            result = promotion.get(
                "panel_f_result" if fig2_schema == fig2_signed_kegg_schema
                else "panel_e_result", {})
            expected_counts = {
                "analysis_status": "POST_HOC_EXPLORATORY",
                "covered_reversals": 4377,
                "resolved": 1279,
                "unresolved": 3098,
                "lost_resolved": 704,
                "lost_rest_sign": 529,
                "lost_48h_sign": 175,
                "gained_resolved": 575,
                "gained_rest_sign": 295,
                "gained_48h_sign": 280,
                "paired_regulators": 75,
                "paired_positive": 40,
                "paired_negative": 14,
                "paired_ties": 21,
                "paired_median_difference_exact": "5/56",
                "wilcoxon_statistic": 417.0,
            }
            for key, expected in expected_counts.items():
                if result.get(key) != expected:
                    failures.append(
                        f"Figure 2{timing_letter} result {key}={result.get(key)!r}; "
                        f"expected {expected!r}")
            timing_p = float(result.get("wilcoxon_two_sided_p", float("nan")))
            if not math.isfinite(timing_p) or abs(
                    timing_p - 0.0049970904) > 5e-11:
                failures.append(f"Figure 2{timing_letter} paired p-value changed")
            if "not an exact switch time" not in result.get("claim_boundary", ""):
                failures.append(
                    f"Figure 2{timing_letter} lost its switch-time claim boundary")
            if result.get("display_unit") != "pooled resolved reversal edges" or \
                    result.get("p_value_unit") != "paired regulator-level fractions":
                failures.append(
                    f"Figure 2{timing_letter} conflates displayed edges with the "
                    "paired-regulator test")
            if result.get("lost_rest_sign", 0) + result.get("lost_48h_sign", 0) != 704 or \
                    result.get("gained_rest_sign", 0) + result.get("gained_48h_sign", 0) != 575:
                failures.append(
                    f"Figure 2{timing_letter} Rest-sign and 48 h-sign shares are not "
                    "complementary")

            display = promotion.get(
                "panel_f_display" if fig2_schema == fig2_signed_kegg_schema
                else "panel_e_display", {})
            expected_display = {
                "lost_colour": "#c1272d",
                "gained_colour": "#2a9d8f",
                "rest_sign_colour": "#2c6fbb",
                "stim48_sign_colour": "#e07b39",
                "ordinary_text_floor_pt": 7.0,
                "compact_key_text_pt": 5.5,
            }
            for key, expected in expected_display.items():
                if display.get(key) != expected:
                    failures.append(
                        f"Figure 2{timing_letter} display {key}={display.get(key)!r}; "
                        f"expected {expected!r}")
            expected_trajectories = {
                "lost": [100.0, 100.0 * 529 / 704, 0.0],
                "gained": [100.0, 100.0 * 295 / 575, 0.0],
            }
            expected_bars = {
                "covered": [1279, 3098],
                "lost_resolved": [529, 175],
                "gained_resolved": [295, 280],
            }
            if display.get("encoding") != \
                    "two Rest-sign-retention lines plus adjacent exact classification bars" or \
                    display.get("sampled_states") != ["Rest", "8 h", "48 h"] or \
                    display.get("classification_bars") != expected_bars or \
                    display.get("trajectories") != expected_trajectories:
                failures.append(
                    f"Figure 2{timing_letter} original-line display contract changed")
            if "fixed by selection" not in display.get("endpoint_markers", "") or \
                    "visual guides" not in display.get("connecting_line_boundary", ""):
                failures.append(
                    f"Figure 2{timing_letter} endpoint/connector claim boundary changed")

            terminal_panel = promotion.get(
                "panel_i_preserved_analysis" if fig2_schema == fig2_signed_kegg_schema
                else "panel_h_preserved_analysis", {})
            terminal_letter = "i" if fig2_schema == fig2_signed_kegg_schema else "h"
            terminal_rho = float(terminal_panel.get(
                "observed_rho", float("nan")))
            if terminal_panel.get("data_artists_unchanged") is not True or \
                    not math.isfinite(terminal_rho) or \
                    abs(terminal_rho - 0.09910995705806) > 5e-12:
                failures.append(
                    f"Figure 2{terminal_letter} no longer preserves the 903-target "
                    "orientation result")

            archive = promotion.get("displaced_canonical_archive", {})
            for label, path_key, hash_key in (
                    ("displaced Figure 2 PDF", "pdf_path", "pdf_sha256"),
                    ("displaced Figure 2 PNG", "png_path", "png_sha256"),
                    ("displaced Figure 2 README", "readme_path", "readme_sha256")):
                _check_fig2_bound_file(archive.get(path_key), archive.get(hash_key), label)
            if archive.get("recoverable") is not True or \
                    archive.get("permanent_deletion_performed") is not False:
                failures.append("Figure 2 archive lost its recoverable no-delete contract")

            retained = promotion.get("retained_predecessors", {})
            for label, path_key, hash_key in (
                    ("retained four-trajectory promotion",
                     "four_trajectory_promotion_manifest_path",
                     "four_trajectory_promotion_manifest_sha256"),
                    ("retained movement result", "negative_movement_result_path",
                     "negative_movement_result_sha256"),
                    ("retained movement table", "negative_movement_table_path",
                     "negative_movement_table_sha256")):
                _check_fig2_bound_file(retained.get(path_key), retained.get(hash_key), label)

            audit = promotion.get("audit", {})
            _check_fig2_bound_file(
                audit.get("path"), audit.get("sha256"), "Figure 2 candidate audit")
            if fig2_schema == fig2_signed_kegg_schema and (
                    audit.get("path") != FIG2_SIGNED_KEGG_V11_EXPECTED[
                        "determinism_receipt"][0] or
                    audit.get("sha256") != FIG2_SIGNED_KEGG_V11_EXPECTED[
                        "determinism_receipt"][1]):
                failures.append("Figure 2 v11 promotion audit receipt binding changed")
            no_delete = promotion.get("no_delete_contract", {})
            if no_delete != {
                    "old_canonical_bytes_archived": True,
                    "old_renderer_and_sources_preserved": True,
                    "old_diagnostics_preserved": True,
                    "files_deleted_or_moved": False}:
                failures.append("Figure 2 promotion no-delete contract changed")

            if check_derived and fig2_schema == fig2_signed_kegg_schema:
                build = promotion.get("build_integration", {})
                if build.get("status") != "COMPLETE" or build.get("sequence") != [
                        "source_guard", "docx_build", "pdf_build",
                        "supplementary_build"]:
                    failures.append("Figure 2 signed-KEGG build integration changed")
                for key in (
                        "source_guard", "docx_build", "pdf_build",
                        "supplementary_build"):
                    item = build.get("logs", {}).get(key, {})
                    _check_fig2_bound_file(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 signed-KEGG {key} log")
                for key, expected_path in (
                        ("docx", "HumanCD4CoDEGNet_signflip.docx"),
                        ("pdf", "HumanCD4CoDEGNet_signflip.pdf"),
                        ("supplementary_pdf",
                         "submission/supplementary/Supplementary_Information.pdf")):
                    item = build.get("outputs", {}).get(key, {})
                    if item.get("path") != expected_path:
                        failures.append(
                            f"Figure 2 signed-KEGG {key} path changed")
                    _check_fig2_shared_record(
                        item.get("path"), item.get("sha256"),
                        f"Figure 2 signed-KEGG built {key}")

            # The separately defined 8-of-131 purity result remains guarded even though it is
            # not a main panel.
            predecessor_path = ROOT / supersedes.get("path", "")
            if predecessor_path.is_file():
                cursor = json.loads(predecessor_path.read_text())
                purity = cursor.get("purity_result_retained_in_text_not_figure", {})
                for _ in range(12):
                    if purity:
                        break
                    base_ref = cursor.get("supersedes", {})
                    base_path = ROOT / base_ref.get("path", "")
                    if not base_path.is_file():
                        break
                    cursor = json.loads(base_path.read_text())
                    purity = cursor.get("purity_result_retained_in_text_not_figure", {})
                _check_fig2_bound_file(
                    purity.get("result_path"), purity.get("result_sha256"),
                    "direction-purity result")
                _check_fig2_bound_file(
                    purity.get("manifest_path"), purity.get("manifest_sha256"),
                    "direction-purity manifest")
                if (purity.get("observed"), purity.get("eligible")) != (8, 131):
                    failures.append("Figure 2 purity result no longer binds 8 of 131")

            # Preserve the historical diagnostic reconciliation as provenance; it no longer
            # represents the active canonical bytes.
            reconciliation_path = ROOT / FIG2_DIAGNOSTIC_RECONCILIATION
            if not reconciliation_path.is_file() or hashlib.sha256(
                    reconciliation_path.read_bytes()).hexdigest() != \
                    FIG2_DIAGNOSTIC_RECONCILIATION_SHA256:
                failures.append("historical Figure 2 diagnostic reconciliation changed")

            legend = (ROOT / "FIGURE_LEGENDS.md").read_text()
            fig2_legend = legend.split("## Figure 2.", 1)[1].split("## Figure 3.", 1)[0]
            if fig2_schema == fig2_signed_kegg_schema:
                required_legend = (
                    "**(d)** Conventional binary KEGG over-representation analysis",
                    "2,577 targets carrying at least one reversal against all 9,554 eligible census targets",
                    "Fourteen of the fixed 113-term KEGG Legacy family",
                    "global *p* = 0.184; 0/113 terms passed single-step maxT correction",
                    "rather than establish pathway enrichment or activity",
                    "12 lost and 43 gained reversal edges",
                    "complementary single-margin tests, not a joint test",
                    "**(e)** Every reversal target of the six regulators inside the "
                    "49-regulator enriched family",
                    "**(f)** Sign state at the unselected 8 h timepoint",
                    "**(g)** Donor AUC",
                    "**(h)** Absolute 0-48 h responses",
                    "**(i)** Across all 903 mapped targets",
                    "`rev/test; rate; fold; L/G`",
                    "529/704 lost effects (75.1%)",
                    "175/704 (24.9%)",
                    "295/575 gained effects (51.3%)",
                    "280/575 (48.7%)",
                    "connecting lines are visual guides rather than continuous trajectories",
                    "not an exact switch time or molecular mechanism",
                )
            else:
                required_legend = (
                    "**(d)** Every reversal target of the six regulators inside the "
                    "49-regulator enriched family",
                    "**(e)** Sign state at the unselected 8 h timepoint",
                    "**(f)** Donor AUC",
                    "**(g)** Absolute 0-48 h responses",
                    "**(h)** Across all 903 mapped targets",
                    "529/704 lost effects (75.1%)",
                    "175/704 (24.9%)",
                    "295/575 gained effects (51.3%)",
                    "280/575 (48.7%)",
                    "Red lost and teal gained lines show the fraction of resolved edges retaining the Rest sign",
                    "the adjacent bars show the complete resolved/unresolved classification",
                    "The main bars retain the linear 0–300 edge-count scale",
                    "Slim colour ticks at the exact midpoint of each one-edge segment",
                    "the underlying bar endpoints remain the count encoding",
                    "`rev/test; rate; fold; L/G`",
                    "connecting lines are visual guides rather than continuous trajectories",
                    "not an exact switch time or molecular mechanism",
                )
            for required in required_legend:
                if required not in fig2_legend:
                    failures.append(f"Figure 2 legend missing required boundary: {required}")
            stale_i = fig2_schema != fig2_signed_kegg_schema and "**(i)**" in fig2_legend
            if stale_i or "p* = 0.572" in fig2_legend or \
                    "The stacked key repeats these proportions" in fig2_legend:
                failures.append("Figure 2 legend retains a stale panel")

            fig2_results = (ROOT / "MANUSCRIPT.md").read_text().split(
                "### Figure 2.", 1)[1].split("### Figure 3.", 1)[0]
            if fig2_schema == fig2_signed_kegg_schema:
                required_results = (
                    "A conventional binary analysis of the 2,577 targets receiving at least one reversal",
                    "14 of 113 filtered KEGG Legacy pathways",
                    "global $p = 0.184$; 0/113 maxT-supported; Figure 2d",
                    "descriptive sensitivity rather than evidence of pathway enrichment or activity",
                    "12 lost and 43 gained reversal edges",
                    "complementary margins and do not constitute a joint test",
                    "Across these six regulators, the 39 reversal targets comprised "
                    "23 lost and 16 gained effects",
                    "common switch time (Figure 2e)",
                    "Figure 2f",
                    "Figure 2g",
                    "Figure 2h",
                    "Figure 2i",
                    "529/704 (75.1%)",
                    "175/704 (24.9%)",
                    "295/575 (51.3%)",
                    "280/575 (48.7%)",
                    "3,098 remained unresolved",
                    "Rest and 48 h defined the reversal set and 8 h classified it afterwards",
                    "Together, panels g–i",
                )
            else:
                required_results = (
                    "Across these six regulators, the 39 reversal targets comprised "
                    "23 lost and 16 gained effects",
                    "common switch time (Figure 2d)",
                    "Figure 2e",
                    "Figure 2f",
                    "Figure 2g",
                    "Figure 2h",
                    "529/704 (75.1%)",
                    "175/704 (24.9%)",
                    "295/575 (51.3%)",
                    "280/575 (48.7%)",
                    "3,098 remained unresolved",
                    "Rest and 48 h defined the reversal set and 8 h classified it afterwards",
                )
            for required in required_results:
                if required not in fig2_results:
                    failures.append(f"Figure 2 Results missing required text: {required}")
            for pattern in (
                    r"0\.337[^\n]{0,100}0\.143",
                    r"1\.9\s*[x×]\s*10\^-?6",
                    r"124/156",
                    r"0\.077[^\n]{0,120}0\.055[^\n]{0,120}Figure 2h"):
                if re.search(pattern, fig2_results, flags=re.I):
                    failures.append(f"Figure 2 Results retains stale panel-h claim /{pattern}/")

    successor_checks = (
        (
            figure3_successor_state,
            ROOT / ACTIVE_FIGURE_SUCCESSORS["active_figure_3"]["promotion_path"],
            audit_candidate_v13,
            check_fig3_unified_zoom_v13_promotion,
        ),
        (
            figure4_successor_state,
            ROOT / ACTIVE_FIGURE_SUCCESSORS["active_figure_4"]["promotion_path"],
            audit_candidate_v7,
            check_fig4_fig3e_crossref_promotion,
        ),
    )
    for state, promotion_path, candidate_audit, promotion_check in successor_checks:
        if legend_concision_active:
            # The legend successor independently hash-locks both sealed parent
            # manifests, their non-text live assets/candidates and archives.
            # Their old whole-file text/guard/build records are intentionally
            # historical, while the new successor owns the live sources and
            # next derived builds.  Do not rewrite the immutable v13/v7 files.
            continue
        if state == "bound":
            promotion_check(
                ROOT, check_derived=check_derived, failures=failures)
        elif state == "pending":
            if promotion_path.is_file():
                promotion_check(ROOT, check_derived=False, failures=failures)
            else:
                candidate_audit(ROOT, failures)

    if check_derived:
        _check_derived_artifacts(failures)

    if failures:
        raise SystemExit("ACTIVE CLAIM GUARD FAILED\n" + "\n".join(f"- {x}" for x in failures))
    scope = "source and derived-artifact" if check_derived else "source"
    print(
        f"ACTIVE CLAIM GUARD PASS ({scope} mode): trans estimand, target-trajectory "
        "provenance, package contract, withdrawn analyses, and reframe NO-GO locks are coherent")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--post-build", action="store_true",
        help="also inspect generated PDF/DOCX/SI files and the packaged S1-S4/Data1-21 contract")
    args = parser.parse_args()
    check_active_sources(check_derived=args.post_build)
