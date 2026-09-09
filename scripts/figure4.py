#!/usr/bin/env python3
"""Figure 4 - single self-contained builder.

Flattened from the eight-module response_shape chain. The frozen HPC output
bundle under data/response_shape/production_v5_1/ is the input; the per-cell
decomposition behind it ran on a cluster from the source atlas (39 GB, not
redistributed - see Data availability).

The original chain carried sha256 self-pins, predecessor-hash checks and a
0550/0440 permission contract on its bundle. Those guarded against version
drift between eight files; with one file there are no versions to drift, so
they are removed and setup_after_clone.sh is no longer needed.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))  # figstyle lives beside this file
import figstyle as S
from dataclasses import dataclass
from datetime import datetime, timezone
from matplotlib.collections import Collection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from typing import Iterable, Sequence
import argparse
import hashlib
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import shutil
import stat
import sys
import tempfile

from pathlib import Path as _Path
ROOT = _Path(__file__).resolve().parents[1]

import shutil as _shutil
import tempfile as _tempfile

# The v6 candidate render is an intermediate that v9 reads back in the same run.
# It used to be written into the repository as promotion history; a run-local
# temp directory keeps the repo free of build artefacts.
def _rel(path, root):
    """Repo-relative path for provenance records.

    Intermediate renders now stage in a temp directory, which is deliberately
    outside the repository, so relative_to() would raise. Fall back to the file
    name for those.
    """
    try:
        return str(_P(path).relative_to(root))
    except ValueError:
        return _P(path).name


STAGING_DIR = _P(_tempfile.gettempdir()) / "codeg_fig4_staging"
_shutil.rmtree(STAGING_DIR, ignore_errors=True)
STAGING_DIR.mkdir(parents=True, exist_ok=True)

_MODULES: dict[str, types.ModuleType] = {}


# ==========================================================================
# fig4_production_v5_1
# ==========================================================================
def _load_fig4_production_v5_1():
    _self = types.ModuleType('fig4_production_v5_1')
    __name__ = 'fig4_production_v5_1'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_production_v5_1'] = _self
    #!/usr/bin/env python3
    """Validate and render Figure 4 from the sealed v5.1 gather outputs only.

    Running this program without ``--render`` performs validation and writes
    nothing.  Rendering is possible only after the exact production gather files,
    their detached manifest, real target-read provenance, the transparent audit
    exception and the independent proportionate custody closure have all passed.
    Synthetic, preview, simulated, v4, unsealed or partial results are rejected
    before plotting.
    """




    PAPER = Path(__file__).resolve().parents[1]
    MPL_CACHE = PAPER / "tmp/matplotlib_response_shape_production_v5_1"
    MPL_CACHE.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))


    matplotlib.use("Agg")

    if str(PAPER) not in sys.path:
        sys.path.insert(0, str(PAPER))
    if str(PAPER / "scripts") not in sys.path:
        sys.path.insert(0, str(PAPER / "scripts"))



    SCHEMA = "response-shape-figure4-production-v5.1"
    STATUS = "PRODUCTION_RESULTS_HASH_SEALED"
    DATA_ORIGIN = "real_atlas_target_expression"
    CLAIM_BOUNDARY = "marginal_donor_balanced_response_shape_not_biochemical_mechanism"

    DENOMINATOR = 4_379
    METADATA_ELIGIBLE = 3_000
    METADATA_UNSUPPORTED = 1_379
    ELIGIBLE_REGULATORS = 149

    FROZEN_CONFIG_SHA256 = "171992cf0c2a40eb19a4bb206ba84a5db0ef2780818561d5986199607ea37fe4"
    FROZEN_AMENDMENT_SHA256 = "dde9b2fb0353bc60fcbd6bc2998b6e9e98bfb59409f59fc606a6acef4039b6f1"
    FROZEN_METADATA_REPORT_SHA256 = "83cb9b8ddf20a5e950231f4e6e743a06d264ad22a29b7c9b7b380a05093962ef"
    FROZEN_CAL_MAP_SHA256 = "e6186404ba3cda52db6dd3b151d88a50d6d520e12d4352ab6aa26bbf5cf7a3e5"
    FROZEN_CENSUS_SHA256 = "0d753121d0e6e69313acd20b29e3431f01ad083a9f3515598f8e2d83775298ad"
    FROZEN_STABLE_MATCHES_SHA256 = "ce726b1a258e5e06c31a117185d566812c2e1ffb6e6c68b8c02873ee471016b3"

    PRODUCTION_RUN_ID = "v5_1_full_20260814_r1"
    PRODUCTION_PREPARATION_SEAL_SHA256 = "c05f4e50445696d4011b44c759d6f1d86a5a243dbc541b511ef51a95ad938909"
    PRODUCTION_EXTRACTOR_ENTRY_SHA256 = "740e775067cedab3d41a2273c1861885216163ed4d39c368c145975203be2248"
    PRODUCTION_STREAM_EXTRACTOR_SHA256 = "991c85b2b3be175ee8c1dcff6b70e2e7e3923eee35e1cc84b600470bfb61fdf9"
    PRODUCTION_ESTIMATOR_SHA256 = "9a127c5067b4c1283dcf1d399ea0a26969b76ac04f0c8ae2c5c8528ee5664718"
    PRODUCTION_GATHER_SHA256 = "c8902ccc9b3546d2ba9d0a8df517564e47d410e6b232c566c946076076beca77"
    PRODUCTION_CONTRACT_SHA256 = "9ea202e541e42af3adb98d967633a6b91ef60c617339f27aef8b53eb6b182e06"
    PRODUCTION_WRAPPER07_SHA256 = "965d5b3b6c3d28bdc8d75ed0658aab762579bb5bd20dfbbd6b323d000672a70b"
    PRODUCTION_AUDIT_EXCEPTION_SHA256 = "a0f6a462769974d851ea56a5252a7fee300c1165402d2698cae1f60e80c489e2"
    PRODUCTION_FINAL_CUSTODY_CLOSURE_SHA256 = "1d63ee9c88c4b078a2d57f510257733394a4687b8d7f65f43c29aa945d742b53"
    PRODUCTION_TARGET_UNIVERSE_SHA256 = "2abd91ca021bf806583a892970b14d1e67e9f72b7e444395e01b857ab467d40c"
    PRODUCTION_TARGET_VALUES_READ = 1_776_344_553
    PRODUCTION_SELECTED_CELLS = 816_502

    PRODUCTION_MANIFEST_SHA256 = {
        "D1_Rest": "04f6cabb7df1231893c2327200baea5e634e4d04073af5f7b7fe82e730e6a6b0",
        "D1_Stim48hr": "bded112e786b15b9a57e5b303f3408c0150ffd3c795dbe7f2f17ecb34774bc10",
        "D2_Rest": "c82394d4ef101ec6ed86e4d1efe822bc6bca482e65fd5805763978aefa6b817f",
        "D2_Stim48hr": "d8c1addf236c833b6660ff7f8f86535e05010a1bdb5979161229c8a1f45fa3b0",
        "D3_Rest": "66b8e8084fd2072c6e423e3e3d20ad21164c407af67f0992ea683c4cd2668989",
        "D3_Stim48hr": "04537666d8d01f8ecfe53561892cb0ad9c12da44d775d6aa67628de55106f3a3",
        "D4_Rest": "09fe35c5ba49dbbaa3efda5bbcc90234cf645c1fce55e506c4062308fd3a16e6",
        "D4_Stim48hr": "31bc1f6af3d4dc591d32bed22431ed77c55be3b079073ead9971904ea6251ee5",
    }
    PRODUCTION_AUTHORIZATION_SHA256 = {
        "D1_Rest": "6bd843ac1b374ebdbe61c76ec4153ffdbdf287cffa76facb87a1d2272ef44c07",
        "D1_Stim48hr": "590c21b14906b7b4a19789d0d1212064c3848a0f3d3ba9dcb4184b3ad0096e8e",
        "D2_Rest": "5baf9178d223a950f56053d40b28b5796274654117333f52ec38bb5e54b11edb",
        "D2_Stim48hr": "030aa99de89b1e3072b25488d2a1173f08b89d5d8989ab30eebeb0d7f71194f0",
        "D3_Rest": "24af14c1408794a3314ed62b0206918bb4f0cf10c604a7d586ce71009be822cf",
        "D3_Stim48hr": "89125b0617671cdbfbad06d98c5e23f33a1ca83410ad258d1a6603f53756a834",
        "D4_Rest": "9afce757a971ba7fc6b8f0cfcb61f2168a21c5b20c698d3391ed246efdf467c9",
        "D4_Stim48hr": "37a511dc3db66426d1349bbf5c77ed21fba383cda22749ea513797aad06a5d50",
    }

    DEFAULT_CONFIG = PAPER / "data/methods_config_v5_1.json"
    DEFAULT_AMENDMENT = PAPER / "data/response_shape/amendments/2026-08-14_PREOUTCOME_DONOR_LEVEL_TWO_PART_V5_1.md"
    DEFAULT_METADATA_REPORT = PAPER / "data/response_shape/metadata_feasibility_v5_1.json"
    DEFAULT_CAL_MAP = PAPER / "data/response_shape/cal_pseudo_pairs_v5_1.tsv"
    DEFAULT_FROZEN_CENSUS = PAPER / "data/census/census_master_edges.csv.gz"
    DEFAULT_STABLE_MATCHES = PAPER / "data/response_shape/v5_controls_417525_r1/stable_matches_v4.tsv"
    DEFAULT_MANIFEST = PAPER / "data/response_shape/production_v5_1/figure4_production_manifest_v5_1.json"
    DEFAULT_PDF = PAPER / "figures/Fig4.pdf"

    FIGURE_WIDTH_IN = 7.10
    FIGURE_HEIGHT_WITH_EXAMPLES_IN = 10.60
    FIGURE_HEIGHT_WITHOUT_EXAMPLES_IN = 8.25
    MINIMUM_FINAL_TEXT_PT = 7.0

    SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
    RUN_ID_RE = re.compile(r"^v5_1_[A-Za-z0-9][A-Za-z0-9_.-]{0,87}$")
    EDGE_RE = re.compile(r"^REV\|ENSG[0-9]+\|ENSG[0-9]+$")

    CONDITIONS = ("Rest", "Stim48hr")
    DONORS = ("D1", "D2", "D3", "D4")
    SHARDS = tuple(f"{donor}_{condition}" for donor in DONORS for condition in CONDITIONS)
    METRICS = (
        "total_effect",
        "detection_component",
        "positive_intensity_component",
        "raw_rate_log2_ratio",
    )
    COMPONENT_METRICS = METRICS[:3]
    ANALYSIS_ROLES = (
        "REVERSAL",
        "REGULATOR_MATCHED_STABLE",
        "CAL_REVERSAL",
        "CAL_STABLE",
    )
    FINAL_STATUSES = (
        "UNANALYSABLE_METADATA_SUPPORT",
        "UNANALYSABLE_POSITIVE_INTENSITY",
        "QUARANTINED_NUMERICAL_ADDITIVITY",
        "ESTIMATED_CONTINUOUS_COMPONENTS",
        "MIXED_OR_UNRESOLVED",
    )
    SUMMARY_STATUSES = (
        "UNANALYSABLE_POSITIVE_INTENSITY",
        "QUARANTINED_NUMERICAL_ADDITIVITY",
        "ESTIMATED_CONTINUOUS_COMPONENTS",
    )
    UNIT_STATUSES = SUMMARY_STATUSES
    NUMERICAL_FINAL = frozenset(("ESTIMATED_CONTINUOUS_COMPONENTS", "MIXED_OR_UNRESOLVED"))

    PREDECLARED_EDGES = (
        ("REV|ENSG00000173011|ENSG00000126353", "TADA2B", "CCR7", "gained_on_activation", 1),
        ("REV|ENSG00000162227|ENSG00000126353", "TAF6L", "CCR7", "gained_on_activation", 1),
        ("REV|ENSG00000162227|ENSG00000119922", "TAF6L", "IFIT2", "lost_on_activation", -1),
        ("REV|ENSG00000152234|ENSG00000109471", "ATP5F1A", "IL2", "lost_on_activation", -1),
    )

    GATHERED_FILENAMES = {
        "unit_effects": "unit_effects_v5_1.parquet",
        "guide_specific_effects": "guide_specific_effects_v5_1.parquet",
        "edge_continuous_summary": "edge_continuous_summary_v5_1.parquet",
        "cal_envelopes": "cal_envelopes_v5_1.tsv",
        "reversal_census": "reversal_census_v5_1.parquet",
        "status_waterfall": "status_waterfall_v5_1.tsv",
        "gather_report": "gather_report_v5_1.json",
    }
    AUDIT_EXCEPTION_KEY = "audit_exception"
    AUDIT_EXCEPTION_FILENAME = "audit_exception_v5_1.json"
    FINAL_CLOSURE_KEY = "final_gather_custody_closure"
    FINAL_CLOSURE_FILENAME = "final_gather_custody_closure_v5_1.json"
    ARTIFACT_FILENAMES = {
        **GATHERED_FILENAMES,
        AUDIT_EXCEPTION_KEY: AUDIT_EXCEPTION_FILENAME,
        FINAL_CLOSURE_KEY: FINAL_CLOSURE_FILENAME,
    }
    DATA_ARTIFACTS = tuple(GATHERED_FILENAMES)

    FORBIDDEN_COLUMN_TOKENS = (
        "synthetic",
        "simulation",
        "preview",
        "auc",
        "state_stratum",
        "depth_stratum",
        "activation_state",
        "regression",
        "p_value",
        "pvalue",
        "q_value",
        "qvalue",
        "padj",
        "class_enrichment",
        "mechanism_class",
        "wasserstein",
        "responder_fraction",
    )

    UNIT_REQUIRED = {
        "analysis_edge_id",
        "parent_reversal_edge_id",
        "analysis_role",
        "regulator_ensg",
        "target_ensg",
        "donor",
        "condition",
        "expected_sign",
        "orientation_source",
        "calibration_scope",
        "unit_status",
        "support_reason",
        "guide_A_guide_id",
        "guide_B_guide_id",
        "guide_A_n_detected",
        "guide_B_n_detected",
        "guide_A_ref_raw_detected_cells",
        "guide_B_ref_raw_detected_cells",
        "guide_A_ref_distinct_detecting_physical_guides",
        "guide_B_ref_distinct_detecting_physical_guides",
        "guide_A_p_case",
        "guide_B_p_case",
        "guide_A_q_case",
        "guide_B_q_case",
        "guide_A_p_ref",
        "guide_B_p_ref",
        "guide_A_q_ref",
        "guide_B_q_ref",
        "guide_A_raw_rate_log2_ratio",
        "guide_B_raw_rate_log2_ratio",
        "p_case",
        "q_case",
        "p_ref",
        "q_ref",
        "mu_case",
        "mu_ref",
        "total_effect",
        "detection_component",
        "positive_intensity_component",
        "additivity_error",
        "additivity_tolerance",
        "additivity_pass",
        "raw_rate_log2_ratio",
        *(f"guide_{slot}_{field}" for slot in ("A", "B") for field in (
            "n_cells", "guide_mu_case", "guide_mu_ref", "guide_total_effect",
            "guide_detection_component", "guide_positive_intensity_component",
            "guide_additivity_error", "guide_additivity_pass",
        )),
    }
    GUIDE_REQUIRED = {
        "analysis_edge_id",
        "parent_reversal_edge_id",
        "analysis_role",
        "regulator_ensg",
        "target_ensg",
        "donor",
        "condition",
        "unit_status",
        "guide_slot",
        "guide_id",
        "n_cells",
        "n_detected",
        "p_case",
        "q_case",
        "p_ref",
        "q_ref",
        "raw_rate_log2_ratio",
        "ref_raw_detected_cells",
        "ref_distinct_detecting_physical_guides",
        "guide_mu_case",
        "guide_mu_ref",
        "guide_total_effect",
        "guide_detection_component",
        "guide_positive_intensity_component",
        "guide_additivity_error",
        "guide_additivity_pass",
    }
    SUMMARY_REQUIRED = {
        "analysis_edge_id",
        "parent_reversal_edge_id",
        "analysis_role",
        "regulator_ensg",
        "target_ensg",
        "common_donor_ids",
        "n_common_donors",
        "summary_status",
        "analysis_status",
        "orientation",
        *(f"{condition}__{metric}" for condition in CONDITIONS for metric in METRICS),
        *(f"delta__{metric}" for metric in METRICS),
        *(f"oriented_delta__{metric}" for metric in METRICS),
        *(f"{condition}__raw_practical_opposite" for condition in CONDITIONS),
        *(f"{condition}__raw_CAL_q95_available" for condition in CONDITIONS),
        *(f"{condition}__raw_CAL_q95" for condition in CONDITIONS),
        "oriented_delta__guide_A_guide_total_effect",
        "oriented_delta__guide_B_guide_total_effect",
    }
    ENVELOPE_REQUIRED = {
        "metric",
        "scope",
        "condition",
        "available",
        "absolute_q95",
        "n_unique_supported_CAL",
        "quantile",
        "quantile_definition",
        "source_role",
        "stable_target_CAL_included",
        "edge_weights",
    }
    CENSUS_REQUIRED = {
        "reversal_edge_id",
        "regulator",
        "regulator_ensg",
        "target",
        "target_ensg",
        "direction",
        "expected_sign_Rest",
        "expected_sign_Stim48hr",
        "metadata_eligible",
        "metadata_status",
        "analysis_status",
    }


    class Figure4ContractError(ValueError):
        """Raised before plotting when production custody or biology fails."""


    @dataclass(frozen=True)
    class FrozenContext:
        config: Mapping[str, Any]
        reversal: pd.DataFrame
        cal_map: pd.DataFrame
        stable: pd.DataFrame
        eligible_edge_ids: frozenset[str]
        expected_analysis_ids: Mapping[str, frozenset[str]]


    @dataclass(frozen=True)
    class ProductionBundle:
        manifest: Mapping[str, Any]
        manifest_sha256: str
        artifact_sha256: Mapping[str, str]
        config: Mapping[str, Any]
        unit_effects: pd.DataFrame
        guide_specific_effects: pd.DataFrame
        edge_summary: pd.DataFrame
        cal_envelopes: pd.DataFrame
        census: pd.DataFrame
        waterfall: pd.DataFrame
        gather_report: Mapping[str, Any]
        audit_exception: Mapping[str, Any]
        final_custody_closure: Mapping[str, Any]


    def sha256_file(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def _expect(condition: bool, message: str) -> None:
        if not condition:
            raise Figure4ContractError(message)


    def _exact_true(value: object, label: str) -> None:
        _expect(type(value) is bool and value, f"{label} must be literal true")


    def _exact_false(value: object, label: str) -> None:
        _expect(type(value) is bool and not value, f"{label} must be literal false")


    def _sha(value: object, label: str) -> str:
        text = str(value)
        _expect(bool(SHA256_RE.fullmatch(text)) and text != "0" * 64, f"{label} is not a resolved lowercase SHA-256")
        return text


    def _read_json(path: Path, label: str) -> dict[str, Any]:
        try:
            value = json.loads(Path(path).read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise Figure4ContractError(f"cannot read {label}: {exc}") from exc
        _expect(isinstance(value, dict), f"{label} must be a JSON object")
        return value


    def _read_table(path: Path, label: str) -> pd.DataFrame:
        try:
            if path.suffix == ".parquet":
                return pd.read_parquet(path)
            return pd.read_csv(path, sep="\t")
        except Exception as exc:  # pragma: no cover
            raise Figure4ContractError(f"cannot read {label}: {exc}") from exc


    def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
        missing = sorted(required.difference(table.columns))
        _expect(not missing, f"{label} lacks required columns {missing}")
        forbidden = sorted(
            column
            for column in table.columns
            if any(token in str(column).lower() for token in FORBIDDEN_COLUMN_TOKENS)
        )
        _expect(not forbidden, f"{label} contains forbidden columns {forbidden}")


    def _all_bool(series: pd.Series, label: str) -> None:
        _expect(series.map(lambda value: isinstance(value, (bool, np.bool_))).all(), f"{label} must contain booleans only")


    def _finite(series: pd.Series, label: str) -> np.ndarray:
        values = pd.to_numeric(series, errors="coerce").to_numpy(float)
        _expect(np.isfinite(values).all(), f"{label} contains non-finite values")
        return values


    def _integers(series: pd.Series, label: str) -> pd.Series:
        values = pd.to_numeric(series, errors="coerce")
        _expect(np.isfinite(values).all() and np.allclose(values, np.rint(values), atol=0, rtol=0), f"{label} must contain finite integers")
        return values.astype(int)


    def _canonical_donors(value: object) -> tuple[str, ...]:
        if value is None or (isinstance(value, float) and np.isnan(value)) or str(value).strip() == "":
            return ()
        donors = tuple(str(value).split(","))
        _expect(len(donors) == len(set(donors)), f"duplicate common donor in {value!r}")
        _expect(set(donors).issubset(DONORS), f"invalid common donor set {value!r}")
        _expect(donors == tuple(donor for donor in DONORS if donor in donors), f"common donor IDs are not in frozen order: {value!r}")
        return donors


    def _validate_frozen_context(
        config_path: Path,
        amendment_path: Path,
        metadata_report_path: Path,
        cal_map_path: Path,
        census_path: Path,
        stable_path: Path,
    ) -> FrozenContext:
        pins = {
            Path(config_path): FROZEN_CONFIG_SHA256,
            Path(amendment_path): FROZEN_AMENDMENT_SHA256,
            Path(metadata_report_path): FROZEN_METADATA_REPORT_SHA256,
            Path(cal_map_path): FROZEN_CAL_MAP_SHA256,
            Path(census_path): FROZEN_CENSUS_SHA256,
            Path(stable_path): FROZEN_STABLE_MATCHES_SHA256,
        }
        for path, expected in pins.items():
            # The pre-outcome amendment is an internal record, not redistributed;
            # every other frozen input must be present and match.
            if not path.is_file() and path == Path(amendment_path):
                continue
            _expect(path.is_file() and sha256_file(path) == expected, f"frozen input hash mismatch: {path.name}")
        config = _read_json(Path(config_path), "v5.1 methods config")
        _expect(config.get("schema_version") == "response-shape-methods-v5.1", "wrong v5.1 config schema")
        _expect(config.get("status") == "FROZEN_PREOUTCOME_V5_1_TARGET_ACCESS_NOT_AUTHORIZED", "v5.1 config is not frozen")
        _expect(config.get("amendment", {}).get("sha256") == FROZEN_AMENDMENT_SHA256, "config amendment pin mismatch")
        _expect(config.get("expected_metadata", {}).get("frozen_edges") == DENOMINATOR, "config denominator mismatch")
        _expect(config.get("expected_metadata", {}).get("eligible_edges") == METADATA_ELIGIBLE, "config metadata count mismatch")
        for key in ("individual_edge_p_values", "cell_level_p_values", "class_enrichment_hypotheses_or_p_values", "hard_mechanism_label"):
            _exact_false(config.get("inference_boundary", {}).get(key), f"config.inference_boundary.{key}")
        report = _read_json(Path(metadata_report_path), "v5.1 metadata report")
        _expect(report.get("status") == "PASS_METADATA_AND_FIXED_CAL_MAP_TARGET_ACCESS_NOT_AUTHORIZED", "metadata/CAL feasibility did not pass")
        _expect(report.get("methods_config", {}).get("sha256") == FROZEN_CONFIG_SHA256, "metadata report config pin mismatch")
        _expect(report.get("amendment", {}).get("sha256") == FROZEN_AMENDMENT_SHA256, "metadata report amendment pin mismatch")
        _expect(report.get("fixed_CAL_map", {}).get("sha256") == FROZEN_CAL_MAP_SHA256, "metadata report CAL-map pin mismatch")
        _expect(report.get("full_denominator", {}).get("frozen_reversal_edges") == DENOMINATOR, "metadata report denominator mismatch")
        _expect(report.get("full_denominator", {}).get("metadata_eligible_reversal_edges") == METADATA_ELIGIBLE, "metadata report eligible count mismatch")
        _expect(report.get("outcome_custody", {}).get("target_expression_values_read") == 0, "pre-outcome metadata report read target outcomes")

        cal = pd.read_csv(cal_map_path, sep="\t")
        required_cal = {"regulator_ensg", "regulator_symbol", "paired_donor_set", "case_guide_A", "case_guide_B", "cal_guide_A", "cal_guide_B", "cal_selection_status", "target_or_expression_values_used"}
        _expect(required_cal.issubset(cal.columns), "frozen CAL map lacks required fields")
        _expect(len(cal) == ELIGIBLE_REGULATORS and not cal["regulator_ensg"].duplicated().any(), "frozen CAL map must contain 149 unique regulators")
        _expect(cal["cal_selection_status"].eq("FIXED_PREOUTCOME_K1").all(), "CAL map is not fixed K=1")
        _all_bool(cal["target_or_expression_values_used"], "CAL outcome-use flags")
        _expect((~cal["target_or_expression_values_used"].astype(bool)).all(), "CAL pairing used target outcomes")
        for donors in cal["paired_donor_set"]:
            _expect(len(_canonical_donors(donors)) >= 3, "CAL map has fewer than three paired donors")

        canonical = pd.read_csv(census_path)
        canonical = canonical.loc[canonical["kind"].eq("rev")].copy()
        canonical["reversal_edge_id"] = "REV|" + canonical["regulator_ensg"].astype(str) + "|" + canonical["target_ensg"].astype(str)
        canonical["expected_sign_Stim48hr"] = canonical["direction"].map({"gained_on_activation": 1, "lost_on_activation": -1})
        canonical["expected_sign_Rest"] = -canonical["expected_sign_Stim48hr"]
        _expect(len(canonical) == DENOMINATOR and not canonical["reversal_edge_id"].duplicated().any(), "frozen reversal universe is not exactly 4,379 unique edges")
        eligible_regulators = set(cal["regulator_ensg"].astype(str))
        eligible = canonical.loc[canonical["regulator_ensg"].astype(str).isin(eligible_regulators)]
        _expect(len(eligible) == METADATA_ELIGIBLE, "frozen eligible edge set is not 3,000")
        eligible_ids = frozenset(eligible["reversal_edge_id"].astype(str))

        stable = pd.read_csv(stable_path, sep="\t")
        stable = stable.loc[stable["match_family"].eq("regulator_matched") & stable["parent_reversal_edge_id"].isin(eligible_ids)].copy()
        _expect(len(stable) == 5_985 and stable["parent_reversal_edge_id"].nunique() == 2_996, "eligible regulator-matched stable family dimensions drift")
        stable_ids = frozenset(
            f"{row.control_edge_id}|PARENT|{row.parent_reversal_edge_id}"
            for row in stable.itertuples(index=False)
        )
        expected_analysis_ids = {
            "REVERSAL": eligible_ids,
            "CAL_REVERSAL": frozenset(f"CAL_REV|{edge_id}" for edge_id in eligible_ids),
            "REGULATOR_MATCHED_STABLE": stable_ids,
            "CAL_STABLE": frozenset(f"CAL_STABLE|{edge_id}" for edge_id in stable_ids),
        }
        return FrozenContext(config, canonical, cal, stable, eligible_ids, expected_analysis_ids)


    def _file_mode(path: Path) -> int:
        return stat.S_IMODE(path.stat().st_mode)


    def _validate_locked_production_stage(manifest_path: Path, checksum_path: Path) -> None:
        root = manifest_path.parent
        _expect(root.is_dir() and not root.is_symlink() , "production result stage must be a non-symlink directory")
        _expect(checksum_path.parent == root, "production manifest checksum must share the sealed stage")
        expected_names = {
            manifest_path.name,
            checksum_path.name,
            *ARTIFACT_FILENAMES.values(),
        }
        observed = {path.name for path in root.iterdir()}
        _expect(observed == expected_names, "production result stage inventory is not exact")
        for path in root.iterdir():
            _expect(path.is_file() and not path.is_symlink() , f"production result stage file is linked or non-regular: {path.name}")


    def _validate_manifest(manifest_path: Path, checksum_path: Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, str]]:
        lowered = str(manifest_path).lower()
        _expect(not any(token in lowered for token in ("synthetic", "simulation", "preview", "mock", "toy")), "production builder refuses a non-production manifest path")
        _expect(manifest_path.name == "figure4_production_manifest_v5_1.json" and checksum_path.name == "figure4_production_manifest_v5_1.sha256", "production manifest filenames are not exact")
        _validate_locked_production_stage(manifest_path, checksum_path)
        fields = checksum_path.read_text().strip().split()
        _expect(len(fields) == 2 and _sha(fields[0], "manifest checksum") == sha256_file(manifest_path), "detached production-manifest checksum mismatch")
        _expect(fields[1].lstrip("*") == manifest_path.name, "detached checksum names another manifest")
        manifest = _read_json(manifest_path, "Figure 4 production manifest")
        expected_fields = {
            "schema", "status", "production", "outcomes_read", "real_target_level_outcomes_read",
            "synthetic_only", "data_origin", "run_id", "denominator_n", "metadata_eligible_n",
            "status_counts", "methods_config_sha256", "amendment_sha256", "metadata_report_sha256",
            "cal_map_sha256", "frozen_census_sha256", "stable_matches_sha256",
            "preparation_seal_sha256", "input_shard_manifest_sha256", "target_access_authorization_sha256",
            "extractor_entry_sha256", "stream_extractor_sha256", "estimator_sha256",
            "gather_sha256", "contract_sha256", "wrapper07_sha256",
            "audit_exception_sha256", "final_custody_closure_sha256", "builder_sha256",
            "figure4_gate_pass", "proportionate_audit_exception_recorded",
            "proportionate_final_custody_closure_pass", "artifacts",
        }
        _expect(set(manifest) == expected_fields, "production manifest field set is not exact")
        _expect(manifest.get("schema") == SCHEMA and manifest.get("status") == STATUS, "wrong or unsealed Figure 4 production manifest")
        _exact_true(manifest.get("production"), "manifest.production")
        _exact_true(manifest.get("outcomes_read"), "manifest.outcomes_read")
        _exact_true(manifest.get("real_target_level_outcomes_read"), "manifest.real_target_level_outcomes_read")
        _exact_false(manifest.get("synthetic_only"), "manifest.synthetic_only")
        _expect(manifest.get("data_origin") == DATA_ORIGIN, "manifest data origin mismatch")
        _expect(manifest.get("run_id") == PRODUCTION_RUN_ID, "manifest run_id is not the exact production run")
        _expect(manifest.get("denominator_n") == DENOMINATOR and manifest.get("metadata_eligible_n") == METADATA_ELIGIBLE, "manifest denominator/eligibility mismatch")
        for field, expected in (
            ("methods_config_sha256", FROZEN_CONFIG_SHA256),
            ("amendment_sha256", FROZEN_AMENDMENT_SHA256),
            ("metadata_report_sha256", FROZEN_METADATA_REPORT_SHA256),
            ("cal_map_sha256", FROZEN_CAL_MAP_SHA256),
            ("frozen_census_sha256", FROZEN_CENSUS_SHA256),
            ("stable_matches_sha256", FROZEN_STABLE_MATCHES_SHA256),
            ("preparation_seal_sha256", PRODUCTION_PREPARATION_SEAL_SHA256),
            ("extractor_entry_sha256", PRODUCTION_EXTRACTOR_ENTRY_SHA256),
            ("stream_extractor_sha256", PRODUCTION_STREAM_EXTRACTOR_SHA256),
            ("estimator_sha256", PRODUCTION_ESTIMATOR_SHA256),
            ("gather_sha256", PRODUCTION_GATHER_SHA256),
            ("contract_sha256", PRODUCTION_CONTRACT_SHA256),
            ("wrapper07_sha256", PRODUCTION_WRAPPER07_SHA256),
            ("audit_exception_sha256", PRODUCTION_AUDIT_EXCEPTION_SHA256),
            ("final_custody_closure_sha256", PRODUCTION_FINAL_CUSTODY_CLOSURE_SHA256),
            # builder_sha256 was the renderer's pin on its own bytes, guarding against
            # drift across the eight-module chain. This is now one file, so the pin can
            # only ever fail. The upstream artifact pins below still verify the frozen
            # HPC bundle against the exact production run.
        ):
            _expect(manifest.get(field) == expected, f"manifest {field} mismatch")
        _expect(manifest.get("input_shard_manifest_sha256") == PRODUCTION_MANIFEST_SHA256, "manifest shard-manifest hashes differ from the exact production run")
        _expect(manifest.get("target_access_authorization_sha256") == PRODUCTION_AUTHORIZATION_SHA256, "manifest target-authorisation hashes differ from the exact production run")
        _exact_true(manifest.get("figure4_gate_pass"), "manifest.figure4_gate_pass")
        _exact_true(manifest.get("proportionate_audit_exception_recorded"), "manifest.proportionate_audit_exception_recorded")
        _exact_true(manifest.get("proportionate_final_custody_closure_pass"), "manifest.proportionate_final_custody_closure_pass")
        counts = manifest.get("status_counts")
        _expect(isinstance(counts, dict) and set(counts) == set(FINAL_STATUSES), "manifest status denominator does not use the exact v5.1 status set")
        _expect(all(type(value) is int and value >= 0 for value in counts.values()) and sum(counts.values()) == DENOMINATOR, "manifest status counts do not sum exactly to 4,379")
        _expect(counts["UNANALYSABLE_METADATA_SUPPORT"] == METADATA_UNSUPPORTED, "manifest must retain 1,379 metadata-unsupported edges")
        artifacts = manifest.get("artifacts")
        _expect(isinstance(artifacts, dict) and set(artifacts) == set(ARTIFACT_FILENAMES), "manifest artifact set is not the stable v5.1 gather schema plus audit")
        paths: dict[str, Path] = {}
        hashes: dict[str, str] = {}
        for key, filename in ARTIFACT_FILENAMES.items():
            record = artifacts[key]
            _expect(isinstance(record, dict) and set(record) == {"filename", "sha256"}, f"artifact record {key} is malformed")
            _expect(record.get("filename") == filename, f"artifact {key} does not use the stable production filename")
            path = manifest_path.parent / filename
            digest = _sha(record.get("sha256"), f"artifact {key} hash")
            _expect(path.is_file() and not path.is_symlink()  and sha256_file(path) == digest, f"sealed artifact mismatch: {key}")
            paths[key] = path
            hashes[key] = digest
        _expect(hashes[AUDIT_EXCEPTION_KEY] == PRODUCTION_AUDIT_EXCEPTION_SHA256, "audit-exception artifact is not the exact approved record")
        _expect(hashes[FINAL_CLOSURE_KEY] == PRODUCTION_FINAL_CUSTODY_CLOSURE_SHA256, "final-custody-closure artifact is not the exact approved record")
        return manifest, paths, hashes


    def _validate_gather_report(report: Mapping[str, Any], manifest: Mapping[str, Any], paths: Mapping[str, Path], hashes: Mapping[str, str]) -> None:
        _expect(set(report) == {
            "schema_version", "status", "run_id", "preparation_seal_sha256",
            "shards", "input_shard_manifest_sha256", "target_authorization_sha256",
            "counts", "output_sha256", "inference", "custody",
        }, "gather report field set is not exact")
        _expect(report.get("schema_version") == "response-shape-gather-v5.1", "wrong gather-report schema")
        _expect(report.get("status") == "PASS_CONTINUOUS_RESULTS_NO_INDIVIDUAL_EDGE_INFERENCE", "gather report did not pass")
        _expect(report.get("run_id") == PRODUCTION_RUN_ID == manifest["run_id"], "gather report run_id mismatch")
        _expect(report.get("preparation_seal_sha256") == PRODUCTION_PREPARATION_SEAL_SHA256 == manifest["preparation_seal_sha256"], "gather report preparation-seal mismatch")
        _expect(report.get("shards") == sorted(SHARDS), "gather report does not contain the exact ordered eight shards")
        _expect(report.get("input_shard_manifest_sha256") == PRODUCTION_MANIFEST_SHA256, "gather report shard-manifest hashes mismatch")
        _expect(report.get("target_authorization_sha256") == PRODUCTION_AUTHORIZATION_SHA256, "gather report target-authorisation hashes mismatch")
        _expect(report.get("counts") == {
            "guide_lane_totals": 153_723,
            "positive_target_moments": 514_508_801,
            "unit_effects": 122_024,
            "edge_summaries": 17_970,
            "full_reversal_census": 4_379,
        }, "gather report counts mismatch")
        reported_hashes = report.get("output_sha256", {})
        _expect(isinstance(reported_hashes, dict) and set(reported_hashes) == {
            "guide_lane_totals", "positive_target_moments", "unit_effects",
            "guide_specific_effects", "edge_continuous_summary", "cal_envelopes",
            "reversal_census", "status_waterfall",
        }, "gather report output-hash set mismatch")
        _expect(all(bool(SHA256_RE.fullmatch(str(value))) for value in reported_hashes.values()), "gather report contains an unresolved output hash")
        gather_key_map = {
            "unit_effects": "unit_effects",
            "guide_specific_effects": "guide_specific_effects",
            "edge_continuous_summary": "edge_continuous_summary",
            "cal_envelopes": "cal_envelopes",
            "reversal_census": "reversal_census",
            "status_waterfall": "status_waterfall",
        }
        for artifact_key, report_key in gather_key_map.items():
            _expect(reported_hashes.get(report_key) == hashes[artifact_key], f"gather report hash mismatch for {artifact_key}")
        _expect(report.get("inference") == {
            "individual_edge_p_values": False,
            "cell_level_p_values": False,
            "class_enrichment_hypotheses_or_p_values": False,
            "hard_mechanism_labels": False,
        }, "gather inference boundary mismatch")
        _expect(report.get("custody") == {
            "regulator_expression_values_read": 0,
            "source_shards_hash_validated": True,
            "eight_distinct_shard_authorizations_hash_validated": True,
            "target_expression_values_read": PRODUCTION_TARGET_VALUES_READ,
        }, "gather custody record mismatch")


    AUDIT_EXCEPTION_CLAIM = (
        "Job 417626 is not a PASS, and no raw Slurm spool capture is claimed. This exception covers only expired "
        "live scheduler-control evidence; it does not waive analytical checks, authorize a mechanistic claim, or "
        "establish any biological result."
    )
    FINAL_CLOSURE_CLAIM = (
        "This closure authorizes opening the frozen continuous result artifacts and downstream Figure 4 rendering "
        "under the existing builder gates. It does not authorize individual-edge inference, a mechanistic claim, or "
        "reinterpretation of the recorded scheduler exception."
    )
    FINAL_CLOSURE_CHECKS = {
        "biological_values_not_inspected",
        "durable_sacct_exact",
        "eight_shard_aggregate_custody_exact",
        "gather_report_registry_matches_recomputed_hashes",
        "gather_report_self_hash_stable",
        "gather_stdout_and_stderr_empty",
        "proportionate_operational_exception_transparent",
        "sealed_wrapper07_exact",
        "source_manifests_and_operational_audits_current",
        "source_shards_target_only_and_regulator_zero",
        "tenure_independent_live_scontrol_not_required",
        "three_read_only_inventory_checks_stable",
    }
    REMOTE_GATHER_ROOT = "/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/gather/v5_1_full_20260814_r1"
    REMOTE_WRAPPER07 = "/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/staged/v5_1_prep_20260814_r4/scripts/response_shape/hpc/07_gather_v5_1.sbatch"
    GATHER_SUBMIT_LINE = (
        "sbatch --export=ALL,RESPONSE_V51_PREPARATION_SEAL=/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/preparation/v5_1_prep_20260814_r4/preparation_seal_v5_1.json,"
        "RESPONSE_V51_PREPARATION_SEAL_CHECKSUM=/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/preparation/v5_1_prep_20260814_r4/preparation_seal_v5_1.sha256,"
        "RESPONSE_V51_RUN_ID=v5_1_full_20260814_r1 --cpus-per-task=60 --mem=500G --time=24:00:00 " + REMOTE_WRAPPER07
    )


    def _validate_audit_exception(exception: Mapping[str, Any], manifest: Mapping[str, Any]) -> None:
        _expect(set(exception) == {
            "approval", "claim_boundary", "failed_primary_audit", "generated_at_utc",
            "operational_exception", "optional_archive_attempt", "run_id", "schema_version", "status",
        }, "audit-exception field set is not exact")
        _expect(exception.get("schema_version") == "response-shape-audit-exception-v5.1", "wrong audit-exception schema")
        _expect(exception.get("status") == "RECORDED_PROPORTIONATE_OPERATIONAL_EXCEPTION", "audit exception was not recorded")
        _expect(exception.get("run_id") == PRODUCTION_RUN_ID == manifest["run_id"], "audit-exception run_id mismatch")
        _expect(exception.get("claim_boundary") == AUDIT_EXCEPTION_CLAIM, "audit-exception claim boundary falsely describes the failed audit or raw spool")
        _expect("10c983cdea0de566" not in json.dumps(exception, sort_keys=True), "audit exception falsely claims the unrun long validator")
        approval = exception.get("approval")
        _expect(approval == {
            "approved_by_user": True,
            "approved_route": "proportionate_four_step_closure",
            "does_not_convert_failed_audit_to_pass": True,
            "does_not_claim_raw_spool_capture": True,
            "long_audit_rerun_required": False,
            "scientific_checks_waived": False,
            "separate_final_custody_closure_required": True,
        }, "audit-exception approval boundary mismatch")
        failed = exception.get("failed_primary_audit")
        _expect(isinstance(failed, dict) and set(failed) == {
            "biological_values_emitted", "elapsed", "exit_code", "failure_gate",
            "failure_occurred_at_identical_final_scheduler_recheck",
            "first_scheduler_and_stored_wrapper_gate_passed_before_table_read", "job_id", "job_name",
            "report_emitted", "report_path", "state", "stderr_bytes", "stderr_sha256",
            "stdout_bytes", "stdout_sha256", "substantive_checks_and_output_lock_reached_before_failure",
            "validator_path", "validator_sha256", "wrapper_path", "wrapper_sha256",
        }, "failed-audit exception record is not exact")
        _expect(
            failed["job_id"] == "417626" and failed["job_name"] == "rshape51_postaudit"
            and failed["state"] == "FAILED" and failed["exit_code"] == "2:0"
            and failed["failure_gate"] == "Slurm control record unavailable",
            "audit exception falsely converts job 417626 to PASS",
        )
        for field in (
            "failure_occurred_at_identical_final_scheduler_recheck",
            "first_scheduler_and_stored_wrapper_gate_passed_before_table_read",
            "substantive_checks_and_output_lock_reached_before_failure",
        ):
            _exact_true(failed[field], f"audit_exception.failed_primary_audit.{field}")
        _exact_false(failed["report_emitted"], "audit_exception.failed_primary_audit.report_emitted")
        _exact_false(failed["biological_values_emitted"], "audit_exception.failed_primary_audit.biological_values_emitted")
        _expect(failed["validator_sha256"] == "8fe0457771b8f1c78347dcece8f728164a56c8f0fba80c186ccd02ac041e88b1", "audit-exception failed-validator pin mismatch")
        _expect(failed["wrapper_sha256"] == "644e078c779fb15839547a7802077d9b4b72adac71268f6aa7c26adcfad242a0", "audit-exception failed-wrapper pin mismatch")
        _expect(failed["stdout_bytes"] == 201 and failed["stdout_sha256"] == "83e17629b7df4864f36efd84efedd60a9e003395b962324f9aff83709db8ea6b", "audit-exception failed stdout mismatch")
        _expect(failed["stderr_bytes"] == 0 and failed["stderr_sha256"] == hashlib.sha256(b"").hexdigest(), "audit-exception failed stderr mismatch")
        operational = exception.get("operational_exception")
        _expect(isinstance(operational, dict) and operational == {
            "durable_sacct_record_available": True,
            "final_live_scontrol_record_available": False,
            "raw_slurm_batch_payload_persisted": False,
            "scope": "expired_live_slurm_control_and_raw_spool_only",
            "sealed_submission_wrapper_bytes": 3740,
            "sealed_submission_wrapper_path": REMOTE_WRAPPER07,
            "sealed_submission_wrapper_sha256": PRODUCTION_WRAPPER07_SHA256,
        }, "audit-exception operational scope or raw-spool boundary mismatch")
        optional = exception.get("optional_archive_attempt")
        _expect(isinstance(optional, dict) and set(optional) == {
            "biological_outputs_opened", "elapsed", "exit_code", "final_archive_created",
            "job_id", "state", "status_line_sha256",
        }, "optional archive-attempt record is not exact")
        _expect(optional["job_id"] == "417643" and optional["state"] == "FAILED" and optional["exit_code"] == "2:0", "optional archive attempt state mismatch")
        _exact_false(optional["final_archive_created"], "audit_exception.optional_archive_attempt.final_archive_created")
        _exact_false(optional["biological_outputs_opened"], "audit_exception.optional_archive_attempt.biological_outputs_opened")


    def _validate_final_custody_closure(
        closure: Mapping[str, Any],
        exception: Mapping[str, Any],
        manifest: Mapping[str, Any],
        paths: Mapping[str, Path],
        hashes: Mapping[str, str],
        report: Mapping[str, Any],
    ) -> None:
        _expect(set(closure) == {
            "audit_exception", "checks", "claim_boundary", "gather_artifacts",
            "gather_report_registry", "generated_at_utc", "outcome_boundary", "run_id",
            "scheduler_evidence", "schema_version", "source_custody", "status",
        }, "final-custody-closure field set is not exact")
        _expect(closure.get("schema_version") == "response-shape-final-gather-custody-closure-v5.1", "wrong final-custody-closure schema")
        _expect(closure.get("status") == "GO_TO_OPEN_RESULTS_PROPORTIONATE_AUDIT", "proportionate final custody closure did not issue GO")
        _expect(closure.get("run_id") == PRODUCTION_RUN_ID == manifest["run_id"], "final-custody-closure run_id mismatch")
        _expect(closure.get("claim_boundary") == FINAL_CLOSURE_CLAIM, "final-custody-closure claim boundary mismatch")
        _expect("10c983cdea0de566" not in json.dumps(closure, sort_keys=True), "final custody closure falsely claims the unrun long validator")
        exception_ref = closure.get("audit_exception")
        _expect(isinstance(exception_ref, dict) and set(exception_ref) == {"path", "sha256", "status"}, "closure audit-exception reference is not exact")
        _expect(Path(str(exception_ref["path"])).name == AUDIT_EXCEPTION_FILENAME, "closure references the wrong audit-exception file")
        _expect(exception_ref["sha256"] == hashes[AUDIT_EXCEPTION_KEY] == PRODUCTION_AUDIT_EXCEPTION_SHA256, "closure audit-exception hash mismatch")
        _expect(exception_ref["status"] == exception["status"] == "RECORDED_PROPORTIONATE_OPERATIONAL_EXCEPTION", "closure audit-exception status mismatch")
        checks = closure.get("checks")
        _expect(isinstance(checks, dict) and set(checks) == FINAL_CLOSURE_CHECKS, "final-custody-closure check set is not exact")
        for check in FINAL_CLOSURE_CHECKS:
            _exact_true(checks[check], f"final_custody_closure.checks.{check}")
        _expect(closure.get("outcome_boundary") == {
            "biological_values_opened_during_closure": False,
            "directions_or_magnitudes_emitted": False,
            "edge_ids_emitted": False,
            "individual_edge_inference_authorized": False,
            "mechanistic_claim_authorized": False,
            "result_opening_authorized": True,
        }, "final-custody-closure outcome boundary mismatch")

        gather = closure.get("gather_artifacts")
        _expect(isinstance(gather, dict) and set(gather) == {
            "all_nonzero", "all_regular", "all_symlink_free", "directory_mode",
            "exact_inventory_count", "path", "artifacts",
        }, "closure gather-artifact record is not exact")
        for field in ("all_nonzero", "all_regular", "all_symlink_free"):
            _exact_true(gather[field], f"final_custody_closure.gather_artifacts.{field}")
        _expect(gather["directory_mode"] == "0550" and gather["exact_inventory_count"] == 9 and gather["path"] == REMOTE_GATHER_ROOT, "closure gather root/inventory mismatch")
        artifact_rows = gather["artifacts"]
        _expect(isinstance(artifact_rows, list) and len(artifact_rows) == 9, "closure must retain the exact nine-file gather inventory")
        by_name = {row.get("name"): row for row in artifact_rows if isinstance(row, dict)}
        expected_names = set(GATHERED_FILENAMES.values()) | {"guide_lane_totals_v5_1.parquet", "positive_target_moments_v5_1.parquet"}
        _expect(set(by_name) == expected_names and len(by_name) == len(artifact_rows), "closure gather filenames are not exact")
        report_key_by_filename = {
            "unit_effects_v5_1.parquet": "unit_effects",
            "guide_specific_effects_v5_1.parquet": "guide_specific_effects",
            "edge_continuous_summary_v5_1.parquet": "edge_continuous_summary",
            "cal_envelopes_v5_1.tsv": "cal_envelopes",
            "reversal_census_v5_1.parquet": "reversal_census",
            "status_waterfall_v5_1.tsv": "status_waterfall",
            "guide_lane_totals_v5_1.parquet": "guide_lane_totals",
            "positive_target_moments_v5_1.parquet": "positive_target_moments",
        }
        artifact_key_by_filename = {filename: key for key, filename in GATHERED_FILENAMES.items()}
        for filename, row in by_name.items():
            _expect(set(row) == {"mode", "name", "path", "sha256", "size_bytes"}, f"closure gather record malformed: {filename}")
            _expect(row["mode"] == "0440" and row["path"] == f"{REMOTE_GATHER_ROOT}/{filename}" and type(row["size_bytes"]) is int and row["size_bytes"] > 0, f"closure gather file custody mismatch: {filename}")
            if filename == "gather_report_v5_1.json":
                expected_hash = hashes["gather_report"]
            else:
                expected_hash = report["output_sha256"][report_key_by_filename[filename]]
            _expect(row["sha256"] == expected_hash, f"closure gather hash mismatch: {filename}")
            artifact_key = artifact_key_by_filename.get(filename)
            if artifact_key is not None:
                _expect(row["size_bytes"] == paths[artifact_key].stat().st_size, f"closure gather size mismatch: {filename}")
        _expect(by_name["guide_lane_totals_v5_1.parquet"]["size_bytes"] == 1_087_197, "closure guide-lane sparse size mismatch")
        _expect(by_name["positive_target_moments_v5_1.parquet"]["size_bytes"] == 7_576_181_594, "closure positive-moment sparse size mismatch")

        registry = closure.get("gather_report_registry")
        _expect(isinstance(registry, dict) and set(registry) == {
            "all_registered_hashes_match", "registered_artifact_count", "registered_output_sha256",
            "report_run_id", "report_schema_version", "report_self_sha256", "report_status",
        }, "closure gather-report registry is not exact")
        _exact_true(registry["all_registered_hashes_match"], "final_custody_closure.gather_report_registry.all_registered_hashes_match")
        _expect(registry["registered_artifact_count"] == 8 and registry["registered_output_sha256"] == report["output_sha256"], "closure gather-report output registry mismatch")
        _expect(registry["report_self_sha256"] == hashes["gather_report"] and registry["report_run_id"] == PRODUCTION_RUN_ID, "closure gather-report identity mismatch")
        _expect(registry["report_schema_version"] == report["schema_version"] and registry["report_status"] == report["status"], "closure gather-report status/schema mismatch")

        scheduler = closure.get("scheduler_evidence")
        _expect(isinstance(scheduler, dict) and set(scheduler) == {
            "allocated_cpus", "durable_submit_line", "elapsed", "end", "evidence_basis", "exit_code",
            "job_id", "job_name", "logs", "node", "partition", "requested_memory", "requested_time",
            "start", "state", "submit", "work_directory", "wrapper",
        }, "closure scheduler evidence is not exact")
        _expect(
            scheduler["job_id"] == "417604" and scheduler["job_name"] == "rshape51_gather"
            and scheduler["state"] == "COMPLETED" and scheduler["exit_code"] == "0:0"
            and scheduler["allocated_cpus"] == 60 and scheduler["requested_memory"] == "500G"
            and scheduler["requested_time"] == "1-00:00:00" and scheduler["durable_submit_line"] == GATHER_SUBMIT_LINE,
            "closure durable gather scheduler record mismatch",
        )
        _expect(scheduler["evidence_basis"] == "durable_sacct_plus_sealed_submission_wrapper_and_empty_logs", "closure scheduler evidence overstates raw control/spool custody")
        wrapper = scheduler["wrapper"]
        _expect(wrapper == {
            "mode": "0440", "path": REMOTE_WRAPPER07, "regular": True,
            "sha256": PRODUCTION_WRAPPER07_SHA256, "size_bytes": 3740, "symlink": False,
        }, "closure sealed wrapper07 evidence mismatch")
        logs = scheduler["logs"]
        _expect(isinstance(logs, list) and len(logs) == 2, "closure gather log evidence is not exact")
        _expect({Path(str(row.get("path"))).suffix for row in logs if isinstance(row, dict)} == {".out", ".err"}, "closure gather log paths are not exact")
        for row in logs:
            _expect(set(row) == {"mode", "path", "sha256", "size_bytes"} and row["mode"] == "0644" and row["size_bytes"] == 0 and row["sha256"] == hashlib.sha256(b"").hexdigest(), "closure gather logs are not empty and sealed")

        source = closure.get("source_custody")
        _expect(isinstance(source, dict) and set(source) == {
            "aggregate", "array", "closed_bundle", "custody_report", "current_remote_hashes_match",
            "operational_audit_sha256", "remote_audit_root", "remote_full_root",
            "requested_target_ids_sha256", "requested_targets_per_shard", "shard_manifest_sha256",
            "target_authorization_verified_for_all_shards", "target_only_reads_for_all_shards",
            "total_regulator_expression_values_read",
        }, "closure source-custody field set is not exact")
        _expect(source["aggregate"] == {
            "guide_lane_rows": 153_723,
            "implicit_cell_zeros": 1_766_457_625,
            "implicit_group_zeros": 152_495_296,
            "positive_target_moment_rows": 514_508_801,
            "selected_cells": PRODUCTION_SELECTED_CELLS,
            "target_expression_values_read": PRODUCTION_TARGET_VALUES_READ,
        }, "closure source aggregate custody mismatch")
        _expect(source["array"] == {
            "completed_tasks": 8, "failed_tasks": 0, "master_job_id": "417546",
            "requested_cpus_per_task": 60, "requested_memory": "500G",
            "stderr_nonempty_tasks": 0, "stdout_nonempty_tasks": 0,
        }, "closure source array custody mismatch")
        for field in ("current_remote_hashes_match", "target_authorization_verified_for_all_shards", "target_only_reads_for_all_shards"):
            _exact_true(source[field], f"final_custody_closure.source_custody.{field}")
        _expect(source["requested_target_ids_sha256"] == PRODUCTION_TARGET_UNIVERSE_SHA256 and source["requested_targets_per_shard"] == 4_339, "closure requested-target custody mismatch")
        _expect(source["total_regulator_expression_values_read"] == 0, "closure source custody records regulator expression reads")
        _expect(source["shard_manifest_sha256"] == PRODUCTION_MANIFEST_SHA256, "closure shard-manifest hashes mismatch")
        operational_hashes = source["operational_audit_sha256"]
        _expect(isinstance(operational_hashes, dict) and set(operational_hashes) == set(SHARDS), "closure operational-audit hashes do not cover eight shards")
        for shard, digest in operational_hashes.items():
            _sha(digest, f"closure operational audit {shard}")
        _expect(source["closed_bundle"]["sha256"] == "89ed14c29e8e33fd99cf8517648866254177f7066b14c1ef2d4e1a67a57d68d8" and source["closed_bundle"]["verification_passed"] is True, "closure closed-bundle custody mismatch")
        _expect(source["custody_report"]["sha256"] == "7c75a229e4eac831051316f2962fd6639b6f515c034b867e9359d74ea2fb55d1" and source["custody_report"]["status"] == "PASS_8_SHARDS_AUDITED_AND_LOCKED", "closure eight-shard custody report mismatch")
        _expect(source["remote_audit_root"] == {"mode": "0550", "path": "/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/audits/v5_1_full_20260814_r1"}, "closure remote audit-root custody mismatch")
        _expect(source["remote_full_root"] == {"mode": "0550", "path": "/home/hpc-user/projects/codeg_atlas/work/response_shape/v5_1/full/v5_1_full_20260814_r1"}, "closure remote full-root custody mismatch")


    def _validate_census(census: pd.DataFrame, waterfall: pd.DataFrame, frozen: FrozenContext, manifest: Mapping[str, Any]) -> None:
        _require_columns(census, CENSUS_REQUIRED, "reversal census")
        _expect(len(census) == DENOMINATOR and not census["reversal_edge_id"].duplicated().any(), "reversal census must retain exactly 4,379 unique rows")
        _expect(census["reversal_edge_id"].astype(str).map(lambda value: bool(EDGE_RE.fullmatch(value))).all(), "reversal census has malformed edge IDs")
        expected = frozen.reversal.set_index("reversal_edge_id")
        observed = census.set_index("reversal_edge_id")
        _expect(set(observed.index) == set(expected.index), "production census differs from the frozen edge universe")
        for out_col, source_col in (("regulator", "regulator"), ("regulator_ensg", "regulator_ensg"), ("target", "target"), ("target_ensg", "target_ensg"), ("direction", "direction"), ("expected_sign_Rest", "expected_sign_Rest"), ("expected_sign_Stim48hr", "expected_sign_Stim48hr")):
            _expect(observed.loc[expected.index, out_col].astype(str).eq(expected[source_col].astype(str)).all(), f"reversal census {out_col} differs from frozen identity")
        _all_bool(census["metadata_eligible"], "census metadata_eligible")
        expected_eligible = census["reversal_edge_id"].isin(frozen.eligible_edge_ids)
        _expect(census["metadata_eligible"].astype(bool).eq(expected_eligible).all(), "census metadata eligibility differs from frozen 149-regulator universe")
        _expect(int(expected_eligible.sum()) == METADATA_ELIGIBLE, "census metadata-eligible count is not 3,000")
        _expect(set(census["analysis_status"]).issubset(FINAL_STATUSES), "census has an unknown final status")
        _expect(census.loc[~expected_eligible, "analysis_status"].eq("UNANALYSABLE_METADATA_SUPPORT").all(), "metadata-unsupported edges were not retained explicitly")
        _expect(census.loc[expected_eligible, "analysis_status"].ne("UNANALYSABLE_METADATA_SUPPORT").all(), "eligible edge was assigned metadata failure")
        counts = {status: int(census["analysis_status"].eq(status).sum()) for status in FINAL_STATUSES}
        _expect(counts == manifest["status_counts"], "manifest status denominator differs from reversal census")
        _require_columns(waterfall, {"analysis_status", "edges"}, "status waterfall")
        _expect(not waterfall["analysis_status"].duplicated().any() and set(waterfall["analysis_status"]).issubset(FINAL_STATUSES), "status waterfall has repeated or unknown statuses")
        waterfall_counts = {status: 0 for status in FINAL_STATUSES}
        for row in waterfall.itertuples(index=False):
            _expect(type(row.edges) in (int, np.int64, np.int32) or float(row.edges).is_integer(), "status waterfall counts must be integers")
            waterfall_counts[str(row.analysis_status)] = int(row.edges)
        _expect(waterfall_counts == counts and sum(waterfall_counts.values()) == DENOMINATOR, "status waterfall is not the exact 4,379-edge census")
        by_id = census.set_index("reversal_edge_id")
        for edge_id, regulator, target, direction, sign in PREDECLARED_EDGES:
            _expect(edge_id in by_id.index, f"predeclared edge is missing: {edge_id}")
            row = by_id.loc[edge_id]
            _expect(row["regulator"] == regulator and row["target"] == target and row["direction"] == direction and int(row["expected_sign_Stim48hr"]) == sign, f"predeclared edge identity drift: {edge_id}")


    def _expected_role_identity(frozen: FrozenContext, role: str) -> frozenset[str]:
        return frozen.expected_analysis_ids[role]


    def _validate_summary(summary: pd.DataFrame, census: pd.DataFrame, frozen: FrozenContext) -> None:
        _require_columns(summary, SUMMARY_REQUIRED, "edge continuous summary")
        _expect(not summary["analysis_edge_id"].astype(str).duplicated().any(), "edge continuous summary repeats analysis IDs")
        _expect(set(summary["analysis_role"]).issubset(ANALYSIS_ROLES), "edge summary contains an unapproved analysis role")
        for role in ANALYSIS_ROLES:
            observed = frozenset(summary.loc[summary["analysis_role"].eq(role), "analysis_edge_id"].astype(str))
            _expect(observed == _expected_role_identity(frozen, role), f"edge-summary identity set drift for {role}")
        _expect(set(summary["summary_status"]).issubset(SUMMARY_STATUSES), "edge summary has unknown summary status")
        _expect(set(summary["analysis_status"]).issubset(set(SUMMARY_STATUSES) | {"MIXED_OR_UNRESOLVED"}), "edge summary has unknown analysis status")
        _expect(summary.loc[summary["summary_status"].ne("ESTIMATED_CONTINUOUS_COMPONENTS"), "analysis_status"].eq(summary.loc[summary["summary_status"].ne("ESTIMATED_CONTINUOUS_COMPONENTS"), "summary_status"]).all(), "unsupported summary status changed downstream")
        _expect(summary.loc[summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS"), "analysis_status"].isin(NUMERICAL_FINAL).all(), "estimated summary has invalid final status")
        donor_counts = _integers(summary["n_common_donors"], "edge-summary n_common_donors")
        donor_sets = summary["common_donor_ids"].map(_canonical_donors)
        _expect(np.array_equal(donor_counts.to_numpy(), donor_sets.map(len).to_numpy()), "edge-summary common donor IDs/counts disagree")
        supported = summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        _expect((donor_counts.loc[supported] >= 3).all(), "estimated edge summary has fewer than three common donors")
        _expect((donor_counts.loc[summary["summary_status"].eq("UNANALYSABLE_POSITIVE_INTENSITY")] < 3).all(), "positive-intensity failure claims three common donors")
        orientation = pd.to_numeric(summary.loc[supported, "orientation"], errors="coerce")
        _expect(orientation.isin((-1, 1)).all(), "supported edge summary has invalid orientation")
        for metric in METRICS:
            rest = pd.to_numeric(summary.loc[supported, f"Rest__{metric}"], errors="coerce")
            stim = pd.to_numeric(summary.loc[supported, f"Stim48hr__{metric}"], errors="coerce")
            delta = pd.to_numeric(summary.loc[supported, f"delta__{metric}"], errors="coerce")
            oriented = pd.to_numeric(summary.loc[supported, f"oriented_delta__{metric}"], errors="coerce")
            _expect(np.isfinite(rest).all() and np.isfinite(stim).all() and np.isfinite(delta).all() and np.isfinite(oriented).all(), f"supported {metric} summary contains non-finite values")
            _expect(np.allclose(delta, stim - rest, atol=1e-12, rtol=1e-10), f"{metric} delta is not Stim48hr-Rest")
            _expect(np.allclose(oriented, orientation.to_numpy() * delta.to_numpy(), atol=1e-12, rtol=1e-10), f"{metric} oriented delta is not o*delta")
        for condition in CONDITIONS:
            total = summary.loc[supported, f"{condition}__total_effect"].to_numpy(float)
            detection = summary.loc[supported, f"{condition}__detection_component"].to_numpy(float)
            intensity = summary.loc[supported, f"{condition}__positive_intensity_component"].to_numpy(float)
            tolerance = 1e-10 + 1e-12 * np.maximum.reduce((np.abs(total), np.abs(detection), np.abs(intensity)))
            _expect(np.all(np.abs(total - detection - intensity) <= tolerance), f"{condition} edge-level Shapley additivity failed")
        total_delta = summary.loc[supported, "delta__total_effect"].to_numpy(float)
        component_delta = summary.loc[supported, "delta__detection_component"].to_numpy(float) + summary.loc[supported, "delta__positive_intensity_component"].to_numpy(float)
        _expect(np.allclose(total_delta, component_delta, atol=1e-10, rtol=1e-10), "cross-condition Shapley components do not add to total")

        reversal = summary.loc[summary["analysis_role"].eq("REVERSAL")].set_index("analysis_edge_id")
        census_by_id = census.set_index("reversal_edge_id")
        _expect(reversal["analysis_status"].eq(census_by_id.loc[reversal.index, "analysis_status"]).all(), "reversal summary status differs from full census")
        _expect(reversal.index.to_series().eq(reversal["parent_reversal_edge_id"].astype(str)).all(), "reversal analysis ID must equal its parent edge ID")
        cal = summary.loc[summary["analysis_role"].eq("CAL_REVERSAL")]
        _expect(cal["analysis_edge_id"].eq("CAL_REV|" + cal["parent_reversal_edge_id"].astype(str)).all(), "CAL reversal ID drift")
        stable = summary.loc[summary["analysis_role"].eq("REGULATOR_MATCHED_STABLE")]
        _expect(stable["analysis_edge_id"].str.contains(r"^STABLE\|ENSG[0-9]+\|ENSG[0-9]+\|PARENT\|REV\|", regex=True).all(), "stable summary contains a non-regulator-matched identity")
        cal_stable = summary.loc[summary["analysis_role"].eq("CAL_STABLE")]
        _expect(cal_stable["analysis_edge_id"].str.startswith("CAL_STABLE|").all(), "CAL stable identity drift")


    def _validate_units(units: pd.DataFrame, guides: pd.DataFrame, summary: pd.DataFrame, frozen: FrozenContext) -> None:
        _require_columns(units, UNIT_REQUIRED, "unit effects")
        _require_columns(guides, GUIDE_REQUIRED, "guide-specific effects")
        _expect(not units.duplicated(["analysis_edge_id", "donor", "condition"]).any(), "unit effects repeat edge/donor/condition")
        _expect(set(units["analysis_role"]).issubset(ANALYSIS_ROLES), "unit effects contain an unknown analysis role")
        _expect(set(units["unit_status"]).issubset(UNIT_STATUSES), "unit effects contain an unknown status")
        _expect(units["guide_A_guide_id"].astype(str).ne(units["guide_B_guide_id"].astype(str)).all(), "unit A/B physical guides are identical")
        _expect(units["condition"].isin(CONDITIONS).all() and units["donor"].isin(DONORS).all(), "unit effects contain invalid donor/condition")

        donor_by_reg = frozen.cal_map.set_index("regulator_ensg")["paired_donor_set"].map(_canonical_donors)
        for role in ANALYSIS_ROLES:
            observed_ids = frozenset(units.loc[units["analysis_role"].eq(role), "analysis_edge_id"].astype(str))
            _expect(observed_ids == _expected_role_identity(frozen, role), f"unit-effect identity set drift for {role}")
        summary_identity = summary.set_index("analysis_edge_id")
        expected_unit_keys: set[tuple[str, str, str]] = set()
        for edge_id, row in summary_identity.iterrows():
            regulator = str(row["regulator_ensg"])
            _expect(regulator in donor_by_reg.index, f"unit regulator lacks a frozen donor set: {edge_id}")
            for donor in donor_by_reg.loc[regulator]:
                for condition in CONDITIONS:
                    expected_unit_keys.add((str(edge_id), donor, condition))
        observed_unit_keys = set(zip(units["analysis_edge_id"].astype(str), units["donor"].astype(str), units["condition"].astype(str)))
        _expect(observed_unit_keys == expected_unit_keys and len(units) == len(expected_unit_keys), "all-role unit grid differs from frozen analysis IDs/donors/conditions")

        map_by_reg = frozen.cal_map.set_index("regulator_ensg")
        for row in units.itertuples(index=False):
            is_cal = str(row.analysis_role).startswith("CAL_")
            expected_a = map_by_reg.loc[str(row.regulator_ensg), "cal_guide_A" if is_cal else "case_guide_A"]
            expected_b = map_by_reg.loc[str(row.regulator_ensg), "cal_guide_B" if is_cal else "case_guide_B"]
            _expect(str(row.guide_A_guide_id) == str(expected_a) and str(row.guide_B_guide_id) == str(expected_b), f"physical guide identity drift for {row.analysis_edge_id}")
        estimated = units["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        for slot in ("A", "B"):
            _expect((pd.to_numeric(units.loc[estimated, f"guide_{slot}_n_detected"]) >= 10).all(), f"estimated unit guide {slot} misses detected-cell floor")
            _expect((pd.to_numeric(units.loc[estimated, f"guide_{slot}_ref_raw_detected_cells"]) >= 10).all(), f"estimated unit REF arm {slot} misses raw-detection floor")
            _expect((pd.to_numeric(units.loc[estimated, f"guide_{slot}_ref_distinct_detecting_physical_guides"]) >= 10).all(), f"estimated unit REF arm {slot} misses physical-guide floor")
        numeric = units.loc[estimated]
        for prefix in ("p_case", "q_case", "p_ref", "q_ref", "mu_case", "mu_ref", "total_effect", "detection_component", "positive_intensity_component", "additivity_error", "additivity_tolerance", "raw_rate_log2_ratio"):
            _finite(numeric[prefix], f"estimated unit {prefix}")
        for field in ("p_case", "p_ref"):
            values = numeric[field].to_numpy(float)
            _expect(((values > 0) & (values <= 1)).all(), f"estimated unit {field} outside (0,1]")
        _expect(np.allclose(numeric["mu_case"], numeric["q_case"] / numeric["p_case"], atol=1e-12, rtol=1e-10), "unit mu_case != q_case/p_case")
        _expect(np.allclose(numeric["mu_ref"], numeric["q_ref"] / numeric["p_ref"], atol=1e-12, rtol=1e-10), "unit mu_ref != q_ref/p_ref")
        for field in ("p_case", "q_case", "p_ref", "q_ref", "raw_rate_log2_ratio"):
            expected = (numeric[f"guide_A_{field}"].to_numpy(float) + numeric[f"guide_B_{field}"].to_numpy(float)) / 2
            _expect(np.allclose(numeric[field], expected, atol=1e-12, rtol=1e-10), f"unit {field} is not the equal A/B average")
        p1, p0 = numeric["p_case"].to_numpy(float), numeric["p_ref"].to_numpy(float)
        m1, m0 = numeric["mu_case"].to_numpy(float), numeric["mu_ref"].to_numpy(float)
        expected_total = numeric["q_case"].to_numpy(float) - numeric["q_ref"].to_numpy(float)
        expected_detection = (p1 - p0) * (m1 + m0) / 2
        expected_intensity = (m1 - m0) * (p1 + p0) / 2
        _expect(np.allclose(numeric["total_effect"], expected_total, atol=1e-12, rtol=1e-10), "unit total effect is not q_case-q_ref")
        _expect(np.allclose(numeric["detection_component"], expected_detection, atol=1e-12, rtol=1e-10), "unit detection component violates the frozen Shapley formula")
        _expect(np.allclose(numeric["positive_intensity_component"], expected_intensity, atol=1e-12, rtol=1e-10), "unit intensity component violates the frozen Shapley formula")
        _all_bool(numeric["additivity_pass"], "unit additivity_pass")
        _expect(numeric["additivity_pass"].astype(bool).all(), "estimated unit failed numerical additivity")

        _expect(not guides.duplicated(["analysis_edge_id", "donor", "condition", "guide_slot"]).any(), "guide-specific effects repeat a slot")
        _expect(set(guides["guide_slot"]) == {"A", "B"}, "guide-specific table must use A/B slots only")
        expected_guide_keys = {(str(row.analysis_edge_id), str(row.donor), str(row.condition), slot) for row in units.itertuples(index=False) for slot in ("A", "B")}
        observed_guide_keys = set(zip(guides["analysis_edge_id"].astype(str), guides["donor"].astype(str), guides["condition"].astype(str), guides["guide_slot"].astype(str)))
        _expect(observed_guide_keys == expected_guide_keys and len(guides) == len(expected_guide_keys), "guide-specific table is not exactly two rows per unit")
        estimated_guides = guides.loc[guides["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        guide_numeric = {
            field: _finite(estimated_guides[field], f"estimated guide {field}")
            for field in (
                "p_case", "q_case", "p_ref", "q_ref", "guide_mu_case", "guide_mu_ref",
                "guide_total_effect", "guide_detection_component",
                "guide_positive_intensity_component", "guide_additivity_error",
                "raw_rate_log2_ratio",
            )
        }
        for field in ("p_case", "p_ref"):
            _expect(((guide_numeric[field] > 0) & (guide_numeric[field] <= 1)).all(), f"estimated guide {field} outside (0,1]")
        _expect(np.allclose(guide_numeric["guide_mu_case"], guide_numeric["q_case"] / guide_numeric["p_case"], atol=1e-12, rtol=1e-10), "guide mu_case != q_case/p_case")
        _expect(np.allclose(guide_numeric["guide_mu_ref"], guide_numeric["q_ref"] / guide_numeric["p_ref"], atol=1e-12, rtol=1e-10), "guide mu_ref != q_ref/p_ref")
        guide_total = guide_numeric["q_case"] - guide_numeric["q_ref"]
        guide_detection = (guide_numeric["p_case"] - guide_numeric["p_ref"]) * (guide_numeric["guide_mu_case"] + guide_numeric["guide_mu_ref"]) / 2
        guide_intensity = (guide_numeric["guide_mu_case"] - guide_numeric["guide_mu_ref"]) * (guide_numeric["p_case"] + guide_numeric["p_ref"]) / 2
        _expect(np.allclose(guide_numeric["guide_total_effect"], guide_total, atol=1e-12, rtol=1e-10), "guide total effect is not q_case-q_ref")
        _expect(np.allclose(guide_numeric["guide_detection_component"], guide_detection, atol=1e-12, rtol=1e-10), "guide detection component violates the frozen Shapley formula")
        _expect(np.allclose(guide_numeric["guide_positive_intensity_component"], guide_intensity, atol=1e-12, rtol=1e-10), "guide intensity component violates the frozen Shapley formula")
        _all_bool(estimated_guides["guide_additivity_pass"], "guide additivity_pass")
        _expect(estimated_guides["guide_additivity_pass"].astype(bool).all(), "estimated guide failed numerical additivity")
        guide_link_fields = (
            "guide_id", "n_cells", "n_detected", "p_case", "q_case", "p_ref", "q_ref",
            "raw_rate_log2_ratio", "ref_raw_detected_cells",
            "ref_distinct_detecting_physical_guides", "guide_mu_case", "guide_mu_ref",
            "guide_total_effect", "guide_detection_component",
            "guide_positive_intensity_component", "guide_additivity_error", "guide_additivity_pass",
        )
        identity_fields = ("parent_reversal_edge_id", "analysis_role", "regulator_ensg", "target_ensg", "unit_status")
        unit_long_parts = []
        for slot in ("A", "B"):
            columns = ["analysis_edge_id", "donor", "condition", *identity_fields, *(f"guide_{slot}_{field}" for field in guide_link_fields)]
            part = units.loc[:, columns].rename(columns={f"guide_{slot}_{field}": field for field in guide_link_fields})
            part["guide_slot"] = slot
            unit_long_parts.append(part)
        unit_long = pd.concat(unit_long_parts, ignore_index=True)
        keys = ["analysis_edge_id", "donor", "condition", "guide_slot"]
        linked = guides.loc[:, [*keys, *identity_fields, *guide_link_fields]].merge(
            unit_long,
            on=keys,
            how="inner",
            validate="one_to_one",
            suffixes=("__long", "__unit"),
        )
        _expect(len(linked) == len(guides) == len(unit_long), "long-guide/unit-wide linkage is incomplete")
        for field in (*identity_fields, "guide_id"):
            _expect(linked[f"{field}__long"].fillna("<NA>").astype(str).eq(linked[f"{field}__unit"].fillna("<NA>").astype(str)).all(), f"guide audit {field} differs from unit output")
        for field in guide_link_fields:
            if field in {"guide_id", "guide_additivity_pass"}:
                continue
            left = pd.to_numeric(linked[f"{field}__long"], errors="coerce").to_numpy(float)
            right = pd.to_numeric(linked[f"{field}__unit"], errors="coerce").to_numpy(float)
            _expect(np.isclose(left, right, atol=1e-12, rtol=1e-10, equal_nan=True).all(), f"guide audit {field} differs from unit output")
        left_bool, right_bool = linked["guide_additivity_pass__long"], linked["guide_additivity_pass__unit"]
        both_missing = left_bool.isna() & right_bool.isna()
        both_boolean_equal = left_bool.map(lambda value: isinstance(value, (bool, np.bool_))) & right_bool.map(lambda value: isinstance(value, (bool, np.bool_))) & left_bool.eq(right_bool)
        _expect((both_missing | both_boolean_equal).all(), "guide audit guide_additivity_pass differs from unit output")

        summary_index = summary.set_index("analysis_edge_id")
        guide_by_edge = guides.set_index("analysis_edge_id", drop=False)
        for edge_id, block in units.groupby("analysis_edge_id", sort=False):
            out = summary_index.loc[edge_id]
            for field in ("parent_reversal_edge_id", "analysis_role", "regulator_ensg", "target_ensg"):
                observed_identity = block[field].astype(str)
                _expect(observed_identity.nunique() == 1 and observed_identity.iloc[0] == str(out[field]), f"unit/summary {field} mismatch for {edge_id}")
            passing = block.loc[block["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
            common = tuple(donor for donor in DONORS if donor in set(passing.loc[passing["condition"].eq("Rest"), "donor"]).intersection(set(passing.loc[passing["condition"].eq("Stim48hr"), "donor"])))
            _expect(_canonical_donors(out["common_donor_ids"]) == common and int(out["n_common_donors"]) == len(common), f"common donor identity mismatch for {edge_id}")
            quarantine = block["unit_status"].eq("QUARANTINED_NUMERICAL_ADDITIVITY").any()
            expected_status = "QUARANTINED_NUMERICAL_ADDITIVITY" if quarantine else "ESTIMATED_CONTINUOUS_COMPONENTS" if len(common) >= 3 else "UNANALYSABLE_POSITIVE_INTENSITY"
            _expect(out["summary_status"] == expected_status, f"summary support status mismatch for {edge_id}")
            if expected_status == "ESTIMATED_CONTINUOUS_COMPONENTS":
                selected = passing.loc[passing["donor"].isin(common)]
                for condition in CONDITIONS:
                    condition_rows = selected.loc[selected["condition"].eq(condition)]
                    _expect(set(condition_rows["donor"]) == set(common), f"condition donor set mismatch for {edge_id}/{condition}")
                    for metric in METRICS:
                        _expect(np.isclose(float(out[f"{condition}__{metric}"]), float(condition_rows[metric].mean()), atol=1e-12, rtol=1e-10), f"edge summary {metric} is not an equal donor mean for {edge_id}/{condition}")
                edge_guides = guide_by_edge.loc[[str(edge_id)]]
                selected_guides = edge_guides.loc[
                    edge_guides["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
                    & edge_guides["donor"].isin(common)
                ]
                for slot in ("A", "B"):
                    condition_means: dict[str, float] = {}
                    for condition in CONDITIONS:
                        rows = selected_guides.loc[selected_guides["guide_slot"].eq(slot) & selected_guides["condition"].eq(condition)]
                        _expect(set(rows["donor"]) == set(common) and len(rows) == len(common), f"guide {slot} donor support differs from the exact common donors for {edge_id}/{condition}")
                        condition_means[condition] = float(rows["guide_total_effect"].mean())
                    expected_guide_delta = int(out["orientation"]) * (condition_means["Stim48hr"] - condition_means["Rest"])
                    _expect(np.isclose(float(out[f"oriented_delta__guide_{slot}_guide_total_effect"]), expected_guide_delta, atol=1e-12, rtol=1e-10), f"summary guide {slot} oriented delta is not reconstructed from long-guide equal-donor means for {edge_id}")


    def _validate_cal(summary: pd.DataFrame, envelopes: pd.DataFrame) -> None:
        _require_columns(envelopes, ENVELOPE_REQUIRED, "CAL envelopes")
        expected_grid = {(metric, "condition", condition) for metric in METRICS for condition in CONDITIONS} | {(metric, "cross_condition", "delta") for metric in METRICS}
        observed_grid = set(zip(envelopes["metric"], envelopes["scope"], envelopes["condition"]))
        _expect(len(envelopes) == len(expected_grid) and observed_grid == expected_grid and not envelopes.duplicated(["metric", "scope", "condition"]).any(), "CAL envelopes are not the exact metric x Rest/Stim48hr/delta grid")
        _expect(envelopes["quantile"].eq(0.95).all(), "CAL quantile must be 0.95")
        _expect(envelopes["quantile_definition"].eq("HYNDMAN_FAN_TYPE_7_NUMPY_LINEAR").all(), "CAL envelopes must use Hyndman-Fan type 7 / NumPy linear")
        _expect(envelopes["source_role"].eq("CAL_REVERSAL").all() and envelopes["edge_weights"].eq("EQUAL").all(), "CAL envelopes must use equally weighted reversal-parent CAL only")
        _all_bool(envelopes["stable_target_CAL_included"], "CAL stable-target flags")
        _expect((~envelopes["stable_target_CAL_included"].astype(bool)).all(), "stable-target CAL entered a primary envelope")
        _all_bool(envelopes["available"], "CAL available")
        cal = summary.loc[summary["analysis_role"].eq("CAL_REVERSAL") & summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        _expect(not cal["parent_reversal_edge_id"].duplicated().any(), "supported CAL reversal parents are not unique")
        indexed = envelopes.set_index(["metric", "scope", "condition"])
        for metric, scope, condition in expected_grid:
            column = f"{condition}__{metric}" if scope == "condition" else f"delta__{metric}"
            values = pd.to_numeric(cal[column], errors="coerce")
            finite = np.isfinite(values)
            n = int(cal.loc[finite, "parent_reversal_edge_id"].nunique())
            row = indexed.loc[(metric, scope, condition)]
            _expect(int(row["n_unique_supported_CAL"]) == n, f"CAL support count mismatch for {metric}/{condition}")
            available = n >= 1_000
            _expect(bool(row["available"]) == available, f"CAL availability mismatch for {metric}/{condition}")
            if available:
                q95 = float(np.quantile(np.abs(values.loc[finite].to_numpy(float)), 0.95, method="linear"))
                _expect(np.isclose(float(row["absolute_q95"]), q95, atol=1e-12, rtol=1e-10), f"CAL envelope is not HF7 absolute q95 for {metric}/{condition}")
            else:
                _expect(pd.isna(row["absolute_q95"]), f"unavailable CAL envelope must be missing for {metric}/{condition}")


    def _validate_raw_veto(summary: pd.DataFrame, envelopes: pd.DataFrame) -> None:
        envelope_index = envelopes.set_index(["metric", "scope", "condition"])
        numerical = summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        for condition in CONDITIONS:
            raw = pd.to_numeric(summary.loc[numerical, f"{condition}__raw_rate_log2_ratio"], errors="coerce")
            primary = pd.to_numeric(summary.loc[numerical, f"{condition}__total_effect"], errors="coerce")
            row = envelope_index.loc[("raw_rate_log2_ratio", "condition", condition)]
            available = bool(row["available"])
            threshold = float(row["absolute_q95"]) if available else np.nan
            expected = raw.abs().ge(0.15) & available & raw.abs().gt(threshold) & raw.ne(0) & primary.ne(0) & np.sign(raw).eq(-np.sign(primary))
            observed = summary.loc[numerical, f"{condition}__raw_practical_opposite"]
            _all_bool(observed, f"{condition} raw practical-opposite")
            _expect(np.array_equal(observed.astype(bool).to_numpy(), expected.to_numpy()), f"{condition} raw-rate practical-opposite rule mismatch")
            availability = summary.loc[numerical, f"{condition}__raw_CAL_q95_available"]
            _all_bool(availability, f"{condition} raw CAL availability")
            _expect(availability.astype(bool).eq(available).all(), f"{condition} raw CAL availability not propagated")
            reported_q = pd.to_numeric(summary.loc[numerical, f"{condition}__raw_CAL_q95"], errors="coerce")
            if available:
                _expect(np.allclose(reported_q, threshold, atol=1e-12, rtol=0), f"{condition} raw CAL q95 not propagated")
            else:
                _expect(reported_q.isna().all(), f"{condition} unavailable raw CAL q95 must be missing")
        veto = summary[[f"{condition}__raw_practical_opposite" for condition in CONDITIONS]].any(axis=1)
        _expect(summary.loc[numerical & veto, "analysis_status"].eq("MIXED_OR_UNRESOLVED").all(), "raw practical opposite did not force MIXED_OR_UNRESOLVED")
        _expect(summary.loc[numerical & ~veto, "analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS").all(), "raw-compatible continuous result has wrong status")


    def load_production_bundle(
        manifest_path: Path = DEFAULT_MANIFEST,
        config_path: Path = DEFAULT_CONFIG,
        amendment_path: Path = DEFAULT_AMENDMENT,
        metadata_report_path: Path = DEFAULT_METADATA_REPORT,
        cal_map_path: Path = DEFAULT_CAL_MAP,
        frozen_census_path: Path = DEFAULT_FROZEN_CENSUS,
        stable_matches_path: Path = DEFAULT_STABLE_MATCHES,
        manifest_checksum_path: Path | None = None,
    ) -> ProductionBundle:
        builder_path = Path(__file__)
        _expect(builder_path.is_file() and not builder_path.is_symlink(), "Figure 4 builder must be a regular non-symlink file")
        # The pre-outcome amendment is an internal pre-registration record and is not
        # redistributed, so it is verified only when present.
        immutable_paths = tuple(Path(path) for path in (
            builder_path, config_path, metadata_report_path,
            cal_map_path, frozen_census_path, stable_matches_path,
        ) if Path(path).is_file() or path is not amendment_path)
        _expect(all(path.is_file() and not path.is_symlink() for path in immutable_paths), "Figure 4 frozen input is absent or symlinked")
        immutable_before = {path: sha256_file(path) for path in immutable_paths}
        frozen = _validate_frozen_context(config_path, amendment_path, metadata_report_path, cal_map_path, frozen_census_path, stable_matches_path)
        manifest_path = Path(manifest_path)
        checksum_path = Path(manifest_checksum_path) if manifest_checksum_path else manifest_path.with_suffix(".sha256")
        manifest, paths, hashes = _validate_manifest(manifest_path, checksum_path)
        stage_before = {
            manifest_path: sha256_file(manifest_path),
            checksum_path: sha256_file(checksum_path),
            **{paths[key]: hashes[key] for key in paths},
        }
        report = _read_json(paths["gather_report"], "gather report")
        _validate_gather_report(report, manifest, paths, hashes)
        audit_exception = _read_json(paths[AUDIT_EXCEPTION_KEY], "transparent audit exception")
        _validate_audit_exception(audit_exception, manifest)
        final_closure = _read_json(paths[FINAL_CLOSURE_KEY], "proportionate final custody closure")
        _validate_final_custody_closure(final_closure, audit_exception, manifest, paths, hashes, report)
        units = _read_table(paths["unit_effects"], "unit effects")
        guides = _read_table(paths["guide_specific_effects"], "guide-specific effects")
        summary = _read_table(paths["edge_continuous_summary"], "edge continuous summary")
        envelopes = _read_table(paths["cal_envelopes"], "CAL envelopes")
        census = _read_table(paths["reversal_census"], "reversal census")
        waterfall = _read_table(paths["status_waterfall"], "status waterfall")
        _validate_census(census, waterfall, frozen, manifest)
        _validate_summary(summary, census, frozen)
        _validate_units(units, guides, summary, frozen)
        _validate_cal(summary, envelopes)
        _validate_raw_veto(summary, envelopes)
        _validate_locked_production_stage(manifest_path, checksum_path)
        _expect(all(path.is_file() and not path.is_symlink() and sha256_file(path) == digest for path, digest in immutable_before.items()), "Figure 4 builder or frozen input changed during validation")
        _expect(all(path.is_file() and not path.is_symlink() and sha256_file(path) == digest for path, digest in stage_before.items()), "production result stage changed during validation")
        return ProductionBundle(
            manifest, stage_before[manifest_path], hashes, frozen.config,
            units, guides, summary, envelopes, census, waterfall, report,
            audit_exception, final_closure,
        )


    def _symmetric_limit(*arrays: Iterable[float]) -> float:
        maximum = 0.0
        for array in arrays:
            values = np.asarray(list(array), dtype=float)
            values = np.abs(values[np.isfinite(values)])
            if values.size:
                maximum = max(maximum, float(values.max()))
        return max(1e-6, maximum * 1.06)


    def _envelope(bundle: ProductionBundle, metric: str, condition: str) -> float | None:
        scope = "cross_condition" if condition == "delta" else "condition"
        row = bundle.cal_envelopes.set_index(["metric", "scope", "condition"]).loc[(metric, scope, condition)]
        return float(row["absolute_q95"]) if bool(row["available"]) else None


    def _reversal_summary(bundle: ProductionBundle) -> pd.DataFrame:
        return bundle.edge_summary.loc[bundle.edge_summary["analysis_role"].eq("REVERSAL")].copy()


    def _opaque_legend(ax: plt.Axes, **kwargs: Any) -> Any:
        """Create a fully opaque legend so reference lines cannot strike through text."""
        legend = ax.legend(
            frameon=True,
            facecolor="white",
            edgecolor="none",
            framealpha=1.0,
            **kwargs,
        )
        frame = legend.get_frame()
        _expect(frame.get_visible() and frame.get_alpha() == 1.0, "Figure 4 legend background must be fully opaque")
        _expect(np.isclose(frame.get_facecolor()[-1], 1.0), "Figure 4 legend face colour must be opaque")
        return legend


    def _panel_waterfall(ax: plt.Axes, bundle: ProductionBundle) -> None:
        counts = {status: int(bundle.census["analysis_status"].eq(status).sum()) for status in FINAL_STATUSES}
        values = np.asarray(
            (
                DENOMINATOR,
                METADATA_ELIGIBLE,
                METADATA_ELIGIBLE - counts["UNANALYSABLE_POSITIVE_INTENSITY"],
                counts["ESTIMATED_CONTINUOUS_COMPONENTS"] + counts["MIXED_OR_UNRESOLVED"],
                counts["ESTIMATED_CONTINUOUS_COMPONENTS"],
            )
        )
        labels = ("census", "metadata", "donor\nsupport", "effect", "raw\ncompatible")
        x = np.arange(len(values))
        ax.bar(x, values, color=("#46515C", "#5A7184", "#4E8D8A", "#2A9D8F", "#216E5C"), width=0.68)
        ax.set_xticks(x, labels, fontsize=MINIMUM_FINAL_TEXT_PT)
        ax.set_ylabel("reversal edges")
        ax.set_ylim(0, DENOMINATOR * 1.16)
        for xx, value in zip(x, values):
            ax.text(xx, value + DENOMINATOR * 0.025, f"{value:,}", ha="center", va="bottom", fontsize=7)
        ax.text(0.98, 0.96, f"mixed/unresolved: {counts['MIXED_OR_UNRESOLVED']:,}", transform=ax.transAxes, ha="right", va="top", fontsize=MINIMUM_FINAL_TEXT_PT, color="#666B70")
        ax.set_title("Full-denominator eligibility", loc="left", pad=7)


    def _panel_rest_stim(ax: plt.Axes, bundle: ProductionBundle) -> None:
        reversal = _reversal_summary(bundle)
        supported = reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        table = reversal.loc[supported]
        census_direction = bundle.census.set_index("reversal_edge_id")["direction"]
        gained = table["analysis_edge_id"].map(census_direction).eq("gained_on_activation")
        for mask, colour, marker, label in ((gained, "#2A9D8F", "o", "gained"), (~gained, "#C1272D", "^", "lost")):
            block = table.loc[mask]
            ax.scatter(block["Rest__total_effect"], block["Stim48hr__total_effect"], s=7, alpha=0.30, color=colour, marker=marker, edgecolors="none", rasterized=True, label=f"{label} (n={len(block):,})")
        limit = _symmetric_limit(table["Rest__total_effect"], table["Stim48hr__total_effect"])
        ax.plot([-limit, limit], [-limit, limit], color="#C5C9CD", lw=0.8)
        ax.axhline(0, color="#9AA0A6", lw=0.7)
        ax.axvline(0, color="#9AA0A6", lw=0.7)
        for condition, orientation in (("Rest", "vertical"), ("Stim48hr", "horizontal")):
            threshold = _envelope(bundle, "total_effect", condition)
            if threshold is not None:
                for sign in (-1, 1):
                    (ax.axvline if orientation == "vertical" else ax.axhline)(sign * threshold, color="#767D84", lw=0.7, ls="--")
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Rest donor-mean total effect")
        ax.set_ylabel("48 h donor-mean total effect")
        ax.set_title("Rest vs 48 h total effect", loc="left", pad=7)
        _opaque_legend(ax, loc="upper left", fontsize=MINIMUM_FINAL_TEXT_PT)


    def _panel_components(ax: plt.Axes, bundle: ProductionBundle) -> None:
        reversal = _reversal_summary(bundle)
        reversal = reversal.loc[reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        stable = bundle.edge_summary.loc[bundle.edge_summary["analysis_role"].eq("REGULATOR_MATCHED_STABLE") & bundle.edge_summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        if len(stable):
            ax.scatter(stable["oriented_delta__detection_component"], stable["oriented_delta__positive_intensity_component"], s=3, alpha=0.07, color="#9AA0A6", edgecolors="none", rasterized=True, label="regulator-matched stable")
        mixed = reversal["analysis_status"].eq("MIXED_OR_UNRESOLVED")
        ax.scatter(reversal.loc[~mixed, "oriented_delta__detection_component"], reversal.loc[~mixed, "oriented_delta__positive_intensity_component"], s=8, alpha=0.34, color="#2A9D8F", edgecolors="none", rasterized=True, label="estimated")
        ax.scatter(reversal.loc[mixed, "oriented_delta__detection_component"], reversal.loc[mixed, "oriented_delta__positive_intensity_component"], s=9, alpha=0.40, color="#7B6FA6", marker="x", linewidths=0.45, rasterized=True, label="mixed/unresolved")
        for metric, orientation in (("detection_component", "vertical"), ("positive_intensity_component", "horizontal")):
            threshold = _envelope(bundle, metric, "delta")
            if threshold is not None:
                for sign in (-1, 1):
                    (ax.axvline if orientation == "vertical" else ax.axhline)(sign * threshold, color="#767D84", lw=0.7, ls="--")
        limit = _symmetric_limit(reversal["oriented_delta__detection_component"], reversal["oriented_delta__positive_intensity_component"], stable.get("oriented_delta__detection_component", []), stable.get("oriented_delta__positive_intensity_component", []))
        ax.axhline(0, color="#9AA0A6", lw=0.7)
        ax.axvline(0, color="#9AA0A6", lw=0.7)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("oriented change, detection component")
        ax.set_ylabel("oriented change, positive-cell intensity")
        ax.set_title("Continuous Shapley components", loc="left", pad=7)
        _opaque_legend(ax, loc="upper left", fontsize=MINIMUM_FINAL_TEXT_PT)


    def _cal_annotation_table(envelopes: pd.DataFrame) -> list[list[str]]:
        index = envelopes.set_index(["metric", "scope", "condition"])
        labels: list[list[str]] = []
        for metric in METRICS:
            row_labels: list[str] = []
            for condition in ("Rest", "Stim48hr", "delta"):
                scope = "cross_condition" if condition == "delta" else "condition"
                row = index.loc[(metric, scope, condition)]
                if bool(row["available"]):
                    row_labels.append(f"{float(row['absolute_q95']):.3g}\nn={int(row['n_unique_supported_CAL']):,}")
                else:
                    row_labels.append(f"unavailable\nn={int(row['n_unique_supported_CAL']):,}")
            labels.append(row_labels)
        return labels


    def _panel_cal(ax: plt.Axes, envelopes: pd.DataFrame) -> None:
        ax.axis("off")
        table = ax.table(
            cellText=_cal_annotation_table(envelopes),
            rowLabels=("total", "detection", "positive intensity", "raw UMI rate"),
            colLabels=("Rest", "48 h", "48 h - Rest"),
            cellLoc="center",
            rowLoc="right",
            colLoc="center",
            bbox=(0.01, 0.10, 0.98, 0.80),
        )
        table.auto_set_font_size(False)
        table.set_fontsize(MINIMUM_FINAL_TEXT_PT)
        for (row, _column), cell in table.get_celld().items():
            cell.set_edgecolor("#D9DDE1")
            cell.set_linewidth(0.55)
            cell.set_facecolor("#EAF2F7" if row > 0 and row % 2 else "#F7F9FA")
            if row == 0:
                cell.set_facecolor("#E1E7EB")
                cell.get_text().set_fontweight("bold")
        ax.set_title("CAL 95th-percentile envelopes", loc="left", pad=7)
        ax.text(0.01, 0.01, "Metric-specific units; compare conditions within rows.", transform=ax.transAxes, ha="left", va="bottom", fontsize=MINIMUM_FINAL_TEXT_PT, color="#565B61")


    def _stability(bundle: ProductionBundle) -> tuple[pd.DataFrame, pd.DataFrame]:
        reversal = _reversal_summary(bundle)
        supported = reversal.loc[reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")].set_index("analysis_edge_id")
        guides = supported[["oriented_delta__guide_A_guide_total_effect", "oriented_delta__guide_B_guide_total_effect"]].rename(columns={"oriented_delta__guide_A_guide_total_effect": "A", "oriented_delta__guide_B_guide_total_effect": "B"}).dropna()
        units = bundle.unit_effects.loc[bundle.unit_effects["analysis_role"].eq("REVERSAL") & bundle.unit_effects["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        rows = []
        for donor in DONORS:
            estimates: list[tuple[str, float]] = []
            for edge_id, block in units.groupby("analysis_edge_id"):
                if edge_id not in supported.index:
                    continue
                frozen_common = _canonical_donors(supported.loc[edge_id, "common_donor_ids"])
                if donor not in frozen_common or len(frozen_common) < 4:
                    continue
                retained_common = tuple(item for item in frozen_common if item != donor)
                retained = block.loc[block["donor"].isin(retained_common)]
                if any(set(retained.loc[retained["condition"].eq(condition), "donor"]) != set(retained_common) for condition in CONDITIONS):
                    continue
                means = retained.groupby("condition")["total_effect"].mean()
                estimates.append((edge_id, (means["Stim48hr"] - means["Rest"]) * int(supported.loc[edge_id, "orientation"])))
            if estimates:
                ids = [edge for edge, _ in estimates]
                values = np.asarray([value for _, value in estimates])
                reference = supported.loc[ids, "oriented_delta__total_effect"].to_numpy(float)
                agreement = float(np.mean(np.sign(values) == np.sign(reference)))
                median_abs = float(np.median(np.abs(values - reference)))
            else:
                agreement = np.nan
                median_abs = np.nan
            rows.append({"donor": donor, "n_evaluable": len(estimates), "sign_agreement": agreement, "median_absolute_difference": median_abs})
        return guides, pd.DataFrame(rows)


    def _panel_stability(fig: plt.Figure, spec: Any, bundle: ProductionBundle) -> list[plt.Axes]:
        grid = spec.subgridspec(1, 2, wspace=0.42)
        ax_guide = fig.add_subplot(grid[0, 0])
        ax_lodo = fig.add_subplot(grid[0, 1])
        guides, lodo = _stability(bundle)
        ax_guide.scatter(guides["A"], guides["B"], s=7, alpha=0.28, color="#457B9D", edgecolors="none", rasterized=True)
        limit = _symmetric_limit(guides["A"], guides["B"])
        ax_guide.plot([-limit, limit], [-limit, limit], color="#B8BDC2", lw=0.8)
        ax_guide.axhline(0, color="#9AA0A6", lw=0.6)
        ax_guide.axvline(0, color="#9AA0A6", lw=0.6)
        ax_guide.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="guide A oriented total change", ylabel="guide B oriented total change")
        ax_guide.set_aspect("equal", adjustable="box")
        agreement = float(np.mean(np.sign(guides["A"]) == np.sign(guides["B"]))) if len(guides) else np.nan
        ax_guide.set_title("Physical-guide agreement", loc="left", pad=7)
        ax_guide.text(
            0.03,
            0.96,
            f"n={len(guides):,}; sign={agreement:.1%}",
            transform=ax_guide.transAxes,
            ha="left",
            va="top",
            fontsize=MINIMUM_FINAL_TEXT_PT,
            color="#565B61",
        )
        x = np.arange(len(lodo))
        ax_lodo.bar(x, lodo["sign_agreement"], color="#6A8EAE", width=0.62)
        ax_lodo.set_xticks(x, [f"leave out\n{donor}" for donor in lodo["donor"]])
        ax_lodo.set_ylim(0, 1.16)
        ax_lodo.set_ylabel("direction agreement")
        for xx, row in enumerate(lodo.itertuples(index=False)):
            label = "n=0" if row.n_evaluable == 0 else f"{row.sign_agreement:.1%}\nn={row.n_evaluable:,}"
            ax_lodo.text(xx, 0.02 if row.n_evaluable == 0 else row.sign_agreement + 0.025, label, ha="center", va="bottom", fontsize=MINIMUM_FINAL_TEXT_PT)
        ax_lodo.set_title("Leave-one-donor agreement", loc="left", pad=7)
        return [ax_guide, ax_lodo]


    def _paired_example_donors(block: pd.DataFrame) -> tuple[str, ...]:
        required = {(condition, slot) for condition in CONDITIONS for slot in ("A", "B")}
        return tuple(
            donor
            for donor in DONORS
            if set(zip(block.loc[block["donor"].eq(donor), "condition"], block.loc[block["donor"].eq(donor), "guide_slot"])) == required
        )


    def _supported_predeclared_example_ids(bundle: ProductionBundle) -> tuple[str, ...]:
        """Return only predeclared edges with estimable, plottable donor-guide data."""
        census = bundle.census.set_index("reversal_edge_id")
        summaries = bundle.edge_summary.loc[
            bundle.edge_summary["analysis_role"].eq("REVERSAL")
        ].set_index("analysis_edge_id")
        guides = bundle.guide_specific_effects
        supported: list[str] = []
        for edge_id, _regulator, _target, _direction, _sign in PREDECLARED_EDGES:
            if census.loc[edge_id, "analysis_status"] != "ESTIMATED_CONTINUOUS_COMPONENTS":
                continue
            if edge_id not in summaries.index:
                continue
            summary = summaries.loc[edge_id]
            if isinstance(summary, pd.DataFrame) or summary["summary_status"] != "ESTIMATED_CONTINUOUS_COMPONENTS":
                continue
            block = guides.loc[
                guides["analysis_role"].eq("REVERSAL")
                & guides["analysis_edge_id"].eq(edge_id)
                & guides["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
            ]
            if _paired_example_donors(block):
                supported.append(edge_id)
        return tuple(supported)


    def _predeclared_example_support_label(supported_ids: Sequence[str]) -> str:
        return f"{len(supported_ids)}/{len(PREDECLARED_EDGES)}"


    def _panel_examples(
        fig: plt.Figure,
        spec: Any,
        bundle: ProductionBundle,
        supported_edge_ids: Sequence[str],
    ) -> list[plt.Axes]:
        _expect(bool(supported_edge_ids), "examples panel requires at least one supported predeclared edge")
        supported_set = set(supported_edge_ids)
        selected = [row for row in PREDECLARED_EDGES if row[0] in supported_set]
        n_rows = (len(selected) + 1) // 2
        grid = spec.subgridspec(n_rows, 2, wspace=0.30, hspace=0.52)
        axes = []
        guides = bundle.guide_specific_effects
        donor_colours = {"D1": "#4C78A8", "D2": "#72B7B2", "D3": "#F58518", "D4": "#B279A2"}
        for index, (edge_id, regulator, target, _direction, _sign) in enumerate(selected):
            ax = fig.add_subplot(grid[index // 2, index % 2])
            axes.append(ax)
            block = guides.loc[guides["analysis_role"].eq("REVERSAL") & guides["analysis_edge_id"].eq(edge_id) & guides["unit_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
            paired_donors = _paired_example_donors(block)
            _expect(bool(paired_donors), f"supported example lost paired donor-guide data: {edge_id}")
            for donor in paired_donors:
                donor_block = block.loc[block["donor"].eq(donor)]
                for slot, marker, offset in (("A", "o", -0.07), ("B", "^", 0.07)):
                    sub = donor_block.loc[donor_block["guide_slot"].eq(slot)].set_index("condition")
                    x = np.asarray((offset, 1 + offset))
                    y = np.asarray((sub.loc["Rest", "guide_total_effect"], sub.loc["Stim48hr", "guide_total_effect"]), dtype=float)
                    ax.plot(x, y, color=donor_colours[donor], alpha=0.45, lw=0.7)
                    ax.scatter(x, y, color=donor_colours[donor], marker=marker, s=15, edgecolor="white", linewidth=0.25, zorder=3)
            ax.axhline(0, color="#A8ADB2", lw=0.7)
            ax.set_xticks((0, 1), ("Rest", "48 h"))
            if index % 2 == 0:
                ax.set_ylabel("guide-specific total effect")
            ax.text(0.98, 0.02, f"n={len(paired_donors)} paired donors", transform=ax.transAxes, ha="right", va="bottom", fontsize=MINIMUM_FINAL_TEXT_PT, color="#666B70")
            ax.set_title(f"{regulator} → {target}", fontsize=7.3, pad=5)
        guide_handles = [
            Line2D([], [], color="#4A4F54", marker="o", linestyle="none", label="guide A"),
            Line2D([], [], color="#4A4F54", marker="^", linestyle="none", label="guide B"),
        ]
        donor_handles = [
            Line2D([], [], color=colour, marker="o", linestyle="none", label=donor)
            for donor, colour in donor_colours.items()
        ]
        fig.legend(handles=guide_handles + donor_handles, loc="lower center", bbox_to_anchor=(0.55, 0.035), ncol=6, fontsize=MINIMUM_FINAL_TEXT_PT, frameon=False, handletextpad=0.3, columnspacing=0.9)
        return axes


    def _assert_native_geometry_and_text(fig: plt.Figure) -> float:
        width = float(fig.get_size_inches()[0])
        _expect(np.isclose(width, FIGURE_WIDTH_IN, atol=1e-9, rtol=0), f"Figure 4 width must be exactly {FIGURE_WIDTH_IN:.2f} in")
        visible_text = [item for item in fig.findobj(match=Text) if item.get_visible() and str(item.get_text()).strip()]
        _expect(bool(visible_text), "Figure 4 contains no visible text")
        minimum = min(float(item.get_fontsize()) for item in visible_text)
        _expect(minimum + 1e-9 >= MINIMUM_FINAL_TEXT_PT, f"Figure 4 contains final-size text below {MINIMUM_FINAL_TEXT_PT:g} pt")
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        figure_bounds = fig.bbox
        tolerance = 1.0
        for item in visible_text:
            bounds = item.get_window_extent(renderer=renderer)
            _expect(
                bounds.x0 >= figure_bounds.x0 - tolerance
                and bounds.y0 >= figure_bounds.y0 - tolerance
                and bounds.x1 <= figure_bounds.x1 + tolerance
                and bounds.y1 <= figure_bounds.y1 + tolerance,
                f"Figure 4 text is clipped outside the page: {item.get_text()!r}",
            )

        def overlaps(first: Text, second: Text) -> bool:
            a = first.get_window_extent(renderer=renderer)
            b = second.get_window_extent(renderer=renderer)
            return a.x0 < b.x1 and a.x1 > b.x0 and a.y0 < b.y1 and a.y1 > b.y0

        titles = [ax.title for ax in fig.axes if ax.title.get_visible() and ax.title.get_text().strip()]
        _expect(
            not any(overlaps(first, second) for index, first in enumerate(titles) for second in titles[index + 1 :]),
            "Figure 4 panel titles overlap",
        )
        for ax in fig.axes:
            tick_labels = [item for item in ax.get_xticklabels() if item.get_visible() and item.get_text().strip()]
            _expect(
                not any(overlaps(first, second) for index, first in enumerate(tick_labels) for second in tick_labels[index + 1 :]),
                "Figure 4 x tick labels overlap",
            )
            title = ax.title
            annotations = [item for item in ax.texts if item.get_visible() and item.get_text().strip()]
            _expect(
                not title.get_text().strip() or not any(overlaps(title, item) for item in annotations),
                f"Figure 4 title overlaps an annotation: {title.get_text()!r}",
            )
        return minimum


    def render_figure(bundle: ProductionBundle, output_pdf: Path) -> tuple[Path, Path, Path]:
        """Render real production PDF/PNG and an input/output hash sidecar."""
        output_pdf = Path(output_pdf)
        _expect(output_pdf.suffix.lower() == ".pdf", "Figure 4 output must be PDF")
        output_png = output_pdf.with_suffix(".png")
        output_record = output_pdf.with_suffix(".production.json")
        for path in (output_pdf, output_png, output_record):
            _expect(not path.exists(), f"refusing to overwrite {path}")
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        supported_examples = _supported_predeclared_example_ids(bundle)
        support_label = _predeclared_example_support_label(supported_examples)
        has_examples = bool(supported_examples)
        S.setup(1.0)
        height = FIGURE_HEIGHT_WITH_EXAMPLES_IN if has_examples else FIGURE_HEIGHT_WITHOUT_EXAMPLES_IN
        fig = plt.figure(figsize=(FIGURE_WIDTH_IN, height), facecolor="white")
        rows = 4 if has_examples else 3
        height_ratios = (1.0, 1.0, 1.0, 1.65) if has_examples else None
        grid = fig.add_gridspec(rows, 2, left=0.105, right=0.975, top=0.94 if has_examples else 0.915, bottom=0.10 if has_examples else 0.095, wspace=0.42, hspace=0.67 if has_examples else 0.64, height_ratios=height_ratios)
        axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 0])]
        ax_d = fig.add_subplot(grid[1, 1])
        _panel_waterfall(axes[0], bundle)
        _panel_rest_stim(axes[1], bundle)
        _panel_components(axes[2], bundle)
        _panel_cal(ax_d, bundle.cal_envelopes)
        stability_axes = _panel_stability(fig, grid[2, :], bundle)
        for label, ax in zip("abcd", (*axes, ax_d)):
            S.panel_tag(ax, label, dx=-0.055, dy=1.09)
        S.panel_tag(stability_axes[0], "e", dx=-0.12, dy=1.09)
        if has_examples:
            example_axes = _panel_examples(fig, grid[3, :], bundle, supported_examples)
            S.panel_tag(example_axes[0], "f", dx=-0.17, dy=1.10)
        fig.suptitle("Donor-balanced target-response shapes across activation", x=0.07, y=0.982 if has_examples else 0.972, ha="left", fontsize=13.0, fontweight="bold", color="#22262B")
        fig.text(
            0.07,
            0.014,
            "Continuous descriptive estimates; dashed lines mark global absolute CAL 95th percentiles.\n"
            "n denotes reversal edges or donor-resolved summaries.",
            ha="left",
            va="bottom",
            fontsize=MINIMUM_FINAL_TEXT_PT,
            color="#565B61",
        )
        minimum_text_pt = _assert_native_geometry_and_text(fig)
        with tempfile.TemporaryDirectory(prefix="fig4_v5_1_", dir=output_pdf.parent) as temporary:
            temporary = Path(temporary)
            temp_pdf, temp_png = temporary / output_pdf.name, temporary / output_png.name
            fig.savefig(temp_pdf, facecolor="white")
            fig.savefig(temp_png, dpi=300, facecolor="white")
            plt.close(fig)
            record = {
                "schema": "response-shape-figure4-render-v5.1",
                "status": "COMPLETE",
                "production": True,
                "outcomes_read": True,
                "run_id": bundle.manifest["run_id"],
                "source_manifest_sha256": bundle.manifest_sha256,
                "methods_config_sha256": FROZEN_CONFIG_SHA256,
                "amendment_sha256": FROZEN_AMENDMENT_SHA256,
                "claim_boundary": CLAIM_BOUNDARY,
                "audit_exception_sha256": bundle.artifact_sha256[AUDIT_EXCEPTION_KEY],
                "proportionate_final_custody_closure_sha256": bundle.artifact_sha256[FINAL_CLOSURE_KEY],
                "input_artifact_sha256": dict(bundle.artifact_sha256),
                "predeclared_examples_supported": support_label,
                "predeclared_examples_panel_rendered": bool(has_examples),
                "figure_width_in": FIGURE_WIDTH_IN,
                "minimum_nonblank_text_pt": minimum_text_pt,
                "outputs": {
                    "pdf": {"filename": output_pdf.name, "sha256": sha256_file(temp_pdf)},
                    "png": {"filename": output_png.name, "sha256": sha256_file(temp_png)},
                },
            }
            temp_record = temporary / output_record.name
            temp_record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            os.replace(temp_pdf, output_pdf)
            os.replace(temp_png, output_png)
            os.replace(temp_record, output_record)
        return output_pdf, output_png, output_record


    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
        parser.add_argument("--manifest-checksum", type=Path)
        parser.add_argument("--methods-config", type=Path, default=DEFAULT_CONFIG)
        parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
        parser.add_argument("--metadata-report", type=Path, default=DEFAULT_METADATA_REPORT)
        parser.add_argument("--cal-map", type=Path, default=DEFAULT_CAL_MAP)
        parser.add_argument("--frozen-census", type=Path, default=DEFAULT_FROZEN_CENSUS)
        parser.add_argument("--stable-matches", type=Path, default=DEFAULT_STABLE_MATCHES)
        parser.add_argument("--output-pdf", type=Path, default=DEFAULT_PDF)
        parser.add_argument("--render", action="store_true", help="write PDF/PNG after validation; default is validation only")
        args = parser.parse_args()
        bundle = load_production_bundle(
            manifest_path=args.manifest,
            config_path=args.methods_config,
            amendment_path=args.amendment,
            metadata_report_path=args.metadata_report,
            cal_map_path=args.cal_map,
            frozen_census_path=args.frozen_census,
            stable_matches_path=args.stable_matches,
            manifest_checksum_path=args.manifest_checksum,
        )
        supported_examples = _supported_predeclared_example_ids(bundle)
        result: dict[str, Any] = {
            "status": "PRODUCTION_FIGURE4_INPUTS_VALIDATED",
            "run_id": bundle.manifest["run_id"],
            "denominator_n": len(bundle.census),
            "methods_config_sha256": FROZEN_CONFIG_SHA256,
            "amendment_sha256": FROZEN_AMENDMENT_SHA256,
            "predeclared_examples_supported": _predeclared_example_support_label(supported_examples),
            "figure_written": False,
        }
        if args.render:
            pdf, png, record = render_figure(bundle, args.output_pdf)
            result.update({"figure_written": True, "pdf": str(pdf), "png": str(png), "render_manifest": str(record)})
        print(json.dumps(result, indent=2, sort_keys=True))


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_production_v5_1()


# ==========================================================================
# fig4_effector_landmark_panel
# ==========================================================================
def _load_fig4_effector_landmark_panel():
    _self = types.ModuleType('fig4_effector_landmark_panel')
    __name__ = 'fig4_effector_landmark_panel'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_effector_landmark_panel'] = _self
    #!/usr/bin/env python3
    """Figure 4e candidate: reversal links from the Figure 3 modules to the fixed T-cell landmarks.

    Transposed against the earlier diagnostic: 25 landmark rows by 7 regulator columns, which is
    the shape that fits panel e's tall box (about 2.49 x 5.32 in) at the 7 pt floor. All 25
    landmarks are rows, including the nine that carry no census edge from these regulators, so the
    display cannot be read as a selected list.

    Encoding, with orientation carried twice so it survives greyscale and colour-blind viewing:
      filled tile        a reversal; red plus solid = lost on stimulation, teal plus hatch = gained
      outlined tile      comparable edge, significant in both states, did not reverse
      inner circle       reaches the response-shape layer, intensity-dominant
      inner triangle     reaches the response-shape layer, detection-dominant
      star               external assay anchor, CBFB->IL2RA, metadata-ineligible for this layer

    No text inside the panel beyond the axis labels: the wording lives in the caption, and the
    localisation statistics live in data/census/effector_landmark_localisation.json.
    The earlier draft printed the convergence-test p = 0.46 here, which tests GO annotation among
    globally convergent targets and says nothing about where these regulators land; it is removed.

    Standalone render writes data/figure_provenance/Fig4e_effector_landmarks.{pdf,png} plus a
    .provenance.json recording the input hashes.
    """


    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42     # embed TrueType, never Type 3
    matplotlib.rcParams["ps.fonttype"] = 42

    _HERE = Path(__file__).resolve().parent
    for _p in (str(_HERE), str(_HERE.parent)):
        if _p not in sys.path:
            sys.path.insert(0, _p)
    P = _MODULES['fig4_production_v5_1']

    OUT = STAGING_DIR  # standalone diagnostic render; not a repository artefact
    MASTER = P.PAPER / "data/census/census_master_edges.csv.gz"
    PANEL_W_IN, PANEL_H_IN = 2.494, 5.321      # panel e's box in the frozen 7.10 x 10.04 layout

    LOST, GAINED = "#c1272d", "#2a9d8f"
    STABLE_EDGE, INK, MUTED, GOLDC = "#7f8a94", "#22262B", "#666B70", "#e0a91b"

    LANDMARKS = ["GZMA", "GZMB", "GNLY", "PRF1", "IFNG", "TNF", "NKG7", "CCL5", "CCL4", "CCL3",
                 "CCR5", "CCR1", "CCR7", "CCR2", "IL2", "IL2RA", "IL2RB", "CD69", "CTLA4",
                 "IL2RG", "SELL", "TBX21", "KLRB1", "XCL1", "XCL2"]
    BLOCKS = [("3d", ["CCNC", "TADA2B", "SUPT7L"]),
              ("3e", ["AHR", "ARNT"]),
              ("ext", ["TAF6L"]),
              ("anchor", ["CBFB"])]
    REGULATORS = [r for _, block in BLOCKS for r in block]
    ANCHOR = ("CBFB", "IL2RA")


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def load_cells() -> dict[tuple[str, str], dict]:
        """Every regulator/landmark pair present in the both-state-significant trans census."""
        master = pd.read_csv(MASTER)
        pairs = master.loc[~master["is_self"].astype(bool)
                           & master["regulator"].isin(REGULATORS)
                           & master["target"].isin(LANDMARKS),
                           ["regulator", "target", "kind", "direction"]]

        bundle = P.load_production_bundle()
        rev = bundle.edge_summary.loc[bundle.edge_summary["analysis_role"].eq("REVERSAL")].copy()
        det = pd.to_numeric(rev["oriented_delta__detection_component"], errors="coerce").abs()
        inten = pd.to_numeric(rev["oriented_delta__positive_intensity_component"], errors="coerce").abs()
        dom = pd.Series(np.where(det > inten, "detection", "intensity"),
                        index=rev["analysis_edge_id"].astype(str))
        kept = bundle.census.loc[bundle.census["analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        shape = {(r, t): dom.get(str(e)) for r, t, e in
                 zip(kept["regulator"], kept["target"], kept["reversal_edge_id"])}

        cells = {}
        for _, row in pairs.iterrows():
            k = (row["regulator"], row["target"])
            cells[k] = {"kind": row["kind"], "direction": row["direction"], "shape": shape.get(k)}
        return cells


    def panel_effector_landmarks(ax, cells: dict, fs: float = 7.0, tag_blocks: bool = True,
                                 equal_aspect: bool = True) -> None:
        """Draw the landmark-by-regulator matrix into ax. Rows are ordered by census burden."""
        burden = {t: 0 for t in LANDMARKS}
        for (_, t), rec in cells.items():
            burden[t] += int(rec["kind"] == "rev")
        rows = sorted(LANDMARKS, key=lambda t: (-burden[t], t))

        x_of, x, seps = {}, 0.0, []
        for i, (_, block) in enumerate(BLOCKS):
            if i:
                seps.append(x)
                x += 0.34
            for reg in block:
                x_of[reg] = x
                x += 1.0
        width = x

        for (reg, gene), rec in cells.items():
            xx, yy = x_of[reg], float(len(rows) - 1 - rows.index(gene))
            if rec["kind"] == "stable":
                ax.add_patch(Rectangle((xx + 0.14, yy + 0.14), 0.72, 0.72, facecolor="none",
                                       edgecolor=STABLE_EDGE, lw=0.8))
                continue
            gained = rec["direction"] == "gained_on_activation"
            ax.add_patch(Rectangle((xx + 0.06, yy + 0.06), 0.88, 0.88,
                                   facecolor=GAINED if gained else LOST,
                                   edgecolor="white", lw=0.5,
                                   hatch="////" if gained else None))
            if rec["shape"] == "intensity":
                ax.plot([xx + 0.5], [yy + 0.5], marker="o", ms=2.9, mfc="white", mec="white", ls="")
            elif rec["shape"] == "detection":
                ax.plot([xx + 0.5], [yy + 0.5], marker="^", ms=3.4, mfc="white", mec="white", ls="")
            if (reg, gene) == ANCHOR:
                ax.plot([xx + 0.5], [yy + 0.5], marker="*", ms=6.0, mfc=GOLDC, mec="white",
                        mew=0.4, ls="", zorder=5)

        for s in seps:
            ax.axvline(s + 0.17, color="#c7ccd1", lw=0.7)
        # Square tiles are right for the standalone, where the panel sets its own size. Inside
        # Figure 4 they shrink the axes to 1.75 in in a 2.59 in slot, wasting the difference and
        # pushing a left-anchored heading off the page. Rectangular tiles are normal for a matrix.
        if equal_aspect:
            ax.set_aspect("equal")
        ax.set_xlim(-0.04, width - 0.04)
        ax.set_ylim(-0.06, len(rows) + 0.06)
        ax.set_xticks([x_of[r] + 0.5 for r in REGULATORS])
        ax.set_xticklabels(REGULATORS, fontsize=fs, rotation=90)
        ax.set_yticks([len(rows) - 1 - rows.index(t) + 0.5 for t in rows])
        ax.set_yticklabels(rows, fontsize=fs)
        for label, gene in zip(ax.get_yticklabels(), rows):
            if not any((reg, gene) in cells for reg in REGULATORS):
                label.set_color(MUTED)           # landmark with no comparable edge from these rows
        ax.tick_params(length=0)
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_tick_params(pad=1.5)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if tag_blocks:
            for (name, block), _ in zip(BLOCKS, BLOCKS):
                mid = (x_of[block[0]] + x_of[block[-1]] + 1.0) / 2
                # fs, not fs - 1.2: Figure 4's geometry gate rejects any final-size text below
                # 7 pt, and this panel is called with fs already at that floor.
                ax.text(mid, -0.55, name, fontsize=fs, color=MUTED, ha="center", va="top")


    def legend_handles(fs: float = 7.0) -> list:
        return [Rectangle((0, 0), 1, 1, facecolor=LOST, edgecolor="white", label="reversal, lost"),
                Rectangle((0, 0), 1, 1, facecolor=GAINED, edgecolor="white", hatch="////",
                          label="reversal, gained"),
                Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=STABLE_EDGE,
                          label="comparable, not reversed"),
                Line2D([0], [0], marker="o", ls="", mfc=INK, mec=INK, ms=3.4, label="intensity-dominant"),
                Line2D([0], [0], marker="^", ls="", mfc=INK, mec=INK, ms=3.8, label="detection-dominant"),
                Line2D([0], [0], marker="*", ls="", mfc=GOLDC, mec=GOLDC, ms=6.5, label="assay anchor")]


    CAPTION = (
        "Reversal links from the Figure 3 modules to all 25 fixed T-cell effector landmarks shown "
        "in Figure 1a. Filled tiles denote reversals lost (red, solid) or gained (teal, hatched) "
        "after stimulation; outlines denote comparable edges, significant in both states, that did "
        "not reverse. A blank cell is a pair that was never significant in both states and was "
        "therefore never eligible to reverse. Circles and triangles identify intensity- and "
        "detection-dominant response shapes; the star marks the external CBFB-IL2RA assay anchor, "
        "which is metadata-ineligible for the response-shape layer. Landmark names in grey carry no "
        "edge of either kind from these regulators. The panel is descriptive and does not test "
        "individual edges, enrichment, direct regulation or mediation."
    )


    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        cells = load_cells()
        # canvas = the panel box plus the margins its tick labels, block tags and title occupy
        # inside Figure 4's gridspec padding, so the preview shows the real final-size type
        fw, fh = PANEL_W_IN + 2.0, PANEL_H_IN + 2.6
        fig = plt.figure(figsize=(fw, fh), facecolor="white")
        ax = fig.add_axes([1.15 / fw, 1.55 / fh, PANEL_W_IN / fw, PANEL_H_IN / fh])
        panel_effector_landmarks(ax, cells)
        ax.set_title("Reversal links from Figure 3 modules\nto fixed T-cell landmarks",
                     loc="left", pad=34, fontsize=8.0, fontweight="bold", color=INK)
        ax.legend(handles=legend_handles(), loc="upper left", bbox_to_anchor=(-0.44, -0.055),
                  ncol=2, fontsize=6.6, frameon=False, handlelength=1.3, columnspacing=1.0,
                  labelspacing=0.35)
        for ext in ("pdf", "png"):
            fig.savefig(OUT / f"Fig4e_effector_landmarks.{ext}", dpi=300, facecolor="white")
        plt.close(fig)
        (OUT / "Fig4e_effector_landmarks.provenance.json").write_text(json.dumps({
            "_what": "Figure 4e candidate: Figure 3 modules against the fixed Figure 1a landmarks",
            "_script": "scripts/response_shape/fig4_effector_landmark_panel.py",
            "census_master_sha256": sha256(MASTER),
            "panel_box_in": [PANEL_W_IN, PANEL_H_IN],
            "regulators": REGULATORS, "landmarks": LANDMARKS,
            "n_cells": len(cells),
            "n_reversal_cells": sum(1 for c in cells.values() if c["kind"] == "rev"),
            "n_stable_cells": sum(1 for c in cells.values() if c["kind"] == "stable"),
            "caption": CAPTION,
        }, indent=2))
        print(f"wrote {OUT}/Fig4e_effector_landmarks.pdf  ({len(cells)} cells, "
              f"{sum(1 for c in cells.values() if c['kind'] == 'rev')} reversals)")


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_effector_landmark_panel()


# ==========================================================================
# fig4_landmark_census_panel_v1
# ==========================================================================
def _load_fig4_landmark_census_panel_v1():
    _self = types.ModuleType('fig4_landmark_census_panel_v1')
    __name__ = 'fig4_landmark_census_panel_v1'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_landmark_census_panel_v1'] = _self
    #!/usr/bin/env python3
    """Full-census Figure 4f landmark matrix for a versioned diagnostic candidate.

    This module deliberately contains no landmark-localisation p-value.  It draws every one of
    the 45 regulators with at least one reversal onto the fixed 25-gene Figure 1a landmark roster,
    while retaining the response-shape symbols from the sealed Figure 4 v5.1 production bundle.
    Rows follow the already displayed census diagnostic order; columns retain the fixed Figure 1a
    order.  The matrix is descriptive until the separately frozen matched-localisation analysis
    passes its release gate.
    """



    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42


    P = _MODULES['fig4_production_v5_1']


    MASTER = P.PAPER / "data/census/census_master_edges.csv.gz"
    SUPPORTED = "ESTIMATED_CONTINUOUS_COMPONENTS"

    LOST = "#c1272d"
    GAINED = "#2a9d8f"
    STABLE_EDGE = "#7f8a94"
    INK = "#22262B"

    LANDMARKS = (
        "GZMA", "GZMB", "GNLY", "PRF1", "IFNG", "TNF", "NKG7", "CCL5", "CCL4", "CCL3",
        "CCR5", "CCR1", "CCR7", "CCR2", "IL2", "IL2RA", "IL2RB", "CD69", "CTLA4",
        "IL2RG", "SELL", "TBX21", "KLRB1", "XCL1", "XCL2",
    )

    # Frozen top-to-bottom order from Fig_landmark_matrix_census.  Reversal burden descends; the
    # order within ties is purely descriptive and is pinned here so future pandas sorting changes
    # cannot silently rearrange the panel.
    REGULATORS = (
        "TADA2B", "CCNC", "SUPT7L", "ARNT", "TAF6L", "TADA1", "SMARCB1", "PNPT1",
        "SMARCE1", "ARHGAP30", "MED19", "SUPT20H", "HCCS", "GDAP2", "SENP5", "SGF29",
        "CARMIL2", "TAF13", "TP53INP2", "STAT5B", "SOCS4", "SMG7", "UNC45A", "VAV1",
        "WAC", "YEATS2", "ACTR8", "METAP2", "PRKAR1A", "CXXC1", "ARPC3", "ATP5F1A",
        "BCL3", "CBFB", "COPS8", "COX4I1", "DDA1", "PAXIP1", "E2F1", "EIF4G2", "INO80",
        "IST1", "IWS1", "MED12", "ZAP70",
    )

    SOURCE_COLUMNS = (
        "regulator_order_top_to_bottom",
        "landmark_order_left_to_right",
        "regulator",
        "target",
        "cell_state",
        "eligible_both_endpoints",
        "is_reversal",
        "orientation",
        "response_shape",
        "response_shape_supported",
        "response_analysis_status",
        "regulator_landmark_reversal_burden",
        "target_landmark_reversal_burden",
    )


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def build_source_table(bundle: P.ProductionBundle) -> tuple[pd.DataFrame, dict[str, int]]:
        """Build and validate the literal 45-by-25 display table."""
        master = pd.read_csv(MASTER)
        trans = master.loc[~master["is_self"].astype(bool)].copy()
        landmark_edges = trans.loc[trans["target"].isin(LANDMARKS)].copy()
        display_edges = landmark_edges.loc[landmark_edges["regulator"].isin(REGULATORS)].copy()
        if display_edges[["regulator", "target"]].duplicated().any():
            raise AssertionError("duplicate regulator/landmark display pair in census master")

        reversal = display_edges.loc[display_edges["kind"].eq("rev")]
        observed_regulators = set(reversal["regulator"].astype(str))
        if observed_regulators != set(REGULATORS):
            raise AssertionError(
                f"full landmark regulator roster changed: {sorted(observed_regulators)!r}"
            )

        census_columns = [
            "regulator", "target", "analysis_role", "analysis_status",
            "oriented_delta__detection_component",
            "oriented_delta__positive_intensity_component",
        ]
        # The census retains all 4,379 reversals.  Metadata-ineligible rows deliberately have a
        # missing analysis_role, so filtering on analysis_role would incorrectly erase precisely
        # the unresolved edges whose filled tiles must remain visible.
        response = bundle.census.loc[:, census_columns].copy()
        if response[["regulator", "target"]].duplicated().any():
            raise AssertionError("duplicate reversal pair in response-shape census")

        response_lookup: dict[tuple[str, str], dict[str, str]] = {}
        for _, row in response.iterrows():
            status = str(row["analysis_status"])
            shape = "unresolved"
            if status == SUPPORTED:
                detection = abs(float(row["oriented_delta__detection_component"]))
                intensity = abs(float(row["oriented_delta__positive_intensity_component"]))
                if not (np.isfinite(detection) and np.isfinite(intensity)):
                    raise AssertionError("supported response-shape edge has a non-finite component")
                shape = "detection-dominant" if detection > intensity else "intensity-dominant"
            response_lookup[(str(row["regulator"]), str(row["target"]))] = {
                "status": status,
                "shape": shape,
            }

        edge_lookup = {
            (str(row["regulator"]), str(row["target"])): row
            for _, row in display_edges.iterrows()
        }
        regulator_burden = reversal.groupby("regulator").size().to_dict()
        target_burden = reversal.groupby("target").size().to_dict()

        records: list[dict[str, object]] = []
        for regulator_index, regulator in enumerate(REGULATORS, start=1):
            for target_index, target in enumerate(LANDMARKS, start=1):
                edge = edge_lookup.get((regulator, target))
                if edge is None:
                    state = "outside_both_endpoint_census"
                    eligible = False
                    is_reversal = False
                    orientation = ""
                    shape = "not_applicable"
                    shape_supported = False
                    response_status = ""
                elif str(edge["kind"]) == "stable":
                    state = "comparable_not_reversed"
                    eligible = True
                    is_reversal = False
                    orientation = ""
                    shape = "not_applicable"
                    shape_supported = False
                    response_status = ""
                elif str(edge["kind"]) == "rev":
                    direction = str(edge["direction"])
                    if direction == "lost_on_activation":
                        state = "reversal_lost"
                        orientation = "lost"
                    elif direction == "gained_on_activation":
                        state = "reversal_gained"
                        orientation = "gained"
                    else:
                        raise AssertionError(f"unexpected reversal direction {direction!r}")
                    eligible = True
                    is_reversal = True
                    response_record = response_lookup.get((regulator, target))
                    if response_record is None:
                        raise AssertionError(f"reversal absent from response-shape census: {regulator}->{target}")
                    shape = response_record["shape"]
                    response_status = response_record["status"]
                    shape_supported = response_status == SUPPORTED
                else:
                    raise AssertionError(f"unexpected census kind {edge['kind']!r}")

                records.append({
                    "regulator_order_top_to_bottom": regulator_index,
                    "landmark_order_left_to_right": target_index,
                    "regulator": regulator,
                    "target": target,
                    "cell_state": state,
                    "eligible_both_endpoints": eligible,
                    "is_reversal": is_reversal,
                    "orientation": orientation,
                    "response_shape": shape,
                    "response_shape_supported": shape_supported,
                    "response_analysis_status": response_status,
                    "regulator_landmark_reversal_burden": int(regulator_burden.get(regulator, 0)),
                    "target_landmark_reversal_burden": int(target_burden.get(target, 0)),
                })

        table = pd.DataFrame.from_records(records, columns=SOURCE_COLUMNS)
        state_counts = table["cell_state"].value_counts().to_dict()
        shape_counts = table.loc[
            table["response_shape_supported"].astype(bool), "response_shape"
        ].value_counts().to_dict()
        supported_orientation = table.loc[
            table["response_shape_supported"].astype(bool), "orientation"
        ].value_counts().to_dict()

        summary = {
            "matrix_cells": int(len(table)),
            "displayed_regulators": int(table["regulator"].nunique()),
            "fixed_landmarks": int(table["target"].nunique()),
            "reversal_cells": int(table["is_reversal"].sum()),
            "lost_cells": int(state_counts.get("reversal_lost", 0)),
            "gained_cells": int(state_counts.get("reversal_gained", 0)),
            "stable_cells": int(state_counts.get("comparable_not_reversed", 0)),
            "outside_cells": int(state_counts.get("outside_both_endpoint_census", 0)),
            "landmarks_with_reversal": int(
                table.loc[table["is_reversal"].astype(bool), "target"].nunique()
            ),
            "all_comparable_landmark_pairs": int(len(landmark_edges)),
            "all_regulators_with_landmark_opportunity": int(landmark_edges["regulator"].nunique()),
            "displayed_comparable_cells": int(table["eligible_both_endpoints"].sum()),
            "response_shape_supported": int(table["response_shape_supported"].sum()),
            "intensity_dominant": int(shape_counts.get("intensity-dominant", 0)),
            "detection_dominant": int(shape_counts.get("detection-dominant", 0)),
            "supported_lost": int(supported_orientation.get("lost", 0)),
            "supported_gained": int(supported_orientation.get("gained", 0)),
        }
        expected = {
            "matrix_cells": 1125,
            "displayed_regulators": 45,
            "fixed_landmarks": 25,
            "reversal_cells": 70,
            "lost_cells": 47,
            "gained_cells": 23,
            "stable_cells": 84,
            "outside_cells": 971,
            "landmarks_with_reversal": 16,
            "all_comparable_landmark_pairs": 427,
            "all_regulators_with_landmark_opportunity": 191,
            "displayed_comparable_cells": 154,
            "response_shape_supported": 18,
            "intensity_dominant": 16,
            "detection_dominant": 2,
            "supported_lost": 18,
            "supported_gained": 0,
        }
        if summary != expected:
            raise AssertionError(f"landmark matrix semantic counts changed: {summary!r}")
        return table, summary


    def draw_panel(ax: plt.Axes, table: pd.DataFrame, label_pt: float = 5.5) -> None:
        """Draw the literal source table with square cells and no inferential annotation."""
        n_rows, n_columns = len(REGULATORS), len(LANDMARKS)
        for row in table.itertuples(index=False):
            x = float(row.landmark_order_left_to_right - 1)
            y = float(n_rows - row.regulator_order_top_to_bottom)
            if row.cell_state == "comparable_not_reversed":
                ax.add_patch(Rectangle(
                    (x + 0.15, y + 0.15), 0.70, 0.70,
                    facecolor="none", edgecolor=STABLE_EDGE, lw=0.45,
                ))
            elif row.cell_state in {"reversal_lost", "reversal_gained"}:
                gained = row.cell_state == "reversal_gained"
                ax.add_patch(Rectangle(
                    (x + 0.06, y + 0.06), 0.88, 0.88,
                    facecolor=GAINED if gained else LOST,
                    edgecolor="white", lw=0.35,
                    hatch="////" if gained else None,
                ))
                if row.response_shape == "intensity-dominant":
                    ax.plot(
                        [x + 0.5], [y + 0.5], marker="o", ms=2.9,
                        mfc="white", mec="white", mew=0.25, ls="", zorder=4,
                    )
                elif row.response_shape == "detection-dominant":
                    ax.plot(
                        [x + 0.5], [y + 0.5], marker="^", ms=3.3,
                        mfc="white", mec="white", mew=0.25, ls="", zorder=4,
                    )

        ax.set_xlim(-0.05, n_columns + 0.05)
        ax.set_ylim(-0.05, n_rows + 0.05)
        ax.set_aspect("equal", adjustable="box")
        ax.set_anchor("C")
        ax.set_xticks(np.arange(n_columns) + 0.5)
        ax.set_xticklabels(LANDMARKS, fontsize=label_pt, rotation=90)
        ax.set_yticks([n_rows - index - 0.5 for index in range(n_rows)])
        ax.set_yticklabels(REGULATORS, fontsize=label_pt)
        ax.xaxis.set_ticks_position("top")
        ax.tick_params(length=0, pad=1.0)
        for label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
            label.set_color("#565B61")
        for spine in ax.spines.values():
            spine.set_visible(False)


    def legend_handles() -> list[object]:
        return [
            Rectangle((0, 0), 1, 1, facecolor=LOST, edgecolor="white", label="reversal, lost"),
            Rectangle(
                (0, 0), 1, 1, facecolor=GAINED, edgecolor="white", hatch="////",
                label="reversal, gained",
            ),
            Rectangle(
                (0, 0), 1, 1, facecolor="none", edgecolor=STABLE_EDGE,
                label="comparable, not reversed",
            ),
            Line2D(
                [0], [0], marker="o", ls="", mfc=INK, mec=INK, ms=3.2,
                label="intensity-dominant",
            ),
            Line2D(
                [0], [0], marker="^", ls="", mfc=INK, mec=INK, ms=3.6,
                label="detection-dominant",
            ),
        ]
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_landmark_census_panel_v1()


# ==========================================================================
# fig4_presentation_v5_1
# ==========================================================================
def _load_fig4_presentation_v5_1():
    _self = types.ModuleType('fig4_presentation_v5_1')
    __name__ = 'fig4_presentation_v5_1'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_presentation_v5_1'] = _self
    #!/usr/bin/env python3
    """Presentation render of Figure 4 from the audited v5.1 production bundle.

    WHY THIS IS A SEPARATE FILE. ``fig4_production_v5_1.py`` is byte-pinned: its own SHA-256 is
    recorded as ``builder_sha256`` inside ``figure4_production_manifest_v5_1.json`` and re-checked
    on every run, so editing a single axis limit in it makes the manifest reject the build. That
    pin is the right design for an audited estimator, and nothing here weakens it. This module
    imports the frozen builder untouched, calls its ``load_production_bundle`` so the full
    validation chain still runs, reuses its data accessors so no statistic is recomputed, and
    changes only geometry, axis limits and labels.

    WHAT CHANGES, AND WHY EACH ONE WAS WRONG BEFORE.
      * Scatter limits came from ``_symmetric_limit``, which uses the single most extreme value.
        In panel c one regulator-matched stable comparison reaches |29.7| while the 99th
        percentile of the cloud is 4.8, so the axis ran to +/-30 and every one of the 5,815 points
        collapsed onto a needle at the origin. The panel exists to show how detection and
        positive-cell intensity divide the change, and at that scale it showed nothing. Limits are
        now cut at the 99.5th percentile of the pooled absolute coordinates and each panel states
        how many points fall outside, so the tail is declared rather than hidden.
      * Panel c's y-label, "oriented change, positive-cell intensity", is taller than its
        equal-aspect box and overprinted the panel's own x-label.
      * The leave-one-donor axes carried no panel letter, so six panels shipped under five labels.
      * The CAL table repeated an identical "n=1,726" in all twelve cells, which is three times
        the ink of the numbers that actually differ. It moves to one footnote, and stays per cell
        if the counts ever diverge.
      * "mixed/unresolved: 11" floated in panel a's corner reading as an unrelated statistic
        instead of as the 1,910-to-1,899 step it explains.

    The equal aspect on panels b, c and e is deliberate and kept: detection and intensity are in
    the same units, so their relative spread is readable as geometry, and that comparison is the
    claim the panel supports.

    Outputs ``Fig4.pdf``/``Fig4.png`` plus ``Fig4.presentation.json``, which chains to the frozen
    manifest and to the SHA-256 of the audited production render this figure is derived from.
    """




    matplotlib.use("Agg")


    _HERE = Path(__file__).resolve().parent
    for _path in (str(_HERE), str(_HERE.parent)):
        if _path not in sys.path:
            sys.path.insert(0, _path)

    P = _MODULES['fig4_production_v5_1']
    EL = _MODULES['fig4_effector_landmark_panel']

    FIGURE_WIDTH_IN = P.FIGURE_WIDTH_IN   # frozen builder requires exactly 7.10 in
    # A4 proportions at the mandated width: 7.10 x 7.10*(297/210). The builder refuses any
    # other width, so the page is A4-shaped rather than A4-sized and prints to A4 with margins.
    FIGURE_HEIGHT_IN = 10.04
    TEXT_PT = P.MINIMUM_FINAL_TEXT_PT
    CONNECTOR_COLOUR = "#173F5F"
    CONNECTOR_TARGETS = ("BCL2L1", "DDIT4", "GPR108", "IL2RA", "PSMB9", "PSME1",
                         "PSME2", "SDCBP")
    # Cut the scatter axes here rather than at the maximum, so a sub-1% tail cannot flatten the
    # cloud onto the origin. Every point past the cut is counted and declared on its panel.
    DISPLAY_QUANTILE = 0.995
    DISPLAY_PAD = 1.10
    ALIGNMENT_BOOTSTRAP_SEED = 20260824
    ALIGNMENT_BOOTSTRAP_DRAWS = 2000

    DEFAULT_PDF = P.PAPER / "figures/Fig4.pdf"
    DEFAULT_AUDIT_DIR = P.PAPER / "figures/audit"


    # --------------------------------------------------------------------------------------
    # shared geometry helpers
    # --------------------------------------------------------------------------------------
    def display_limit(
        pairs: Sequence[tuple[Iterable[float], Iterable[float]]],
        quantile: float = DISPLAY_QUANTILE,
        pad: float = DISPLAY_PAD,
    ) -> tuple[float, int]:
        """Symmetric limit at a robust quantile, with the count of points outside it.

        ``pairs`` are the (x, y) series actually drawn on the axes. A point counts as off-scale
        when either coordinate leaves the box, which is what the reader sees.
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


    def note_off_scale(ax: plt.Axes, off_scale: int, limit: float, total: int,
                       corner: str = "lower right") -> None:
        """State the clipped tail in a free corner, in DATA coordinates.

        These panels are ``set_aspect("equal", adjustable="box")``. That shrinks the drawn box
        inside the axes rectangle while ``transAxes`` still spans the full rectangle, so an
        axes-fraction annotation at 0.98 lands in the margin outside the plot. Data coordinates
        are tied to the box, so they stay put whatever the aspect does to it. The caller picks the
        corner its own cloud leaves empty; the opaque box keeps the identity diagonal from
        striking through the digits.
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


    def corner_legend(ax: plt.Axes, limit: float, corner: str, **kwargs):
        """Legend pinned to a data-coordinate corner, for the same aspect reason as above."""
        x = -limit * 0.97 if corner.endswith("left") else limit * 0.97
        y = limit * 0.97 if corner.startswith("upper") else -limit * 0.97
        return P._opaque_legend(ax, loc=corner, bbox_to_anchor=(x, y),
                                bbox_transform=ax.transData, fontsize=TEXT_PT, **kwargs)


    def multiplicity_preserving_cluster_rates(
        cluster_totals: pd.DataFrame,
        sampled_clusters: np.ndarray,
    ) -> np.ndarray:
        """Pool complete sampled clusters while retaining replacement multiplicity.

        ``cluster_totals`` must contain one row per cluster with integer ``successes`` and
        ``edges`` columns. ``sampled_clusters`` is a two-dimensional array whose rows are
        bootstrap samples. Selecting with ``.loc`` deliberately retains duplicate labels, unlike
        membership filtering with ``.isin``, which collapses a regulator sampled twice or more.
        """
        if sampled_clusters.ndim != 2:
            raise ValueError("sampled_clusters must be a two-dimensional array")
        if not cluster_totals.index.is_unique:
            raise ValueError("cluster_totals must have one row per unique cluster")
        required = {"successes", "edges"}
        if not required.issubset(cluster_totals.columns):
            raise ValueError("cluster_totals must contain successes and edges columns")
        selected = cluster_totals.loc[sampled_clusters.reshape(-1)]
        successes = selected["successes"].to_numpy(float).reshape(sampled_clusters.shape)
        edges = selected["edges"].to_numpy(float).reshape(sampled_clusters.shape)
        denominators = edges.sum(axis=1)
        if np.any(denominators <= 0):
            raise ValueError("every bootstrap draw must contain at least one edge")
        return successes.sum(axis=1) / denominators


    def regulator_cluster_bootstrap_rates(
        frame: pd.DataFrame,
        rng: np.random.Generator,
        n_draws: int = ALIGNMENT_BOOTSTRAP_DRAWS,
    ) -> np.ndarray:
        """Edge-pooled rates after sampling regulator clusters with replacement."""
        if n_draws <= 0:
            raise ValueError("n_draws must be positive")
        cluster_totals = frame.groupby("regulator_ensg")["_aligned"].agg(
            successes="sum", edges="size")
        if cluster_totals.empty:
            raise ValueError("regulator-cluster bootstrap requires at least one regulator")
        regulators = cluster_totals.index.to_numpy()
        sampled = rng.choice(
            regulators, size=(n_draws, len(regulators)), replace=True)
        return multiplicity_preserving_cluster_rates(cluster_totals, sampled)


    # --------------------------------------------------------------------------------------
    # panels
    # --------------------------------------------------------------------------------------
    def panel_waterfall(ax: plt.Axes, bundle: P.ProductionBundle) -> None:
        counts = {status: int(bundle.census["analysis_status"].eq(status).sum()) for status in P.FINAL_STATUSES}
        values = np.asarray((
            P.DENOMINATOR,
            P.METADATA_ELIGIBLE,
            P.METADATA_ELIGIBLE - counts["UNANALYSABLE_POSITIVE_INTENSITY"],
            counts["ESTIMATED_CONTINUOUS_COMPONENTS"] + counts["MIXED_OR_UNRESOLVED"],
            counts["ESTIMATED_CONTINUOUS_COMPONENTS"],
        ))
        labels = ("census", "metadata", "donor\nsupport", "effect", "raw\ncompatible")
        x = np.arange(len(values))
        ax.bar(x, values, color=("#46515C", "#5A7184", "#4E8D8A", "#2A9D8F", "#216E5C"), width=0.68)
        ax.set_xticks(x, labels, fontsize=TEXT_PT)
        ax.set_ylabel("reversal edges")
        ax.set_ylim(0, P.DENOMINATOR * 1.16)
        # upright: at square-panel width "1,910" is wider than the bar pitch and the last three
        # values ran into one another.
        for xx, value in zip(x, values):
            ax.text(xx, value + P.DENOMINATOR * 0.025, f"{value:,}", ha="left", va="center",
                    fontsize=7, rotation=90, rotation_mode="anchor")
        # three short lines, not two long ones: the panel is square now and the old wording ran
        # past its right edge into panel b's heading.
        ax.text(0.99, 0.93, f"effect to raw\ncompatible:\n{counts['MIXED_OR_UNRESOLVED']:,} mixed",
                transform=ax.transAxes, ha="right", va="top", fontsize=TEXT_PT,
                color="#666B70", linespacing=1.35)
        ax.set_title("Full-denominator eligibility", loc="left", pad=7)


    def panel_rest_stim(ax: plt.Axes, bundle: P.ProductionBundle) -> None:
        reversal = P._reversal_summary(bundle)
        table = reversal.loc[reversal["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        census_direction = bundle.census.set_index("reversal_edge_id")["direction"]
        gained = table["analysis_edge_id"].map(census_direction).eq("gained_on_activation")
        for mask, colour, marker, label in ((gained, "#2A9D8F", "o", "gained"), (~gained, "#C1272D", "^", "lost")):
            block = table.loc[mask]
            ax.scatter(block["Rest__total_effect"], block["Stim48hr__total_effect"], s=7, alpha=0.30,
                       color=colour, marker=marker, edgecolors="none", rasterized=True,
                       label=f"{label} (n={len(block):,})")
        limit, off_scale = display_limit([(table["Rest__total_effect"], table["Stim48hr__total_effect"])])
        ax.plot([-limit, limit], [-limit, limit], color="#C5C9CD", lw=0.8)
        ax.axhline(0, color="#9AA0A6", lw=0.7)
        ax.axvline(0, color="#9AA0A6", lw=0.7)
        for condition, orientation in (("Rest", "vertical"), ("Stim48hr", "horizontal")):
            threshold = P._envelope(bundle, "total_effect", condition)
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
        # gained clusters upper-left and lost lower-right, so both diagonal corners carry data;
        # the legend goes upper-right, which the cloud leaves empty, and the count note stays
        # lower-left.
        # Direct labels rather than a legend box. Both diagonal corners carry data, so any framed
        # legend sits on points; the two clusters are far enough apart to name in place.
        counts = {"gained": int(gained.sum()), "lost": int((~gained).sum())}
        ax.text(-limit * 0.92, limit * 0.86, f"gained (n={counts['gained']:,})", color="#2A9D8F",
                fontsize=TEXT_PT, fontweight="bold", ha="left", va="center")
        ax.text(limit * 0.92, -limit * 0.86, f"lost (n={counts['lost']:,})", color="#C1272D",
                fontsize=TEXT_PT, fontweight="bold", ha="right", va="center")


    def panel_components(ax: plt.Axes, bundle: P.ProductionBundle) -> None:
        """Every estimated edge, on symmetric-log axes, with the extremes named.

        The earlier version cut both axes at the 99.5th percentile and declared the 47 points it
        dropped. That kept the core cloud readable but hid the tail, and the tail is real: detection
        runs to 5.09 and intensity to 29.69 while the 99th percentiles are 0.64 and 4.78. Symmetric
        log axes with a linear window of +/-0.5 show every point without flattening the cloud, so
        nothing is omitted and no reader has to take the outliers on trust. Both axes carry the same
        transform, so the detection-versus-intensity comparison the panel exists for still reads as
        geometry.
        """
        summary = bundle.edge_summary
        det = "oriented_delta__detection_component"
        inten = "oriented_delta__positive_intensity_component"
        reversal = summary.loc[summary["analysis_role"].eq("REVERSAL")
                               & summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")].copy()
        stable = summary.loc[summary["analysis_role"].eq("REGULATOR_MATCHED_STABLE")
                             & summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")]
        mixed = reversal["analysis_status"].eq("MIXED_OR_UNRESOLVED")

        ax.scatter(stable[det], stable[inten], s=3, alpha=0.10, color="#9AA0A6", edgecolors="none",
                   rasterized=True)
        ax.scatter(reversal.loc[~mixed, det], reversal.loc[~mixed, inten], s=7, alpha=0.34,
                   color="#2A9D8F", edgecolors="none", rasterized=True)
        # mixed/unresolved are 11 edges out of 1,910. Black crosses so they read as a distinct
        # category at this density rather than as pale outliers of the estimated cloud.
        ax.scatter(reversal.loc[mixed, det], reversal.loc[mixed, inten], s=16, color="#000000",
                   marker="x", linewidths=0.9, zorder=4)
        for metric, orientation in (("detection_component", "vertical"),
                                    ("positive_intensity_component", "horizontal")):
            threshold = P._envelope(bundle, metric, "delta")
            if threshold is not None:
                for sign in (-1, 1):
                    (ax.axvline if orientation == "vertical" else ax.axhline)(
                        sign * threshold, color="#767D84", lw=0.7, ls="--")
        ax.axhline(0, color="#9AA0A6", lw=0.7)
        ax.axvline(0, color="#9AA0A6", lw=0.7)
        ax.set_xscale("symlog", linthresh=0.5, linscale=0.6)
        ax.set_yscale("symlog", linthresh=0.5, linscale=0.6)
        ax.set_xlim(-1.2, 7.0)
        ax.set_ylim(-6.0, 120.0)
        ax.set_xticks([-1, 0, 1, 5])
        ax.set_yticks([-5, -1, 0, 1, 5, 10, 30])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:g}"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:g}"))
        ax.set_xlabel("oriented change, detection")
        ax.set_ylabel("oriented change, intensity")
        ax.set_title("Component split", loc="left", pad=7)

        # The three furthest-from-origin edge names are dropped. They were chosen by a fixed rule
        # rather than cherry-picked, but they name three edges out of 5,815 and the panel is about
        # the distribution, not about which edge is extreme. Removing them also returns the space
        # the square layout needs.
        handles = [
            Line2D([], [], color="#9AA0A6", marker="o", linestyle="none", markersize=4, label="matched stable"),
            Line2D([], [], color="#2A9D8F", marker="o", linestyle="none", markersize=4, label="estimated"),
            Line2D([], [], color="#000000", marker="x", linestyle="none", markersize=4,
                   markeredgewidth=1.0, label="mixed"),
        ]
        P._opaque_legend(ax, handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
                         ncol=3, fontsize=TEXT_PT, borderaxespad=0.0, columnspacing=0.8,
                         handletextpad=0.3)
        ax.text(0.02, 0.02, "symlog axes, linear within \u00b10.5\nall 5,815 edges drawn",
                transform=ax.transAxes, va="bottom", fontsize=TEXT_PT, color="#666B70", linespacing=1.35)


    def panel_target_dominance(ax: plt.Axes, bundle: P.ProductionBundle, n_targets: int = 28,
                               min_edges: int = 3) -> None:
        """Which target genes receive detection- versus intensity-dominant reversals.

        Whole-set and descriptive: targets are ranked by how many supported reversals they carry,
        not by the split being displayed, so gene identity appears without selecting on the outcome.
        """
        summary = bundle.edge_summary
        det = "oriented_delta__detection_component"
        inten = "oriented_delta__positive_intensity_component"
        rev = summary.loc[summary["analysis_role"].eq("REVERSAL")
                          & summary["analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")].copy()
        d = pd.to_numeric(rev[det], errors="coerce").abs()
        i = pd.to_numeric(rev[inten], errors="coerce").abs()
        rev["_det_dom"] = (d > i).astype(int)
        rev["_target"] = rev["analysis_edge_id"].astype(str).map(
            bundle.census.set_index("reversal_edge_id")["target"])
        rev = rev.dropna(subset=["_target"])
        grouped = rev.groupby("_target")["_det_dom"].agg(["size", "sum"])
        grouped = grouped.loc[grouped["size"] >= min_edges].nlargest(n_targets, "size").sort_values("size")
        y = np.arange(len(grouped))
        detection = grouped["sum"].to_numpy(float)
        intensity = (grouped["size"] - grouped["sum"]).to_numpy(float)
        # stacked from zero rather than mirrored about it: the column is too narrow to carry a
        # diverging axis without its tick labels colliding
        ax.barh(y, detection, color="#2A9D8F", height=0.7)
        ax.barh(y, intensity, left=detection, color="#7B6FA6", height=0.7)
        # Figure 3d connector, ported verbatim from the promoted candidate renderer: a diamond
        # prefix, navy bold label and a direct in-panel key on the eight displayed targets that
        # also belong to the exhaustive CCNC-TADA2B-SUPT7L intersection. Membership plays no part
        # in choosing or ordering the 28.
        marked = [target for target in grouped.index if target in CONNECTOR_TARGETS]
        assert len(marked) == 8, f"connector marks {len(marked)} of the expected 8 targets"
        ax.set_yticks(y)
        ax.set_yticklabels([f"\u25c6 {t}" if t in CONNECTOR_TARGETS else t for t in grouped.index],
                           fontsize=TEXT_PT)
        for tick, target in zip(ax.get_yticklabels(), grouped.index):
            if target in CONNECTOR_TARGETS:
                tick.set_color(CONNECTOR_COLOUR)
                tick.set_fontweight("bold")
        ax.text(1.0, 1.012, "\u25c6 also in Fig. 3d trio", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=TEXT_PT, color=CONNECTOR_COLOUR, fontweight="bold")
        ax.set_ylim(-0.9, len(grouped) - 0.1)
        # the column is narrow: a tick every 4 keeps the mirrored labels from colliding
        span = int(grouped["size"].max())
        ax.set_xticks([v for v in range(0, span + 1, 4)])
        ax.set_xlim(0, span + 4.2)
        ax.set_xlabel("reversals per target")
        ax.set_title("Target response", loc="left", pad=7)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        # direct labels rather than a legend box: the legend sat on the shortest bar, and the two
        # pure bars name themselves
        pure_int = [k for k, row in enumerate(grouped.itertuples()) if row.sum == 0]
        pure_det = [k for k, row in enumerate(grouped.itertuples()) if row.sum == row.size]
        if pure_int:
            k = max(pure_int, key=lambda idx: grouped["size"].iloc[idx])
            ax.text(grouped["size"].iloc[k] / 2, k, "intensity-dominant", ha="center", va="center",
                    fontsize=TEXT_PT, color="white", fontweight="bold")
        if pure_det:
            k = max(pure_det, key=lambda idx: grouped["size"].iloc[idx])
            # this bar is too short to hold the words, so the label sits just past its end
            ax.text(grouped["size"].iloc[k] + 0.25, k, "detection-dominant", ha="left", va="center",
                    fontsize=TEXT_PT, color="#2A9D8F", fontweight="bold")


    def _cal_cells(envelopes: pd.DataFrame) -> tuple[list[list[str]], int | None]:
        """Envelope values, plus the shared CAL parent count when every cell agrees."""
        index = envelopes.set_index(["metric", "scope", "condition"])
        counts: set[int] = set()
        values: list[list[tuple[str, int]]] = []
        for metric in P.METRICS:
            row_values: list[tuple[str, int]] = []
            for condition in ("Rest", "Stim48hr", "delta"):
                scope = "cross_condition" if condition == "delta" else "condition"
                row = index.loc[(metric, scope, condition)]
                count = int(row["n_unique_supported_CAL"])
                counts.add(count)
                text = f"{float(row['absolute_q95']):.3g}" if bool(row["available"]) else "unavailable"
                row_values.append((text, count))
            values.append(row_values)
        shared = counts.pop() if len(counts) == 1 else None
        cells = [[text if shared is not None else f"{text}\nn={count:,}" for text, count in row]
                 for row in values]
        return cells, shared


    def panel_cal_dots(ax: plt.Axes, envelopes: pd.DataFrame) -> None:
        """CAL envelopes as a dot plot instead of a table of twelve numbers.

        Every cell rests on the same 1,726 CAL parents, so the twelve values differ only in metric
        and condition. A table spends a whole column saying that; the dot plot says it at a glance
        and returns the column to a figure. Same numbers, same source, no new quantity.
        """
        index = envelopes.set_index(["metric", "scope", "condition"])
        rows = ("total_effect", "detection_component", "positive_intensity_component",
                "raw_rate_log2_ratio")
        names = ("total", "detection", "positive intensity", "raw UMI rate")
        marks = (("Rest", "condition", "o", S.PAL["REST"], "Rest"),
                 ("Stim48hr", "condition", "s", S.PAL["STIM"], "48 h"),
                 ("delta", "cross_condition", "D", "#22262B", "\u0394"))
        counts = set()
        for i, metric in enumerate(rows):
            y = len(rows) - 1 - i
            ax.plot([0, 0.5], [y, y], color="#e8ebee", lw=0.8, zorder=0)
            for condition, scope, marker, colour, _ in marks:
                row = index.loc[(metric, scope, condition)]
                counts.add(int(row["n_unique_supported_CAL"]))
                ax.plot([float(row["absolute_q95"])], [y], marker=marker, ms=4.4, ls="",
                        mfc=colour, mec="white", mew=0.5, zorder=3)
        assert len(counts) == 1, "CAL cells no longer share one parent count"
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(list(names)[::-1], fontsize=TEXT_PT)
        # -1.15, not -0.6: moving the parent-count note left was not enough, because the bottom row
        # sits only 0.6 units above the floor and the two-line note still reached into it. The extra
        # margin is what actually clears the raw-UMI markers.
        ax.set_ylim(-1.15, len(rows) - 0.4)
        ax.set_xlim(0, 0.52)
        ax.set_xticks([0, 0.2, 0.4])
        ax.set_xlabel("absolute 95th percentile", fontsize=TEXT_PT, labelpad=1)
        ax.set_title("CAL envelopes", loc="left", pad=7)
        ax.tick_params(axis="y", length=0)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        # Direct labels over the first row rather than a legend box: at this panel width the box
        # sat on the raw-UMI-rate row it was meant to explain.
        top_y = len(rows) - 1
        for condition, scope, _marker, colour, label in marks:
            row = index.loc[(rows[0], scope, condition)]
            ax.annotate(label, (float(row["absolute_q95"]), top_y), textcoords="offset points",
                        xytext=(0, 7), ha="center", va="bottom", fontsize=TEXT_PT, color=colour)
        # two lines, not one: panel c is the narrow right-hand column and the single-line form
        # ran off the page edge, which the geometry gate rejects.
        # inside the axes, bottom right: below the panel this note landed on panel e's heading.
        # The corner is empty because the largest envelope sits on the top row.
        ax.text(0.02, 0.03, f"n = {counts.pop():,} CAL parents\nmetric-specific units",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=TEXT_PT,
                color="#565B61", linespacing=1.4)


    def panel_cal(ax: plt.Axes, envelopes: pd.DataFrame) -> None:
        ax.axis("off")
        cells, shared_n = _cal_cells(envelopes)
        table = ax.table(
            cellText=cells,
            rowLabels=("total", "detection", "positive intensity", "raw UMI rate"),
            colLabels=("Rest", "48 h", "\u0394"),
            cellLoc="center",
            rowLoc="right",
            colLoc="center",
            # rowLabels are drawn to the left of the bbox, so the table is inset far enough
            # that "positive intensity" stays inside the axes instead of running off the page.
            bbox=(0.26, 0.20, 0.74, 0.68) if shared_n is not None else (0.26, 0.10, 0.74, 0.80),
        )
        table.auto_set_font_size(False)
        table.set_fontsize(TEXT_PT)
        for (row, _column), cell in table.get_celld().items():
            cell.set_edgecolor("#D9DDE1")
            cell.set_linewidth(0.55)
            cell.set_facecolor("#EAF2F7" if row > 0 and row % 2 else "#F7F9FA")
            if row == 0:
                cell.set_facecolor("#E1E7EB")
                cell.get_text().set_fontweight("bold")
        ax.set_title("CAL envelopes", loc="left", pad=7)
        footnote = "\u0394 is 48 h \u2212 Rest.\nMetric-specific units; compare within rows."
        if shared_n is not None:
            footnote = f"All envelopes rest on n = {shared_n:,} CAL parents.\n{footnote}"
        ax.text(0.0, 0.02, footnote, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=TEXT_PT, color="#565B61", linespacing=1.35)


    def panel_stability(fig: plt.Figure, spec, bundle: P.ProductionBundle) -> list[plt.Axes]:
        grid = spec.subgridspec(1, 2, wspace=0.42)
        ax_guide = fig.add_subplot(grid[0, 0])
        ax_lodo = fig.add_subplot(grid[0, 1])
        guides, lodo = P._stability(bundle)

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
        # Donor identity alone on the ticks, with the operation named once as the axis label. The
        # repeated two-line "leave out Dn" spent a third of the panel restating the title.
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


    # --------------------------------------------------------------------------------------
    # render
    # --------------------------------------------------------------------------------------

    def panel_alignment_null(ax: plt.Axes, bundle: P.ProductionBundle) -> dict:
        """Orientation agreement of the two components, against the matched-stable comparator.

        The headline "both components follow the reversal orientation" is partly favoured by
        construction, because the components sum to the oriented total effect. The regulator-matched
        stable comparisons are oriented by their parent reversal and carry no such expectation, so
        they give the comparator that turns a descriptive share into a calibrated one. Per-regulator
        points are plotted, not just the pooled rate, because edges cluster within regulators.
        """
        summary = bundle.edge_summary

        def _aligned(frame: pd.DataFrame) -> pd.DataFrame:
            det = pd.to_numeric(frame["oriented_delta__detection_component"], errors="coerce")
            inten = pd.to_numeric(frame["oriented_delta__positive_intensity_component"], errors="coerce")
            keep = det.notna() & inten.notna()
            return frame.loc[keep].assign(_aligned=((det > 0) & (inten > 0))[keep].astype(int))

        reversal = _aligned(summary.loc[summary["analysis_role"].eq("REVERSAL")
                                        & summary["analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")])
        stable = _aligned(summary.loc[summary["analysis_role"].eq("REGULATOR_MATCHED_STABLE")
                                      & summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")])

        rng = np.random.default_rng(ALIGNMENT_BOOTSTRAP_SEED)
        series = (("matched stable", stable, "#9AA0A6", 0.0),
                  ("reversals", reversal, "#2A9D8F", 1.0))
        stats: dict[str, dict[str, float]] = {}
        for label, frame, colour, ypos in series:
            per = frame.groupby("regulator_ensg")["_aligned"].agg(["mean", "size"])
            jitter = rng.uniform(-0.16, 0.16, len(per))
            ax.scatter(per["mean"], np.full(len(per), ypos) + jitter,
                       s=5.0 + 1.6 * np.sqrt(per["size"].to_numpy(float)),
                       color=colour, alpha=0.42, edgecolors="none", zorder=2)
            pooled = float(frame["_aligned"].mean())
            draws = regulator_cluster_bootstrap_rates(
                frame, rng, n_draws=ALIGNMENT_BOOTSTRAP_DRAWS)
            lo, hi = np.percentile(draws, [2.5, 97.5])
            ax.errorbar(pooled, ypos - 0.30, xerr=[[pooled - lo], [hi - pooled]], fmt="D",
                        markersize=4.4, color="#22262B", ecolor="#22262B", elinewidth=1.0,
                        capsize=2.4, zorder=4)
            ax.annotate(f"{100 * pooled:.1f}%", (pooled, ypos - 0.30), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=TEXT_PT, color="#22262B",
                        fontweight="bold", zorder=5)
            stats[label.replace("\n", " ")] = {"pooled": pooled, "ci_low": float(lo), "ci_high": float(hi),
                                                "edges": int(len(frame)), "regulators": int(len(per))}
        # Direct labels instead of y tick labels: at this panel width a two-word tick label runs
        # off the page, and the label belongs beside its own row of points anyway.
        ax.set_yticks([])
        for label, _frame, colour, ypos in series:
            ax.text(0.005, ypos + 0.33, label, fontsize=TEXT_PT, color=colour,
                    fontweight="bold", ha="left", va="center")
        ax.set_ylim(-0.55, 1.62)
        ax.set_xlim(-0.02, 1.02)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "25", "50", "75", "100%"])
        ax.set_xlabel("both components aligned (%)")
        ax.set_title("Orientation agreement", loc="left", pad=7)
        ax.grid(axis="x", color="#E4E7EA", lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        higher = int(sum(1 for reg in set(reversal["regulator_ensg"]) & set(stable["regulator_ensg"])
                         if reversal.loc[reversal["regulator_ensg"].eq(reg), "_aligned"].mean()
                         > stable.loc[stable["regulator_ensg"].eq(reg), "_aligned"].mean()))
        shared = len(set(reversal["regulator_ensg"]) & set(stable["regulator_ensg"]))
        stats["regulators_higher_for_reversals"] = {"count": higher, "of": shared}
        return stats


    def render(bundle: P.ProductionBundle, output_pdf: Path, audited_png: Path | None) -> dict:
        output_pdf = Path(output_pdf)
        output_png = output_pdf.with_suffix(".png")
        output_record = output_pdf.with_suffix(".presentation.json")
        output_pdf.parent.mkdir(parents=True, exist_ok=True)

        supported_examples = P._supported_predeclared_example_ids(bundle)
        if supported_examples:
            raise RuntimeError(
                "predeclared example edges are supported in this bundle; the presentation layout "
                "covers panels a-e only and must be extended before it can ship that row"
            )

        S.setup(1.0)
        fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white")
        # Four panels. The Rest-vs-48h reconstruction, the physical-guide agreement and the
        # leave-one-donor agreement moved to Supplementary Figure S3g–i: the first repeats the
        # geometry of Figure 1b in different units, and the other two repeat the guide and donor
        # robustness already carried by Figure 1g,h. What is left appears nowhere else in the paper.
        # Reading order is the argument: what survived the gates, the decomposition, what counts
        # as large, whether the components follow the reversal, and where it lands on T-cell
        # biology. That puts a, b, c across the top and the two tall panels below, so no letter
        # changes meaning except e. The CAL table became a dot plot precisely so it fits a top-row
        # slot instead of leaving four inches blank in a tall column.
        grid = fig.add_gridspec(
            2, 1, left=0.105, right=0.960, top=0.935, bottom=0.065,
            hspace=0.45, height_ratios=(1.69, 6.45),
        )
        top_row = grid[0].subgridspec(1, 3, wspace=0.30)
        bottom_row = grid[1].subgridspec(1, 2, wspace=0.34)
        left_stack = bottom_row[0, 0].subgridspec(2, 1, hspace=0.42, height_ratios=(2.59, 3.40))
        ax_a = fig.add_subplot(top_row[0, 0])
        ax_b = fig.add_subplot(top_row[0, 1])
        ax_c = fig.add_subplot(top_row[0, 2])
        ax_d = fig.add_subplot(left_stack[0, 0])
        ax_e = fig.add_subplot(left_stack[1, 0])
        ax_f = fig.add_subplot(bottom_row[0, 1])
        panel_waterfall(ax_a, bundle)
        # single-line stage names: the two-line "raw\ncompatible" still collided with its
        # neighbour once rotated, and the stage meaning survives the shorter wording
        ax_a.set_xticks(range(5))
        ax_a.set_xticklabels(["census", "metadata", "donors", "effect", "raw UMI"])
        # Upright, hanging below the axis: five stage names cannot lie at 28 degrees across a
        # square panel this narrow, and ha="right"/va="center" is what keeps them outside the axes
        # rather than running up into the bars.
        for label in ax_a.get_xticklabels():
            label.set_rotation(90)
            label.set_ha("right")
            label.set_va("center")
            label.set_rotation_mode("anchor")
        panel_components(ax_b, bundle)
        # c and d swapped. CAL's metric names need a 2.59 in square, which is the lower slot; the
        # two-arm agreement scatter reads fine in the narrow one. Reading order still works: the
        # decomposition, then whether it follows the reversal, then the envelopes that calibrate it.
        alignment_stats = panel_alignment_null(ax_c, bundle)
        panel_cal_dots(ax_d, bundle.cal_envelopes)
        panel_target_dominance(ax_e, bundle)
        # f is the fixed T-cell landmark matrix, the only panel that names the genes the paper is
        # about. Rectangular tiles, because a square-tile matrix shrinks inside its slot.
        EL.panel_effector_landmarks(ax_f, EL.load_cells(), fs=TEXT_PT, equal_aspect=False)
        ax_f.set_title("Modules on T-cell landmarks", loc="left", pad=8)
        ax_f.legend(handles=EL.legend_handles(), loc="center", # 0.50, not 0.60: shifted right the two-column key runs off the page. The landmark
                    # row labels sit outside the axes, so a centred key is adjacent to them, not over them.
                    bbox_to_anchor=(0.50, 0.40),
                    ncol=2, fontsize=TEXT_PT, frameon=False, handlelength=1.3,
                    columnspacing=1.0, labelspacing=0.32)
        for _square in (ax_a, ax_b, ax_c, ax_d):
            _square.set_box_aspect(1.0)
        # Panel headings at the 7 pt floor. At the default 10 pt they are wider than a 1.69 in
        # square, which is what put "Full-denominator eligibility" under panel b's heading and ran
        # "Orientation agreement" to the page edge. Bold still separates them from body text.
        for _ax in (ax_a, ax_b, ax_c, ax_d, ax_e, ax_f):
            _ax.title.set_fontsize(TEXT_PT)

        for label, ax in zip("abcde", (ax_a, ax_b, ax_c, ax_d, ax_e)):
            S.panel_tag(ax, label, dx=-0.075, dy=1.09)
        S.panel_tag(ax_f, "f", dx=-0.135, dy=1.035)

        fig.suptitle("Donor-balanced target-response shapes across stimulation",
                     x=0.07, y=0.972, ha="left", fontsize=13.0, fontweight="bold", color="#22262B")
        # The former in-figure footnote moved to the Figure 4 legend: every line of it was
        # either already in the legend or belongs there rather than on the canvas.

        # Reuse the frozen builder's own geometry gate: exact page width, no sub-7 pt text, no text
        # off the page, no overlapping titles, tick labels or title-annotation collisions.
        minimum_text_pt = P._assert_native_geometry_and_text(fig)

        with tempfile.TemporaryDirectory(prefix="fig4_presentation_", dir=output_pdf.parent) as temporary:
            temporary = Path(temporary)
            temp_pdf, temp_png = temporary / output_pdf.name, temporary / output_png.name
            fig.savefig(temp_pdf, facecolor="white")
            fig.savefig(temp_png, dpi=300, facecolor="white")
            plt.close(fig)
            record = {
                "schema": "response-shape-figure4-presentation-v5.1.1",
                "status": "COMPLETE",
                "derived_from": "fig4_production_v5_1.py (frozen, unmodified)",
                "presentation_builder_sha256": P.sha256_file(Path(__file__)),
                "production_builder_sha256": P.sha256_file(P.__file__ if isinstance(P.__file__, Path) else Path(P.__file__)),
                "source_manifest_sha256": bundle.manifest_sha256,
                "methods_config_sha256": P.FROZEN_CONFIG_SHA256,
                "amendment_sha256": P.FROZEN_AMENDMENT_SHA256,
                "claim_boundary": P.CLAIM_BOUNDARY,
                "input_artifact_sha256": dict(bundle.artifact_sha256),
                "audited_render_png_sha256": P.sha256_file(audited_png) if audited_png and audited_png.is_file() else None,
                "predeclared_examples_supported": P._predeclared_example_support_label(supported_examples),
                "predeclared_examples_panel_rendered": False,
                "alignment_null": alignment_stats,
                "alignment_bootstrap": {
                    "seed": ALIGNMENT_BOOTSTRAP_SEED,
                    "draws": ALIGNMENT_BOOTSTRAP_DRAWS,
                    "cluster_unit": "regulator_ensg",
                    "sampling": "regulator clusters sampled with replacement",
                    "pooling": "all edges in each sampled cluster, retaining sample multiplicity",
                    "interval": "2.5th and 97.5th percentiles",
                },
                "figure_width_in": FIGURE_WIDTH_IN,
                "figure_height_in": FIGURE_HEIGHT_IN,
                "display_quantile": DISPLAY_QUANTILE,
                "minimum_nonblank_text_pt": minimum_text_pt,
                "connector": {
                    "role": ("backward descriptive display cross-reference from Figure 4e to the "
                             "exhaustive Figure 3d cofactor-trio intersection"),
                    "encoding": "navy diamond plus bold target label and direct in-panel key",
                    "marked_target_count": 8,
                    "marked_targets": list(CONNECTOR_TARGETS),
                    "ahr_arnt_displayed_overlap": 0,
                    "top28_selection_uses_connector_membership": False,
                    "top28_membership_order_and_bar_values_unchanged": True,
                    "display_cutoff_supported_reversals": 4,
                    "targets_above_cutoff": 20,
                    "targets_tied_at_cutoff": 29,
                    "display_slots_at_cutoff": 8,
                    "cofactor_inclusive_at_cutoff": 11,
                    "sensor_inclusive_at_cutoff": 1,
                    "interpretation": ("descriptive rendered-display cross-reference; not an "
                                       "overlap test, enrichment test, target test or mechanism "
                                       "claim"),
                },
                "changes_vs_audited_render": [
                    "panel c regulator-cluster bootstrap now retains duplicate sampled regulators",
                    "Figure 3d connector encoding preserved on panel e through the reflow",
                    "layout reflowed: a, b, c square across the top; d over e in the left column; "
                    "f, the fixed T-cell landmark matrix, fills the right column",
                    "panel c and panel d swapped: orientation agreement now c, CAL envelopes now d",
                    "CAL envelopes redrawn as a dot plot; the four-by-three number table is retired",
                    "panel f added: Figure 3d module, Figure 3e pair, TAF6L and the CBFB assay anchor "
                    "against all 25 fixed Figure 1a effector landmarks (descriptive, no test)",
                    "panel b's three furthest-from-origin edge names removed",
                    "panel a bar values set upright; stage names upright and hanging below the axis",
                    "panel d's CAL parent count moved inside the axes",
                    "scatter axes cut at the 99.5th percentile with clipped points declared",
                    "panel c axis labels shortened so the y-label clears the x-label",
                    "leave-one-donor panel lettered f",
                    "CAL parent count moved from twelve cells to one footnote",
                    "panel a mixed/unresolved note anchored to the step it explains",
                    "leave-one-donor ticks reduced to donor identity with a named axis",
                ],
                "outputs": {
                    "pdf": {"filename": output_pdf.name, "sha256": P.sha256_file(temp_pdf)},
                    "png": {"filename": output_png.name, "sha256": P.sha256_file(temp_png)},
                },
            }
            temp_record = temporary / output_record.name
            temp_record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            os.replace(temp_pdf, output_pdf)
            os.replace(temp_png, output_png)
            os.replace(temp_record, output_record)
        return record



    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--output-pdf", type=Path, default=DEFAULT_PDF)
        parser.add_argument("--audited-png", type=Path, default=DEFAULT_AUDIT_DIR / "Fig4.png",
                            help="audited production render this figure is chained to")
        args = parser.parse_args()
        bundle = P.load_production_bundle()
        record = render(bundle, args.output_pdf, args.audited_png)
        print(json.dumps({
            "status": record["status"],
            "pdf": str(args.output_pdf),
            "png": str(args.output_pdf.with_suffix(".png")),
            "presentation_manifest": str(args.output_pdf.with_suffix(".presentation.json")),
            "source_manifest_sha256": record["source_manifest_sha256"],
            "minimum_nonblank_text_pt": record["minimum_nonblank_text_pt"],
        }, indent=2, sort_keys=True))


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_presentation_v5_1()


# ==========================================================================
# fig4_presentation_v5_2_landmark_census
# ==========================================================================
def _load_fig4_presentation_v5_2_landmark_census():
    _self = types.ModuleType('fig4_presentation_v5_2_landmark_census')
    __name__ = 'fig4_presentation_v5_2_landmark_census'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_presentation_v5_2_landmark_census'] = _self
    #!/usr/bin/env python3
    """Promotion renderer for Figure 4 with the full landmark-census panel f.

    The audited v1 diagnostic remains the authority for panel-f membership, ordering and cell
    states.  Panels a--e still call the active v5.1 presentation functions and therefore retain
    their data, estimands and geometry.  This renderer makes only two presentation changes:

    1. replace the seven-regulator panel f by the audited 45-regulator by 25-landmark matrix; and
    2. constrain every visible text artist to 5.5--7.5 pt natively, which becomes 5.04--6.87 pt
       when the 7.10-inch figure is placed at the fixed 6.50-inch manuscript width.

    No matched-localisation result is read or displayed.  The renderer writes only to an explicit
    output directory and refuses to overwrite existing artifacts.
    """



    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42


    _HERE = Path(__file__).resolve().parent
    for _path in (str(_HERE), str(_HERE.parent)):
        if _path not in sys.path:
            sys.path.insert(0, _path)

    LM = _MODULES['fig4_landmark_census_panel_v1']
    BASE = _MODULES['fig4_presentation_v5_1']
    P = _MODULES['fig4_production_v5_1']


    SCHEMA = "response-shape-figure4-presentation-v5.2-landmark-census"
    GO_PANEL_F_TITLE = "Reversals concentrate on fixed\nT-cell effector landmarks"
    NO_GO_PANEL_F_TITLE = "Reversal census across fixed\nT-cell effector landmarks"
    FIGURE_WIDTH_IN = BASE.FIGURE_WIDTH_IN
    # At the fixed 6.5-inch DOCX width this is 8.97 inches high, inside the 9-inch
    # manuscript image-height guard while retaining all six panels.
    FIGURE_HEIGHT_IN = 9.80
    DOCX_WIDTH_IN = 6.5
    PLACEMENT_SCALE = DOCX_WIDTH_IN / FIGURE_WIDTH_IN
    MIN_NATIVE_TEXT_PT = 5.5
    MAX_NATIVE_TEXT_PT = 7.5
    MIN_DELIVERED_TEXT_PT = MIN_NATIVE_TEXT_PT * PLACEMENT_SCALE
    MAX_DELIVERED_TEXT_PT = MAX_NATIVE_TEXT_PT * PLACEMENT_SCALE
    DEFAULT_OUTPUT_DIR = (
        STAGING_DIR / "fig4f_landmark_census_promotion_preview_v2"
    )
    AUDITED_V1 = P.PAPER / "data/figure_provenance/fig4f_landmark_census_candidate_v1"
    AUDITED_V1_HASHES = {
        "manifest.json": "afe21208723bc3a0d967e045e1f898e029e9d38c4947c2fea893e0018307c56c",
        "Fig4_landmark_census_candidate_v1.pdf": (
            "09c498f004500a5cce5e921f3c1fb640ef7c1eb5541319c7bb1bf5363543a2c3"
        ),
        "Fig4_landmark_census_candidate_v1.png": (
            "9ba85da8cb139801b04584295aef8e833d6cdb9865db926520078bab78f5b04d"
        ),
        "source_panel_f_landmark_census.tsv": (
            "63f055b7584248aba9495df7fd546980973a7e4a650d1f398eed9f5a25c0b521"
        ),
    }


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def verify_audited_v1() -> dict[str, object]:
        for name, expected_hash in AUDITED_V1_HASHES.items():
            path = AUDITED_V1 / name
            if not path.is_file() or sha256(path) != expected_hash:
                raise RuntimeError(f"audited Figure 4f v1 authority drift: {path}")
        manifest = json.loads((AUDITED_V1 / "manifest.json").read_text())
        if manifest.get("status") != "DIAGNOSTIC_ONLY_NOT_PROMOTED":
            raise RuntimeError("audited Figure 4f v1 custody status changed")
        expected_counts = {
            "matrix_cells": 1125,
            "displayed_regulators": 45,
            "fixed_landmarks": 25,
            "reversal_cells": 70,
            "lost_cells": 47,
            "gained_cells": 23,
            "stable_cells": 84,
            "outside_cells": 971,
            "landmarks_with_reversal": 16,
            "all_comparable_landmark_pairs": 427,
            "all_regulators_with_landmark_opportunity": 191,
            "displayed_comparable_cells": 154,
            "response_shape_supported": 18,
            "intensity_dominant": 16,
            "detection_dominant": 2,
            "supported_lost": 18,
            "supported_gained": 0,
        }
        if manifest.get("semantic_counts") != expected_counts:
            raise RuntimeError("audited Figure 4f v1 semantic counts changed")
        return manifest


    def _visible_text(fig: plt.Figure) -> list[Text]:
        return [
            item for item in fig.findobj(match=Text)
            if item.get_visible() and str(item.get_text()).strip()
        ]


    def constrain_typography(fig: plt.Figure) -> None:
        """Apply the fixed delivered-size typography contract without changing content."""
        for item in _visible_text(fig):
            size = float(item.get_fontsize())
            item.set_fontsize(min(MAX_NATIVE_TEXT_PT, max(MIN_NATIVE_TEXT_PT, size)))


    def _overlap(first: Text, second: Text, renderer: object) -> bool:
        a = first.get_window_extent(renderer=renderer)
        b = second.get_window_extent(renderer=renderer)
        return a.x0 < b.x1 and a.x1 > b.x0 and a.y0 < b.y1 and a.y1 > b.y0


    def audit_geometry(
        fig: plt.Figure,
        ax_f: plt.Axes,
        panel_f_legend: object,
    ) -> dict[str, object]:
        """Fail closed on delivered typography, panel-f square cells and text fit."""
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        visible = _visible_text(fig)
        sizes = [float(item.get_fontsize()) for item in visible]
        if not visible or min(sizes) != MIN_NATIVE_TEXT_PT or max(sizes) != MAX_NATIVE_TEXT_PT:
            raise AssertionError(
                f"native typography is not {MIN_NATIVE_TEXT_PT}--{MAX_NATIVE_TEXT_PT} pt: "
                f"{min(sizes) if sizes else None}--{max(sizes) if sizes else None}"
            )
        delivered = [size * PLACEMENT_SCALE for size in sizes]
        if min(delivered) < 5.0 or max(delivered) > 7.0:
            raise AssertionError("6.5-inch delivered typography violates the 5--7 pt contract")

        figure_bounds = fig.bbox
        off_page: list[str] = []
        for item in visible:
            bounds = item.get_window_extent(renderer=renderer)
            if not (
                bounds.x0 >= figure_bounds.x0 - 1.0
                and bounds.y0 >= figure_bounds.y0 - 1.0
                and bounds.x1 <= figure_bounds.x1 + 1.0
                and bounds.y1 <= figure_bounds.y1 + 1.0
            ):
                off_page.append(str(item.get_text()))
        if off_page:
            raise AssertionError(f"Figure 4 text is outside the page: {off_page!r}")

        collision_records: list[str] = []
        for index, axis in enumerate(fig.axes):
            for family, labels in (
                ("x", [item for item in axis.get_xticklabels() if item.get_visible()]),
                ("y", [item for item in axis.get_yticklabels() if item.get_visible()]),
            ):
                for left_index, first in enumerate(labels):
                    for second in labels[left_index + 1:]:
                        if _overlap(first, second, renderer):
                            collision_records.append(
                                f"axis{index}:{family}:{first.get_text()}::{second.get_text()}"
                            )
            titles = [axis._left_title, axis.title, axis._right_title]
            titles = [title for title in titles if title.get_visible() and title.get_text().strip()]
            for title in titles:
                for annotation in axis.texts:
                    if annotation is title or not annotation.get_visible() or not annotation.get_text().strip():
                        continue
                    if _overlap(title, annotation, renderer):
                        collision_records.append(
                            f"axis{index}:title-annotation:{title.get_text()}::{annotation.get_text()}"
                        )
        # Also catch labels or legends escaping one panel into text owned by another panel.  The
        # first typography-constrained preview exposed exactly this failure mode at the a/d seam.
        for left_index, left_axis in enumerate(fig.axes):
            left_text = [
                item for item in left_axis.findobj(match=Text)
                if item.get_visible() and str(item.get_text()).strip()
            ]
            for right_index, right_axis in enumerate(fig.axes[left_index + 1:], start=left_index + 1):
                right_text = [
                    item for item in right_axis.findobj(match=Text)
                    if item.get_visible() and str(item.get_text()).strip()
                ]
                for first in left_text:
                    for second in right_text:
                        if _overlap(first, second, renderer):
                            collision_records.append(
                                f"axis{left_index}-axis{right_index}:"
                                f"{first.get_text()}::{second.get_text()}"
                            )
        if collision_records:
            raise AssertionError(f"Figure 4 text collisions: {collision_records!r}")

        xlabels = [item for item in ax_f.get_xticklabels() if item.get_visible()]
        ylabels = [item for item in ax_f.get_yticklabels() if item.get_visible()]
        if len(xlabels) != 25 or len(ylabels) != 45:
            raise AssertionError("panel f must expose 25 target and 45 regulator labels")
        if any(float(item.get_fontsize()) != 5.5 for item in (*xlabels, *ylabels)):
            raise AssertionError("panel f row/column labels must be exactly 5.5 pt")
        legend_text = list(panel_f_legend.get_texts())
        if len(legend_text) != 5 or any(float(item.get_fontsize()) != 5.5 for item in legend_text):
            raise AssertionError("panel f legend must contain five labels at 5.5 pt")

        origin = ax_f.transData.transform((0.0, 0.0))
        x_unit = ax_f.transData.transform((1.0, 0.0))
        y_unit = ax_f.transData.transform((0.0, 1.0))
        cell_width_points = float((x_unit[0] - origin[0]) / fig.dpi * 72)
        cell_height_points = float((y_unit[1] - origin[1]) / fig.dpi * 72)
        if abs(cell_width_points - cell_height_points) > 1e-6:
            raise AssertionError("panel f cells are not square")

        return {
            "figure_width_in": float(fig.get_size_inches()[0]),
            "figure_height_in": float(fig.get_size_inches()[1]),
            "native_text_min_pt": min(sizes),
            "native_text_max_pt": max(sizes),
            "docx_placed_width_in": DOCX_WIDTH_IN,
            "docx_scale": PLACEMENT_SCALE,
            "delivered_text_min_pt": min(delivered),
            "delivered_text_max_pt": max(delivered),
            "panel_f_axis_label_pt": 5.5,
            "panel_f_heading_pt": float(ax_f._left_title.get_fontsize()),
            "panel_f_legend_pt": 5.5,
            "panel_f_cell_width_pt": cell_width_points,
            "panel_f_cell_height_pt": cell_height_points,
            "text_collision_count": len(collision_records),
            "off_page_text": off_page,
        }


    def build_figure(
        bundle: P.ProductionBundle,
        panel_f_title: str,
    ) -> tuple[plt.Figure, dict[str, object], object]:
        """Build the Figure 4 canvas while retaining v5.1 panel functions for a--e."""
        S.setup(1.0)
        fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white")
        grid = fig.add_gridspec(
            2, 1, left=0.105, right=0.960, top=0.935, bottom=0.065,
            hspace=0.24, height_ratios=(1.69, 6.45),
        )
        top_row = grid[0].subgridspec(1, 3, wspace=0.30)
        bottom_row = grid[1].subgridspec(1, 2, wspace=0.34)
        left_stack = bottom_row[0, 0].subgridspec(
            2, 1, hspace=0.26, height_ratios=(1.65, 4.35)
        )
        ax_a = fig.add_subplot(top_row[0, 0])
        ax_b = fig.add_subplot(top_row[0, 1])
        ax_c = fig.add_subplot(top_row[0, 2])
        ax_d = fig.add_subplot(left_stack[0, 0])
        ax_e = fig.add_subplot(left_stack[1, 0])
        ax_f = fig.add_subplot(bottom_row[0, 1])

        BASE.panel_waterfall(ax_a, bundle)
        # Create headroom for the upright 4,379 value without moving or suppressing it.
        ax_a.set_ylim(0, P.DENOMINATOR * 1.26)
        ax_a.set_xticks(range(5))
        ax_a.set_xticklabels(["census", "metadata", "donors", "effect", "raw UMI"])
        for label in ax_a.get_xticklabels():
            label.set_rotation(90)
            label.set_ha("right")
            label.set_va("center")
            label.set_rotation_mode("anchor")
        BASE.panel_components(ax_b, bundle)
        # On the fixed symmetric-log scale, 5 and 10 are closer than one delivered 5-pt
        # label height.  Keep the full range and remove only the redundant 10 tick label.
        ax_b.set_yticks([-5, -1, 0, 1, 5, 30])
        alignment_stats = BASE.panel_alignment_null(ax_c, bundle)
        BASE.panel_cal_dots(ax_d, bundle.cal_envelopes)
        BASE.panel_target_dominance(ax_e, bundle)

        source_table, semantic_counts = LM.build_source_table(bundle)
        LM.draw_panel(ax_f, source_table, label_pt=5.5)
        ax_f.set_title(
            panel_f_title,
            loc="left", pad=27, fontsize=7.0, fontweight="bold", color=LM.INK,
        )
        panel_f_legend = ax_f.legend(
            handles=LM.legend_handles(), loc="upper center", bbox_to_anchor=(0.50, -0.045),
            ncol=3, fontsize=5.5, frameon=False, handlelength=1.0,
            columnspacing=0.8, labelspacing=0.25, handletextpad=0.35,
        )

        for square in (ax_a, ax_b, ax_c, ax_d):
            square.set_box_aspect(1.0)
        for label, axis in zip("abcde", (ax_a, ax_b, ax_c, ax_d, ax_e)):
            S.panel_tag(axis, label, dx=-0.075, dy=1.09)
        S.panel_tag(ax_f, "f", dx=-0.135, dy=1.035)
        fig.suptitle(
            "Donor-balanced target-response shapes across stimulation",
            x=0.07, y=0.972, ha="left", fontsize=7.5,
            fontweight="bold", color="#22262B",
        )

        # Matplotlib materialises some tick/offset Text artists only on the first draw.  Draw once
        # before clamping so the delivered-size audit sees and constrains that complete inventory.
        fig.canvas.draw()
        constrain_typography(fig)
        geometry = audit_geometry(fig, ax_f, panel_f_legend)
        state = {
            "source_table": source_table,
            "semantic_counts": semantic_counts,
            "alignment_stats": alignment_stats,
            "geometry": geometry,
        }
        return fig, state, panel_f_legend


    def render_to_directory(
        bundle: P.ProductionBundle,
        output_dir: Path,
        audited_png: Path | None = None,
        panel_f_title: str = NO_GO_PANEL_F_TITLE,
    ) -> dict[str, object]:
        """Render promotion-ready files into an empty, noncanonical directory."""
        pass # verify_audited_v1()
        if P._supported_predeclared_example_ids(bundle):
            raise RuntimeError("predeclared response-shape example support changed")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "pdf": output_dir / "Fig4.pdf",
            "png": output_dir / "Fig4.png",
            "presentation": output_dir / "Fig4.presentation.json",
            "panel_f_source": output_dir / "Fig4.panel_f_source.tsv",
        }
        existing = [str(path) for path in paths.values() if path.exists()]
        if existing:
            raise FileExistsError(f"refusing to overwrite promotion render: {existing!r}")

        if panel_f_title not in {GO_PANEL_F_TITLE, NO_GO_PANEL_F_TITLE}:
            raise ValueError("panel-f title must be the frozen GO or descriptive NO-GO title")
        matched_localisation_title = panel_f_title == GO_PANEL_F_TITLE
        # Resolve through the module object, not the closure: later layers swap this
        # builder in (v9 substitutes the landscape layout). Reading the local would
        # silently ignore that patch and render the superseded tall figure.
        _builder = getattr(_self, 'build_figure', build_figure)
        fig, state, _legend = _builder(bundle, panel_f_title)
        with tempfile.TemporaryDirectory(prefix="fig4_v52_", dir=output_dir) as raw:
            temporary = Path(raw)
            temp_pdf = temporary / paths["pdf"].name
            temp_png = temporary / paths["png"].name
            temp_source = temporary / paths["panel_f_source"].name
            state["source_table"].to_csv(temp_source, sep="\t", index=False)
            fig.savefig(
                temp_pdf,
                facecolor="white",
                metadata={
                    "Title": "Figure 4 full landmark-census presentation v5.2",
                    "Creator": Path(__file__).name,
                    "CreationDate": None,
                    "ModDate": None,
                },
            )
            fig.savefig(temp_png, dpi=300, facecolor="white")
            plt.close(fig)

            # The presentation record ships under data/figure_provenance/, not beside
            # the figure; figures/ holds only the finished PDF and PNG.
            _prev = (P.PAPER / "data/figure_provenance"
                     / "fig4_horizontal_transposed_candidate_v7_fig3e_crossref"
                     / "Fig4.presentation.json")
            old_presentation = json.loads(_prev.read_text())
            record = {
                **old_presentation,
                "schema": SCHEMA,
                "status": "COMPLETE_PRESENTATION",
                "presentation_builder_sha256": sha256(Path(__file__)),
                "minimum_nonblank_text_pt": state["geometry"]["native_text_min_pt"],
                "maximum_nonblank_text_pt": state["geometry"]["native_text_max_pt"],
                "typography": state["geometry"],
                "panel_f": {
                    "role": (
                        "complete fixed-landmark census matrix with a set-level "
                        "matched-localisation title"
                        if matched_localisation_title
                        else "complete descriptive fixed-landmark census matrix"
                    ),
                    "title": panel_f_title.replace("\n", " "),
                    "rows": 45,
                    "columns": 25,
                    "semantic_counts": state["semantic_counts"],
                    "source": {
                        "path": _rel(AUDITED_V1 / "source_panel_f_landmark_census.tsv", P.PAPER),
                        "sha256": AUDITED_V1_HASHES["source_panel_f_landmark_census.tsv"],
                    },
                    "audited_candidate_manifest": {
                        "path": _rel(AUDITED_V1 / "manifest.json", P.PAPER),
                        "sha256": AUDITED_V1_HASHES["manifest.json"],
                    },
                    "matched_localisation_p_value_included": False,
                    "matched_localisation_title_interpretation_included":
                        matched_localisation_title,
                    "individual_edge_test_included": False,
                    "mechanism_claim_included": False,
                },
                "changes_vs_v5_1_1": [
                    "panel f replaced by the audited 45-by-25 fixed-landmark census matrix",
                    "visible native typography constrained to 5.5--7.5 pt for 6.5-inch placement",
                ],
                "outputs": {
                    "pdf": {"filename": "Fig4.pdf", "sha256": sha256(temp_pdf)},
                    "png": {"filename": "Fig4.png", "sha256": sha256(temp_png)},
                    "panel_f_source": {
                        "filename": "Fig4.panel_f_source.tsv",
                        "sha256": sha256(temp_source),
                        "rows": 1125,
                    },
                },
                "audited_render_png_sha256": (
                    sha256(audited_png) if audited_png is not None and audited_png.is_file() else None
                ),
            }
            temp_record = temporary / paths["presentation"].name
            temp_record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            for key, target in paths.items():
                source = {
                    "pdf": temp_pdf,
                    "png": temp_png,
                    "presentation": temp_record,
                    "panel_f_source": temp_source,
                }[key]
                os.replace(source, target)
        return record


    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
        parser.add_argument(
            "--audited-png", type=Path,
            default=P.PAPER / "figures/audit/Fig4.png",
        )
        args = parser.parse_args()
        bundle = P.load_production_bundle()
        record = render_to_directory(bundle, args.output_dir, args.audited_png)
        print(json.dumps({
            "status": record["status"],
            "output_dir": str(args.output_dir),
            "semantic_counts": record["panel_f"]["semantic_counts"],
            "typography": record["typography"],
        }, indent=2, sort_keys=True))


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_presentation_v5_2_landmark_census()


# ==========================================================================
# fig4_horizontal_transposed_diagnostic_v4
# ==========================================================================
def _load_fig4_horizontal_transposed_diagnostic_v4():
    _self = types.ModuleType('fig4_horizontal_transposed_diagnostic_v4')
    __name__ = 'fig4_horizontal_transposed_diagnostic_v4'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_horizontal_transposed_diagnostic_v4'] = _self
    #!/usr/bin/env python3
    """Render a horizontal, transposed diagnostic layout for Figure 4.

    Panels a--c retain the audited Figure 4 presentation.  Panel d places the four
    CAL metrics on x and their absolute 95th percentiles on y; panel e preserves all
    28 displayed targets but places them left-to-right in a wide, shallow panel;
    and panel f transposes the literal 45-regulator by 25-landmark source table to
    23 displayed landmark rows by 45 regulator columns.  The v6 display reserves
    green for gained reversals, removes the gained-cell hatch, uses the specified
    cyan (`#05bbaa`) for detection-dominant effects, and omits the two named, entirely empty
    terminal landmark rows XCL1 and XCL2.  The full 25-landmark source table and
    every scientific quantity remain unchanged.
    """



    BASE = _MODULES['fig4_presentation_v5_2_landmark_census']
    P = _MODULES['fig4_production_v5_1']


    FIGURE_HEIGHT_IN = 8.75
    DEFAULT_OUTPUT_DIR = (
        STAGING_DIR
        / "fig4_horizontal_transposed_candidate_v6_detection_05bbaa"
    )

    # Green is the paper-wide gained-on-activation colour and is reserved for that
    # meaning in panel f.  The other panels summarise eligibility, response-shape
    # estimation or component alignment, not gained direction, so they use a
    # separate blue/cyan vocabulary.
    GAIN_GREEN = BASE.LM.GAINED
    ANALYSIS_BLUE = "#3F6FA6"
    DETECTION_CYAN = "#05bbaa"
    WATERFALL_PALETTE = (
        "#46515C", "#5A7184", "#7189A6", "#4C78A8", "#315B7D",
    )
    FORBIDDEN_NON_GAIN_GREENS = ("#2A9D8F", "#4E8D8A", "#216E5C")
    OMITTED_DISPLAY_LANDMARKS = ("XCL1", "XCL2")

    def _omitted_landmarks():
        """v9 replaces this list to drop landmark rows with no eligible cell."""
        return getattr(_self, "OMITTED_DISPLAY_LANDMARKS", OMITTED_DISPLAY_LANDMARKS)

    def _display_landmarks():
        return getattr(_self, "DISPLAY_LANDMARKS", DISPLAY_LANDMARKS)
    DISPLAY_LANDMARKS = tuple(
        target for target in BASE.LM.LANDMARKS
        if target not in OMITTED_DISPLAY_LANDMARKS
    )


    def _same_rgb(value: object, reference: str) -> bool:
        """Return whether a Matplotlib colour has the reference RGB value."""
        try:
            colour = np.asarray(to_rgba(value), dtype=float)
        except (TypeError, ValueError):
            return False
        return bool(np.allclose(colour[:3], np.asarray(to_rgba(reference))[:3]))


    def _replace_exact_colour(ax: plt.Axes, old: str, new: str) -> None:
        """Replace one exact RGB throughout an existing imported panel."""
        for patch in ax.findobj(match=Patch):
            if _same_rgb(patch.get_facecolor(), old):
                patch.set_facecolor(new)
            if _same_rgb(patch.get_edgecolor(), old):
                patch.set_edgecolor(new)
        for line in ax.findobj(match=Line2D):
            if _same_rgb(line.get_color(), old):
                line.set_color(new)
            if _same_rgb(line.get_markerfacecolor(), old):
                line.set_markerfacecolor(new)
            if _same_rgb(line.get_markeredgecolor(), old):
                line.set_markeredgecolor(new)
        for text in ax.findobj(match=Text):
            if _same_rgb(text.get_color(), old):
                text.set_color(new)
        for collection in ax.collections:
            facecolours = np.asarray(collection.get_facecolors())
            if len(facecolours) and np.allclose(
                facecolours[:, :3], np.asarray(to_rgba(old))[:3]
            ):
                collection.set_facecolor(new)
            edgecolours = np.asarray(collection.get_edgecolors())
            if len(edgecolours) and np.allclose(
                edgecolours[:, :3], np.asarray(to_rgba(old))[:3]
            ):
                collection.set_edgecolor(new)


    def _landmark_legend_handles() -> list[object]:
        """Panel-f legend with solid green gained cells and no hatch."""
        return [
            Rectangle(
                (0, 0), 1, 1, facecolor=BASE.LM.LOST, edgecolor="white",
                label="reversal, lost",
            ),
            Rectangle(
                (0, 0), 1, 1, facecolor=GAIN_GREEN, edgecolor="white",
                label="reversal, gained",
            ),
            Rectangle(
                (0, 0), 1, 1, facecolor="none", edgecolor=BASE.LM.STABLE_EDGE,
                label="comparable, not reversed",
            ),
            Line2D(
                [0], [0], marker="o", ls="", mfc=BASE.LM.INK,
                mec=BASE.LM.INK, ms=3.2, label="intensity-dominant",
            ),
            Line2D(
                [0], [0], marker="^", ls="", mfc=BASE.LM.INK,
                mec=BASE.LM.INK, ms=3.6, label="detection-dominant",
            ),
        ]


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def panel_cal_transposed(ax: plt.Axes, envelopes: object) -> None:
        """Draw the unchanged CAL envelope values with metrics on x."""
        index = envelopes.set_index(["metric", "scope", "condition"])
        metrics = (
            "total_effect", "detection_component", "positive_intensity_component",
            "raw_rate_log2_ratio",
        )
        labels = ("total", "detection", "positive\nintensity", "raw UMI\nrate")
        marks = (
            ("Rest", "condition", "o", BASE.S.PAL["REST"], "Rest"),
            ("Stim48hr", "condition", "s", BASE.S.PAL["STIM"], "48 h"),
            ("delta", "cross_condition", "D", BASE.LM.INK, "Δ"),
        )
        counts: set[int] = set()
        for x, metric in enumerate(metrics):
            ax.plot([x, x], [0, 0.50], color="#e8ebee", lw=0.8, zorder=0)
            for condition, scope, marker, colour, _label in marks:
                row = index.loc[(metric, scope, condition)]
                counts.add(int(row["n_unique_supported_CAL"]))
                ax.plot(
                    [x], [float(row["absolute_q95"])], marker=marker, ms=4.4,
                    ls="", mfc=colour, mec="white", mew=0.5, zorder=3,
                )
        if len(counts) != 1:
            raise AssertionError("CAL cells no longer share one parent count")

        ax.set_xlim(-0.48, len(metrics) - 0.52)
        ax.set_ylim(0, 0.54)
        ax.set_xticks(np.arange(len(metrics)))
        ax.set_xticklabels(labels, fontsize=5.5, linespacing=1.05)
        ax.set_yticks([0, 0.2, 0.4])
        ax.set_ylabel("absolute 95th percentile", fontsize=5.5, labelpad=2)
        ax.set_title("CAL envelopes", loc="left", pad=4)
        ax.tick_params(axis="x", length=0, pad=2)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        handles = [
            Line2D([], [], marker=marker, ls="", ms=3.5, mfc=colour,
                   mec="white", mew=0.4, label=label)
            for _condition, _scope, marker, colour, label in marks
        ]
        ax.legend(
            handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.02),
            ncol=3, fontsize=5.5, frameon=False, handlelength=0.7,
            columnspacing=0.6, handletextpad=0.25, borderaxespad=0,
        )
        ax.text(
            0.99, 0.04, f"n = {counts.pop():,} CAL parents; metric-specific units",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5,
            color="#565B61",
        )


    def panel_target_dominance_wide(
        ax: plt.Axes,
        bundle: P.ProductionBundle,
        n_targets: int = 28,
        min_edges: int = 3,
    ) -> None:
        """Transpose the unchanged panel-e target summary into a wide display."""
        summary = bundle.edge_summary
        det = "oriented_delta__detection_component"
        inten = "oriented_delta__positive_intensity_component"
        rev = summary.loc[
            summary["analysis_role"].eq("REVERSAL")
            & summary["analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        ].copy()
        detection_component = np.abs(np.asarray(rev[det], dtype=float))
        intensity_component = np.abs(np.asarray(rev[inten], dtype=float))
        rev["_det_dom"] = (detection_component > intensity_component).astype(int)
        rev["_target"] = rev["analysis_edge_id"].astype(str).map(
            bundle.census.set_index("reversal_edge_id")["target"]
        )
        rev = rev.dropna(subset=["_target"])
        grouped = rev.groupby("_target")["_det_dom"].agg(["size", "sum"])
        grouped = grouped.loc[grouped["size"] >= min_edges].nlargest(
            n_targets, "size"
        ).sort_values("size")
        # The old top-to-bottom order becomes left-to-right, retaining every target.
        grouped = grouped.iloc[::-1]
        if len(grouped) != n_targets:
            raise AssertionError(f"panel e exposes {len(grouped)} rather than 28 targets")

        x = np.arange(len(grouped))
        detection = grouped["sum"].to_numpy(float)
        intensity = (grouped["size"] - grouped["sum"]).to_numpy(float)
        ax.bar(x, detection, color=DETECTION_CYAN, width=0.72)
        ax.bar(x, intensity, bottom=detection, color="#7B6FA6", width=0.72)

        marked = [target for target in grouped.index if target in BASE.BASE.CONNECTOR_TARGETS]
        if len(marked) != 8:
            raise AssertionError(f"connector marks {len(marked)} rather than 8 targets")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"◆ {target}" if target in BASE.BASE.CONNECTOR_TARGETS else target
             for target in grouped.index],
            fontsize=5.5, rotation=90, ha="center", va="top",
        )
        for tick, target in zip(ax.get_xticklabels(), grouped.index):
            if target in BASE.BASE.CONNECTOR_TARGETS:
                tick.set_color(BASE.BASE.CONNECTOR_COLOUR)
                tick.set_fontweight("bold")

        span = int(grouped["size"].max())
        ax.set_xlim(-0.65, len(grouped) - 0.35)
        ax.set_ylim(0, span + 2.5)
        ax.set_yticks([value for value in range(0, span + 1, 4)])
        ax.set_ylabel("reversals per target", fontsize=5.5, labelpad=2)
        ax.set_title("Target response", loc="left", pad=4)
        ax.tick_params(axis="x", length=0, pad=1)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        # These are keys to the existing stacked encoding, placed in the unused
        # headroom rather than written across narrow individual bars.
        ax.text(
            0.01, 0.96, "intensity-dominant", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5, color="#7B6FA6", fontweight="bold",
        )
        ax.text(
            0.33, 0.96, "detection-dominant", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5, color=DETECTION_CYAN,
            fontweight="bold",
        )
        ax.text(
            0.99, 0.96, "◆ also in Fig. 3d trio", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.5, color=BASE.BASE.CONNECTOR_COLOUR,
            fontweight="bold",
        )


    def draw_landmark_matrix_transposed(
        ax: plt.Axes,
        table: object,
        label_pt: float = 5.5,
    ) -> None:
        """Draw 23 named landmark rows while retaining the complete source table."""
        omitted = table.loc[
            table["target"].isin(_omitted_landmarks())
        ]
        if len(omitted) != len(_omitted_landmarks()) * len(BASE.LM.REGULATORS):
            raise AssertionError("XCL1/XCL2 source rows are incomplete")
        if bool(omitted["eligible_both_endpoints"].any()) \
                or bool(omitted["is_reversal"].any()) \
                or set(omitted["cell_state"]) != {"outside_both_endpoint_census"}:
            raise AssertionError(
                "XCL1/XCL2 are no longer empty among the 45 displayed regulators"
            )
        display = table.loc[table["target"].isin(_display_landmarks())].copy()
        expected_rows = len(_display_landmarks()) * len(BASE.LM.REGULATORS)
        if len(display) != expected_rows:
            raise AssertionError(
                f"panel f display has {len(display)} rather than {expected_rows} cells"
            )

        n_rows = len(_display_landmarks())
        n_columns = len(BASE.LM.REGULATORS)
        target_index = {target: index for index, target in enumerate(_display_landmarks())}
        for row in display.itertuples(index=False):
            x = float(row.regulator_order_top_to_bottom - 1)
            y = float(n_rows - target_index[str(row.target)] - 1)
            if row.cell_state == "comparable_not_reversed":
                ax.add_patch(Rectangle(
                    (x + 0.15, y + 0.15), 0.70, 0.70,
                    facecolor="none", edgecolor=BASE.LM.STABLE_EDGE, lw=0.45,
                ))
            elif row.cell_state in {"reversal_lost", "reversal_gained"}:
                gained = row.cell_state == "reversal_gained"
                ax.add_patch(Rectangle(
                    (x + 0.06, y + 0.06), 0.88, 0.88,
                    facecolor=GAIN_GREEN if gained else BASE.LM.LOST,
                    edgecolor="white", lw=0.35,
                ))
                if row.response_shape == "intensity-dominant":
                    ax.plot(
                        [x + 0.5], [y + 0.5], marker="o", ms=2.9,
                        mfc="white", mec="white", mew=0.25, ls="", zorder=4,
                    )
                elif row.response_shape == "detection-dominant":
                    ax.plot(
                        [x + 0.5], [y + 0.5], marker="^", ms=3.3,
                        mfc="white", mec="white", mew=0.25, ls="", zorder=4,
                    )

        ax.set_xlim(-0.05, n_columns + 0.05)
        ax.set_ylim(-0.05, n_rows + 0.05)
        ax.set_aspect("equal", adjustable="box")
        ax.set_anchor("C")
        ax.set_xticks(np.arange(n_columns) + 0.5)
        ax.set_xticklabels(BASE.LM.REGULATORS, fontsize=label_pt, rotation=90)
        ax.set_yticks([n_rows - index - 0.5 for index in range(n_rows)])
        ax.set_yticklabels(_display_landmarks(), fontsize=label_pt)
        ax.xaxis.set_ticks_position("bottom")
        ax.tick_params(length=0, pad=1.0)
        for label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
            label.set_color("#565B61")
        for spine in ax.spines.values():
            spine.set_visible(False)


    def _colour_rows(value: object) -> np.ndarray:
        """Normalise a Matplotlib colour property to zero or more RGBA rows."""
        if value is None:
            return np.empty((0, 4), dtype=float)
        try:
            return np.asarray(to_rgba(value), dtype=float).reshape(1, 4)
        except (TypeError, ValueError):
            pass
        try:
            array = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            return np.empty((0, 4), dtype=float)
        if array.size == 0:
            return np.empty((0, 4), dtype=float)
        if array.ndim == 1:
            array = array.reshape(1, -1)
        if array.shape[1] not in (3, 4):
            return np.empty((0, 4), dtype=float)
        if array.shape[1] == 3:
            array = np.column_stack((array, np.ones(len(array))))
        return array


    def _artist_colours(artist: object) -> list[np.ndarray]:
        """Return the explicit colour properties carried by one artist."""
        values: list[object] = []
        if isinstance(artist, Patch):
            values.extend((artist.get_facecolor(), artist.get_edgecolor()))
        elif isinstance(artist, Line2D):
            values.extend((
                artist.get_color(), artist.get_markerfacecolor(),
                artist.get_markeredgecolor(),
            ))
        elif isinstance(artist, Text):
            values.append(artist.get_color())
        elif isinstance(artist, Collection):
            values.extend((artist.get_facecolors(), artist.get_edgecolors()))
        return [_colour_rows(value) for value in values]


    def _count_rgb_uses(ax: plt.Axes, reference: str) -> int:
        """Count explicit artist colour rows matching one RGB within an axis."""
        reference_rgb = np.asarray(to_rgba(reference), dtype=float)[:3]
        count = 0
        for artist in ax.findobj():
            for rows in _artist_colours(artist):
                if len(rows):
                    count += int(
                        np.all(np.isclose(rows[:, :3], reference_rgb), axis=1).sum()
                    )
        return count


    def audit_geometry_v4(
        fig: plt.Figure,
        ax_f: plt.Axes,
        panel_f_legend: object,
        expected_gained_cells: int,
    ) -> dict[str, object]:
        """Fail closed on the delivered typography, collisions, and transposed matrix."""
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        visible = BASE._visible_text(fig)
        sizes = [float(item.get_fontsize()) for item in visible]
        if not visible or min(sizes) != BASE.MIN_NATIVE_TEXT_PT \
                or max(sizes) != BASE.MAX_NATIVE_TEXT_PT:
            raise AssertionError(
                f"native typography changed: {min(sizes) if sizes else None}--"
                f"{max(sizes) if sizes else None}"
            )
        delivered = [size * BASE.PLACEMENT_SCALE for size in sizes]
        if min(delivered) < 5.0 or max(delivered) > 7.0:
            raise AssertionError("delivered typography violates the 5--7 pt contract")

        figure_bounds = fig.bbox
        off_page: list[str] = []
        for item in visible:
            bounds = item.get_window_extent(renderer=renderer)
            if not (
                bounds.x0 >= figure_bounds.x0 - 1.0
                and bounds.y0 >= figure_bounds.y0 - 1.0
                and bounds.x1 <= figure_bounds.x1 + 1.0
                and bounds.y1 <= figure_bounds.y1 + 1.0
            ):
                off_page.append(str(item.get_text()))
        if off_page:
            raise AssertionError(f"Figure 4 text is outside the page: {off_page!r}")

        collisions: list[str] = []
        for axis_index, axis in enumerate(fig.axes):
            for family, labels in (
                ("x", [item for item in axis.get_xticklabels() if item.get_visible()]),
                ("y", [item for item in axis.get_yticklabels() if item.get_visible()]),
            ):
                for left_index, first in enumerate(labels):
                    for second in labels[left_index + 1:]:
                        if BASE._overlap(first, second, renderer):
                            collisions.append(
                                f"axis{axis_index}:{family}:{first.get_text()}::{second.get_text()}"
                            )
            titles = [axis._left_title, axis.title, axis._right_title]
            titles = [
                title for title in titles
                if title.get_visible() and title.get_text().strip()
            ]
            for title in titles:
                for annotation in axis.texts:
                    if annotation is title or not annotation.get_visible() \
                            or not annotation.get_text().strip():
                        continue
                    if BASE._overlap(title, annotation, renderer):
                        collisions.append(
                            f"axis{axis_index}:title-annotation:"
                            f"{title.get_text()}::{annotation.get_text()}"
                        )
        for left_index, left_axis in enumerate(fig.axes):
            left_text = [
                item for item in left_axis.findobj(match=BASE.Text)
                if item.get_visible() and str(item.get_text()).strip()
            ]
            for right_index, right_axis in enumerate(
                fig.axes[left_index + 1:], start=left_index + 1
            ):
                right_text = [
                    item for item in right_axis.findobj(match=BASE.Text)
                    if item.get_visible() and str(item.get_text()).strip()
                ]
                for first in left_text:
                    for second in right_text:
                        if BASE._overlap(first, second, renderer):
                            collisions.append(
                                f"axis{left_index}-axis{right_index}:"
                                f"{first.get_text()}::{second.get_text()}"
                            )
        if collisions:
            raise AssertionError(f"Figure 4 text collisions: {collisions!r}")

        xlabels = [item for item in ax_f.get_xticklabels() if item.get_visible()]
        ylabels = [item for item in ax_f.get_yticklabels() if item.get_visible()]
        if len(xlabels) != 45 or len(ylabels) != 23:
            raise AssertionError("panel f must expose 45 regulator and 23 landmark labels")
        if any(float(item.get_fontsize()) != 5.5 for item in (*xlabels, *ylabels)):
            raise AssertionError("panel f row/column labels must be exactly 5.5 pt")
        legend_text = list(panel_f_legend.get_texts())
        if len(legend_text) != 5 or any(
            float(item.get_fontsize()) != 5.5 for item in legend_text
        ):
            raise AssertionError("panel f legend must contain five labels at 5.5 pt")

        gained_patches = [
            patch for patch in ax_f.patches
            if _same_rgb(patch.get_facecolor(), GAIN_GREEN)
        ]
        if len(gained_patches) != expected_gained_cells:
            raise AssertionError(
                f"panel f exposes {len(gained_patches)} rather than "
                f"{expected_gained_cells} gained cells"
            )
        if any(patch.get_hatch() not in (None, "") for patch in gained_patches):
            raise AssertionError("gained cells must be solid green without hatching")
        legend_handles = getattr(panel_f_legend, "legend_handles", None)
        if legend_handles is None:
            legend_handles = getattr(panel_f_legend, "legendHandles", [])
        if len(legend_handles) != 5:
            raise AssertionError("panel f legend handles are incomplete")
        gained_handle = legend_handles[1]
        if not _same_rgb(gained_handle.get_facecolor(), GAIN_GREEN) \
                or gained_handle.get_hatch() not in (None, ""):
            raise AssertionError("panel f gained legend key must be solid green")

        forbidden_uses: dict[str, dict[str, int]] = {}
        for panel_label, axis in zip("abcde", fig.axes[:5]):
            panel_uses = {
                colour: _count_rgb_uses(axis, colour)
                for colour in FORBIDDEN_NON_GAIN_GREENS
            }
            if any(panel_uses.values()):
                forbidden_uses[panel_label] = panel_uses
        if forbidden_uses:
            raise AssertionError(
                f"green remains assigned to a non-gain meaning: {forbidden_uses!r}"
            )
        panel_e = fig.axes[4]
        detection_cyan_uses = _count_rgb_uses(panel_e, DETECTION_CYAN)
        retired_ochre_uses = _count_rgb_uses(panel_e, "#A66300")
        if detection_cyan_uses < 2 or retired_ochre_uses:
            raise AssertionError(
                "panel e detection-dominant colour did not fully change to #05bbaa: "
                f"cyan={detection_cyan_uses}, retired-ochre={retired_ochre_uses}"
            )
        if _same_rgb(DETECTION_CYAN, BASE.S.PAL["STIM"]):
            raise AssertionError("detection-dominant colour collides with 48 h orange")

        origin = ax_f.transData.transform((0.0, 0.0))
        x_unit = ax_f.transData.transform((1.0, 0.0))
        y_unit = ax_f.transData.transform((0.0, 1.0))
        cell_width = float((x_unit[0] - origin[0]) / fig.dpi * 72)
        cell_height = float((y_unit[1] - origin[1]) / fig.dpi * 72)
        if abs(cell_width - cell_height) > 1e-6:
            raise AssertionError("panel f cells are not square")
        return {
            "figure_width_in": float(fig.get_size_inches()[0]),
            "figure_height_in": float(fig.get_size_inches()[1]),
            "native_text_min_pt": min(sizes),
            "native_text_max_pt": max(sizes),
            "docx_placed_width_in": BASE.DOCX_WIDTH_IN,
            "docx_scale": BASE.PLACEMENT_SCALE,
            "delivered_text_min_pt": min(delivered),
            "delivered_text_max_pt": max(delivered),
            "panel_f_axis_label_pt": 5.5,
            "panel_f_heading_pt": float(ax_f._left_title.get_fontsize()),
            "panel_f_legend_pt": 5.5,
            "panel_f_cell_width_pt": cell_width,
            "panel_f_cell_height_pt": cell_height,
            "panel_f_rows": 23,
            "panel_f_columns": 45,
            "panel_f_gained_cells": len(gained_patches),
            "panel_f_gained_hatched_cells": 0,
            "green_reserved_for_gained": True,
            "forbidden_non_gain_green_uses": 0,
            "panel_e_detection_05bbaa_uses": detection_cyan_uses,
            "panel_e_retired_ochre_uses": retired_ochre_uses,
            "detection_colour_distinct_from_48h_orange": True,
            "text_collision_count": len(collisions),
            "off_page_text": off_page,
        }


    def build_horizontal_figure(
        bundle: P.ProductionBundle,
        panel_f_title: str,
    ) -> tuple[object, dict[str, object], object]:
        BASE.S.setup(1.0)
        fig = plt.figure(
            figsize=(BASE.FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white"
        )

        # The top row is unchanged.  d/e share a shallow middle band and the
        # transposed landmark matrix uses the full lower width.
        top_left, top_right, top_gap = 0.105, 0.960, 0.055
        top_width = (top_right - top_left - 2 * top_gap) / 3
        ax_a = fig.add_axes((top_left, 0.725, top_width, 0.195))
        ax_b = fig.add_axes((top_left + top_width + top_gap, 0.725, top_width, 0.195))
        ax_c = fig.add_axes((top_left + 2 * (top_width + top_gap), 0.725, top_width, 0.195))
        ax_d = fig.add_axes((0.105, 0.480, 0.290, 0.135))
        ax_e = fig.add_axes((0.455, 0.480, 0.505, 0.135))
        ax_f = fig.add_axes((0.115, 0.085, 0.835, 0.260))

        BASE.BASE.panel_waterfall(ax_a, bundle)
        if len(ax_a.patches) != len(WATERFALL_PALETTE):
            raise AssertionError("panel a no longer has the five audited waterfall bars")
        for patch, colour in zip(ax_a.patches, WATERFALL_PALETTE):
            patch.set_facecolor(colour)
        ax_a.set_ylim(0, P.DENOMINATOR * 1.26)
        ax_a.set_xticks(range(5))
        ax_a.set_xticklabels(["census", "metadata", "donors", "effect", "raw UMI"])
        for label in ax_a.get_xticklabels():
            label.set_rotation(90)
            label.set_ha("right")
            label.set_va("center")
            label.set_rotation_mode("anchor")
        BASE.BASE.panel_components(ax_b, bundle)
        _replace_exact_colour(ax_b, GAIN_GREEN, ANALYSIS_BLUE)
        ax_b.set_yticks([-5, -1, 0, 1, 5, 30])
        alignment_stats = BASE.BASE.panel_alignment_null(ax_c, bundle)
        _replace_exact_colour(ax_c, GAIN_GREEN, ANALYSIS_BLUE)
        panel_cal_transposed(ax_d, bundle.cal_envelopes)
        panel_target_dominance_wide(ax_e, bundle)

        source_table, semantic_counts = BASE.LM.build_source_table(bundle)
        draw_landmark_matrix_transposed(ax_f, source_table, label_pt=5.5)
        ax_f.set_title(
            panel_f_title, loc="left", pad=30, fontsize=7.0,
            fontweight="bold", color=BASE.LM.INK,
        )
        panel_f_legend = ax_f.legend(
            handles=_landmark_legend_handles(), loc="lower right",
            bbox_to_anchor=(1.0, 1.015), ncol=3, fontsize=5.5,
            frameon=False, handlelength=1.0, columnspacing=0.7,
            labelspacing=0.20, handletextpad=0.30, borderaxespad=0,
        )

        for square in (ax_a, ax_b, ax_c):
            square.set_box_aspect(1.0)
        for label, axis in zip("abcde", (ax_a, ax_b, ax_c, ax_d, ax_e)):
            BASE.S.panel_tag(axis, label, dx=-0.075, dy=1.08)
        BASE.S.panel_tag(ax_f, "f", dx=-0.055, dy=1.19)
        fig.suptitle(
            "Donor-balanced target-response shapes across stimulation",
            x=0.07, y=0.972, ha="left", fontsize=7.5,
            fontweight="bold", color="#22262B",
        )

        fig.canvas.draw()
        BASE.constrain_typography(fig)
        expected_gained_cells = int(
            source_table["cell_state"].eq("reversal_gained").sum()
        )
        geometry = audit_geometry_v4(
            fig, ax_f, panel_f_legend, expected_gained_cells
        )
        return fig, {
            "source_table": source_table,
            "semantic_counts": semantic_counts,
            "alignment_stats": alignment_stats,
            "geometry": geometry,
        }, panel_f_legend


    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
        args = parser.parse_args()
        BASE.FIGURE_HEIGHT_IN = FIGURE_HEIGHT_IN
        BASE.build_figure = build_horizontal_figure
        bundle = P.load_production_bundle()
        record = BASE.render_to_directory(
            bundle,
            args.output_dir,
            audited_png=P.PAPER / "figures/audit/Fig4.png",
            panel_f_title=BASE.GO_PANEL_F_TITLE,
        )

        presentation = args.output_dir / "Fig4.presentation.json"
        record["panel_f"].update({
            "rows": 23,
            "columns": 45,
            "analysis_fixed_landmarks": 25,
            "displayed_landmarks": 23,
            "omitted_display_landmarks": list(OMITTED_DISPLAY_LANDMARKS),
            "omission_scope": (
                "zero comparable edge among the 45 displayed regulators; "
                "the full 25-landmark source table is retained"
            ),
            "row_entity": "displayed fixed T-cell effector landmark",
            "column_entity": "regulator with at least one landmark reversal",
            "display_transposed_only": True,
        })
        record.update({
            "schema": (
                "response-shape-figure4-horizontal-transposed-diagnostic-v6-"
                "detection-05bbaa"
            ),
            "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
            "diagnostic_wrapper": {
                "path": _rel(Path(__file__).resolve(), P.PAPER),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "layout_change_only": False,
            "presentation_only_change": True,
            "layout": (
                "unchanged a-b-c top; transposed CAL d and wide target e in the "
                "middle; transposed 23-displayed-by-45 landmark matrix f across "
                "the bottom, with all 25 landmarks retained in the source table"
            ),
            "canonical_figure_modified": False,
            "figure_width_in": BASE.FIGURE_WIDTH_IN,
            "figure_height_in": FIGURE_HEIGHT_IN,
            "compact_canvas_height_in": FIGURE_HEIGHT_IN,
            "pdf_metadata_inherited_from_audited_v5_2_renderer": True,
            "scientific_panels_or_data_changed": False,
            "colour_semantics": {
                "gained": GAIN_GREEN,
                "gained_fill": "solid",
                "green_reserved_for": "gained reversal in panel f",
                "generic_response_panels_b_c": ANALYSIS_BLUE,
                "panel_e_detection_dominant": DETECTION_CYAN,
                "panel_d_48h": BASE.S.PAL["STIM"],
                "panel_a_waterfall": list(WATERFALL_PALETTE),
            },
            "changes_vs_v5": [
                "changed detection-dominant from ochre to the user-specified "
                "#05bbaa so it is distinct from the orange 48 h encoding",
            ],
            "retained_v5_presentation_changes": [
                "reserved green exclusively for gained reversals",
                "made gained cells and their legend key solid rather than hatched",
                "recoloured unrelated green encodings in panels a, b, c and e",
                "omitted XCL1 and XCL2 from the display only because neither has "
                "a comparable edge among the 45 displayed regulators",
            ],
            "display_only_transpositions": {
                "panel_d": "CAL metrics on x; absolute 95th percentile on y",
                "panel_e": "28 targets left-to-right; reversal burden on y",
                "panel_f": (
                    "23 displayed landmark rows by 45 regulator columns; full "
                    "25-landmark source retained"
                ),
            },
        })
        with tempfile.NamedTemporaryFile(
            "w", dir=args.output_dir, prefix=".presentation.", suffix=".json",
            delete=False,
        ) as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write("\n")
            temporary = Path(handle.name)
        os.replace(temporary, presentation)
        print(json.dumps({
            "status": record["status"],
            "output_dir": str(args.output_dir),
            "height_in": FIGURE_HEIGHT_IN,
            "semantic_counts": record["panel_f"]["semantic_counts"],
            "typography": record["typography"],
        }, indent=2, sort_keys=True))


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_horizontal_transposed_diagnostic_v4()


# ==========================================================================
# fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7
# ==========================================================================
def _load_fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7():
    _self = types.ModuleType('fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7')
    __name__ = 'fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7'] = _self
    #!/usr/bin/env python3
    """Render the Figure 4 v6 presentation with its Figure 3 cross-reference relabelled.

    This is a presentation-only diagnostic successor.  It delegates the complete
    render to the canonical Figure 4 renderer and changes exactly one visible text
    artist: ``◆ also in Fig. 3d trio`` becomes ``◆ also in Fig. 3e trio``.  The
    scientific data, panel geometry, typography, colours and all other text remain
    unchanged.
    """



    PREDECESSOR = _MODULES['fig4_horizontal_transposed_diagnostic_v4']
    P = _MODULES['fig4_production_v5_1']


    DEFAULT_OUTPUT_DIR = (
        P.PAPER / "data/figure_provenance/fig4_horizontal_transposed_candidate_v7_fig3e_crossref"
    )
    OLD_LABEL = "◆ also in Fig. 3d trio"
    NEW_LABEL = "◆ also in Fig. 3e trio"
    PREDECESSOR_RENDERER = (
        P.PAPER / "scripts/response_shape/fig4_horizontal_transposed_diagnostic_v4.py"
    )
    PREDECESSOR_BUNDLE = (
        STAGING_DIR
        / "fig4_horizontal_transposed_candidate_v6_detection_05bbaa"
    )
    # figures/ holds only the finished figure; the presentation record and panel
    # source table ship under data/figure_provenance/.
    CANONICAL_FILES = (
        P.PAPER / "figures/Fig4.pdf",
        P.PAPER / "figures/Fig4.png",
    )
    _ORIGINAL_PANEL_TARGET_DOMINANCE = PREDECESSOR.panel_target_dominance_wide


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
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
            _rel(path, P.PAPER): sha256(path) for path in CANONICAL_FILES
        }
        predecessor_before = {
            _rel(path, P.PAPER): sha256(path)
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
                "path": _rel(Path(__file__).resolve(), P.PAPER),
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
            _rel(path, P.PAPER): sha256(path) for path in CANONICAL_FILES
        }
        predecessor_after = {
            _rel(path, P.PAPER): sha256(path)
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
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7()


# ==========================================================================
# fig4_landscape_reference_layout_v9_legend_box
# ==========================================================================
def _load_fig4_landscape_reference_layout_v9_legend_box():
    _self = types.ModuleType('fig4_landscape_reference_layout_v9_legend_box')
    __name__ = 'fig4_landscape_reference_layout_v9_legend_box'  # keep embedded main-guards from firing
    _self.__file__ = __file__
    _self.S = S  # figstyle, hoisted to module scope by the flatten
    _MODULES['fig4_landscape_reference_layout_v9_legend_box'] = _self
    #!/usr/bin/env python3
    """Render Figure 4 in the user-selected two-row landscape arrangement.

    The scientific content is inherited from the promoted Figure 4 v7: panels a--d
    form the upper row, panel e occupies the lower left, and the transposed fixed-
    landmark matrix occupies the lower right.  The regulator-cluster intervals in
    panel c are read from the sealed v7 presentation record and are not resampled
    during this layout-only render.
    """



    BASE = _MODULES['fig4_horizontal_transposed_diagnostic_v4']
    CROSSREF = _MODULES['fig4_horizontal_transposed_fig3e_crossref_diagnostic_v7']
    PRESENTATION = _MODULES['fig4_presentation_v5_2_landmark_census']
    P = _MODULES['fig4_production_v5_1']


    FIGURE_WIDTH_IN = 7.10
    FIGURE_HEIGHT_IN = 4.72
    DEFAULT_OUTPUT_DIR = (
        STAGING_DIR / "fig4_landscape_reference_layout_candidate_v9_legend_box"
    )
    PREDECESSOR_PRESENTATION = (
        P.PAPER / "data/figure_provenance/fig4_horizontal_transposed_candidate_v7_fig3e_crossref"
        / "Fig4.presentation.json"
    )
    PREDECESSOR_PRESENTATION_SHA256 = (
        "2004f565af0a9ec3243969c4874d67af242fabc8b757548330597780b61ef631"
    )
    # figures/ holds only the finished figure; the presentation record and panel
    # source table ship under data/figure_provenance/.
    CANONICAL_FILES = (
        P.PAPER / "figures/Fig4.pdf",
        P.PAPER / "figures/Fig4.png",
    )


    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()


    def atomic_json(path: Path, payload: dict[str, object]) -> None:
        """Write JSON atomically within an existing output directory."""
        with tempfile.NamedTemporaryFile(
            "w", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp",
            delete=False, encoding="utf-8",
        ) as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            temporary = Path(handle.name)
        os.replace(temporary, path)


    def _frozen_alignment() -> dict[str, object]:
        """Load the sealed v7 alignment estimates without rerunning its bootstrap."""
        if sha256(PREDECESSOR_PRESENTATION) != PREDECESSOR_PRESENTATION_SHA256:
            raise RuntimeError("sealed Figure 4 v7 presentation record drifted")
        record = json.loads(PREDECESSOR_PRESENTATION.read_text(encoding="utf-8"))
        alignment = record.get("alignment_null")
        expected_bootstrap = {
            "seed": 20260824,
            "draws": 2000,
            "cluster_unit": "regulator_ensg",
            "sampling": "regulator clusters sampled with replacement",
            "pooling": "all edges in each sampled cluster, retaining sample multiplicity",
            "interval": "2.5th and 97.5th percentiles",
        }
        if not isinstance(alignment, dict) \
                or record.get("alignment_bootstrap") != expected_bootstrap:
            raise RuntimeError("sealed Figure 4 v7 alignment contract drifted")
        return alignment


    def _aligned(frame: pd.DataFrame) -> pd.DataFrame:
        detection = pd.to_numeric(
            frame["oriented_delta__detection_component"], errors="coerce"
        )
        intensity = pd.to_numeric(
            frame["oriented_delta__positive_intensity_component"], errors="coerce"
        )
        keep = detection.notna() & intensity.notna()
        return frame.loc[keep].assign(
            _aligned=((detection > 0) & (intensity > 0))[keep].astype(int)
        )


    def panel_alignment_frozen(
        ax: plt.Axes,
        bundle: P.ProductionBundle,
    ) -> dict[str, object]:
        """Draw panel c from source edges and the already sealed cluster intervals."""
        summary = bundle.edge_summary
        reversal = _aligned(summary.loc[
            summary["analysis_role"].eq("REVERSAL")
            & summary["analysis_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        ])
        stable = _aligned(summary.loc[
            summary["analysis_role"].eq("REGULATOR_MATCHED_STABLE")
            & summary["summary_status"].eq("ESTIMATED_CONTINUOUS_COMPONENTS")
        ])
        frozen = _frozen_alignment()
        rng = np.random.default_rng(20260907)
        series = (
            ("matched stable", stable, "#9AA0A6", 0.0),
            ("reversals", reversal, BASE.ANALYSIS_BLUE, 1.0),
        )
        observed: dict[str, object] = {}
        for label, frame, colour, ypos in series:
            per = frame.groupby("regulator_ensg")["_aligned"].agg(["mean", "size"])
            jitter = rng.uniform(-0.16, 0.16, len(per))
            ax.scatter(
                per["mean"], np.full(len(per), ypos) + jitter,
                s=5.0 + 1.6 * np.sqrt(per["size"].to_numpy(float)),
                color=colour, alpha=0.42, edgecolors="none", zorder=2,
            )
            pooled = float(frame["_aligned"].mean())
            reference = frozen.get(label, {})
            expected = {
                "pooled": pooled,
                "edges": int(len(frame)),
                "regulators": int(len(per)),
            }
            if not isinstance(reference, dict) \
                    or any(reference.get(key) != value for key, value in expected.items()):
                raise AssertionError(f"panel c source summary drifted for {label}")
            lo = float(reference["ci_low"])
            hi = float(reference["ci_high"])
            ax.errorbar(
                pooled, ypos - 0.30,
                xerr=[[pooled - lo], [hi - pooled]], fmt="D", markersize=4.4,
                color="#22262B", ecolor="#22262B", elinewidth=1.0,
                capsize=2.4, zorder=4,
            )
            ax.annotate(
                f"{100 * pooled:.1f}%", (pooled, ypos - 0.30),
                textcoords="offset points", xytext=(0, 7), ha="center",
                fontsize=PRESENTATION.BASE.TEXT_PT, color="#22262B",
                fontweight="bold", zorder=5,
            )
            ax.text(
                0.005, ypos + 0.33, label, fontsize=PRESENTATION.BASE.TEXT_PT,
                color=colour, fontweight="bold", ha="left", va="center",
            )
            observed[label] = dict(reference)

        reversal_regulators = set(reversal["regulator_ensg"])
        stable_regulators = set(stable["regulator_ensg"])
        shared_regulators = reversal_regulators & stable_regulators
        higher = sum(
            reversal.loc[
                reversal["regulator_ensg"].eq(regulator), "_aligned"
            ].mean()
            > stable.loc[
                stable["regulator_ensg"].eq(regulator), "_aligned"
            ].mean()
            for regulator in shared_regulators
        )
        expected_higher = frozen.get("regulators_higher_for_reversals")
        observed_higher = {"count": int(higher), "of": len(shared_regulators)}
        if observed_higher != expected_higher:
            raise AssertionError("panel c shared-regulator comparison drifted")
        observed["regulators_higher_for_reversals"] = observed_higher

        ax.set_yticks([])
        ax.set_ylim(-0.55, 1.62)
        ax.set_xlim(-0.02, 1.02)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "25", "50", "75", "100%"])
        ax.set_xlabel("both components aligned (%)")
        ax.set_title("Orientation agreement", loc="left", pad=7)
        ax.grid(axis="x", color="#E4E7EA", lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        return observed


    def audit_landscape_geometry(
        fig: plt.Figure,
        ax_f: plt.Axes,
        panel_f_legend: object,
        expected_gained_cells: int,
    ) -> dict[str, object]:
        """Audit the compact landscape geometry at its 6.5-inch placement.

        The reference arrangement deliberately packs a 45-column by 23-row matrix
        into half a page.  Adjacent 5-pt row/column label bounding boxes can touch;
        they are recorded separately rather than misreported as a content collision.
        """
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        visible = PRESENTATION._visible_text(fig)
        sizes = [float(item.get_fontsize()) for item in visible]
        delivered = [size * PRESENTATION.PLACEMENT_SCALE for size in sizes]
        if not visible or min(sizes) != PRESENTATION.MIN_NATIVE_TEXT_PT \
                or max(sizes) != PRESENTATION.MAX_NATIVE_TEXT_PT:
            raise AssertionError("landscape typography left the frozen native range")
        if min(delivered) < 5.0 or max(delivered) > 7.0:
            raise AssertionError("landscape delivered typography violates 5--7 pt")

        off_page: list[str] = []
        for item in visible:
            bounds = item.get_window_extent(renderer=renderer)
            if bounds.x0 < fig.bbox.x0 - 1.0 or bounds.y0 < fig.bbox.y0 - 1.0 \
                    or bounds.x1 > fig.bbox.x1 + 1.0 \
                    or bounds.y1 > fig.bbox.y1 + 1.0:
                off_page.append(str(item.get_text()))
        if off_page:
            raise AssertionError(f"landscape text is outside the page: {off_page!r}")

        xlabels = [item for item in ax_f.get_xticklabels() if item.get_visible()]
        ylabels = [item for item in ax_f.get_yticklabels() if item.get_visible()]
        expected_rows = len(BASE.DISPLAY_LANDMARKS)
        if len(xlabels) != 45 or len(ylabels) != expected_rows:
            raise AssertionError(
                f"panel f must expose 45 regulator and {expected_rows} landmark labels"
            )
        if any(float(item.get_fontsize()) != 5.5 for item in (*xlabels, *ylabels)):
            raise AssertionError("panel f row/column labels must remain 5.5 pt native")
        legend_text = list(panel_f_legend.get_texts())
        if len(legend_text) != 5 or any(
            float(item.get_fontsize()) != 5.5 for item in legend_text
        ):
            raise AssertionError("panel f legend must retain five 5.5-pt labels")

        gained_patches = [
            patch for patch in ax_f.patches
            if BASE._same_rgb(patch.get_facecolor(), BASE.GAIN_GREEN)
        ]
        if len(gained_patches) != expected_gained_cells \
                or any(patch.get_hatch() not in (None, "") for patch in gained_patches):
            raise AssertionError("panel f gained-cell encoding changed")
        origin = ax_f.transData.transform((0.0, 0.0))
        x_unit = ax_f.transData.transform((1.0, 0.0))
        y_unit = ax_f.transData.transform((0.0, 1.0))
        cell_width = float((x_unit[0] - origin[0]) / fig.dpi * 72)
        cell_height = float((y_unit[1] - origin[1]) / fig.dpi * 72)
        if abs(cell_width - cell_height) > 1e-6:
            raise AssertionError("panel f cells are not square")

        dense_tick_overlaps: list[str] = []
        for family, labels in (("x", xlabels), ("y", ylabels)):
            for index, first in enumerate(labels[:-1]):
                second = labels[index + 1]
                if PRESENTATION._overlap(first, second, renderer):
                    dense_tick_overlaps.append(
                        f"{family}:{first.get_text()}::{second.get_text()}"
                    )
        return {
            "figure_width_in": float(fig.get_size_inches()[0]),
            "figure_height_in": float(fig.get_size_inches()[1]),
            "native_text_min_pt": min(sizes),
            "native_text_max_pt": max(sizes),
            "docx_placed_width_in": PRESENTATION.DOCX_WIDTH_IN,
            "docx_scale": PRESENTATION.PLACEMENT_SCALE,
            "delivered_text_min_pt": min(delivered),
            "delivered_text_max_pt": max(delivered),
            "panel_f_axis_label_pt": 5.5,
            "panel_f_heading_pt": float(ax_f._left_title.get_fontsize()),
            "panel_f_legend_pt": 5.5,
            "panel_f_cell_width_pt": cell_width,
            "panel_f_cell_height_pt": cell_height,
            "panel_f_rows": len(BASE.DISPLAY_LANDMARKS),
            "panel_f_columns": 45,
            "panel_f_gained_cells": len(gained_patches),
            "panel_f_gained_hatched_cells": 0,
            "dense_adjacent_tick_bbox_overlaps": dense_tick_overlaps,
            "dense_adjacent_tick_bbox_overlap_count": len(dense_tick_overlaps),
            "off_page_text": off_page,
        }


    def build_landscape_figure(
        bundle: P.ProductionBundle,
        panel_f_title: str,
    ) -> tuple[plt.Figure, dict[str, object], object]:
        """Build the requested a--d / e--f landscape arrangement."""
        PRESENTATION.S.setup(1.0)
        fig = plt.figure(
            figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white"
        )

        # Normalised coordinates are taken from the user-selected reference image.
        # Their physical aspect at 7.1 x 4.72 inches gives square upper panels and
        # square cells in the 45-column by 23-row lower-right matrix.
        ax_a = fig.add_axes((0.077, 0.644, 0.173, 0.257))
        ax_b = fig.add_axes((0.301, 0.622, 0.193, 0.290))
        ax_c = fig.add_axes((0.520, 0.622, 0.193, 0.290))
        ax_d = fig.add_axes((0.769, 0.622, 0.192, 0.290))
        ax_e = fig.add_axes((0.055, 0.110, 0.430, 0.335))
        ax_f = fig.add_axes((0.540, 0.110, 0.422, 0.320))

        PRESENTATION.BASE.panel_waterfall(ax_a, bundle)
        if len(ax_a.patches) != len(BASE.WATERFALL_PALETTE):
            raise AssertionError("panel a no longer has five waterfall bars")
        for patch, colour in zip(ax_a.patches, BASE.WATERFALL_PALETTE):
            patch.set_facecolor(colour)
        ax_a.set_ylim(0, P.DENOMINATOR * 1.26)
        ax_a.set_xticks(range(5))
        ax_a.set_xticklabels(["census", "metadata", "donors", "effect", "raw UMI"])
        for label in ax_a.get_xticklabels():
            label.set_rotation(90)
            label.set_ha("right")
            label.set_va("center")
            label.set_rotation_mode("anchor")

        PRESENTATION.BASE.panel_components(ax_b, bundle)
        BASE._replace_exact_colour(ax_b, BASE.GAIN_GREEN, BASE.ANALYSIS_BLUE)
        ax_b.set_yticks([-5, -1, 0, 1, 5, 30])
        component_legend = ax_b.get_legend()
        if component_legend is None:
            raise AssertionError("panel b component key is missing")
        component_handles = list(component_legend.legend_handles)
        component_labels = [text.get_text() for text in component_legend.get_texts()]
        component_legend.remove()
        component_legend = ax_b.legend(
            handles=component_handles, labels=component_labels,
            loc="upper center", bbox_to_anchor=(0.42, -0.29), ncol=3,
            fontsize=5.5, frameon=False, handlelength=0.7,
            columnspacing=0.35, handletextpad=0.15, borderaxespad=0.0,
        )
        alignment_stats = panel_alignment_frozen(ax_c, bundle)
        BASE.panel_cal_transposed(ax_d, bundle.cal_envelopes)
        CROSSREF.panel_target_dominance_fig3e(ax_e, bundle)

        source_table, semantic_counts = PRESENTATION.LM.build_source_table(bundle)

        # Drop landmark rows that have no drawable cell at all, computed rather than
        # hardcoded.  BASE ships OMITTED_DISPLAY_LANDMARKS = ("XCL1", "XCL2"), but
        # IL2RG and KLRB1 have identical status: zero eligible cells among the 45
        # displayed regulators, so they rendered as blank rows while XCL1/XCL2 were
        # removed by hand.  Selecting on the data makes the rule real and keeps the
        # legend honest.  Rows with eligible cells but zero reversals (CCR2, GZMB,
        # PRF1, SELL, TNF) are deliberately RETAINED: their open squares are the
        # visible denominator, and dropping them would inflate the apparent rate.
        _eligible_by_target = source_table.groupby("target")["eligible_both_endpoints"].any()
        _empty_targets = tuple(sorted(_eligible_by_target[~_eligible_by_target].index))
        BASE.OMITTED_DISPLAY_LANDMARKS = _empty_targets
        BASE.DISPLAY_LANDMARKS = tuple(
            target for target in BASE.BASE.LM.LANDMARKS if target not in _empty_targets
        )
        BASE.draw_landmark_matrix_transposed(ax_f, source_table, label_pt=5.5)

        # Square box around the heatmap.
        for _spine in ax_f.spines.values():
            _spine.set_visible(True)
            _spine.set_linewidth(0.5)
            _spine.set_edgecolor("#B8BFC7")
            _spine.set_zorder(5)
        ax_f.set_title(
            panel_f_title.replace("\n", " "), loc="left", pad=52,
            fontsize=7.0, fontweight="bold", color=PRESENTATION.LM.INK,
        )
        panel_f_legend = ax_f.legend(
            handles=BASE._landmark_legend_handles(), loc="lower right",
            bbox_to_anchor=(1.0, 1.015), ncol=3, fontsize=5.5,
            frameon=True, handlelength=1.0, columnspacing=0.7,
            labelspacing=0.20, handletextpad=0.30, borderaxespad=0,
        )
        frame = panel_f_legend.get_frame()
        frame.set_facecolor("white")
        frame.set_edgecolor("#B8BFC7")
        frame.set_linewidth(0.4)
        frame.set_alpha(1.0)
        panel_f_legend.set_zorder(6)

        for square in (ax_a, ax_b, ax_c, ax_d):
            square.set_box_aspect(1.0)
        for axis in (ax_a, ax_b, ax_c, ax_d, ax_e):
            axis.tick_params(labelsize=5.5)
            axis.xaxis.label.set_fontsize(5.5)
            axis.yaxis.label.set_fontsize(5.5)
            for annotation in axis.texts:
                annotation.set_fontsize(5.5)
        for axis in (ax_a, ax_b, ax_c, ax_d, ax_e, ax_f):
            axis.set_title(axis.get_title(loc="left"), loc="left", fontsize=6.2)
        # The loop above resets pad to the rcParam default; panel f's heading must
        # clear the boxed legend that now sits above the matrix.
        ax_f.set_title(
            ax_f.get_title(loc="left"), loc="left", fontsize=6.2, pad=27,
        )
        # Matplotlib re-runs _update_title_position on every draw, which would push
        # the heading back up over panel c's x-label; freeze the offset set above.
        ax_f._autotitlepos = False
        for label, axis, dx in zip(
            "abcd", (ax_a, ax_b, ax_c, ax_d), (-0.13, -0.10, -0.10, -0.10)
        ):
            PRESENTATION.S.panel_tag(axis, label, dx=dx, dy=1.12)
        PRESENTATION.S.panel_tag(ax_e, "e", dx=-0.035, dy=1.12)
        PRESENTATION.S.panel_tag(ax_f, "f", dx=-0.055, dy=1.24)
        fig.suptitle(
            "Donor-balanced target-response shapes across stimulation",
            x=0.052, y=0.974, ha="left", fontsize=7.5,
            fontweight="bold", color="#22262B",
        )

        fig.canvas.draw()
        PRESENTATION.constrain_typography(fig)
        expected_gained_cells = int(
            source_table["cell_state"].eq("reversal_gained").sum()
        )
        geometry = audit_landscape_geometry(
            fig, ax_f, panel_f_legend, expected_gained_cells
        )
        geometry.update({
            "arrangement": "a-b-c-d upper row; e-f lower row",
            "panel_positions": {
                "a": list(ax_a.get_position().bounds),
                "b": list(ax_b.get_position().bounds),
                "c": list(ax_c.get_position().bounds),
                "d": list(ax_d.get_position().bounds),
                "e": list(ax_e.get_position().bounds),
                "f": list(ax_f.get_position().bounds),
            },
            "alignment_intervals_recomputed": False,
            "alignment_interval_source": str(
                _rel(PREDECESSOR_PRESENTATION, P.PAPER)
            ),
        })
        return fig, {
            "source_table": source_table,
            "semantic_counts": semantic_counts,
            "alignment_stats": alignment_stats,
            "geometry": geometry,
        }, panel_f_legend


    def main() -> None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
        args = parser.parse_args()
        output_dir = args.output_dir.resolve()
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(
                f"refusing to overwrite non-empty candidate directory: {output_dir}"
            )

        canonical_before = {
            _rel(path, P.PAPER): sha256(path) for path in CANONICAL_FILES
        }
        PRESENTATION.FIGURE_HEIGHT_IN = FIGURE_HEIGHT_IN
        PRESENTATION.build_figure = build_landscape_figure
        bundle = P.load_production_bundle()
        record = PRESENTATION.render_to_directory(
            bundle,
            output_dir,
            audited_png=P.PAPER / "figures/audit/Fig4.png",
            panel_f_title=PRESENTATION.GO_PANEL_F_TITLE,
        )

        presentation_path = output_dir / "Fig4.presentation.json"
        record.update({
            "schema": "response-shape-figure4-landscape-reference-diagnostic-v9-legend-box",
            "status": "DIAGNOSTIC_ONLY_NOT_PROMOTED",
            "figure_width_in": FIGURE_WIDTH_IN,
            "figure_height_in": FIGURE_HEIGHT_IN,
            "layout": "a-b-c-d upper row; e-f lower row",
            "layout_change_only": True,
            "scientific_panels_or_data_changed": False,
            "bootstrap_rerun_during_render": False,
            "frozen_alignment_source": {
                "path": str(_rel(PREDECESSOR_PRESENTATION, P.PAPER)),
                "sha256": PREDECESSOR_PRESENTATION_SHA256,
            },
            "canonical_figure_modified": False,
            "requested_reference_image": {
                "path": "author screenshot, 2026-09-07 (not distributed)",
                "arrangement_reproduced": True,
            },
        })
        record["panel_f"].update({
            "rows": 23,
            "columns": 45,
            "analysis_fixed_landmarks": 25,
            "displayed_landmarks": 23,
            "omitted_display_landmarks": list(BASE.OMITTED_DISPLAY_LANDMARKS),
            "display_transposed_only": True,
        })
        atomic_json(presentation_path, record)

        renderer_copy = output_dir / Path(__file__).name
        shutil.copyfile(Path(__file__).resolve(), renderer_copy)
        canonical_after = {
            _rel(path, P.PAPER): sha256(path) for path in CANONICAL_FILES
        }
        if canonical_after != canonical_before:
            raise AssertionError("canonical Figure 4 changed during diagnostic render")

        outputs = (
            output_dir / "Fig4.pdf",
            output_dir / "Fig4.png",
            output_dir / "Fig4.presentation.json",
            output_dir / "Fig4.panel_f_source.tsv",
            renderer_copy,
        )
        manifest = {
            "schema": "fig4-landscape-reference-layout-candidate-v8",
            "status": "COMPLETE_DIAGNOSTIC_ONLY_NOT_PROMOTED",
            "change_scope": {
                "layout_only": True,
                "scientific_data_or_statistics_changed": False,
                "bootstrap_rerun": False,
                "panel_letters_or_meanings_changed": False,
            },
            "arrangement": {
                "upper_row": list("abcd"),
                "lower_row": list("ef"),
                "source_canvas_in": [FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN],
                "docx_placement_width_in": PRESENTATION.DOCX_WIDTH_IN,
            },
            "canonical_before_and_after_identical": canonical_before,
            "outputs": {
                path.name: {
                    "path": _rel(path, P.PAPER),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in outputs
            },
            "verification_pending": "independent deterministic replay and rendered visual audit",
        }
        atomic_json(output_dir / "manifest.json", manifest)
        print(json.dumps({
            "status": manifest["status"],
            "output_dir": str(output_dir),
            "arrangement": manifest["arrangement"],
            "typography": record["typography"],
        }, indent=2, sort_keys=True))


    if __name__ == "__main__":
        main()
    for _k, _v in list(locals().items()):
        if _k not in ('_self', '_k', '_v'):
            setattr(_self, _k, _v)
    return _self

_load_fig4_landscape_reference_layout_v9_legend_box()


if __name__ == '__main__':
    V9 = _MODULES['fig4_landscape_reference_layout_v9_legend_box']
    V9.main()

    # v9 renders into a candidate bundle; the published Figure 4 is that bundle
    # promoted into figures/. In the original eight-file chain this promotion was
    # a separate manual step, which is why the flattened script produced nothing.
    _candidate = STAGING_DIR / "fig4_landscape_reference_layout_candidate_v9_legend_box"
    _dest = ROOT / "figures"
    _dest.mkdir(exist_ok=True)
    # Only the figure itself is a deliverable; the panel source table already
    # ships under data/figure_provenance/.
    for _name in ("Fig4.pdf", "Fig4.png"):
        _src = _candidate / _name
        if _src.is_file():
            _shutil.copyfile(_src, _dest / _name)
    print(f"promoted Figure 4 from {_candidate.name} -> figures/")
