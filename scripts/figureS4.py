#!/usr/bin/env python3
"""Render Supplementary Figure S4, the broad all-target KEGG/GO rank enrichment.

The script is deliberately fail closed. It does not accept the exploratory v1 result and writes
nothing unless the v2 production checksum ledger is bound to a passing independent full replay.

Inputs. The run directory given by CODEG_ENRICHMENT_DIR, or by default
data/census/broad_alltarget_kegg_go_rank_enrichment_v2/production_v2, must hold the
complete v2 production output. Only 01_terms.tsv is redistributed here: the two
permutation null banks are 803 MB and are regenerated, with their checksum ledger, by
scripts/hpc/enrichment_permutation_banks.sbatch. Run that first, then point this script
at its output directory.

The figure ships as figures/FigS4.pdf; this renderer writes the same bytes under its own
candidate name, which the manuscript carries renumbered as S4.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import textwrap
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_ROOT = Path(
    os.environ.get(
        "CODEG_ENRICHMENT_DIR",
        ROOT / "data/census/broad_alltarget_kegg_go_rank_enrichment_v2",
    )
)
INPUT_DIR = ANALYSIS_ROOT / "production_v2"
AUDIT_PATH = ANALYSIS_ROOT / "independent_audit_v2.json"
OUTPUT_DIR = ROOT / "figures"

DECISION_SCHEMA = "broad-alltarget-kegg-go-rank-enrichment-v2-go-no-go"
AUDIT_SCHEMA = "broad-alltarget-kegg-go-rank-enrichment-v2-independent-audit"
ANALYSIS_SCHEMA = "broad-alltarget-kegg-go-rank-enrichment-v2"
ANALYSIS_STATUS = "POSTHOC_EXPLORATORY_METHOD_CORRECTION"
EXPECTED_DRAWS_PER_BANK = 10_000
EXPECTED_TERMS_PER_DRAW = 2_593
ALLOWED_DECISIONS = {
    "BROAD_ALLTARGET_RANK_ENRICHMENT_SUPPORTED",
    "POOLED_HIT_BANK_UNSTABLE_NO_STATEMENT",
    "NO_BROAD_ALLTARGET_RANK_ENRICHMENT_SUPPORTED",
}
KEGG_LIBRARY = "KEGG_Legacy_v2026.1"
GO_LIBRARY = "GO_Biological_Process_2023"
LIBRARIES = {KEGG_LIBRARY, GO_LIBRARY}

CONSUMED_FILES = (
    "00_run_manifest.json",
    "01_terms.tsv",
    "06_analysis_summary.json",
    "07_go_no_go.json",
    "08_leading_edges.tsv.gz",
)

TERM_COLUMNS = {
    "library",
    "term_id",
    "term_name",
    "member_targets",
    "observed_es",
    "observed_nes",
    "pooled_marginal_p",
    "bh_q_within_library",
    "leading_edge_targets",
    "leading_edge_targets_tie_expanded",
    "n_leading_edge",
    "n_leading_edge_tie_expanded",
    "statement_eligible",
    "global_familywise_supported",
    "global_familywise_stable",
}

LEADING_COLUMNS = {
    "library",
    "term_id",
    "term_name",
    "target",
    "target_ensg",
    "target_rank",
    "target_score",
    "observed_reversals",
    "expected_reversals_under_network_null",
    "term_member",
    "leading_edge",
    "leading_edge_tie_expanded",
    "leading_edge_order",
}

GO_OVERLAP_THRESHOLD = 0.50
MAX_GO_REPRESENTATIVES = 12
MAX_KEGG_PATHWAYS = 24
MAX_SHARED_TARGETS = 28

INK = "#25272a"
MUTED = "#666b70"
GRID = "#d8dce0"
HIGH = "#6f3c99"
MATRIX = "#6f3c99"
ABSENT = "#f4f5f6"
HEAD_PT = 7.0
BODY_PT = 6.0
SMALL_PT = 5.8
FIG_W = 7.20
FIG_H = 8.60
# The Supplement builder embeds figures at 100% of a conservatively estimated 6.25-inch
# text block.  Source text therefore remains at least 5 pt after that deterministic scaling.
SUPPLEMENT_TEXT_WIDTH_IN = 6.25
PDF_METADATA_DATE = datetime(2026, 9, 4, tzinfo=timezone.utc)
PDF_METADATA = {
    "Title": "Network-calibrated broad KEGG and GO rank enrichment",
    "Author": "CoDEG T-cell project",
    "Subject": "Post-hoc exploratory Supplementary diagnostic candidate",
    "Keywords": "CRISPR; reversal burden; KEGG; GO Biological Process; permutation",
    "Creator": "scripts/fig_supp_broad_gsea_v2.py",
    "CreationDate": PDF_METADATA_DATE,
    "ModDate": PDF_METADATA_DATE,
}
PNG_METADATA = {
    "Software": "scripts/fig_supp_broad_gsea_v2.py",
    "Creation Time": "2026-09-04T00:00:00+00:00",
}


class FigureInputError(RuntimeError):
    """Raised when an audited production input contract is not satisfied."""


@dataclass(frozen=True)
class AuditedBundle:
    input_dir: Path
    audit_path: Path
    checksums_path: Path
    terms: pd.DataFrame
    memberships: pd.DataFrame
    manifest: dict
    summary: dict
    decision: dict
    audit: dict
    input_hashes: dict[str, str]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FigureInputError(f"cannot read valid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise FigureInputError(f"expected a JSON object: {path}")
    return value


def parse_sha256s(path: Path) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise FigureInputError(f"missing checksum ledger: {path}") from exc
    parsed: dict[str, str] = {}
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64:
            raise FigureInputError(f"malformed SHA256SUMS line {line_number}")
        digest, name = parts[0].lower(), parts[1].strip()
        if name.startswith("*"):
            name = name[1:]
        candidate = Path(name)
        if candidate.is_absolute() or ".." in candidate.parts or name in parsed:
            raise FigureInputError(f"unsafe or duplicate SHA256SUMS entry: {name}")
        if any(ch not in "0123456789abcdef" for ch in digest):
            raise FigureInputError(f"invalid SHA-256 digest for {name}")
        parsed[name] = digest
    return parsed


def strict_bool(value: object, *, column: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == "true":
            return True
        if normalized == "false":
            return False
    raise FigureInputError(f"{column} must contain only literal True/False values")


def require_json_bool(payload: dict, key: str, expected: bool, *, label: str) -> None:
    value = payload.get(key)
    if type(value) is not bool or value is not expected:
        raise FigureInputError(f"{label}.{key} must be exactly {str(expected).lower()}")


def require_json_int(payload: dict, key: str, expected: int, *, label: str) -> None:
    value = payload.get(key)
    if type(value) is not int or value != expected:
        raise FigureInputError(f"{label}.{key} must be exactly {expected}")


def split_symbols(value: object) -> set[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    return {item for item in str(value).split("|") if item}


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise FigureInputError(f"{label} missing required columns: {', '.join(missing)}")


def load_audited_bundle(
    input_dir: Path = INPUT_DIR,
    audit_path: Path = AUDIT_PATH,
) -> AuditedBundle:
    """Load only a checksum-bound, independently replayed v2 production bundle."""
    input_dir = input_dir.resolve()
    audit_path = audit_path.resolve()
    checksums_path = input_dir / "SHA256SUMS.txt"
    decision_path = input_dir / "07_go_no_go.json"

    decision = read_json(decision_path)
    if decision.get("schema_version") != DECISION_SCHEMA:
        raise FigureInputError("07_go_no_go.json has the wrong schema_version")
    if decision.get("analysis_status") != ANALYSIS_STATUS:
        raise FigureInputError(f"analysis must remain {ANALYSIS_STATUS}")
    if decision.get("decision") not in ALLOWED_DECISIONS:
        raise FigureInputError("07_go_no_go.json has an unrecognised decision")
    require_json_bool(
        decision, "promotion_authorized", False, label="07_go_no_go.json"
    )

    audit = read_json(audit_path)
    if audit.get("schema_version") != AUDIT_SCHEMA:
        raise FigureInputError("independent_audit_v2.json has the wrong schema_version")
    if audit.get("audit_status") != "PASS":
        raise FigureInputError("independent audit has not passed")
    if audit.get("analysis_status") != ANALYSIS_STATUS:
        raise FigureInputError("independent audit has the wrong analysis_status")
    require_json_bool(audit, "all_checks_passed", True, label="independent audit")
    require_json_bool(
        audit, "full_permutation_replay", True, label="independent audit"
    )
    require_json_bool(
        audit, "production_module_imported", False, label="independent audit"
    )
    require_json_int(
        audit,
        "draws_replayed_per_bank",
        EXPECTED_DRAWS_PER_BANK,
        label="independent audit",
    )
    require_json_int(
        audit,
        "terms_replayed_per_draw",
        EXPECTED_TERMS_PER_DRAW,
        label="independent audit",
    )
    require_json_bool(
        audit,
        "scientific_promotion_authorized",
        False,
        label="independent audit",
    )
    require_json_bool(
        audit, "canonical_artifacts_modified", False, label="independent audit"
    )
    require_json_bool(
        audit,
        "cross_bank_calibration_evaluable",
        True,
        label="independent audit",
    )
    calibration = audit.get("cross_bank_null_calibration")
    if not isinstance(calibration, dict):
        raise FigureInputError("independent audit lacks cross_bank_null_calibration")
    require_json_bool(
        calibration,
        "all_reference_checks_passed",
        True,
        label="cross_bank_null_calibration",
    )
    require_json_bool(
        calibration,
        "same_bank_leave_one_out_used",
        False,
        label="cross_bank_null_calibration",
    )
    require_json_bool(
        calibration,
        "production_technical_calibration_gate",
        True,
        label="cross_bank_null_calibration",
    )
    require_json_bool(
        calibration,
        "biological_statement_gate",
        False,
        label="cross_bank_null_calibration",
    )
    if calibration.get("evaluation_status") != "EVALUABLE":
        raise FigureInputError(
            "cross_bank_null_calibration.evaluation_status must be EVALUABLE"
        )
    try:
        ledger_hash = sha256(checksums_path)
    except OSError as exc:
        raise FigureInputError(f"missing checksum ledger: {checksums_path}") from exc
    if audit.get("production_sha256s_sha256") != ledger_hash:
        raise FigureInputError("independent audit is not bound to this SHA256SUMS.txt")

    ledger = parse_sha256s(checksums_path)
    input_hashes: dict[str, str] = {}
    for name in CONSUMED_FILES:
        if name not in ledger:
            raise FigureInputError(f"consumed file is absent from SHA256SUMS.txt: {name}")
        path = input_dir / name
        try:
            observed = sha256(path)
        except OSError as exc:
            raise FigureInputError(f"missing consumed file: {name}") from exc
        if observed != ledger[name]:
            raise FigureInputError(f"checksum mismatch for {name}")
        input_hashes[name] = observed
    input_hashes["SHA256SUMS.txt"] = ledger_hash
    input_hashes["independent_audit_v2.json"] = sha256(audit_path)

    manifest = read_json(input_dir / "00_run_manifest.json")
    summary = read_json(input_dir / "06_analysis_summary.json")
    if manifest.get("schema_version") != ANALYSIS_SCHEMA:
        raise FigureInputError("00_run_manifest.json has the wrong schema_version")
    if manifest.get("analysis_status") != ANALYSIS_STATUS:
        raise FigureInputError("00_run_manifest.json has the wrong analysis_status")
    configuration = manifest.get("configuration")
    if not isinstance(configuration, dict):
        raise FigureInputError("00_run_manifest.json lacks configuration")
    require_json_int(
        configuration,
        "draws_per_bank",
        EXPECTED_DRAWS_PER_BANK,
        label="00_run_manifest.json.configuration",
    )
    if summary.get("analysis_status") != ANALYSIS_STATUS:
        raise FigureInputError("06_analysis_summary.json has the wrong analysis_status")
    try:
        terms = pd.read_csv(input_dir / "01_terms.tsv", sep="\t", dtype={"term_id": str})
        memberships = pd.read_csv(
            input_dir / "08_leading_edges.tsv.gz",
            sep="\t",
            dtype={"term_id": str, "target": str, "target_ensg": str},
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise FigureInputError("could not parse v2 result tables") from exc

    require_columns(terms, TERM_COLUMNS, "01_terms.tsv")
    require_columns(memberships, LEADING_COLUMNS, "08_leading_edges.tsv.gz")
    if set(terms["library"].dropna().unique()) != LIBRARIES:
        raise FigureInputError("01_terms.tsv must contain exactly the declared KEGG and GO libraries")
    if set(memberships["library"].dropna().unique()) != LIBRARIES:
        raise FigureInputError("08_leading_edges.tsv.gz has unexpected library values")
    if terms.duplicated(["library", "term_id"]).any():
        raise FigureInputError("duplicate library/term_id rows in 01_terms.tsv")
    if memberships.duplicated(["library", "term_id", "target"]).any():
        raise FigureInputError("duplicate library/term_id/target rows in 08_leading_edges.tsv.gz")

    for column in (
        "statement_eligible",
        "global_familywise_supported",
        "global_familywise_stable",
    ):
        terms[column] = terms[column].map(lambda value: strict_bool(value, column=column))
    if (terms["global_familywise_stable"] & ~terms["global_familywise_supported"]).any():
        raise FigureInputError("global_familywise_stable is true without pooled support")
    for column in ("term_member", "leading_edge", "leading_edge_tie_expanded"):
        memberships[column] = memberships[column].map(
            lambda value: strict_bool(value, column=column)
        )
    if not memberships["term_member"].all():
        raise FigureInputError("08_leading_edges.tsv.gz must contain only true term-member rows")
    if (memberships["leading_edge"] & ~memberships["term_member"]).any():
        raise FigureInputError("a leading-edge row is not a term member")
    if (memberships["leading_edge_tie_expanded"] & ~memberships["term_member"]).any():
        raise FigureInputError("a tie-expanded leading-edge row is not a term member")
    if (memberships["leading_edge"] & ~memberships["leading_edge_tie_expanded"]).any():
        raise FigureInputError("the tie-expanded leading edge omits an original leading-edge row")

    n_eligible = int(terms["statement_eligible"].sum())
    supported_decision = decision["decision"] == "BROAD_ALLTARGET_RANK_ENRICHMENT_SUPPORTED"
    if supported_decision != (n_eligible > 0):
        raise FigureInputError("decision and statement_eligible rows disagree")
    require_json_int(
        audit,
        "statement_eligible_terms",
        n_eligible,
        label="independent audit",
    )
    require_json_int(
        audit,
        "global_familywise_supported_terms",
        int(terms["global_familywise_supported"].sum()),
        label="independent audit",
    )

    numeric_term_columns = (
        "observed_es",
        "observed_nes",
        "pooled_marginal_p",
        "bh_q_within_library",
    )
    for column in numeric_term_columns:
        terms[column] = pd.to_numeric(terms[column], errors="coerce")
        if terms[column].isna().any() or ~np.isfinite(terms[column]).all():
            raise FigureInputError(f"non-finite {column} in 01_terms.tsv")
    for column in ("target_rank", "target_score", "leading_edge_order"):
        memberships[column] = pd.to_numeric(memberships[column], errors="coerce")
    if memberships["target_rank"].isna().any() or memberships["target_score"].isna().any():
        raise FigureInputError("target rank and score must be numeric for every membership")

    known_terms = set(zip(terms["library"], terms["term_id"]))
    membership_terms = set(zip(memberships["library"], memberships["term_id"]))
    if known_terms != membership_terms:
        raise FigureInputError("term IDs disagree between term and membership tables")
    grouped = memberships.groupby(["library", "term_id"], sort=False)
    for _, row in terms.iterrows():
        key = (row["library"], row["term_id"])
        group = grouped.get_group(key)
        stored_members = split_symbols(row["member_targets"])
        table_members = set(group["target"])
        if stored_members != table_members:
            raise FigureInputError(f"member_targets disagree for {key}")
        stored_leading = split_symbols(row["leading_edge_targets"])
        table_leading = set(group.loc[group["leading_edge"], "target"])
        if stored_leading != table_leading:
            raise FigureInputError(f"leading_edge_targets disagree for {key}")
        stored_expanded = split_symbols(row["leading_edge_targets_tie_expanded"])
        table_expanded = set(group.loc[group["leading_edge_tie_expanded"], "target"])
        if stored_expanded != table_expanded:
            raise FigureInputError(f"leading_edge_targets_tie_expanded disagree for {key}")
        if int(row["n_leading_edge"]) != len(table_leading):
            raise FigureInputError(f"n_leading_edge disagrees for {key}")
        if int(row["n_leading_edge_tie_expanded"]) != len(table_expanded):
            raise FigureInputError(f"n_leading_edge_tie_expanded disagrees for {key}")
        if row["statement_eligible"] and (row["observed_es"] <= 0 or not table_expanded):
            raise FigureInputError(f"eligible term lacks a positive ES/leading edge: {key}")

    target_consistency = memberships.groupby("target", sort=False).agg(
        n_rank=("target_rank", "nunique"), n_score=("target_score", "nunique")
    )
    if (target_consistency["n_rank"] != 1).any() or (target_consistency["n_score"] != 1).any():
        raise FigureInputError("a target has inconsistent rank or score across term rows")

    return AuditedBundle(
        input_dir=input_dir,
        audit_path=audit_path,
        checksums_path=checksums_path,
        terms=terms,
        memberships=memberships,
        manifest=manifest,
        summary=summary,
        decision=decision,
        audit=audit,
        input_hashes=input_hashes,
    )


def overlap_coefficient(left: set[str], right: set[str]) -> float:
    denominator = min(len(left), len(right))
    return 0.0 if denominator == 0 else len(left & right) / denominator


def select_go_representatives(
    eligible_go: pd.DataFrame,
    *,
    threshold: float = GO_OVERLAP_THRESHOLD,
    maximum: int = MAX_GO_REPRESENTATIVES,
) -> pd.DataFrame:
    """Greedily select reproducible GO representatives by actual member-set overlap."""
    ordered = eligible_go.sort_values(
        ["bh_q_within_library", "observed_nes", "term_name", "term_id"],
        ascending=[True, False, True, True],
        kind="mergesort",
    ).copy()
    representatives: list[tuple[str, set[str]]] = []
    records: list[dict[str, object]] = []
    for _, row in ordered.iterrows():
        members = split_symbols(row["member_targets"])
        best_id: str | None = None
        best_overlap = -1.0
        for representative_id, representative_members in representatives:
            overlap = overlap_coefficient(members, representative_members)
            if overlap > best_overlap:
                best_id, best_overlap = representative_id, overlap
        is_representative = best_id is None or best_overlap < threshold
        if is_representative:
            representative_id = str(row["term_id"])
            representatives.append((representative_id, members))
            assigned_overlap = 1.0
        else:
            representative_id = str(best_id)
            assigned_overlap = best_overlap
        record = row.to_dict()
        record.update(
            {
                "term_id": str(row["term_id"]),
                "representative_term_id": representative_id,
                "overlap_with_representative": assigned_overlap,
                "is_representative": is_representative,
            }
        )
        records.append(record)
    selected_ids = [term_id for term_id, _ in representatives[:maximum]]
    table = pd.DataFrame.from_records(records)
    if table.empty:
        return pd.DataFrame(
            columns=[
                "term_id",
                "term_name",
                "representative_term_id",
                "overlap_with_representative",
                "is_representative",
                "bh_q_within_library",
                "observed_nes",
                "representative_term_name",
                "cluster_size",
                "displayed_representative",
                "display_order",
                "display_label_raw",
                "display_label_wrapped",
            ]
        )
    sizes = table.groupby("representative_term_id")["term_id"].size().to_dict()
    representative_names = dict(zip(table["term_id"], table["term_name"]))
    table["representative_term_name"] = table["representative_term_id"].map(
        representative_names
    )
    table["cluster_size"] = table["representative_term_id"].map(sizes).astype(int)
    table["displayed_representative"] = (
        table["is_representative"] & table["term_id"].isin(selected_ids)
    )
    order_by_id = {term_id: index + 1 for index, term_id in enumerate(selected_ids)}
    table["display_order"] = pd.array(
        [order_by_id.get(term_id, pd.NA) for term_id in table["term_id"]],
        dtype="Int64",
    )
    table["display_label_raw"] = ""
    table["display_label_wrapped"] = ""
    return table


def ordered_kegg(terms: pd.DataFrame) -> pd.DataFrame:
    rows = terms.loc[(terms["library"] == KEGG_LIBRARY) & terms["statement_eligible"]].copy()
    rows = rows.sort_values(
        ["bh_q_within_library", "observed_nes", "term_id"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    if len(rows) > MAX_KEGG_PATHWAYS:
        raise FigureInputError(
            f"{len(rows)} eligible KEGG pathways exceed the readable {MAX_KEGG_PATHWAYS}-row design"
        )
    rows["pathway_number"] = np.arange(1, len(rows) + 1)
    return rows


def build_leading_edge_matrix(
    kegg: pd.DataFrame,
    memberships: pd.DataFrame,
    *,
    maximum_targets: int = MAX_SHARED_TARGETS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the complete target-by-term matrix and metadata marking displayed rows."""
    term_ids = list(kegg["term_id"])
    if not term_ids:
        return pd.DataFrame(), pd.DataFrame(
            columns=["target", "n_pathways", "target_rank", "target_score", "displayed", "display_order"]
        )
    subset = memberships.loc[
        (memberships["library"] == KEGG_LIBRARY)
        & memberships["term_id"].isin(term_ids)
        & memberships["leading_edge_tie_expanded"]
    ].copy()
    if subset.empty:
        return pd.DataFrame(columns=term_ids), pd.DataFrame(
            columns=["target", "n_pathways", "target_rank", "target_score", "displayed", "display_order"]
        )
    counts = subset.groupby("target")["term_id"].nunique()
    shared = counts.loc[counts >= 2]
    if shared.empty:
        return pd.DataFrame(columns=term_ids), pd.DataFrame(
            columns=["target", "n_pathways", "target_rank", "target_score", "displayed", "display_order"]
        )
    ranks = subset.groupby("target")["target_rank"].first()
    scores = subset.groupby("target")["target_score"].first()
    meta = pd.DataFrame(
        {
            "target": shared.index,
            "n_pathways": shared.values,
            "target_rank": [ranks[target] for target in shared.index],
            "target_score": [scores[target] for target in shared.index],
        }
    ).sort_values(
        ["n_pathways", "target_rank", "target"],
        ascending=[False, True, True],
        kind="mergesort",
    )
    meta = meta.reset_index(drop=True)
    meta["displayed"] = np.arange(len(meta)) < maximum_targets
    meta["display_order"] = pd.array(
        [index + 1 if index < maximum_targets else pd.NA for index in range(len(meta))],
        dtype="Int64",
    )
    matrix = pd.DataFrame(0, index=list(meta["target"]), columns=term_ids, dtype=int)
    for _, row in subset.iterrows():
        matrix.loc[row["target"], row["term_id"]] = 1
    return matrix, meta


def compact_label(value: object, maximum: int = 34) -> str:
    label = " ".join(str(value).replace("_", " ").split())
    return label if len(label) <= maximum else label[: maximum - 1].rstrip() + "…"


def unique_compact_labels(
    names: Iterable[object], term_ids: Iterable[object], *, maximum: int
) -> list[str]:
    """Compact labels deterministically and disambiguate any resulting collision."""
    labels = [compact_label(value, maximum) for value in names]
    counts = pd.Series(labels, dtype="object").value_counts().to_dict()
    output = []
    for label, term_id in zip(labels, term_ids):
        output.append(f"{label} ({term_id})" if counts[label] > 1 else label)
    if len(output) != len(set(output)):
        raise FigureInputError("term labels remain duplicated after stable-ID disambiguation")
    return output


def wrapped_full_labels(
    rows: pd.DataFrame,
    *,
    cluster_sizes: dict[str, int],
    width: int = 40,
) -> tuple[list[str], list[str]]:
    """Return collision-safe full GO labels and their deterministic wrapped forms."""
    names = [" ".join(str(value).replace("_", " ").split()) for value in rows["term_name"]]
    duplicate_names = pd.Series(names, dtype="object").value_counts().to_dict()
    raw: list[str] = []
    wrapped: list[str] = []
    for name, term_id, n_leading in zip(
        names, rows["term_id"], rows["n_leading_edge_tie_expanded"]
    ):
        identity = f" ({term_id})" if duplicate_names[name] > 1 else ""
        grouped = cluster_sizes[str(term_id)] - 1
        suffix = (
            f"(+{grouped} related; LE {int(n_leading)})"
            if grouped
            else f"[LE {int(n_leading)}]"
        )
        full_name = f"{name}{identity}"
        raw.append(f"{full_name} {suffix}")
        lines = textwrap.wrap(
            full_name,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
        if lines and len(lines[-1]) + 1 + len(suffix) <= width:
            lines[-1] = f"{lines[-1]} {suffix}"
        else:
            lines.append(suffix)
        wrapped.append("\n".join(lines))
    if len(wrapped) != len(set(wrapped)):
        raise FigureInputError("wrapped GO display labels are not unique")
    return raw, wrapped


def shared_scope_label(total: int, displayed: int) -> str:
    if total < 0 or displayed < 0 or displayed > total:
        raise FigureInputError("invalid shared-target display counts")
    return f"top {displayed} of {total}" if total > displayed else f"all {total}"


def style_axis(ax) -> None:
    ax.tick_params(axis="both", labelsize=SMALL_PT, length=0, pad=1.5, colors=INK)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["bottom"].set_linewidth(0.5)


def lollipop_panel(ax, rows: pd.DataFrame, labels: list[str], *, empty_message: str) -> None:
    if rows.empty:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        ax.text(0.5, 0.5, empty_message, ha="center", va="center", fontsize=BODY_PT, color=MUTED)
        return
    values = np.asarray(rows["observed_nes"], dtype=float)
    y = np.arange(len(rows))[::-1]
    if np.any(values <= 0):
        raise FigureInputError("the one-sided high-burden display received a non-positive NES")
    for yy, (_, row), value in zip(y, rows.iterrows(), values):
        ax.plot([0, value], [yy, yy], color=GRID, lw=0.8, solid_capstyle="round", zorder=1)
        ax.scatter([value], [yy], s=15, color=HIGH, edgecolor="white", linewidth=0.35, zorder=2)
        if row["global_familywise_stable"]:
            ax.scatter(
                [value], [yy], s=25, facecolor="none", edgecolor=INK, linewidth=0.55, zorder=3
            )
    ax.axvline(0, color=MUTED, lw=0.5, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=SMALL_PT, color=INK)
    for label in ax.get_yticklabels():
        label.set_linespacing(0.95)
        label.set_multialignment("right")
    ax.set_ylim(-0.7, len(rows) - 0.3)
    bound = max(1.0, float(np.max(np.abs(values))) * 1.16)
    ax.set_xlim(min(-0.08, float(np.min(values)) * 1.12), bound)
    ax.set_xlabel("network-null normalized ES", fontsize=BODY_PT, color=INK, labelpad=2)
    style_axis(ax)


def panel_leading_matrix(ax, matrix: pd.DataFrame, kegg: pd.DataFrame) -> None:
    if matrix.empty:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No target belongs to the leading edge of two eligible KEGG pathways.",
            ha="center",
            va="center",
            fontsize=BODY_PT,
            color=MUTED,
        )
        return
    # Pathways are rows and targets are columns so equal data units yield genuinely square cells
    # for the expected roughly 14-by-28 display instead of a stretched wide-row heatmap.
    values = matrix.transpose().to_numpy(dtype=int)
    ax.imshow(
        values,
        cmap=ListedColormap([ABSENT, MATRIX]),
        vmin=0,
        vmax=1,
        aspect="equal",
        origin="upper",
        interpolation="none",
    )
    ax.set_gid("panel_c_square_matrix")
    number_of = dict(zip(kegg["term_id"], kegg["pathway_number"]))
    ax.set_yticks(np.arange(len(matrix.columns)))
    ax.set_yticklabels(
        [str(number_of[term_id]) for term_id in matrix.columns], fontsize=SMALL_PT
    )
    ax.set_xticks(np.arange(len(matrix.index)))
    ax.set_xticklabels(
        list(matrix.index), fontsize=SMALL_PT, rotation=90, ha="center", va="top"
    )
    ax.tick_params(axis="both", length=0, pad=1.4)
    ax.set_ylabel("KEGG pathway number (panel a)", fontsize=BODY_PT, color=INK, labelpad=3)
    ax.set_xlabel("shared leading-edge target", fontsize=BODY_PT, color=INK, labelpad=3)
    ax.set_xticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.columns), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.55)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _visible_text(fig) -> list:
    items = list(fig.texts)
    for legend in fig.legends:
        items.extend(legend.get_texts())
    for ax in fig.axes:
        items.extend(ax.texts)
        xlo, xhi = sorted(ax.get_xlim())
        ylo, yhi = sorted(ax.get_ylim())
        items.extend(
            label
            for label in ax.get_xticklabels()
            if xlo - 1e-12 <= float(label.get_position()[0]) <= xhi + 1e-12
        )
        items.extend(
            label
            for label in ax.get_yticklabels()
            if ylo - 1e-12 <= float(label.get_position()[1]) <= yhi + 1e-12
        )
        items.extend((ax.xaxis.label, ax.yaxis.label, ax.title))
        legend = ax.get_legend()
        if legend is not None:
            items.extend(legend.get_texts())
    return [item for item in items if item.get_visible() and item.get_text()]


def _collision_pairs(items: list, renderer, *, minimum_area_px2: float = 0.5) -> list[dict[str, str]]:
    boxes = [
        (item.get_text(), item.get_window_extent(renderer=renderer))
        for item in items
        if item.get_visible() and item.get_text()
    ]
    collisions: list[dict[str, str]] = []
    for index, (left_text, left_box) in enumerate(boxes):
        for right_text, right_box in boxes[index + 1 :]:
            width = max(
                0.0,
                min(left_box.x1, right_box.x1) - max(left_box.x0, right_box.x0),
            )
            height = max(
                0.0,
                min(left_box.y1, right_box.y1) - max(left_box.y0, right_box.y0),
            )
            if width * height > minimum_area_px2:
                collisions.append({"left": left_text[:80], "right": right_text[:80]})
    return collisions


def audit_figure(fig) -> dict:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = _visible_text(fig)
    sizes = [float(item.get_fontsize()) for item in texts]
    if not sizes:
        raise AssertionError("figure contains no visible text")
    if min(sizes) < 5.0 or max(sizes) > 7.0:
        raise AssertionError(f"visible font sizes leave the 5--7 pt range: {min(sizes)}--{max(sizes)}")
    supplement_scale = min(1.0, SUPPLEMENT_TEXT_WIDTH_IN / FIG_W)
    delivered_sizes = [size * supplement_scale for size in sizes]
    if min(delivered_sizes) < 5.0 or max(delivered_sizes) > 7.0:
        raise AssertionError(
            "visible font sizes leave the 5--7 pt range after Supplement scaling: "
            f"{min(delivered_sizes):.3f}--{max(delivered_sizes):.3f}"
        )
    figure_box = fig.get_window_extent(renderer=renderer)
    off_page: list[str] = []
    for item in texts:
        box = item.get_window_extent(renderer=renderer)
        if box.width <= 0 or box.height <= 0:
            continue
        if (
            box.x0 < figure_box.x0 - 0.5
            or box.y0 < figure_box.y0 - 0.5
            or box.x1 > figure_box.x1 + 0.5
            or box.y1 > figure_box.y1 + 0.5
        ):
            off_page.append(item.get_text()[:60])
    if off_page:
        raise AssertionError(f"off-page text: {off_page}")

    tick_collisions: list[dict[str, str]] = []
    duplicate_tick_labels: list[dict[str, object]] = []
    for axis_index, ax in enumerate(fig.axes):
        for axis_name, labels in (
            ("x", list(ax.get_xticklabels())),
            ("y", list(ax.get_yticklabels())),
        ):
            visible = [item for item in labels if item.get_visible() and item.get_text()]
            collisions = _collision_pairs(visible, renderer)
            tick_collisions.extend(
                {
                    "axis": str(axis_index),
                    "dimension": axis_name,
                    **collision,
                }
                for collision in collisions
            )
            strings = [item.get_text() for item in visible]
            duplicates = sorted({value for value in strings if strings.count(value) > 1})
            if duplicates:
                duplicate_tick_labels.append(
                    {"axis": axis_index, "dimension": axis_name, "labels": duplicates}
                )
    if tick_collisions:
        raise AssertionError(f"tick-label collisions: {tick_collisions}")
    if duplicate_tick_labels:
        raise AssertionError(f"duplicate tick labels: {duplicate_tick_labels}")

    figure_text_collisions = _collision_pairs(list(fig.texts), renderer)
    if figure_text_collisions:
        raise AssertionError(f"figure-text collisions: {figure_text_collisions}")

    square_cell_audit: dict[str, float] | None = None
    matrix_axes = [ax for ax in fig.axes if ax.get_gid() == "panel_c_square_matrix"]
    if len(matrix_axes) > 1:
        raise AssertionError("multiple panel-c square matrices found")
    if matrix_axes:
        matrix_ax = matrix_axes[0]
        if len(matrix_ax.images) != 1:
            raise AssertionError("panel c must contain exactly one matrix image")
        n_rows, n_columns = matrix_ax.images[0].get_array().shape
        axes_box = matrix_ax.get_window_extent(renderer=renderer)
        cell_width = axes_box.width / n_columns
        cell_height = axes_box.height / n_rows
        ratio = cell_width / cell_height
        if not 0.995 <= ratio <= 1.005:
            raise AssertionError(f"panel-c cells are not square: width/height={ratio:.6f}")
        square_cell_audit = {
            "cell_width_px": cell_width,
            "cell_height_px": cell_height,
            "width_to_height_ratio": ratio,
        }
    return {
        "visible_text_items": len(texts),
        "min_font_pt": min(sizes),
        "max_font_pt": max(sizes),
        "supplement_text_width_in": SUPPLEMENT_TEXT_WIDTH_IN,
        "supplement_scale": supplement_scale,
        "min_font_at_supplement_width_pt": min(delivered_sizes),
        "max_font_at_supplement_width_pt": max(delivered_sizes),
        "off_page_text_items": 0,
        "tick_label_collisions": 0,
        "duplicate_tick_labels": 0,
        "figure_text_collisions": 0,
        "panel_c_square_cells": square_cell_audit,
        "figure_width_in": FIG_W,
        "figure_height_in": FIG_H,
    }


def build_figure(
    bundle: AuditedBundle,
) -> tuple[plt.Figure, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    kegg = ordered_kegg(bundle.terms)
    eligible_go = bundle.terms.loc[
        (bundle.terms["library"] == GO_LIBRARY) & bundle.terms["statement_eligible"]
    ].copy()
    go_assignment = select_go_representatives(eligible_go)
    displayed_assignments = go_assignment.loc[
        go_assignment["displayed_representative"]
    ].sort_values("display_order", kind="stable")
    go_display = displayed_assignments.copy()
    cluster_size = dict(
        zip(displayed_assignments["term_id"], displayed_assignments["cluster_size"])
    )
    if not go_display.empty:
        go_raw_labels, go_labels = wrapped_full_labels(
            go_display, cluster_sizes=cluster_size
        )
        label_by_id = dict(zip(go_display["term_id"], go_labels))
        raw_by_id = dict(zip(go_display["term_id"], go_raw_labels))
        go_assignment["display_label_raw"] = go_assignment["term_id"].map(
            raw_by_id
        ).fillna("")
        go_assignment["display_label_wrapped"] = go_assignment["term_id"].map(
            label_by_id
        ).fillna("")
        go_display = go_assignment.loc[
            go_assignment["displayed_representative"]
        ].sort_values("display_order", kind="stable")
    else:
        go_labels = []

    matrix, matrix_meta = build_leading_edge_matrix(kegg, bundle.memberships)
    displayed_targets = (
        list(matrix_meta.loc[matrix_meta["displayed"], "target"])
        if not matrix_meta.empty
        else []
    )
    displayed_matrix = matrix.loc[displayed_targets] if displayed_targets else matrix.iloc[0:0]

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")
    ax_a = fig.add_axes([0.245, 0.565, 0.240, 0.360])
    ax_b = fig.add_axes([0.745, 0.565, 0.225, 0.360])
    ax_c = fig.add_axes([0.105, 0.115, 0.855, 0.330])

    compact_kegg = unique_compact_labels(
        kegg["term_name"], kegg["term_id"], maximum=23
    )
    kegg_labels = [
        f"{number}. {name}  [LE {int(n_le)}]"
        for number, name, n_le in zip(
            kegg["pathway_number"], compact_kegg, kegg["n_leading_edge_tie_expanded"]
        )
    ]
    kegg["display_label"] = kegg_labels
    lollipop_panel(ax_a, kegg, kegg_labels, empty_message="No statement-eligible KEGG pathway")
    lollipop_panel(ax_b, go_display, go_labels, empty_message="No statement-eligible GO term")
    panel_leading_matrix(ax_c, displayed_matrix, kegg)

    fig.text(0.018, 0.975, "a", fontsize=HEAD_PT, fontweight="bold", color=INK, va="top")
    fig.text(
        0.045,
        0.975,
        "KEGG pathways",
        fontsize=HEAD_PT,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(0.515, 0.975, "b", fontsize=HEAD_PT, fontweight="bold", color=INK, va="top")
    fig.text(
        0.542,
        0.975,
        f"GO-BP representatives (up to {MAX_GO_REPRESENTATIVES})",
        fontsize=HEAD_PT,
        fontweight="bold",
        color=INK,
        va="top",
    )
    total_shared = len(matrix_meta)
    displayed_shared = len(displayed_matrix)
    shared_scope = shared_scope_label(total_shared, displayed_shared)
    fig.text(0.018, 0.510, "c", fontsize=HEAD_PT, fontweight="bold", color=INK, va="top")
    fig.text(
        0.045,
        0.510,
        f"Shared KEGG leading-edge targets ({shared_scope})",
        fontsize=HEAD_PT,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.045,
        0.487,
        "purple = stored leading edge, including peak-score ties; targets shown in ≥2 pathways",
        fontsize=SMALL_PT,
        color=MUTED,
        va="top",
    )
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=HIGH,
                   markeredgecolor="none", label="high excess-reversal-burden enrichment"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
                   markeredgecolor=INK, markeredgewidth=0.55,
                   label="also stable under global maxT"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.046),
        ncol=2,
        frameon=False,
        fontsize=SMALL_PT,
        handletextpad=0.3,
        columnspacing=1.2,
        borderpad=0,
    )
    fig.text(
        0.045,
        0.030,
        "network-null normalized ES = observed ES / mean permutation ES",
        fontsize=SMALL_PT,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.985,
        0.014,
        "Post-hoc exploratory; eligibility is defined by the audited analysis, not by this display.",
        fontsize=SMALL_PT,
        color=MUTED,
        ha="right",
        va="bottom",
    )
    return fig, kegg, go_assignment, matrix, matrix_meta


def records_sha256(paths: Iterable[Path]) -> dict[str, str]:
    return {path.name: sha256(path) for path in paths}


def render_figure_bytes(fig: plt.Figure) -> tuple[bytes, bytes]:
    """Serialize with fixed metadata so independent rerenders are byte-identical."""
    pdf_buffer = io.BytesIO()
    png_buffer = io.BytesIO()
    fig.savefig(
        pdf_buffer,
        format="pdf",
        facecolor="white",
        metadata=PDF_METADATA,
    )
    fig.savefig(
        png_buffer,
        format="png",
        dpi=300,
        facecolor="white",
        metadata=PNG_METADATA,
    )
    return pdf_buffer.getvalue(), png_buffer.getvalue()


def audit_saved_artifacts(pdf_path: Path, png_path: Path) -> dict[str, object]:
    """Validate saved bytes, page geometry, raster geometry and actual PDF font operators."""
    try:
        from PIL import Image
        from pypdf import PdfReader
        from pypdf.generic import ContentStream
    except ImportError as exc:
        raise FigureInputError("saved-artifact QA requires Pillow and pypdf") from exc

    if not pdf_path.read_bytes().startswith(b"%PDF-"):
        raise FigureInputError("saved PDF lacks a PDF header")
    reader = PdfReader(str(pdf_path))
    if len(reader.pages) != 1:
        raise FigureInputError("saved PDF must contain exactly one page")
    page = reader.pages[0]
    width_pt = float(page.mediabox.width)
    height_pt = float(page.mediabox.height)
    if abs(width_pt - FIG_W * 72.0) > 0.05 or abs(height_pt - FIG_H * 72.0) > 0.05:
        raise FigureInputError(
            f"saved PDF geometry drift: {width_pt:.3f} x {height_pt:.3f} pt"
        )
    content = ContentStream(page.get_contents(), reader)
    pdf_font_sizes = [
        float(operands[1])
        for operands, operator in content.operations
        if operator == b"Tf"
    ]
    if not pdf_font_sizes:
        raise FigureInputError("saved PDF contains no font-size operators")
    supplement_scale = min(1.0, SUPPLEMENT_TEXT_WIDTH_IN / FIG_W)
    delivered_font_sizes = [size * supplement_scale for size in pdf_font_sizes]
    if min(pdf_font_sizes) < 5.0 or max(pdf_font_sizes) > 7.0:
        raise FigureInputError(
            "saved PDF font operators leave the 5--7 pt source range: "
            f"{min(pdf_font_sizes):.3f}--{max(pdf_font_sizes):.3f}"
        )
    if min(delivered_font_sizes) < 5.0 or max(delivered_font_sizes) > 7.0:
        raise FigureInputError(
            "saved PDF font operators leave the 5--7 pt range after Supplement scaling: "
            f"{min(delivered_font_sizes):.3f}--{max(delivered_font_sizes):.3f}"
        )

    try:
        with Image.open(png_path) as image:
            image.load()
            png_format = image.format
            png_pixels = list(image.size)
    except OSError as exc:
        raise FigureInputError("saved PNG is invalid") from exc
    expected_pixels = [round(FIG_W * 300), round(FIG_H * 300)]
    if png_format != "PNG" or png_pixels != expected_pixels:
        raise FigureInputError(
            f"saved PNG geometry drift: format={png_format}, pixels={png_pixels}"
        )
    return {
        "status": "PASS",
        "pdf_pages": 1,
        "pdf_page_points": [width_pt, height_pt],
        "pdf_font_operator_count": len(pdf_font_sizes),
        "pdf_min_font_operator_pt": min(pdf_font_sizes),
        "pdf_max_font_operator_pt": max(pdf_font_sizes),
        "pdf_min_font_at_supplement_width_pt": min(delivered_font_sizes),
        "pdf_max_font_at_supplement_width_pt": max(delivered_font_sizes),
        "png_format": png_format,
        "png_pixels": png_pixels,
    }


def write_candidate(bundle: AuditedBundle, output_dir: Path = OUTPUT_DIR) -> None:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FigureInputError(f"refusing to overwrite existing candidate directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage_dir = output_dir.parent / f".{output_dir.name}.incomplete-{os.getpid()}"
    if stage_dir.exists():
        raise FigureInputError(f"refusing to overwrite incomplete candidate directory: {stage_dir}")

    fig, kegg, go_assignment, matrix, matrix_meta = build_figure(bundle)
    render_audit = audit_figure(fig)
    pdf_bytes, png_bytes = render_figure_bytes(fig)
    plt.close(fig)

    # An independent construction and serialization must reproduce every displayed/source byte.
    replay_fig, replay_kegg, replay_go, replay_matrix, replay_meta = build_figure(bundle)
    replay_audit = audit_figure(replay_fig)
    replay_pdf, replay_png = render_figure_bytes(replay_fig)
    plt.close(replay_fig)
    pd.testing.assert_frame_equal(kegg, replay_kegg, check_exact=True)
    pd.testing.assert_frame_equal(go_assignment, replay_go, check_exact=True)
    pd.testing.assert_frame_equal(matrix, replay_matrix, check_exact=True)
    pd.testing.assert_frame_equal(matrix_meta, replay_meta, check_exact=True)
    if render_audit != replay_audit or pdf_bytes != replay_pdf or png_bytes != replay_png:
        raise FigureInputError("independent deterministic figure rerender did not match")
    render_audit["deterministic_rerender"] = {
        "status": "PASS",
        "pdf_byte_identical": True,
        "png_byte_identical": True,
        "pdf_sha256": hashlib.sha256(pdf_bytes).hexdigest(),
        "png_sha256": hashlib.sha256(png_bytes).hexdigest(),
    }

    stage_dir.mkdir(exist_ok=False)
    pdf_path = stage_dir / "FigS_broad_alltarget_kegg_go_rank_enrichment_v2.pdf"
    png_path = stage_dir / "FigS_broad_alltarget_kegg_go_rank_enrichment_v2.png"
    kegg_path = stage_dir / "source_panel_a_kegg.tsv"
    go_path = stage_dir / "source_panel_b_go_representatives.tsv"
    matrix_path = stage_dir / "source_panel_c_leading_edge_matrix.tsv"
    meta_path = stage_dir / "source_panel_c_targets.tsv"

    pdf_path.write_bytes(pdf_bytes)
    png_path.write_bytes(png_bytes)
    kegg.to_csv(kegg_path, sep="\t", index=False, lineterminator="\n")
    go_assignment.to_csv(go_path, sep="\t", index=False, lineterminator="\n")
    matrix_frame = matrix.rename_axis("target").reset_index()
    matrix_source = matrix_meta.merge(
        matrix_frame, on="target", how="left", validate="one_to_one"
    )
    matrix_source.to_csv(matrix_path, sep="\t", index=False, lineterminator="\n")
    matrix_meta.to_csv(meta_path, sep="\t", index=False, lineterminator="\n")
    saved_artifact_audit = audit_saved_artifacts(pdf_path, png_path)

    output_paths = [pdf_path, png_path, kegg_path, go_path, matrix_path, meta_path]
    manifest = {
        "schema_version": "broad-alltarget-kegg-go-rank-enrichment-v2-figure",
        "status": "SUPPLEMENTARY_DIAGNOSTIC_CANDIDATE_NOT_PROMOTED",
        "analysis_status": ANALYSIS_STATUS,
        "script": "scripts/fig_supp_broad_gsea_v2.py",
        "design": "docs/SUPPLEMENTARY_BROAD_GSEA_V2_FIGURE_DESIGN_2026-09-04.md",
        "analysis_decision": bundle.decision["decision"],
        "input_hashes": bundle.input_hashes,
        "render_audit": render_audit,
        "display": {
            "eligible_kegg_pathways": int(len(kegg)),
            "eligible_go_terms": int((bundle.terms["library"].eq(GO_LIBRARY) & bundle.terms["statement_eligible"]).sum()),
            "displayed_go_representatives": int(go_assignment["displayed_representative"].sum()) if not go_assignment.empty else 0,
            "go_overlap_coefficient_threshold": GO_OVERLAP_THRESHOLD,
            "shared_leading_edge_targets_total": int(len(matrix_meta)),
            "shared_leading_edge_targets_displayed": int(matrix_meta["displayed"].sum()) if not matrix_meta.empty else 0,
            "shared_leading_edge_targets_truncated": bool(
                len(matrix_meta) > int(matrix_meta["displayed"].sum())
            ) if not matrix_meta.empty else False,
            "shared_leading_edge_min_pathways": 2,
            "panel_c_orientation": "KEGG pathways as rows; shared targets as columns; square cells",
        },
        "saved_artifact_audit": saved_artifact_audit,
        "output_hashes": records_sha256(output_paths),
    }
    (stage_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(stage_dir, output_dir)


def main() -> None:
    bundle = load_audited_bundle()
    write_candidate(bundle)
    print(f"wrote audited Supplementary diagnostic candidate to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
