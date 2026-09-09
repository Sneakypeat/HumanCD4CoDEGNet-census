#!/usr/bin/env python3
"""Build the CX2 network files deposited in NDEx, straight from the census table.

Writes two networks into results/ndex/:

  codeg_cd4_trans_regulatory_census_full.cx2   9,581 nodes / 83,489 edges
  codeg_cd4_reversal_subnetwork.cx2            2,802 nodes /  4,379 edges

The full network is every trans regulator->target edge significant in both the
resting and 48 h stimulated states. The reversal subnetwork is the subset whose
sign differs between the two states.

A previous deposit carried 83,409 edges, 80 short of the census, which would have
contradicted the edge count reported in the manuscript. The assertions below fail
loudly rather than let a short network reach NDEx again: the counts are checked
against the census table, and every edge is required to survive into the CX2.

Upload is deliberately not part of this script; it needs NDEx credentials and is
run separately.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT / "data/census/census_master_edges.csv.gz"
OUTDIR = ROOT / "results/ndex"

# Authoritative totals, identical to the values reported in the manuscript.
EXPECT_FULL = {"nodes": 9581, "edges": 83489}
EXPECT_REV = {"nodes": 2802, "edges": 4379}

EDGE_ATTRS = {
    "lfc_Rest": "double",
    "lfc_Stim48": "double",
    "se_Rest": "double",
    "se_Stim48": "double",
    "adjp_Rest": "double",
    "adjp_Stim48": "double",
    "kind": "string",
    "direction": "string",
}


def build_cx2(edges: pd.DataFrame, name: str, description: str) -> list:
    symbols = sorted(set(edges.regulator) | set(edges.target))
    index = {symbol: i for i, symbol in enumerate(symbols)}

    nodes = [{"id": index[s], "v": {"name": s}} for s in symbols]
    rows = []
    for i, row in enumerate(edges.itertuples(index=False)):
        values = {}
        for attr in EDGE_ATTRS:
            value = getattr(row, attr)
            if pd.isna(value):
                continue
            values[attr] = float(value) if EDGE_ATTRS[attr] == "double" else str(value)
        rows.append({"id": i, "s": index[row.regulator], "t": index[row.target], "v": values})

    # every edge must survive into the CX2; this is the check the lost 80 needed
    assert len(rows) == len(edges), f"{name}: {len(edges) - len(rows)} edges dropped during export"

    return [
        {"CXVersion": "2.0", "hasFragments": False},
        {"metaData": [{"name": "nodes", "elementCount": len(nodes)},
                      {"name": "edges", "elementCount": len(rows)}]},
        {"attributeDeclarations": [{
            "nodes": {"name": {"d": "string", "a": "name"}},
            "edges": {k: {"d": v, "a": k} for k, v in EDGE_ATTRS.items()},
        }]},
        {"networkAttributes": [{"name": name, "description": description}]},
        {"nodes": nodes},
        {"edges": rows},
        {"status": [{"success": True}]},
    ]


def main() -> int:
    census = pd.read_csv(CENSUS)
    trans = census[~census.is_self]

    full = trans
    reversal = trans[trans.kind == "rev"]

    for frame, expect, label in ((full, EXPECT_FULL, "full census"),
                                 (reversal, EXPECT_REV, "reversal subnetwork")):
        n_nodes = len(set(frame.regulator) | set(frame.target))
        got = {"nodes": n_nodes, "edges": len(frame)}
        assert got == expect, f"{label} does not match the published totals: {got} != {expect}"

    OUTDIR.mkdir(parents=True, exist_ok=True)
    targets = (
        (full, "CoDEG CD4 T-cell trans-regulatory census (full)",
         "Every trans regulator-target edge significant in both resting and 48 h "
         "stimulated primary human CD4+ T cells, carrying the per-state log2 fold "
         "change, standard error and adjusted p-value.",
         "codeg_cd4_trans_regulatory_census_full.cx2"),
        (reversal, "CoDEG CD4 T-cell sign-reversal subnetwork",
         "The subset of the trans census whose regulatory effect reverses sign "
         "between the resting and 48 h stimulated states.",
         "codeg_cd4_reversal_subnetwork.cx2"),
    )
    for frame, name, description, filename in targets:
        cx2 = build_cx2(frame, name, description)
        out = OUTDIR / filename
        out.write_text(json.dumps(cx2))
        n_nodes = sum(a["elementCount"] for a in cx2[1]["metaData"] if a["name"] == "nodes")
        n_edges = sum(a["elementCount"] for a in cx2[1]["metaData"] if a["name"] == "edges")
        print(f"[ndex] {filename}: {n_nodes:,} nodes / {n_edges:,} edges")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
