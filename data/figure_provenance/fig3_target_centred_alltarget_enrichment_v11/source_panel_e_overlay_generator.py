#!/usr/bin/env python
"""Overlay the actual reversal pattern onto the CORRECTED, promoted excess-burden enrichment
result, rather than reporting NES/q numbers alone.

Version 2: sourced from the authoritative, independently audited, promoted analysis
(results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/production_v2/01_terms.tsv,
Supplementary Results R7, Supplementary Figure S5), NOT the exploratory v1 GSEA
(scripts/gsea_reversal_burden.py), which is superseded. v2 uses
a standardised score (exact conditional mean and variance under a regulator-by-Rest-sign null,
not the plain excess-burden of v1), two independent 10,000-draw permutation banks, and a global
maxT correction across all 2,593 pooled terms in addition to the within-library BH rule.

Scope, matching what R7 itself calls out as the headline sets: the 9 KEGG Legacy pathways meeting
the repeatable within-library-and-both-banks rule (statement_eligible), and the 10 GO Biological
Process terms that additionally survive the stricter global maxT correction in both permutation
banks separately (global_familywise_stable). GO's other 151 statement-eligible-but-not-globally-
stable terms are not shown here, matching S5's own choice to display only representative subsets.

R7 states explicitly that the enrichment result "does not identify whether their reversals were
lost or gained" -- that is exactly what this script adds, descriptively, computing no new
enrichment statistic.

Member gene lists are read directly from 01_terms.tsv's member_targets column, not re-derived from
the GMT files, so this cannot drift from the authoritative analysis's own membership.

Writes results_data/census/gsea_pathway_reversal_overlay.json
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TERMS = ROOT / "results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/production_v2/01_terms.tsv"
FINAL_AUDIT = (ROOT / "results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2"
              / "supplement_package_final_audit.json")
CENSUS = ROOT / "results_data/census/census_master_edges.csv.gz"
OUT = ROOT / "results_data/census/gsea_pathway_reversal_overlay.json"

SHARED_MIN_PATHWAYS = 4

# KEGG_PRION_DISEASES shares ubiquitin-proteasome genes with KEGG_PROTEASOME; declared here, not
# silently dropped, so the figure can label rather than hide the overlap R7 itself flags.
LIKELY_ANNOTATION_OVERLAP = {
    "Prion Diseases": "shares ubiquitin-proteasome genes with Proteasome (R7)",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_terms() -> list[dict]:
    with open(TERMS) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    kegg = [r for r in rows if r["library"] == "KEGG_Legacy_v2026.1"
            and r["statement_eligible"] == "True"]
    go = [r for r in rows if r["library"] == "GO_Biological_Process_2023"
          and r["global_familywise_stable"] == "True"]
    audit = json.loads(FINAL_AUDIT.read_text())["authoritative_results"]
    assert len(kegg) == audit["kegg_statement_eligible"] == 9
    assert len(go) == audit["go_bp_global_maxT_stable"] == 10
    return kegg + go


def main() -> None:
    rows = load_terms()

    m = pd.read_csv(CENSUS)
    m = m.loc[~m["is_self"].astype(bool)]
    universe = set(m["target"].unique())
    rev = m.loc[m["kind"].eq("rev")]

    orient = {}
    for t, g in rev.groupby("target"):
        dirs = set(g["direction"])
        orient[t] = next(iter(dirs)) if len(dirs) == 1 else "mixed"

    pathway_rows, membership_by_pathway = [], {}
    for r in rows:
        name, members = r["term_name"], set(r["member_targets"].split("|")) & universe
        assert len(members) == int(r["member_targets_n"]), name
        rc = rev.loc[rev["target"].isin(members)]
        rt = sorted(rc["target"].unique())
        lost = sum(1 for t in rt if orient[t] == "lost_on_activation")
        gained = sum(1 for t in rt if orient[t] == "gained_on_activation")
        mixed = sum(1 for t in rt if orient[t] == "mixed")
        pathway_rows.append({
            "term": name, "library": r["library"].replace("_v2026.1", "").replace("_2023", ""),
            "NES": round(float(r["observed_nes"]), 4), "q": round(float(r["bh_q_within_library"]), 4),
            "global_familywise_stable": r["global_familywise_stable"] == "True",
            "n_members_in_universe": len(members), "n_members_with_reversal": len(rt),
            "n_reversal_edges": int(len(rc)), "n_lost": lost, "n_gained": gained, "n_mixed": mixed,
            "annotation_overlap_flag": LIKELY_ANNOTATION_OVERLAP.get(name),
            "members_with_reversal": rt,
        })
        membership_by_pathway[name] = set(rt)

    from collections import Counter
    cnt = Counter()
    for name, members in membership_by_pathway.items():
        for t in members:
            cnt[t] += 1
    shared = sorted((t for t, c in cnt.items() if c >= SHARED_MIN_PATHWAYS), key=lambda t: -cnt[t])
    matrix = {t: {"n_pathways": cnt[t], "orientation": orient[t],
                  "pathways": sorted(n for n in membership_by_pathway if t in membership_by_pathway[n])}
              for t in shared}

    res = {"_what": ("reversal-pattern overlay on the CORRECTED, promoted excess-burden "
                     "enrichment result (v2)"),
           "_script": "scripts/gsea_pathway_reversal_overlay.py",
           "_source": "Supplementary Results R7 / Supplementary Figure S5, "
                      "results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2",
           "_supersedes": ("v1 used the exploratory scripts/gsea_reversal_burden.py result "
                           "(14 KEGG, weaker method); it is not used for this overlay"),
           "_note": ("descriptive only; no new enrichment statistic. A target's orientation is "
                     "its own dataset-wide direction, not restricted to pathway-member regulators. "
                     "R7 explicitly does not test lost-versus-gained direction; this does."),
           "inputs": {"01_terms.tsv": sha256(TERMS),
                      "supplement_package_final_audit.json": sha256(FINAL_AUDIT),
                      "census_master_edges.csv.gz": sha256(CENSUS)},
           "shared_min_pathways": SHARED_MIN_PATHWAYS,
           "n_kegg_statement_eligible": 9, "n_go_global_stable": 10,
           "pathways": pathway_rows,
           "cross_pathway_drivers": matrix}
    OUT.write_text(json.dumps(res, indent=2) + "\n")

    print(f"{len(pathway_rows)} pathways (9 KEGG statement-eligible + 10 GO globally stable)")
    for r in pathway_rows:
        flag = f"  [{r['annotation_overlap_flag']}]" if r["annotation_overlap_flag"] else ""
        print(f"  {r['library']:7s} {r['term'][:42]:42s} rev {r['n_members_with_reversal']:3d}/"
              f"{r['n_members_in_universe']:3d}  L{r['n_lost']:3d} G{r['n_gained']:3d} "
              f"M{r['n_mixed']:2d}  NES {r['NES']:+.2f}{flag}")
    print(f"\n{len(shared)} targets reversal-positive in >= {SHARED_MIN_PATHWAYS} pathways")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
