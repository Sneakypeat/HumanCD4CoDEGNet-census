#!/usr/bin/env python
"""Emit the canonical, machine-readable Arce validation artefacts.

  paper/data/arce/arce_validation_edges.csv.gz    one row per covered edge, both assets
  paper/data/arce/arce_validation_results.json    headline stats

Everything downstream (Fig 7) is built from the CSV, never from numbers
retyped out of ARCE_VALIDATION.md.

Sign convention: our lfc and Arce's effects are both DESeq2 log2 fold changes
of the perturbation vs its control, so they are on the same scale and signs are
directly comparable. (This is NOT true of the pooled screen, whose LFC is
log2(IL2Ra-low/high) and therefore sign-inverted; the screen is not in this
table.)
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
import arce_stats as S
import os
from pathlib import Path
_HOME = os.path.expanduser("~")   # was hardcoded to a specific home directory

PAPER = str(Path(__file__).resolve().parents[2])
DER = f"{_HOME}/CoDEG_Tcell/data/arce/derived"

# ------------------------------------------------------------ donor-reproduced
h = pd.read_csv(f"{PAPER}/data/holdout/holdout_donors_instances.csv.gz")
hr = h[h.kind == "rev"].groupby(["regulator", "target"])["repro"].mean().reset_index()
repro_set = set(map(tuple, hr[hr.repro >= 0.5][["regulator", "target"]].values))

COLS = ["regulator", "target", "target_ensg", "our_lfc_Rest", "our_lfc_Stim48hr",
        "our_kind", "minmag", "arce_dataset", "arce_effect_rest", "arce_effect_stim",
        "arce_se_rest", "arce_se_stim", "arce_padj_rest", "arce_padj_stim",
        "arce_basemean_rest", "arce_flipped", "arce_dir_match", "donor_reproduced"]


def norm(path, dataset):
    d = pd.read_csv(path)
    o = pd.DataFrame({
        "regulator": d.regulator, "target": d.target, "target_ensg": d.target_ensg,
        "our_lfc_Rest": d.lfc_Rest, "our_lfc_Stim48hr": d.lfc_Stim48hr,
        "our_kind": d.kind, "minmag": d.minmag,
        "arce_dataset": dataset,
        "arce_effect_rest": d.a_rest, "arce_effect_stim": d.a_stim,
        "arce_se_rest": d.a_rest_se, "arce_se_stim": d.a_stim_se,
        "arce_padj_rest": d.a_rest_padj, "arce_padj_stim": d.a_stim_padj,
        "arce_basemean_rest": d.a_rest_bm,
        "arce_flipped": d.arce_flip.astype(int),
        "arce_dir_match": d.dir_match.astype(int),
    })
    o["donor_reproduced"] = [int(t in repro_set) for t in
                             zip(o.regulator, o.target)]
    # recompute the flip flag from the effects themselves, as a self-check that
    # the exported columns are internally consistent
    chk = (np.sign(o.arce_effect_rest) != np.sign(o.arce_effect_stim)).astype(int)
    assert (chk == o.arce_flipped).all(), f"{dataset}: arce_flipped inconsistent"
    return o[COLS]


b = norm(f"{DER}/assetB_edges.csv", "perturbseq")
c = norm(f"{DER}/assetC_bulk_edges.csv", "bulk")
edges = pd.concat([b, c], ignore_index=True)
edges.to_csv(f"{PAPER}/data/arce/arce_validation_edges.csv.gz", index=False)
print(f"wrote data/arce/arce_validation_edges.csv.gz: {len(edges)} rows")
for ds, g in edges.groupby("arce_dataset"):
    print(f"  {ds:10s} rev={int((g.our_kind=='rev').sum()):4d} "
          f"stable={int((g.our_kind=='stable').sum()):4d} "
          f"regulators={g.regulator.nunique()}")
print(f"  donor_reproduced rows: {int(edges.donor_reproduced.sum())} "
       f"({edges[edges.donor_reproduced==1].groupby('regulator').target.unique().to_dict()})")

# ---------------------------------------------------------------- headline stats
# Recomputed FROM THE EXPORTED TABLE, so the JSON and the CSV cannot drift.
E = pd.read_csv(f"{PAPER}/data/arce/arce_validation_edges.csv.gz")
E["kind"] = E.our_kind
E["arce_flip"] = E.arce_flipped.astype(bool)

R = {"_note": ("Cluster on REGULATOR, never on edges. All CIs are cluster "
               "bootstraps over regulators. No genome-wide rate is estimable "
               "from these data: Arce's 28-gene library was selected for "
               "state-specific effects.")}

for ds in ["perturbseq", "bulk"]:
    d = E[E.arce_dataset == ds]
    tabs = S.strata_tables(d)
    or_mh = S.mh_odds_ratio(tabs)
    chi2, p = S.cmh_test(tabs)
    rev, st = d[d.kind == "rev"], d[d.kind == "stable"]
    entry = {
        "n_rev": len(rev), "n_stable": len(st),
        "n_regulators": int(d.regulator.nunique()),
        "rev_flip_rate": float(rev.arce_flip.mean()),
        "stable_flip_rate": float(st.arce_flip.mean()),
        "or_mh": float(or_mh), "cmh_p_edge_level": float(p),
        "per_regulator": S.per_regulator_rates(d),
    }
    if d.regulator.nunique() >= 3:
        ors = S.cluster_bootstrap_or(d, n_boot=10000, seed=1)
        lo, hi = np.percentile(ors, [2.5, 97.5])
        entry["cluster_boot_ci"] = [float(lo), float(hi)]
        entry["cluster_boot_p_le1"] = float((ors <= 1).mean())
        entry["leave_one_out"] = S.leave_one_out(d)
        entry["loo_or_min"] = float(min(x["or_mh"] for x in entry["leave_one_out"]))
        entry["loo_or_max"] = float(max(x["or_mh"] for x in entry["leave_one_out"]))
        pr = pd.DataFrame(entry["per_regulator"])
        nb = int((pr.rev_flip_rate > pr.stable_flip_rate).sum())
        entry["sign_test_n_better"] = nb
        entry["sign_test_n_reg"] = len(pr)
        entry["sign_test_p"] = float(stats.binomtest(nb, len(pr), 0.5).pvalue)
        entry["wilcoxon_p"] = float(stats.wilcoxon(pr.rev_flip_rate,
                                                   pr.stable_flip_rate).pvalue)
    # magnitude-matched (1:1 on minmag, within regulator)
    m = S.match_within_regulator(d, on="minmag", k=1)
    tabs_m = S.strata_tables(m)
    mr, ms = m[m.kind == "rev"], m[m.kind == "stable"]
    entry["matched"] = {
        "n_rev": len(mr), "n_stable": len(ms),
        "or_mh": float(S.mh_odds_ratio(tabs_m)),
        "cmh_p_edge_level": float(S.cmh_test(tabs_m)[1]),
        "rev_flip_rate": float(mr.arce_flip.mean()),
        "stable_flip_rate": float(ms.arce_flip.mean()),
        "balance_mwu_p": float(stats.mannwhitneyu(mr.minmag, ms.minmag).pvalue),
    }
    if m.regulator.nunique() >= 3:
        ors_m = S.cluster_bootstrap_or(m, n_boot=10000, seed=2)
        entry["matched"]["cluster_boot_ci"] = [float(x) for x in
                                               np.percentile(ors_m, [2.5, 97.5])]
    # MED12
    m12 = d[d.regulator == "MED12"]
    r12, s12 = m12[m12.kind == "rev"], m12[m12.kind == "stable"]
    a, bb = int(r12.arce_flip.sum()), len(r12) - int(r12.arce_flip.sum())
    cc, dd = int(s12.arce_flip.sum()), len(s12) - int(s12.arce_flip.sum())
    orr, pf = stats.fisher_exact([[a, bb], [cc, dd]])
    entry["med12"] = {
        "n_rev": len(r12), "rev_flip": a, "rev_flip_rate": float(r12.arce_flip.mean()),
        "n_stable": len(s12), "stable_flip": cc,
        "stable_flip_rate": float(s12.arce_flip.mean()),
        "fisher_or": float(orr), "fisher_p": float(pf),
    }
    R[ds] = entry

cb = E[(E.regulator == "CBFB") & (E.target == "IL2RA")]
R["cbfb_il2ra"] = {r.arce_dataset: {
    "our_rest": float(r.our_lfc_Rest), "our_stim": float(r.our_lfc_Stim48hr),
    "arce_rest": float(r.arce_effect_rest), "arce_stim": float(r.arce_effect_stim),
    "arce_rest_padj": float(r.arce_padj_rest), "arce_stim_padj": float(r.arce_padj_stim),
    "flipped": int(r.arce_flipped), "dir_match": int(r.arce_dir_match),
    "donor_reproduced": int(r.donor_reproduced)} for _, r in cb.iterrows()}
# The pooled CRISPR-KO screen (orthogonal modality AND readout). Read straight
# from Arce's table -- never retyped. Sign is INVERTED relative to a knockdown
# effect: LFC = log2(IL2Ra-low bin / IL2Ra-high bin), so a POSITIVE screen LFC
# means the gene is a POSITIVE regulator (its loss lowers IL-2Ra), which
# corresponds to a NEGATIVE knockdown effect on the transcript.
S1 = pd.read_excel(f"{_HOME}/CoDEG_Tcell/data/arce/zenodo/"
                   "data_tables/S1_all_screens_gene_summary.xlsx")
row = S1[S1.id == "CBFB"].iloc[0]
gata3 = float(S1[S1.id == "GATA3"].iloc[0]["neg|lfc.Stimulated_Teff"])
R["cbfb_il2ra"]["screen_ko_surface_protein"] = {
    "lfc_rest_lowhigh": float(row["neg|lfc.Resting_Teff"]),
    "fdr_rest_neg_direction": float(row["neg|fdr.Resting_Teff"]),
    "lfc_stim_lowhigh": float(row["neg|lfc.Stimulated_Teff"]),
    "fdr_stim_pos_direction": float(row["pos|fdr.Stimulated_Teff"]),
    "sign_relative_to_knockdown_effect": -1,
    "interpretation": ("negative regulator at rest, positive when stimulated; "
                       "both FDR-significant, opposite directions; matches ours "
                       "after the low/high sign inversion"),
    "polarity_anchor": {
        "gene": "GATA3", "published_stim_teff_lowhigh_lfc": 2.47,
        "table_value": gata3,
        "source": "Arce et al. 2024 main text vs Zenodo S1"},
    "note": "72 h restimulation, not 48 h; source Zenodo S1_all_screens_gene_summary",
}
assert abs(gata3 - 2.47) < 0.01, "screen polarity anchor failed"
json.dump(R, open(f"{PAPER}/data/arce/arce_validation_results.json", "w"), indent=2)
print("wrote data/arce/arce_validation_results.json")

print("\n--- headline, recomputed from the exported CSV ---")
for ds in ["perturbseq", "bulk"]:
    e = R[ds]
    ci = e.get("cluster_boot_ci")
    cis = f"[{ci[0]:.2f}, {ci[1]:.2f}]" if ci else "n/a (<3 clusters)"
    print(f"  {ds:10s} rev {e['rev_flip_rate']:.1%} vs stable {e['stable_flip_rate']:.1%}  "
          f"OR={e['or_mh']:.2f}  boot95 {cis}  matched OR={e['matched']['or_mh']:.2f}")
