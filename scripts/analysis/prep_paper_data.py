#!/usr/bin/env python
"""
prep_paper_data.py  — canonical number/array lock for the sign-reversal paper.
Reads only guaranteed-local artifacts (no S3 layer streaming):
  - paper/data/census/sign_inverting_edges.csv     (4,379 inverting edges + lfc/se both states)
  - paper/data/census/census_master_edges.csv.gz   (authoritative trans totals)
Emits paper/data/census/figure_data.json — the single source of truth for all figures + text.

The former ``stage2_replication.json`` pooled-score stratification is deliberately not read:
it was not an edge-level donor/guide holdout. Current holdout evidence lives under
``data/holdout/`` and must not be copied into this descriptive figure cache.
"""
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = str(Path(__file__).resolve().parents[2])
edges = pd.read_csv(os.path.join(ROOT, "data/census/sign_inverting_edges.csv"))
feas   = json.load(open(os.path.join(ROOT, "data/census/signflip_feasibility.json")))
master = pd.read_csv(os.path.join(ROOT, "data/census/census_master_edges.csv.gz"),
                     usecols=["regulator", "target", "kind", "is_self"])

print("edges columns:", list(edges.columns))
print("sample rows:\n", edges.head(3).to_string())
print("n inverting edges:", len(edges))

lfcR = edges["lfc_Rest"].to_numpy(float)
lfcS = edges["lfc_Stim48"].to_numpy(float)
min_abs = np.minimum(np.abs(lfcR), np.abs(lfcS))
delta   = np.abs(lfcR - lfcS)

# --- census composition: the paper estimand is trans, so self edges are masked here ---
trans = master[~master.is_self]
n_self = int(master.is_self.sum())
n_both = len(trans)                                      # 83,489
n_same = int(trans.kind.eq("stable").sum())              # 79,110
n_inv = int(trans.kind.eq("rev").sum())                  # 4,379
n_reg = int(trans.loc[trans.kind.eq("rev"), "regulator"].nunique())  # 324
rate     = n_inv / n_both
assert (n_self, n_both, n_same, n_inv, n_reg) == (555, 83_489, 79_110, 4_379, 324)

pos_to_neg = int(np.sum((lfcR > 0) & (lfcS < 0)))   # up in Rest -> down in Stim48
neg_to_pos = int(np.sum((lfcR < 0) & (lfcS > 0)))   # down in Rest -> up in Stim48
assert pos_to_neg + neg_to_pos == n_inv, (pos_to_neg, neg_to_pos, n_inv)

# --- magnitude survival (threshold on the WEAKER of the two states) ---
mag = {f">{t}": int(np.sum(min_abs > t)) for t in (0.0, 0.25, 0.5, 1.0)}
print("magnitude survival (min|lfc| over both states):", mag)

# --- per-regulator inverted-target counts (ALL 324, recomputed from edges) ---
reg_counts = (edges.groupby("regulator").size()
              .sort_values(ascending=False).astype(int))
n_reg_edges = len(reg_counts)
print("n regulators in edge file:", n_reg_edges)

# functional-class annotation for the top inverting regulators (curated)
CLASS = {}
for g in ["TADA2B","SGF29","SUPT7L","TAF6L","SUPT20H","TADA1","TAF5L","ATXN7L3",
          "USP22","ENY2","CCDC101","KAT2A","KAT2B","TRRAP"]: CLASS[g]="SAGA/co-activator (HAT)"
for g in ["CCNC","MED12","MED24","MED13","MED13L","CDK8","CDK19","MED1","MED23"]: CLASS[g]="Mediator"
for g in ["CREBBP","EP300"]: CLASS[g]="CBP/p300 (HAT)"
for g in ["LEO1","PAF1","CTR9","CDC73","RTF1","WDR61"]: CLASS[g]="PAF1C (elongation)"
for g in ["KDM1A","RCOR1","RCOR3","HDAC1","HDAC2"]: CLASS[g]="CoREST/LSD1 (co-repressor)"
for g in ["NCOR1","NCOR2","TBL1X","GPS2"]: CLASS[g]="NCoR (co-repressor)"
for g in ["STAT5B","STAT3","STAT5A","BATF","BATF3","IRF4","REL","RELA","NFKB1","FOXP3"]: CLASS[g]="signalling / activation TF"
for g in ["AHR","ARNT","HIF1A","EPAS1"]: CLASS[g]="bHLH-PAS TF"
for g in ["NCKAP1L","DOCK2","ARHGAP30","ARHGDIB","WAS","RAC2"]: CLASS[g]="cytoskeletal / immune-synapse"
for g in ["ATAD5","ELOF1","SENP5","EIF4G2","DHX36","ATP2A2","STIP1","POLR2","GTF2"]: CLASS[g]="other (replication/translation/misc)"

def klass(g):
    if g in CLASS: return CLASS[g]
    for pre,lab in [("MED","Mediator"),("TAF","SAGA/co-activator (HAT)"),
                    ("SUPT","SAGA/co-activator (HAT)"),("TADA","SAGA/co-activator (HAT)"),
                    ("STAT","signalling / activation TF"),("POLR","other (replication/translation/misc)")]:
        if g.startswith(pre): return lab
    return "other / unclassified"

top_regs = []
for g, c in reg_counts.head(25).items():
    top_regs.append({"gene": str(g), "n_inverted_targets": int(c), "class": klass(str(g))})

# Descriptive positive-control rates are regrouped from the current identity-safe master.
# This preserves Figure S1b without importing the retired pooled-score replication object.
positive_control_genes = ["MED12", "MED24", "CBFB", "PRDM1", "STAT5B", "STAT3",
                          "BATF", "CCNC", "TADA2B"]
positive_controls = {}
for gene in positive_control_genes:
    sub = trans[trans.regulator.eq(gene)]
    if sub.empty:
        continue
    n_gene_inv = int(sub.kind.eq("rev").sum())
    positive_controls[gene] = {
        "n_inv": n_gene_inv,
        "n_both": int(len(sub)),
        "rate": round(n_gene_inv / len(sub), 3),
    }

# exemplar edges for slope plot: strongest reversals (both |lfc|>1), most extreme by delta.
# The census admits edges at adjusted p < 0.10, but an exemplar shown as a concrete case must
# not sit near that boundary: ranking on magnitude alone previously surfaced edges whose weaker
# side was only q ~ 0.05-0.08. Exemplars are therefore additionally required to clear q < 0.01
# in BOTH states. Adjusted p-values are not carried in sign_inverting_edges.csv, so they are
# joined from the frozen edge table.
_q = pd.read_csv(os.path.join(ROOT, "data/paper1_flips/frozen_edges.tsv"), sep="\t",
                 usecols=["regulator", "target", "adjp_Rest", "adjp_Stim48"])
strong = edges.merge(_q, on=["regulator", "target"], how="left")
strong = strong[(np.abs(strong.lfc_Rest) > 1) & (np.abs(strong.lfc_Stim48) > 1)].copy()
EXEMPLAR_Q_MAX = 0.01
strong = strong[strong[["adjp_Rest", "adjp_Stim48"]].max(axis=1) < EXEMPLAR_Q_MAX]
strong["delta"] = (strong.lfc_Rest - strong.lfc_Stim48).abs()
strong = strong.sort_values("delta", ascending=False)
# keep a readable handful spanning both flip directions
ex_pn = strong[(strong.lfc_Rest > 0)].head(4)
ex_np = strong[(strong.lfc_Rest < 0)].head(4)
exemplars = pd.concat([ex_pn, ex_np]).to_dict("records")
exemplars = [{k: (float(v) if isinstance(v,(int,float,np.floating)) else str(v))
              for k,v in r.items()} for r in exemplars]

out = {
    "census": {
        "n_eligibility_gate_passed": int(feas["n_candidate_perturbations"]),
        "n_census_regulators": int(trans["regulator"].nunique()),
        "n_gate_passed_without_shared_edge": int(
            feas["n_candidate_perturbations"] - trans["regulator"].nunique()),
        # Backward-compatible alias for the represented census family. New code should use
        # n_census_regulators; this is not the upstream gate-passing count.
        "n_candidate_regulators": int(trans["regulator"].nunique()),
        "min_outdeg": feas["min_outdeg"],
        "self_effects_excluded": n_self,
        "trans_edge_definition": "regulator_ensg != target_ensg (master is_self == false)",
        "edges_sig_both": n_both, "same_sign": n_same, "inverting": n_inv,
        "inversion_rate": round(rate, 4), "n_inverting_regulators": n_reg,
        "pos_to_neg": pos_to_neg, "neg_to_pos": neg_to_pos,
        "median_min_abs_lfc": round(float(np.median(min_abs)), 4),
    },
    "magnitude_survival": mag,
    "positive_controls": positive_controls,
    "top_regulators": top_regs,
    "exemplars": exemplars,
    # arrays for the scatter / hist (embedded; 4,379 pts is small)
    "arr": {
        "lfc_Rest": [round(float(x),4) for x in lfcR],
        "lfc_Stim48": [round(float(x),4) for x in lfcS],
        "min_abs": [round(float(x),4) for x in min_abs],
    },
    "all_regulator_counts": {str(k): int(v) for k, v in reg_counts.items()},
}
json.dump(out, open(os.path.join(ROOT, "data/census/figure_data.json"), "w"), indent=1)
print("\nWROTE data/census/figure_data.json")
print(f"  census: {n_inv}/{n_both} = {rate:.2%} inverting, {n_reg} regulators")
print(f"  split : {pos_to_neg} up->down, {neg_to_pos} down->up")
print(f"  core  : {mag['>0.5']} at |lfc|>0.5, {mag['>1.0']} at |lfc|>1.0")
print(f"  top regs: {[ (r['gene'],r['n_inverted_targets']) for r in top_regs[:6] ]}")
print(f"  exemplars: {[(e['regulator'],e['target']) for e in exemplars]}")
