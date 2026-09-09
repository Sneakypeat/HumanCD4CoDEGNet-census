#!/usr/bin/env python
"""
Is the convergence of reversing edges on shared targets REAL, or an artefact of which genes are
heavily regulated in the first place?

THE CLAIM SPLITS IN TWO, AND ONLY ONE IS SAFELY TESTABLE:
  A. CONVERGENCE  -- do reversing edges concentrate on fewer targets than chance? No annotation
     needed, so no annotation-density trap. This is the testable half.
  B. EFFECTOR IDENTITY -- are the convergent targets *the effector programme*? This needs a gene
     list, and effector genes are both heavily regulated AND heavily studied, which is exactly the
     confound that has already produced three artefacts in this project (GO annotation density
     rho=-0.448; the dark-matter split reporting input composition; AF3 n_high tracking input size).

THE NULL (this is the whole design). Naive nulls fail here:
  * shuffling targets genome-wide ignores that a regulator can only reverse a target it regulates;
  * any null ignoring target in-degree "discovers" that popular genes are hit often -- by
    construction, since a gene regulated by 300 regulators will be reversed by some of them.
So: DEGREE-PRESERVING BIPARTITE null. For each regulator independently, resample its OBSERVED
number of reversals uniformly from ITS OWN both-significant targets. This preserves
  (i)  every regulator's out-degree (number of reversals it contributes), and
  (ii) target availability -- a target can only be drawn by regulators that actually regulate it,
       so widely-regulated genes are frequently drawn IN THE NULL TOO.
"Effector genes are popular" therefore cannot generate a positive result: the null is already
allowed to exploit it. Only concentration BEYOND that is counted.

For B we reuse the same null: gene-set membership is external (GO BP 2023 via Enrichr), and the
comparison is observed convergent targets vs NULL convergent targets drawn from the same pools --
so annotation density cancels, because the null's targets are drawn from the same regulated genes.

The analysis universe is the full Rest--Stim48 trans census. The master remains inclusive for
provenance, so self effects are excluded by its Ensembl-ID-derived `is_self` flag. The 8 h
subset is not admissible because 8 h availability is unrelated to this Rest--Stim48 hypothesis.

log_fc is log2 (DESeq2). Usage: convergence_test.py
"""
import os
from pathlib import Path
import json
import numpy as np
import pandas as pd
from collections import Counter

ROOT = str(Path(__file__).resolve().parents[2])
GS = str(Path(__file__).resolve().parents[2] / "data/genesets")
NPERM = 5000
SEED = 20260717
rng = np.random.default_rng(SEED)

SOURCE = "data/census/census_master_edges.csv.gz"
df = pd.read_csv(f"{ROOT}/{SOURCE}",
                 usecols=["regulator", "target", "kind", "is_self"])
df = df.loc[~df.is_self].copy()
rev = df[df.kind == "rev"]
assert len(df) == 83_489
assert int((df.kind == "stable").sum()) == 79_110
assert len(rev) == 4_379
assert rev.regulator.nunique() == 324
print(f"edges: {len(df):,}  reversing: {len(rev):,}  regulators with reversals: {rev.regulator.nunique()}")

# per-regulator candidate pools (both-significant targets) and observed reversal counts
pools, kobs = {}, {}
for r, g in df.groupby("regulator"):
    n = int((g.kind == "rev").sum())
    if n == 0: continue
    pools[r] = g.target.values          # only targets this regulator actually regulates
    kobs[r] = n
print(f"regulators contributing: {len(pools)}   total reversals: {sum(kobs.values()):,}")

def indeg(counter): return np.array(list(counter.values()))
def stats(counter):
    v = indeg(counter)
    return {"n_targets": len(v), "ge2": int((v >= 2).sum()), "ge3": int((v >= 3).sum()),
            "ge5": int((v >= 5).sum()), "max": int(v.max()) if len(v) else 0}

obs_c = Counter(rev.target)
obs = stats(obs_c)
print(f"\nOBSERVED: {obs['n_targets']:,} distinct targets receive reversals; "
      f">=2 regulators: {obs['ge2']:,}   >=3: {obs['ge3']:,}   >=5: {obs['ge5']:,}   max: {obs['max']}")

# ---- A. degree-preserving bipartite null ----
print(f"\nrunning {NPERM} degree-preserving permutations ...", flush=True)
null = {k: np.empty(NPERM) for k in ("n_targets", "ge2", "ge3", "ge5", "max")}
null_ge5_sets = []
for i in range(NPERM):
    c = Counter()
    for r, pool in pools.items():
        c.update(rng.choice(pool, kobs[r], replace=False))
    s = stats(c)
    for k in null: null[k][i] = s[k]
    if i < 200:  # keep some null convergent sets for part B
        null_ge5_sets.append({t for t, n in c.items() if n >= 5})

res = {
    "source": SOURCE,
    "trans_definition": "master is_self == false; is_self is endpoint Ensembl-ID equality",
    "universe": {"edges": int(len(df)), "stable": int((df.kind == "stable").sum()),
                 "reversing": int(len(rev)),
                 "reversal_bearing_regulators": int(rev.regulator.nunique())},
    "seed": SEED,
    "observed": obs,
    "nperm": NPERM,
    "null": {},
}
print(f"\n{'statistic':12s} {'observed':>10s} {'null mean':>10s} {'null 95%':>18s} {'p':>10s}")
for k in ("n_targets", "ge2", "ge3", "ge5", "max"):
    d = null[k]; o = obs[k]
    lo, hi = np.percentile(d, [2.5, 97.5])
    # two-sided-ish: fewer distinct targets / more shared targets both = more convergent
    p_hi = (np.sum(d >= o) + 1) / (NPERM + 1); p_lo = (np.sum(d <= o) + 1) / (NPERM + 1)
    p = min(p_hi, p_lo) * 2
    res["null"][k] = {"observed": o, "null_mean": round(float(d.mean()), 1),
                      "null_95": [round(float(lo), 1), round(float(hi), 1)], "p_two_sided": round(float(min(p, 1)), 5)}
    print(f"{k:12s} {o:10,d} {d.mean():10.1f} {'['+format(lo,'.0f')+', '+format(hi,'.0f')+']':>18s} {min(p,1):10.4g}")

# ---- B. are the convergent targets the effector programme? same null, so density cancels ----
def read_gmt(p):
    d = {}
    for l in open(p):
        f = l.rstrip("\n").split("\t"); g = {x.split(",")[0].strip() for x in f[2:] if x.strip()}
        if g: d[f[0]] = g
    return d
GO = read_gmt(f"{GS}/GO_Biological_Process_2023.gmt")
EFF = set()
for name, genes in GO.items():
    n = name.lower()
    if any(k in n for k in ("t cell activation", "t-cell activation", "cytokine production",
                            "immune effector process", "lymphocyte activation")):
        EFF |= genes
obs_ge5 = {t for t, n in obs_c.items() if n >= 5}
o_frac = len(obs_ge5 & EFF) / max(len(obs_ge5), 1)
n_fracs = np.array([len(s & EFF) / max(len(s), 1) for s in null_ge5_sets if len(s)])
p_eff = (np.sum(n_fracs >= o_frac) + 1) / (len(n_fracs) + 1) if len(n_fracs) else float("nan")
print(f"\nB. effector identity of the CONVERGENT (>=5) targets, vs the SAME null:")
print(f"   external gene set: GO BP activation/cytokine/effector terms, {len(EFF):,} genes")
print(f"   observed convergent targets: {len(obs_ge5)}   annotated: {len(obs_ge5 & EFF)} = {o_frac*100:.1f}%")
print(f"   null convergent targets    : {np.mean([len(s) for s in null_ge5_sets]):.0f}   annotated: {n_fracs.mean()*100:.1f}%"
      f"  (95% [{np.percentile(n_fracs,2.5)*100:.1f}, {np.percentile(n_fracs,97.5)*100:.1f}])")
verdict = ("EXCESS over null: annotation density cannot explain it, since the null draws from the\n      same regulated targets" if p_eff < 0.05 else
           "NO EXCESS over null: the convergent targets are annotated at the null rate, so the\n      effector-programme framing is NOT supported -- convergence is real, its identity is not")
print(f"   p = {p_eff:.4g}   -> {verdict}")
res["effector"] = {"n_geneset": len(EFF), "n_convergent": len(obs_ge5),
                   "obs_frac": round(float(o_frac), 4), "null_frac_mean": round(float(n_fracs.mean()), 4),
                   "p": round(float(p_eff), 5)}
json.dump(res, open(f"{ROOT}/data/census/convergence_test_results.json", "w"), indent=2)
print(f"\nwrote data/census/convergence_test_results.json")
