#!/usr/bin/env python3
"""Try to kill B2 and B3 from regulator_target_structure.py.

B3 (69% of multi-reversal targets reverse unanimously, against a 37% null) is large enough to
be suspicious, so the job here is to find the mundane explanation before it is written up.
Five attacks, each targeting a different way the result could be an artifact:

  1. PAIRS-ONLY.  Targets with exactly 2 reversing edges make unanimity cheap. Restrict to
     >=3, where chance unanimity collapses. If the effect lives only at n=2, it is a
     small-count artifact.
  2. LEAVE-ONE-REGULATOR-OUT.  Five regulators carry most of the census. If dropping one
     collapses the effect, it is that regulator's direction bias leaking through, not a
     property of targets.
  3. DROP THE DOMINANT BLOC.  CCNC, AHR and ARNT together carry a large share of reversals and
     are functionally linked. Remove all three at once, which is a harder test than dropping
     any one.
  4. BASELINE-EXPRESSION STRATIFICATION (B2).  If highly expressed genes simply move the same
     way under most perturbations, sign concordance is a expression artifact rather than
     target-specific regulation. Recompute within baseMean deciles.
  5. IN-DEGREE STRATIFICATION (B2).  Same logic for connectivity.

Everything reuses the within-regulator null from the parent script, so a failure here is a
failure of the result, not of a changed null.

Reads   data/census/census_master_edges.csv.gz
Writes  data/census/target_theory_falsification.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
RD = PAPER / "data"

SEED = 20260812
N_PERM = 2000


def within_group_shuffle(values, group_codes, rng):
    r = rng.random(values.size)
    order = np.lexsort((r, group_codes))
    home = np.argsort(group_codes, kind="stable")
    out = np.empty_like(values)
    out[home] = values[order]
    return out


def unanimity(sub: pd.DataFrame, rng, min_rev: int = 2, n_perm: int = N_PERM) -> dict:
    """Fraction of targets with >=min_rev reversals whose reversals all run the same way."""
    if len(sub) == 0:
        return {"insufficient": True}
    reg = pd.Categorical(sub.regulator).codes
    tgt = pd.Categorical(sub.target).codes
    n_t = tgt.max() + 1
    up = (sub.sign_rest > 0).to_numpy().astype(float)
    n_r = np.bincount(tgt, minlength=n_t)
    keep = n_r >= min_rev
    if keep.sum() < 20:
        return {"insufficient": True, "n_targets": int(keep.sum())}

    def stat(lab):
        k = np.bincount(tgt, weights=lab, minlength=n_t)
        return float(np.mean((k[keep] == 0) | (k[keep] == n_r[keep])))

    obs = stat(up)
    null = np.array([stat(within_group_shuffle(up, reg, rng)) for _ in range(n_perm)])
    return {"n_targets": int(keep.sum()), "n_edges": int(len(sub)),
            "observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "excess": round(obs - float(null.mean()), 4),
            "z": round(float((obs - null.mean()) / null.std(ddof=1)), 2),
            "p": round(float((np.sum(null >= obs) + 1) / (n_perm + 1)), 5)}


def concordance(sub: pd.DataFrame, rng, n_perm: int = N_PERM) -> dict:
    """Same-sign edge-pair fraction within target, resting state."""
    reg = pd.Categorical(sub.regulator).codes
    tgt = pd.Categorical(sub.target).codes
    n_t = tgt.max() + 1
    pos = (sub.sign_rest > 0).to_numpy().astype(float)
    n = np.bincount(tgt, minlength=n_t).astype(float)
    keep = n >= 2
    if keep.sum() < 20:
        return {"insufficient": True, "n_targets": int(keep.sum())}

    def stat(lab):
        k = np.bincount(tgt, weights=lab, minlength=n_t)
        same = (k[keep] * (k[keep] - 1) / 2
                + (n[keep] - k[keep]) * (n[keep] - k[keep] - 1) / 2)
        return float(same.sum() / (n[keep] * (n[keep] - 1) / 2).sum())

    obs = stat(pos)
    null = np.array([stat(within_group_shuffle(pos, reg, rng)) for _ in range(n_perm)])
    return {"n_targets": int(keep.sum()), "n_edges": int(len(sub)),
            "observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "excess": round(obs - float(null.mean()), 4),
            "z": round(float((obs - null.mean()) / null.std(ddof=1)), 2),
            "p": round(float((np.sum(null >= obs) + 1) / (n_perm + 1)), 5)}


def main() -> None:
    rng = np.random.default_rng(SEED)
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "lfc_Rest", "baseMean_Rest",
                             "is_self"])
    e = e[~e.is_self].copy()
    assert len(e) == 83_489 and int((e.kind == "rev").sum()) == 4_379
    e["sign_rest"] = np.sign(e.lfc_Rest)
    e = e[e.sign_rest != 0].dropna(subset=["sign_rest"])
    rev = e[e.kind == "rev"].copy()

    res = {"seed": SEED, "n_perm": N_PERM}

    # --- attack 1: is it only the n=2 targets? -------------------------------
    res["A1_min_reversals_per_target"] = {
        str(m): unanimity(rev, rng, min_rev=m) for m in (2, 3, 4, 5)
    }

    # --- attack 2: leave one dominant regulator out --------------------------
    top = rev.regulator.value_counts().head(5)
    res["A2_leave_one_regulator_out"] = {
        "_top5_by_reversal_count": {k: int(v) for k, v in top.items()},
        **{g: unanimity(rev[rev.regulator != g], rng) for g in top.index}
    }

    # --- attack 3: drop the CCNC / AHR / ARNT bloc together ------------------
    bloc = ["CCNC", "AHR", "ARNT"]
    kept = rev[~rev.regulator.isin(bloc)]
    res["A3_drop_bloc"] = {
        "_dropped": bloc,
        "n_reversals_dropped": int(len(rev) - len(kept)),
        "frac_of_census_dropped": round(1 - len(kept) / len(rev), 3),
        "result": unanimity(kept, rng),
        "result_min3": unanimity(kept, rng, min_rev=3),
    }

    # --- attack 4/5: does B2 survive stratification? -------------------------
    e = e.copy()
    # Bin TARGETS, then carry those bins back to their edges. Ranking edges directly would let
    # high-in-degree targets occupy several bins and would not be a target-level stratification.
    target_meta = e.groupby("target").agg(
        baseMean_Rest=("baseMean_Rest", "median"),
        indeg=("target", "size"),
    )
    target_meta["bm_decile"] = pd.qcut(
        target_meta.baseMean_Rest.rank(method="first"), 10,
        labels=False, duplicates="drop")
    target_meta["indeg_decile"] = pd.qcut(
        target_meta.indeg.rank(method="first"), 10,
        labels=False, duplicates="drop")
    e = e.join(target_meta[["bm_decile", "indeg_decile"]], on="target")

    res["A4_B2_within_baseMean_decile"] = {
        "_design": ("if concordance is an expression artifact it should vanish inside a "
                    "decile of target baseline expression"),
        **{str(d): concordance(sub, rng) for d, sub in e.groupby("bm_decile")}
    }
    res["A5_B2_within_indegree_decile"] = {
        str(d): concordance(sub, rng) for d, sub in e.groupby("indeg_decile")
    }
    res["A0_B2_unstratified"] = concordance(e, rng)

    (RD / "census/target_theory_falsification.json").write_text(json.dumps(res, indent=2) + "\n")

    print("ATTACK 1  minimum reversals per target")
    for m, v in res["A1_min_reversals_per_target"].items():
        if v.get("insufficient"):
            print(f"   >={m}: insufficient targets ({v.get('n_targets','?')})"); continue
        print(f"   >={m}: {v['observed']:.3f} vs null {v['null_mean']:.3f}  "
              f"excess {v['excess']:+.3f}  z={v['z']:.1f}  p={v['p']}  n={v['n_targets']}")

    print("\nATTACK 2  leave-one-regulator-out")
    for g, v in res["A2_leave_one_regulator_out"].items():
        if g.startswith("_") or v.get("insufficient"):
            continue
        print(f"   drop {g:9s}: {v['observed']:.3f} vs {v['null_mean']:.3f}  "
              f"excess {v['excess']:+.3f}  n={v['n_targets']}")

    print("\nATTACK 3  drop CCNC+AHR+ARNT together "
          f"({res['A3_drop_bloc']['frac_of_census_dropped']:.0%} of reversals)")
    for k in ("result", "result_min3"):
        v = res["A3_drop_bloc"][k]
        if v.get("insufficient"):
            print(f"   {k}: insufficient ({v.get('n_targets','?')} targets)"); continue
        print(f"   {k}: {v['observed']:.3f} vs {v['null_mean']:.3f}  "
              f"excess {v['excess']:+.3f}  z={v['z']:.1f}  p={v['p']}  n={v['n_targets']}")

    print("\nATTACK 4/5  B2 concordance, unstratified then within deciles")
    u = res["A0_B2_unstratified"]
    print(f"   all edges          : {u['observed']:.4f} vs {u['null_mean']:.4f}  "
          f"excess {u['excess']:+.4f}")
    for tag, key in [("baseMean", "A4_B2_within_baseMean_decile"),
                     ("in-degree", "A5_B2_within_indegree_decile")]:
        ex = [v["excess"] for v in res[key].values()
              if isinstance(v, dict) and not v.get("insufficient")]
        ps = [v["p"] for v in res[key].values()
              if isinstance(v, dict) and not v.get("insufficient")]
        print(f"   within {tag:10s}: excess min {min(ex):+.4f} median "
              f"{np.median(ex):+.4f} max {max(ex):+.4f}; "
              f"{sum(p < 0.05 for p in ps)}/{len(ps)} deciles p<0.05")


if __name__ == "__main__":
    main()
