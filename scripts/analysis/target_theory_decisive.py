#!/usr/bin/env python3
"""The tests that decide whether B3 is a real target property or a restatement of B2.

THE CIRCULARITY. B2 established that edges landing on one target share resting sign more often
than chance. B3 then reported that reversals onto one target run the same way. But direction of
reversal IS the resting sign. So if a target's edges are mostly positive at rest, its reversing
edges are automatically mostly positive-to-negative, and B3 follows from B2 with no extra
biology. The within-regulator null used in B3 cannot see this, because permuting inside a
regulator destroys exactly the target-level sign bias that B2 says exists.

D1 provides the primary target-conditioned null. Within each target, permute WHICH of its
edges are the reversing ones,
holding every sign attached to its edge. This preserves the target's entire sign composition
and its reversal count exactly, so the only thing that can move the statistic is whether
reversal picks out sign-coherent edges beyond the target's own baseline coherence. If B3 is a
restatement of B2, D1 returns nothing. It does not simultaneously preserve regulator reversal
margins, so it supports target-associated concordance rather than an independent causal share.

D2 audits whether the concordance can be tested in held-out donors. The canonical
holdout_donors_edges.csv.gz contains stable Ensembl IDs, symbols, directional fold identifiers,
signs and effects in one row. Its explicit regulator_id-target_id-split-fold key replaces the
retired cross-file join on a floating-point magnitude.

D3 audits the independent Arce dataset, using only edges that actually reverse in Arce and
keeping its two assays separate.

D2 and D3 are too small for a confirmatory target-direction replication after one row per
regulator-target edge is enforced. They are retained as a power audit, not evidence.

Reads   data/census/census_master_edges.csv.gz
        data/holdout/holdout_donors_edges.csv.gz
        data/arce/arce_validation_edges.csv.gz
Writes  data/census/target_theory_decisive.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
RD = PAPER / "data"

SEED = 20260812
N_PERM = 5000


def within_group_shuffle(values: np.ndarray, group_codes: np.ndarray,
                         rng: np.random.Generator) -> np.ndarray:
    """Shuffle values independently within each integer-coded group."""
    out = np.empty_like(values)
    for code in np.unique(group_codes):
        idx = np.flatnonzero(group_codes == code)
        out[idx] = rng.permutation(values[idx])
    return out


def unanimity_conditional(df: pd.DataFrame, rng, min_rev: int = 2,
                          n_perm: int = N_PERM) -> dict:
    """Unanimity of reversing-edge direction, against a within-TARGET null.

    df needs columns: target, is_rev (0/1), pos (0/1 resting sign positive).
    Null reassigns the is_rev flag among a target's own edges, so the target's sign
    composition and reversal count are both held fixed.
    """
    d = df.dropna(subset=["pos"]).copy()
    # only targets that carry >=min_rev reversals AND have a spare non-reversing edge can move
    n_rev = d.groupby("target").is_rev.transform("sum")
    n_all = d.groupby("target").target.transform("size")
    d = d[(n_rev >= min_rev) & (n_all > n_rev)]
    if d.target.nunique() < 15:
        return {"insufficient": True, "n_targets": int(d.target.nunique())}

    codes = pd.Categorical(d.target).codes
    n_t = codes.max() + 1
    pos = d.pos.to_numpy().astype(float)
    isrev = d.is_rev.to_numpy().astype(float)
    k_rev = np.bincount(codes, weights=isrev, minlength=n_t)

    def stat(flag):
        k = np.bincount(codes, weights=flag * pos, minlength=n_t)
        return float(np.mean((k == 0) | (k == k_rev)))

    def shuffle(flag):
        r = rng.random(flag.size)
        order = np.lexsort((r, codes))
        home = np.argsort(codes, kind="stable")
        out = np.empty_like(flag)
        out[home] = flag[order]
        return out

    obs = stat(isrev)
    n_unanimous_obs = int(round(obs * n_t))
    null = np.array([stat(shuffle(isrev)) for _ in range(n_perm)])
    # Regulator-margin sensitivity on the SAME eligible targets. This is not nested within the
    # target-conditional null: it preserves regulator direction composition instead of target
    # sign composition, so the two nulls are shown as complementary checks, never subtracted.
    rev_mask = isrev == 1
    rev_tgt = codes[rev_mask]
    rev_reg = pd.Categorical(d.loc[rev_mask, "regulator"]).codes
    rev_pos = pos[rev_mask]

    def weak_stat(lab):
        k = np.bincount(rev_tgt, weights=lab, minlength=n_t)
        return float(np.mean((k == 0) | (k == k_rev)))

    weak_null = np.array([
        weak_stat(within_group_shuffle(rev_pos, rev_reg, rng))
        for _ in range(n_perm)
    ])
    return {"n_targets": int(n_t), "n_unanimous_observed": n_unanimous_obs,
            "n_edges": int(len(d)),
            "n_reversing": int(isrev.sum()),
            "observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "null_q025": round(float(np.quantile(null, 0.025)), 4),
            "null_q975": round(float(np.quantile(null, 0.975)), 4),
            "weak_null_mean_same_targets": round(float(weak_null.mean()), 4),
            "excess": round(obs - float(null.mean()), 4),
            "z": round(float((obs - null.mean()) / max(null.std(ddof=1), 1e-12)), 2),
            "p": round(float((np.sum(null >= obs) + 1) / (n_perm + 1)), 5)}


def unanimity_simple(df: pd.DataFrame, rng, min_rev: int = 2, n_perm: int = N_PERM) -> dict:
    """Unanimity among reversing edges only, null = shuffle direction within regulator.

    Used for D2/D3, where the non-reversing edges needed by the conditional null are not
    all present. Weaker null, and labelled as such wherever it is reported.
    """
    d = df[df.is_rev == 1].dropna(subset=["pos"]).copy()
    n_r = d.groupby("target").target.transform("size")
    d = d[n_r >= min_rev]
    if d.target.nunique() < 15:
        return {"insufficient": True, "n_targets": int(d.target.nunique()),
                "n_reversing_edges": int(len(d))}
    tgt = pd.Categorical(d.target).codes
    reg = pd.Categorical(d.regulator).codes
    n_t = tgt.max() + 1
    cnt = np.bincount(tgt, minlength=n_t)
    pos = d.pos.to_numpy().astype(float)

    def stat(lab):
        k = np.bincount(tgt, weights=lab, minlength=n_t)
        return float(np.mean((k == 0) | (k == cnt)))

    def shuffle(lab):
        r = rng.random(lab.size)
        order = np.lexsort((r, reg))
        home = np.argsort(reg, kind="stable")
        out = np.empty_like(lab)
        out[home] = lab[order]
        return out

    obs = stat(pos)
    null = np.array([stat(shuffle(pos)) for _ in range(n_perm)])
    return {"n_targets": int(n_t), "n_reversing_edges": int(len(d)),
            "observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "excess": round(obs - float(null.mean()), 4),
            "z": round(float((obs - null.mean()) / max(null.std(ddof=1), 1e-12)), 2),
            "p": round(float((np.sum(null >= obs) + 1) / (n_perm + 1)), 5)}


def main() -> None:
    rng = np.random.default_rng(SEED)
    res = {"seed": SEED, "n_perm": N_PERM}

    # ---------------- D1: conditional on the target's own sign composition ----
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "lfc_Rest", "is_self"])
    e = e[~e.is_self].copy()
    assert len(e) == 83_489 and int((e.kind == "rev").sum()) == 4_379
    e = e[np.sign(e.lfc_Rest) != 0].dropna(subset=["lfc_Rest"])
    e["is_rev"] = (e.kind == "rev").astype(int)
    e["pos"] = (e.lfc_Rest > 0).astype(float)

    res["D1_conditional_on_target_sign_composition"] = {
        "_design": ("within-TARGET permutation of which edges are reversing; the target's "
                    "resting knockdown-effect sign composition among its both-significant trans "
                    "edges and its reversal count are both held fixed. regulator "
                    "reversal margins are not simultaneously fixed, so this tests "
                    "target-associated concordance rather than a causal target contribution."),
        **{f"min_rev_{m}": unanimity_conditional(e, rng, min_rev=m) for m in (2, 3, 4)}
    }

    # ---------------- D2: held-out donor pair --------------------------------
    h = pd.read_csv(RD / "holdout/holdout_donors_edges.csv.gz")
    key = ["regulator_id", "target_id", "split", "fold"]
    required = set(key + ["regulator", "target", "kind", "repro", "sRV"])
    assert required <= set(h.columns), f"held-out edge table lacks {sorted(required-set(h.columns))}"
    assert not h.duplicated(key).any(), "held-out edge key is not unique"
    join_ok = {
        "canonical_file": "data/holdout/holdout_donors_edges.csv.gz",
        "n_rows": int(len(h)), "key": key, "key_unique": True,
        "identity_and_signs_colocated": True,
    }
    # Reciprocal train -> held-out folds must not be pooled: doing so repeats the same
    # regulator-target pair and manufactures multi-edge targets. Require held-out reproduction
    # and keep one row per regulator-target pair within each directional fold.
    d2 = {"_join_check": join_ok,
          "_design": ("one record per regulator-target edge within each directional donor fold; "
                      "training reversal must reproduce in the held-out pair; direction is the "
                      "held-out resting sign; folds are not pooled"),
          "_verdict": ("underpowered: only two of six directional folds retain at least 15 "
                       "targets with >=2 distinct reproduced reversal edges; no external "
                       "replication claim is made")}
    for i, (fold, sub) in enumerate(sorted(h.groupby("fold"), key=lambda z: z[0])):
        d = sub[(sub.kind == "rev") & (sub.repro == 1)].copy()
        d = d.drop_duplicates(["regulator", "target"])
        d["is_rev"] = 1
        d["pos"] = np.where(d.sRV > 0, 1.0, np.where(d.sRV < 0, 0.0, np.nan))
        v = unanimity_simple(d, rng)
        v.update({"fold": fold, "n_unique_reproduced_edges": int(len(d))})
        d2[f"fold_{i + 1}"] = v
    res["D2_held_out_donors"] = d2

    # ---------------- D3: independent Arce dataset ---------------------------
    a = pd.read_csv(RD / "arce/arce_validation_edges.csv.gz")
    d3 = {
        "_design": ("one record per regulator-target edge within each Arce assay; restricted "
                    "to our reversing edges that also reverse in Arce; direction is read from "
                    "Arce's own resting effect; assays are not pooled"),
        "_verdict": ("underpowered: neither assay retains 15 targets with >=2 distinct Arce-"
                     "reversing edges; no independent target-direction replication claim is made")
    }
    for dataset, sub in a.groupby("arce_dataset"):
        d = sub[(sub.our_kind == "rev") & (sub.arce_flipped == 1)].copy()
        d = d.drop_duplicates(["regulator", "target"])
        d["is_rev"] = 1
        d["pos"] = np.where(d.arce_effect_rest > 0, 1.0,
                            np.where(d.arce_effect_rest < 0, 0.0, np.nan))
        v = unanimity_simple(d, rng)
        v.update({"n_unique_arce_reversing_edges": int(len(d))})
        d3[str(dataset)] = v
    res["D3_arce_independent"] = d3

    (RD / "census/target_theory_decisive.json").write_text(json.dumps(res, indent=2) + "\n")

    print("D1  CONDITIONAL ON THE TARGET'S OWN SIGN COMPOSITION  (primary null)")
    for k, v in res["D1_conditional_on_target_sign_composition"].items():
        if k.startswith("_"):
            continue
        if v.get("insufficient"):
            print(f"   {k}: insufficient ({v['n_targets']} targets)"); continue
        print(f"   {k}: {v['observed']:.4f} vs null {v['null_mean']:.4f}  "
              f"excess {v['excess']:+.4f}  z={v['z']:.1f}  p={v['p']}  "
              f"n={v['n_targets']} targets")

    print("\nD2  HELD-OUT DONOR PAIR")
    print(f"   canonical keyed input: {join_ok['n_rows']} rows; "
          f"key unique {join_ok['key_unique']}; identities and signs colocated")
    for k, v in res["D2_held_out_donors"].items():
        if k.startswith("_"):
            continue
        if v.get("insufficient"):
            print(f"   {k}: insufficient ({v.get('n_targets','?')} targets, "
                  f"{v.get('n_reversing_edges','?')} reversing edges)"); continue
        print(f"   {k}: {v['observed']:.4f} vs {v['null_mean']:.4f}  "
              f"excess {v['excess']:+.4f}  p={v['p']}  n={v['n_targets']}")

    print("\nD3  INDEPENDENT ARCE DATASET")
    for k, v in res["D3_arce_independent"].items():
        if k.startswith("_"):
            continue
        if v.get("insufficient"):
            print(f"   {k}: insufficient ({v['n_targets']} targets with >=2, "
                  f"{v['n_reversing_edges']} edges)")
        else:
            print(f"   {k}: {v['observed']:.4f} vs {v['null_mean']:.4f}  "
                  f"excess {v['excess']:+.4f} p={v['p']} n={v['n_targets']}")


if __name__ == "__main__":
    main()
