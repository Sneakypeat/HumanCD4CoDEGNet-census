#!/usr/bin/env python
"""Shared statistics for the Arce external validation.

Design logic (model B):
  A genuinely STABLE edge carrying ONE ERRONEOUS ENDPOINT in our atlas looks
  exactly like a reversal in our data. An endpoint error in OUR data cannot
  reproduce in an INDEPENDENT experiment. So:
    - if our reversals are real state-dependent sign changes, they should
      re-flip in Arce more often than our stable edges do  -> OR > 1
    - if model B holds, our "reversals" are stable edges + noise, and they
      should flip in Arce no more often than our stable edges -> OR = 1

  Stratifying by regulator is what cancels Arce's library selection bias
  (28 genes chosen "prioritizing trans-factors with state-specific effects"):
  both arms are drawn from the SAME perturbations within a stratum.

  Inference clusters on REGULATOR, never on edges: edges sharing a regulator
  share the perturbation, the cells and the knockdown efficiency, so they are
  not independent.
"""
import numpy as np
import pandas as pd
from scipy import stats


def strata_tables(df, reg_col="regulator", kind_col="kind", flip_col="arce_flip"):
    """Per-regulator 2x2: rows = (rev, stable), cols = (flip, no-flip)."""
    out = []
    for reg, d in df.groupby(reg_col):
        r, s = d[d[kind_col] == "rev"], d[d[kind_col] == "stable"]
        a = int(r[flip_col].sum()); b = len(r) - a
        c = int(s[flip_col].sum()); dd = len(s) - c
        out.append((reg, a, b, c, dd))
    return out


def mh_odds_ratio(tabs):
    """Mantel-Haenszel pooled odds ratio. tabs = [(reg,a,b,c,d), ...]."""
    num = den = 0.0
    for _, a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
    if den == 0:
        return np.inf if num > 0 else np.nan
    return num / den


def cmh_test(tabs, correct=True):
    """Cochran-Mantel-Haenszel chi-square (1 df). Edge-level; assumes
    independence WITHIN stratum -- reported only alongside the cluster
    bootstrap, which does not."""
    A = E = V = 0.0
    for _, a, b, c, d in tabs:
        n = a + b + c + d
        if n < 2:
            continue
        r1, r2 = a + b, c + d
        c1, c2 = a + c, b + d
        A += a
        E += r1 * c1 / n
        V += r1 * r2 * c1 * c2 / (n * n * (n - 1))
    if V == 0:
        return np.nan, np.nan
    num = abs(A - E) - (0.5 if correct else 0.0)
    num = max(num, 0.0)
    chi2 = num ** 2 / V
    return chi2, stats.chi2.sf(chi2, 1)


def rbg_ci(tabs, alpha=0.05):
    """Robins-Breslow-Greenland CI for the MH OR (edge-level, for reference)."""
    R = S = 0.0
    PR = PS_QR = QS = 0.0
    for _, a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        Rk, Sk = a * d / n, b * c / n
        R += Rk; S += Sk
        P, Q = (a + d) / n, (b + c) / n
        PR += P * Rk
        PS_QR += P * Sk + Q * Rk
        QS += Q * Sk
    if R == 0 or S == 0:
        return np.nan, np.nan
    var = PR / (2 * R * R) + PS_QR / (2 * R * S) + QS / (2 * S * S)
    or_mh = R / S
    z = stats.norm.ppf(1 - alpha / 2)
    lo = np.exp(np.log(or_mh) - z * np.sqrt(var))
    hi = np.exp(np.log(or_mh) + z * np.sqrt(var))
    return lo, hi


def cluster_bootstrap_or(df, n_boot=10000, seed=0, reg_col="regulator",
                         flip_col="arce_flip"):
    """Resample REGULATORS (clusters) with replacement; recompute the MH OR.

    This is the primary inference: it propagates the fact that we have 11
    independent units, not 3,467 independent edges.
    """
    rng = np.random.default_rng(seed)
    regs = df[reg_col].unique()
    by = {r: d for r, d in df.groupby(reg_col)}
    ors = []
    for _ in range(n_boot):
        pick = rng.choice(regs, size=len(regs), replace=True)
        tabs = []
        for i, r in enumerate(pick):
            d = by[r]
            rv, st = d[d.kind == "rev"], d[d.kind == "stable"]
            a = int(rv[flip_col].sum()); b = len(rv) - a
            c = int(st[flip_col].sum()); dd = len(st) - c
            tabs.append((f"{r}_{i}", a, b, c, dd))
        o = mh_odds_ratio(tabs)
        if np.isfinite(o):
            ors.append(o)
    ors = np.array(ors)
    return ors


def leave_one_out(df, reg_col="regulator", flip_col="arce_flip"):
    """Drop each regulator in turn; recompute the MH OR and CMH p."""
    rows = []
    for r in sorted(df[reg_col].unique()):
        sub = df[df[reg_col] != r]
        tabs = strata_tables(sub, reg_col=reg_col, flip_col=flip_col)
        o = mh_odds_ratio(tabs)
        _, p = cmh_test(tabs)
        n_rev = int((sub.kind == "rev").sum())
        rows.append({"dropped": r, "or_mh": o, "cmh_p": p, "n_rev_left": n_rev})
    return rows


def match_within_regulator(df, on="minmag", reg_col="regulator", k=1,
                           caliper=None, seed=0):
    """1:k nearest-neighbour match of stable edges to each rev edge, WITHIN
    regulator, on `on` (without replacement).

    Why this is mandatory. Our reversal edges carry systematically smaller
    effect magnitudes than our stable edges (selection: an edge is called a
    reversal partly because one endpoint is small). A genuinely STABLE edge
    with a SMALL true effect flips in Arce by measurement noise alone, and its
    flip probability rises towards 0.5 as the effect shrinks. So an unmatched
    stable arm understates the model-B flip baseline and inflates the odds
    ratio. Matching on magnitude makes the stable arm a like-for-like model-B
    null.
    """
    keep_idx = []          # collect INDEX LABELS, never row-frames:
    for reg, d in df.groupby(reg_col):          # .to_frame().T casts everything
        rev = d[d.kind == "rev"]                # to object dtype and silently
        stab = d[d.kind == "stable"]            # breaks downstream numerics.
        if len(rev) == 0 or len(stab) == 0:
            continue
        used = set()
        for ridx, row in rev.iterrows():
            avail = stab[~stab.index.isin(used)]
            if len(avail) == 0:
                break
            dist = (avail[on] - row[on]).abs()
            if caliper is not None:
                dist = dist[dist <= caliper]
                if len(dist) == 0:
                    continue
            take = list(dist.nsmallest(min(k, len(dist))).index)
            used.update(take)
            keep_idx.append(ridx)
            keep_idx.extend(take)
    if not keep_idx:
        return df.iloc[0:0]
    return df.loc[keep_idx].reset_index(drop=True)


def quantile_strata(df, on="minmag", reg_col="regulator", q=5):
    """Strata = regulator x within-regulator quantile bin of `on`.
    Keeps the CMH framework while controlling magnitude."""
    out = []
    for reg, d in df.groupby(reg_col):
        d = d.copy()
        try:
            d["_bin"] = pd.qcut(d[on], q=min(q, d[on].nunique()),
                                labels=False, duplicates="drop")
        except (ValueError, IndexError):
            d["_bin"] = 0
        d["_stratum"] = d[reg_col].astype(str) + "|q" + d["_bin"].astype(str)
        out.append(d)
    return pd.concat(out, ignore_index=True)


def per_regulator_rates(df, reg_col="regulator", flip_col="arce_flip"):
    rows = []
    for reg, d in df.groupby(reg_col):
        r, s = d[d.kind == "rev"], d[d.kind == "stable"]
        a = int(r[flip_col].sum()); b = len(r) - a
        c = int(s[flip_col].sum()); dd = len(s) - c
        # Haldane-Anscombe correction for display only
        or_k = ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({
            "regulator": reg, "n_rev": len(r), "n_stable": len(s),
            "rev_flip": a, "rev_flip_rate": a / len(r) if len(r) else np.nan,
            "stable_flip": c,
            "stable_flip_rate": c / len(s) if len(s) else np.nan,
            "or_haldane": or_k,
        })
    return rows
