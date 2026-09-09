#!/usr/bin/env python3
"""Who carries reversal: the regulator, the target, or both?

Three questions, answered on the SAME edge table so the parts are comparable.

A. REGULATORS.
   Figure 2a shows all 324 regulators with at least one Rest--Stim48 reversal and labels the
   20 strongest two-sided deviations from their out-degree-matched expectation. To avoid selecting
   on the outcome, the test is run first over all 617 regulators, including the 293 with zero
   reversals: within each out-degree decile, a conditional Fisher/hypergeometric test compares one
   regulator with all other edges in its decile, followed by BH-FDR over all 617. The original
   one-sided enrichment result is retained in the output as a sensitivity analysis.
   Out-degree predicts reversal rate (Spearman 0.235, p = 1.2e-4), so a per-regulator test
   against the pooled baseline rediscovers connectivity. Each regulator is therefore tested
   against the pooled rate of the OTHER regulators in its own out-degree decile (leave-one-out,
   exact binomial), with BH-FDR over candidates. The naive pooled-baseline test is reported
   alongside purely to show how many calls it inflates.

B. Target-associated direction: is reversal direction more concordant among edges landing on
   the same target? Three predictions are explored, each with its own null. Because CRISPRi
   coefficients are knockdown-versus-control effects and may be indirect, this does not identify
   a target as causally deciding whether the intact regulator activates or inhibits it.

   B1 concentration. Does reversal pile onto specific targets beyond what each regulator's own
      reversal load explains? Null permutes the rev label WITHIN each regulator, so every
      regulator keeps its exact reversal count and only the choice of target is randomised.
      The mirror test (permute within target, ask about regulators) uses the same statistic, so
      the two sides are directly comparable.
   B2 sign concordance. If direction is associated with target identity, edges landing on one target should
      share sign more than chance. Measured on ALL edges in the resting state, as the fraction
      of same-sign pairs within a target. Null permutes sign within regulator, which preserves
      each regulator's own positive/negative knockdown-effect balance.
   B3 direction concordance. Among targets carrying >=2 reversing edges, do those reversals go
      the same way (positive-to-negative vs negative-to-positive)? Same within-regulator null.

   B1 and B3 differ, and the difference is what makes the theory falsifiable: B1 can pass on
   target connectivity alone while B3 fails, which would mean targets attract reversal without
   carrying additional directional information.

C. ORA over targets. Per-target Fisher exact for over-representation among reversing edges,
   BH-FDR, run twice: against all other edges, and against other edges in the same in-degree
   decile. In-degree also predicts reversal rate (Spearman 0.106, p = 8.3e-17), so the
   unstratified version is the one to distrust.

Reads   data/census/census_master_edges.csv.gz (validated Rest--Stim48 census)
Writes  data/census/regulator_target_structure.json
        data/census/per_regulator_enrichment_all.csv
        data/census/per_regulator_tested_all.csv
        data/census/per_target_ora.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import binomtest, fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

PAPER = Path(__file__).resolve().parents[2]
RD = PAPER / "data"

SEED = 20260812
N_PERM = 2000
MIN_REG_EDGES = 50     # census candidate criterion, same as perreg_rates_and_enrichment.py
MIN_TGT_EDGES = 5
N_DECILES = 10


# ---------------------------------------------------------------- edge table

def edge_table() -> pd.DataFrame:
    """All trans edges with reversal label and resting sign from the master census.

    Regulator membership and rate denominators are defined only from Rest and Stim48hr.
    The intermediate 8 h subset is therefore not an admissible source for this analysis.
    """
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "lfc_Rest", "is_self"])
    e["regulator"] = e.regulator.astype(str)
    e["target"] = e.target.astype(str)
    e = e[~e.is_self].copy()
    e["rev"] = (e.kind == "rev").astype(int)
    e["sign_rest"] = np.sign(e.lfc_Rest)
    assert len(e) == 83_489
    assert e.regulator.nunique() == 617
    assert int(e.rev.sum()) == 4_379
    assert e.loc[e.rev == 1, "regulator"].nunique() == 324
    return e


# ------------------------------------------------- within-group permutation

def within_group_shuffle(values: np.ndarray, group_codes: np.ndarray,
                         rng: np.random.Generator) -> np.ndarray:
    """Permute `values` inside each group. Group sums are preserved exactly."""
    r = rng.random(values.size)
    order = np.lexsort((r, group_codes))          # sort by group, random inside group
    home = np.argsort(group_codes, kind="stable")  # positions grouped, original order inside
    out = np.empty_like(values)
    out[home] = values[order]
    return out


def concentration(labels: np.ndarray, codes: np.ndarray, n_groups: int, p: float) -> float:
    """Chi-square-like concentration of `labels` over groups: sum (k - n p)^2 / (n p (1-p)).

    Scale-free in the number of groups, so the regulator side and the target side of B1 are
    read off the same statistic.
    """
    n = np.bincount(codes, minlength=n_groups).astype(float)
    k = np.bincount(codes, weights=labels, minlength=n_groups)
    keep = n > 0
    return float(np.sum((k[keep] - n[keep] * p) ** 2 / (n[keep] * p * (1 - p))))


def perm_test(stat_fn, labels, shuffle_codes, rng, n_perm=N_PERM) -> dict:
    obs = stat_fn(labels)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stat_fn(within_group_shuffle(labels, shuffle_codes, rng))
    return {"observed": round(float(obs), 2),
            "null_mean": round(float(null.mean()), 2),
            "null_sd": round(float(null.std(ddof=1)), 2),
            "z": round(float((obs - null.mean()) / null.std(ddof=1)), 2),
            "fold": round(float(obs / null.mean()), 3),
            "p_one_sided": round(float((np.sum(null >= obs) + 1) / (n_perm + 1)), 5),
            "n_perm": n_perm}


# ------------------------------------------------------------ A: regulators

def test_figure2a_enrichment(e: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Exact degree-stratified regulator tests; display filtering happens afterwards.

    Fisher's exact test is the conditional hypergeometric test for one regulator's reversal
    count given the total reversals and edges in its out-degree decile. It remains exact for
    small denominators, so the arbitrary >=50-edge gate is not needed. The main Figure 2 test is
    two-sided, allowing both enrichment and depletion; the earlier one-sided enrichment result is
    retained in parallel. BH is applied over all 617 regulators before the 324 with n_inv > 0 are
    selected for plotting; testing only those 324 would select on the observed outcome.
    """
    g = (e.groupby("regulator").agg(n_both=("rev", "size"), n_inv=("rev", "sum"))
         .reset_index().rename(columns={"regulator": "gene"}))
    g["rate"] = g.n_inv / g.n_both
    g["decile"] = pd.qcut(g.n_both, N_DECILES, labels=False, duplicates="drop")

    dec_n = g.groupby("decile").n_both.transform("sum")
    dec_k = g.groupby("decile").n_inv.transform("sum")
    other_n = dec_n - g.n_both
    other_k = dec_k - g.n_inv
    g["p_expected"] = other_k / other_n
    g["enrichment_vs_decile"] = g.rate / g.p_expected

    odds, pvals_enrichment, pvals_two_sided = [], [], []
    for k, n, ok, on in zip(g.n_inv, g.n_both, other_k, other_n):
        table = [[int(k), int(n - k)], [int(ok), int(on - ok)]]
        result_enrichment = fisher_exact(table, alternative="greater")
        result_two_sided = fisher_exact(table, alternative="two-sided")
        odds.append(float(result_enrichment.statistic))
        pvals_enrichment.append(float(result_enrichment.pvalue))
        pvals_two_sided.append(float(result_two_sided.pvalue))
    g["odds_ratio"] = odds
    g["p_enrichment"] = pvals_enrichment
    g["q_enrichment"] = multipletests(g.p_enrichment, method="fdr_bh")[1]
    g["sig_enriched"] = g.q_enrichment < 0.05
    g["p_two_sided"] = pvals_two_sided
    g["q_two_sided"] = multipletests(g.p_two_sided, method="fdr_bh")[1]
    g["sig_two_sided"] = g.q_two_sided < 0.05
    g["direction"] = np.where(g.enrichment_vs_decile >= 1, "up", "down")
    g["shown_in_fig2a"] = g.n_inv > 0
    g["label_rank_fig2a"] = 0
    g["label_rank_fig2a_two_sided"] = 0

    label_enrichment = (g[g.shown_in_fig2a]
                        .sort_values(["p_enrichment", "enrichment_vs_decile", "n_inv", "gene"],
                                     ascending=[True, False, False, True])
                        .head(20))
    assert label_enrichment.sig_enriched.all(), (
        "a non-significant regulator entered the one-sided Figure 2a sensitivity top 20")
    g.loc[label_enrichment.index, "label_rank_fig2a"] = np.arange(
        1, len(label_enrichment) + 1)

    # Rank a two-sided result first by its exact p-value, then by absolute rate-ratio departure.
    # The latter is only a deterministic tie-break and does not change the inferential quantity.
    shown = g[g.shown_in_fig2a].copy()
    shown["_abs_log2_effect"] = np.abs(np.log2(shown.enrichment_vs_decile))
    label_two_sided = (shown
                       .sort_values(["p_two_sided", "_abs_log2_effect", "n_inv", "gene"],
                                    ascending=[True, False, False, True])
                       .head(20))
    assert label_two_sided.sig_two_sided.all(), (
        "a non-significant regulator entered the two-sided Figure 2a top 20")
    g.loc[label_two_sided.index, "label_rank_fig2a_two_sided"] = np.arange(
        1, len(label_two_sided) + 1)

    g = g.sort_values(["p_enrichment", "enrichment_vs_decile", "n_inv", "gene"],
                      ascending=[True, False, False, True])
    shown = g[g.shown_in_fig2a]
    labelled = g[g.label_rank_fig2a > 0].sort_values("label_rank_fig2a")
    labelled_two_sided = g[g.label_rank_fig2a_two_sided > 0].sort_values(
        "label_rank_fig2a_two_sided")
    summary = {
        "_design": ("two-sided Fisher/hypergeometric per regulator against all other edges in "
                    "its out-degree decile; BH-FDR over all 617 regulators before selecting "
                    "the 324 reversal-bearing regulators for display; one-sided enrichment "
                    "columns retained as a sensitivity analysis"),
        "n_regulators_tested": int(len(g)),
        "n_reversal_bearing_shown": int(g.shown_in_fig2a.sum()),
        "n_zero_reversal_in_multiplicity_family": int((g.n_inv == 0).sum()),
        "q_enrichment_lt_0.05_all": int(g.sig_enriched.sum()),
        "q_enrichment_lt_0.05_shown": int(shown.sig_enriched.sum()),
        "q_two_sided_lt_0.05_all": int(g.sig_two_sided.sum()),
        "q_two_sided_lt_0.05_all_enriched": int(
            (g.sig_two_sided & (g.direction == "up")).sum()),
        "q_two_sided_lt_0.05_all_depleted": int(
            (g.sig_two_sided & (g.direction == "down")).sum()),
        "q_two_sided_lt_0.05_shown": int(shown.sig_two_sided.sum()),
        "q_two_sided_lt_0.05_shown_enriched": int(
            (shown.sig_two_sided & (shown.direction == "up")).sum()),
        "q_two_sided_lt_0.05_shown_depleted": int(
            (shown.sig_two_sided & (shown.direction == "down")).sum()),
        "top20_labelled": labelled[
            ["gene", "n_both", "n_inv", "rate", "p_expected", "enrichment_vs_decile",
             "odds_ratio", "p_enrichment", "q_enrichment", "label_rank_fig2a"]
        ].round(10).to_dict("records"),
        "top20_labelled_two_sided": labelled_two_sided[
            ["gene", "n_both", "n_inv", "rate", "p_expected", "enrichment_vs_decile",
             "odds_ratio", "p_two_sided", "q_two_sided", "direction",
             "label_rank_fig2a_two_sided"]
        ].round(10).to_dict("records"),
    }
    return g, summary

def test_all_regulators(e: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    g = (e.groupby("regulator").agg(n_both=("rev", "size"), n_inv=("rev", "sum"))
         .reset_index().rename(columns={"regulator": "gene"}))
    baseline = g.n_inv.sum() / g.n_both.sum()
    g["rate"] = g.n_inv / g.n_both

    cand = g[g.n_both >= MIN_REG_EDGES].copy()
    cand["decile"] = pd.qcut(cand.n_both, N_DECILES, labels=False, duplicates="drop")

    # leave-one-out decile rate: the expectation this regulator is tested against
    dec_n = cand.groupby("decile").n_both.transform("sum")
    dec_k = cand.groupby("decile").n_inv.transform("sum")
    cand["p_expected"] = (dec_k - cand.n_inv) / (dec_n - cand.n_both)

    cand["p_naive"] = [binomtest(int(k), int(n), baseline, "two-sided").pvalue
                       for k, n in zip(cand.n_inv, cand.n_both)]
    cand["p_degree_matched"] = [binomtest(int(k), int(n), float(pe), "two-sided").pvalue
                                for k, n, pe in zip(cand.n_inv, cand.n_both, cand.p_expected)]
    for src, dst in [("p_naive", "q_naive"), ("p_degree_matched", "q_degree_matched")]:
        cand[dst] = multipletests(cand[src], method="fdr_bh")[1]
    cand["enrichment_vs_decile"] = cand.rate / cand.p_expected

    # A logistic regression with a per-regulator indicator and log10 out-degree, cluster-robust
    # by regulator, was tried as a second correction. It is KEPT ON RECORD BUT NOT USED, because
    # enrichment_test_calibration.py shows it is unusable at this level: `member` is constant
    # within regulator and there is exactly one treated cluster, the known failure case for
    # cluster-robust standard errors. On a null census with every regulator effect removed it
    # returns 239 and 232 false positives out of 264 at FDR. The degree-matched binomial returns
    # 0. Figure marks therefore come from the binomial ALONE; adding a broken second test would
    # have manufactured agreement rather than evidence.
    d = e[["regulator", "rev"]].copy()
    d["log_od"] = np.log10(d.groupby("regulator").regulator.transform("size"))
    rows = []
    for gene in cand.gene:          # NOT `g`: that is the regulator table in this scope
        d["member"] = (d.regulator == gene).astype(int)
        m = sm.GLM(d.rev, sm.add_constant(d[["member", "log_od"]]),
                   family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d.regulator})
        rows.append({"gene": gene, "OR_logistic": float(np.exp(m.params["member"])),
                     "p_logistic": float(m.pvalues["member"])})
    lg = pd.DataFrame(rows)
    lg["q_logistic"] = multipletests(lg.p_logistic, method="fdr_bh")[1]
    cand = cand.merge(lg, on="gene", how="left")

    # the annotation rule: the calibrated test only, and direction is carried because three of
    # the top 20 by count are significantly DEPLETED, where an unsigned mark would misread
    cand["sig"] = cand.q_degree_matched < 0.05
    cand["direction"] = np.where(cand.enrichment_vs_decile > 1, "up", "down")
    cand["logistic_uncalibrated_do_not_use"] = cand.pop("q_logistic")

    cand = cand.sort_values("p_degree_matched")
    summary = {
        "n_regulators_total": int(len(g)),
        "n_candidates_ge_%d_edges" % MIN_REG_EDGES: int(len(cand)),
        "n_reversal_bearing_shown_in_fig2a": int((g.n_inv > 0).sum()),
        "pooled_baseline_rate": round(float(baseline), 5),
        "spearman_outdegree_vs_rate": [round(float(x), 4) for x in
                                       spearmanr(cand.n_both, cand.rate)],
        "naive_q_lt_0.05": int((cand.q_naive < 0.05).sum()),
        "degree_matched_q_lt_0.05": int((cand.q_degree_matched < 0.05).sum()),
        "degree_matched_q_lt_0.05_enriched": int(((cand.q_degree_matched < 0.05) &
                                                  (cand.enrichment_vs_decile > 1)).sum()),
        "degree_matched_q_lt_0.05_depleted": int(((cand.q_degree_matched < 0.05) &
                                                  (cand.enrichment_vs_decile < 1)).sum()),
        "lost_by_degree_matching": int(((cand.q_naive < 0.05) &
                                        (cand.q_degree_matched >= 0.05)).sum()),
        "_logistic_withdrawn": ("a per-regulator cluster-robust logistic was computed and is "
                                "retained in the CSV as logistic_uncalibrated_do_not_use; it "
                                "has a ~90% false-positive rate on a null census (see "
                                "enrichment_test_calibration.json) and annotates nothing"),
        "sig_up": int((cand.sig & (cand.direction == "up")).sum()),
        "sig_down": int((cand.sig & (cand.direction == "down")).sum()),
        "top15_degree_matched": cand.head(15)[
            ["gene", "n_both", "n_inv", "rate", "p_expected", "enrichment_vs_decile",
             "q_degree_matched", "OR_logistic", "sig", "direction"]
        ].round(6).to_dict("records"),
        "fig2a_top20_test_ranked_enriched": cand[
            cand.sig & (cand.direction == "up")
        ].head(20)[
            ["gene", "n_both", "n_inv", "rate", "p_expected", "enrichment_vs_decile",
             "q_degree_matched"]
        ].round(8).to_dict("records"),
    }
    return cand, summary


# ------------------------------------------- B: is reversal target-determined

def target_theory(e: pd.DataFrame, rng: np.random.Generator) -> dict:
    reg_codes = pd.Categorical(e.regulator).codes
    tgt_cat = pd.Categorical(e.target)
    tgt_codes = tgt_cat.codes
    n_reg, n_tgt = reg_codes.max() + 1, tgt_codes.max() + 1
    rev = e.rev.to_numpy().astype(float)
    p = rev.mean()

    out = {"n_edges": int(len(e)), "n_regulators": int(n_reg), "n_targets": int(n_tgt),
           "pooled_rev_rate": round(float(p), 5)}

    # --- B1 concentration, both directions, same statistic -------------------
    out["B1_concentration"] = {
        "_design": ("target side: rev permuted WITHIN regulator, so each regulator keeps its "
                    "reversal count and only the choice of target moves. regulator side is "
                    "the mirror. same statistic, so fold values are comparable."),
        "targets_given_regulator": perm_test(
            lambda lab: concentration(lab, tgt_codes, n_tgt, p), rev, reg_codes, rng),
        "regulators_given_target": perm_test(
            lambda lab: concentration(lab, reg_codes, n_reg, p), rev, tgt_codes, rng),
    }

    # --- B2 sign concordance on ALL edges, resting state ----------------------
    s = e.dropna(subset=["sign_rest"])
    s = s[s.sign_rest != 0]
    s_reg = pd.Categorical(s.regulator).codes
    s_tgt_cat = pd.Categorical(s.target)
    s_tgt = s_tgt_cat.codes
    n_s_tgt = s_tgt.max() + 1
    pos = (s.sign_rest > 0).to_numpy().astype(float)

    def same_sign_frac(lab: np.ndarray) -> float:
        """Fraction of within-target edge PAIRS that share sign. Targets with n<2 drop out."""
        n = np.bincount(s_tgt, minlength=n_s_tgt).astype(float)
        k = np.bincount(s_tgt, weights=lab, minlength=n_s_tgt)
        keep = n >= 2
        same = k[keep] * (k[keep] - 1) / 2 + (n[keep] - k[keep]) * (n[keep] - k[keep] - 1) / 2
        tot = n[keep] * (n[keep] - 1) / 2
        return float(same.sum() / tot.sum())

    out["B2_sign_concordance_all_edges"] = {
        "_design": ("fraction of same-sign edge pairs within a target, resting state. null "
                    "permutes sign WITHIN regulator, preserving each regulator's own "
                    "positive/negative knockdown-effect balance."),
        "n_edges_with_sign": int(len(s)),
        "n_targets_ge2_edges": int((np.bincount(s_tgt, minlength=n_s_tgt) >= 2).sum()),
        "global_frac_positive": round(float(pos.mean()), 4),
        **perm_test(same_sign_frac, pos, s_reg, rng),
    }

    # --- B3 direction concordance among reversing edges ----------------------
    r = s[s.rev == 1].copy()
    r_reg = pd.Categorical(r.regulator).codes
    r_tgt_cat = pd.Categorical(r.target)
    r_tgt = r_tgt_cat.codes
    n_r_tgt = r_tgt.max() + 1
    up = (r.sign_rest > 0).to_numpy().astype(float)   # positive at rest, negative when activated
    n_r = np.bincount(r_tgt, minlength=n_r_tgt)

    def unanimous_frac(lab: np.ndarray) -> float:
        """Fraction of multi-reversal targets whose reversals all go the same way."""
        k = np.bincount(r_tgt, weights=lab, minlength=n_r_tgt)
        keep = n_r >= 2
        return float(np.mean((k[keep] == 0) | (k[keep] == n_r[keep])))

    out["B3_direction_concordance_reversing"] = {
        "_design": ("among targets with >=2 reversing edges, fraction where every reversal runs "
                    "the same way. this is what separates 'targets attract reversal' (B1) from "
                    "'direction is associated with target identity'. this weak null alone does "
                    "not condition on the target's baseline sign composition; see "
                    "target_theory_decisive.py for the primary conditional test."),
        "n_reversing_edges_with_sign": int(len(r)),
        "n_targets_ge2_reversals": int((n_r >= 2).sum()),
        "frac_pos_to_neg": round(float(up.mean()), 4),
        **perm_test(unanimous_frac, up, r_reg, rng),
    }
    return out


# --------------------------------------------------------------- C: ORA

def target_ora(e: pd.DataFrame, baseline: float) -> tuple[pd.DataFrame, dict]:
    t = (e.groupby("target").agg(n_both=("rev", "size"), n_inv=("rev", "sum"))
         .reset_index().rename(columns={"target": "gene"}))
    t = t[t.n_both >= MIN_TGT_EDGES].copy()
    N, K = int(e.rev.size), int(e.rev.sum())
    t["rate"] = t.n_inv / t.n_both
    t["decile"] = pd.qcut(t.n_both, N_DECILES, labels=False, duplicates="drop")

    dec_n = t.groupby("decile").n_both.transform("sum")
    dec_k = t.groupby("decile").n_inv.transform("sum")
    t["p_expected"] = (dec_k - t.n_inv) / (dec_n - t.n_both)
    t["enrichment_vs_decile"] = t.rate / t.p_expected

    # BOTH columns are the same test with the same tail; only the expected rate differs.
    # An earlier version ran the unstratified arm as Fisher and the stratified arm as a
    # binomial, so the difference between them mixed a change of expectation with a change of
    # test form and the comparison meant nothing. Fisher is kept as a third column so the
    # classical ORA is still on record, but the "lost to stratification" count is now computed
    # between the two binomials, which is the only pair that isolates the confound.
    t["p_ora_pooled"] = [binomtest(int(k), int(n), baseline, "greater").pvalue
                         for k, n in zip(t.n_inv, t.n_both)]
    t["p_ora_indeg_matched"] = [binomtest(int(k), int(n), float(pe), "greater").pvalue
                                for k, n, pe in zip(t.n_inv, t.n_both, t.p_expected)]
    t["p_ora_fisher"] = [fisher_exact([[k, n - k], [K - k, (N - n) - (K - k)]], "greater")[1]
                         for k, n in zip(t.n_inv, t.n_both)]
    for src, dst in [("p_ora_pooled", "q_ora_pooled"),
                     ("p_ora_indeg_matched", "q_ora_indeg_matched"),
                     ("p_ora_fisher", "q_ora_fisher")]:
        t[dst] = multipletests(t[src], method="fdr_bh")[1]

    t = t.sort_values("p_ora_indeg_matched")
    summary = {
        "_design": ("one-sided exact binomial per target for over-representation among "
                    "reversing edges, BH-FDR. pooled vs in-degree-matched differ ONLY in the "
                    "expected rate, so their difference isolates the in-degree confound. the "
                    "in-degree-matched column is the one to read."),
        "n_targets_tested_ge_%d_edges" % MIN_TGT_EDGES: int(len(t)),
        "q_ora_pooled_lt_0.05": int((t.q_ora_pooled < 0.05).sum()),
        "q_ora_indeg_matched_lt_0.05": int((t.q_ora_indeg_matched < 0.05).sum()),
        "q_ora_fisher_lt_0.05": int((t.q_ora_fisher < 0.05).sum()),
        "lost_by_indegree_matching": int(((t.q_ora_pooled < 0.05) &
                                          (t.q_ora_indeg_matched >= 0.05)).sum()),
        "gained_by_indegree_matching": int(((t.q_ora_pooled >= 0.05) &
                                            (t.q_ora_indeg_matched < 0.05)).sum()),
        "top15": t.head(15)[["gene", "n_both", "n_inv", "rate", "p_expected",
                             "enrichment_vs_decile", "p_ora_pooled", "q_ora_pooled",
                             "p_ora_indeg_matched", "q_ora_indeg_matched"]]
                  .round(8).to_dict("records"),
    }
    return t, summary


def main() -> None:
    rng = np.random.default_rng(SEED)
    e = edge_table()

    reg_fig2a, reg_fig2a_sum = test_figure2a_enrichment(e)
    reg_fig2a.to_csv(RD / "census/per_regulator_enrichment_all.csv", index=False)

    reg, reg_sum = test_all_regulators(e)
    reg.to_csv(RD / "census/per_regulator_tested_all.csv", index=False)

    theory = target_theory(e, rng)

    tgt, tgt_sum = target_ora(e, float(e.rev.mean()))
    tgt.to_csv(RD / "census/per_target_ora.csv", index=False)

    res = {"seed": SEED, "n_perm": N_PERM,
           "A0_figure2a_enrichment_all": reg_fig2a_sum,
           "A_all_regulators": reg_sum,
           "B_target_determines_direction": theory,
           "C_target_ora": tgt_sum}
    (RD / "census/regulator_target_structure.json").write_text(json.dumps(res, indent=2) + "\n")

    print("A0. FIGURE 2A TWO-SIDED TEST (ALL REGULATORS TESTED BEFORE DISPLAY FILTER)")
    print(f"   {reg_fig2a_sum['n_regulators_tested']} tested; "
          f"{reg_fig2a_sum['n_reversal_bearing_shown']} reversal-bearing shown; "
          f"among those shown, {reg_fig2a_sum['q_two_sided_lt_0.05_shown']} differ at q<0.05 "
          f"({reg_fig2a_sum['q_two_sided_lt_0.05_shown_enriched']} enriched, "
          f"{reg_fig2a_sum['q_two_sided_lt_0.05_shown_depleted']} depleted); "
          f"{reg_fig2a_sum['q_two_sided_lt_0.05_all']} calls in the full family")
    print("   top 20 labels: " + ", ".join(
        r["gene"] for r in reg_fig2a_sum["top20_labelled_two_sided"]))

    print("\nA1. >=50-EDGE TWO-SIDED REGULATOR TEST")
    print(f"   {reg_sum['n_regulators_total']} regulators total, "
          f"{reg_sum['n_candidates_ge_50_edges']} candidates, "
          f"{reg_sum['n_reversal_bearing_shown_in_fig2a']} reversal-bearing regulators "
          f"plotted in Fig 2a")
    print(f"   naive q<0.05:          {reg_sum['naive_q_lt_0.05']}")
    print(f"   degree-matched q<0.05: {reg_sum['degree_matched_q_lt_0.05']} "
          f"({reg_sum['degree_matched_q_lt_0.05_enriched']} up, "
          f"{reg_sum['degree_matched_q_lt_0.05_depleted']} down); "
          f"{reg_sum['lost_by_degree_matching']} calls lost to degree matching")
    print(f"   FIGURE MARKS (degree-matched binomial only): "
          f"{reg_sum['sig_up']} up, {reg_sum['sig_down']} down")
    for r in reg_sum["top15_degree_matched"][:8]:
        print(f"     {r['gene']:10s} {r['n_inv']:5d}/{r['n_both']:5d}  "
              f"{r['enrichment_vs_decile']:.2f}x vs decile  q={r['q_degree_matched']:.3g}  "
              f"{'SIG' if r['sig'] else '-':3s} {r['direction']}")

    print("\nB. IS REVERSAL DIRECTION ASSOCIATED WITH TARGET IDENTITY?")
    b1 = theory["B1_concentration"]
    print(f"   B1 targets | regulator : {b1['targets_given_regulator']['fold']:.3f}x  "
          f"z={b1['targets_given_regulator']['z']:.1f}  "
          f"p={b1['targets_given_regulator']['p_one_sided']}")
    print(f"   B1 regulators | target : {b1['regulators_given_target']['fold']:.3f}x  "
          f"z={b1['regulators_given_target']['z']:.1f}  "
          f"p={b1['regulators_given_target']['p_one_sided']}")
    b2 = theory["B2_sign_concordance_all_edges"]
    print(f"   B2 same-sign pairs     : {b2['observed']:.4f} vs null {b2['null_mean']:.4f}  "
          f"z={b2['z']:.1f}  p={b2['p_one_sided']}")
    b3 = theory["B3_direction_concordance_reversing"]
    print(f"   B3 unanimous direction : {b3['observed']:.4f} vs null {b3['null_mean']:.4f}  "
          f"z={b3['z']:.1f}  p={b3['p_one_sided']}  "
          f"(n={b3['n_targets_ge2_reversals']} targets)")

    print("\nC. TARGET ORA")
    print(f"   {tgt_sum['n_targets_tested_ge_5_edges']} targets tested (same test, "
          f"only the expected rate differs)")
    print(f"     pooled baseline   q<0.05: {tgt_sum['q_ora_pooled_lt_0.05']}")
    print(f"     in-degree matched q<0.05: {tgt_sum['q_ora_indeg_matched_lt_0.05']} "
          f"({tgt_sum['lost_by_indegree_matching']} lost, "
          f"{tgt_sum['gained_by_indegree_matching']} gained)")
    print(f"     Fisher (on record)  q<0.05: {tgt_sum['q_ora_fisher_lt_0.05']}")
    for r in tgt_sum["top15"][:8]:
        print(f"     {r['gene']:10s} {r['n_inv']:3d}/{r['n_both']:3d}  "
              f"{r['enrichment_vs_decile']:.2f}x  p={r['p_ora_indeg_matched']:.3g}  "
              f"q={r['q_ora_indeg_matched']:.3g}")


if __name__ == "__main__":
    main()
