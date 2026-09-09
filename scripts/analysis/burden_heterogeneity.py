#!/usr/bin/env python3
"""Reversal-burden heterogeneity across the 617 census regulators.

Replaces the power-law/topology framing in `degree_distribution_test.py`, which the census
cannot identify: the upstream gate admitted a regulator only with at least 50 significant targets
in each state, so the observed out-degree distribution is a property of that gate rather than of
the network, and sampled subnets do not inherit a parent network's degree form in any case
(Stumpf, Wiuf & May, PNAS 2005).

The estimable question is whether the per-edge reversal probability is constant across
regulators. Modelling regulator i as k_i reversals out of its own n_i both-significant edges
conditions on that regulator's denominator, so the gate decides which n_i are observed rather
than the k-given-n relationship being tested.

Writes data/burden_heterogeneity/burden_heterogeneity.json.
Figure: scripts/fig_burden_heterogeneity.py.

CAVEAT CARRIED IN THE OUTPUT: the companion covariate closure found that measured edge
precision, on-target knockdown efficacy, nonlinear degree and available cell support did not
collapse rho. This excludes those measured explanations, not unmeasured technical variation,
and does not make the residual a regulator-intrinsic or biological mechanism.
"""
from __future__ import annotations
import hashlib, json, os, pathlib
from pathlib import Path

import numpy as np, pandas as pd
from scipy import stats, optimize

# Repo-relative so the scripts run from any checkout, not only this machine.
ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data/census/per_regulator_enrichment_all.csv"
OUT = ROOT / "data/burden_heterogeneity/burden_heterogeneity.json"


def load():
    d = pd.read_csv(SRC)
    assert len(d) == 617, f"expected the 617-regulator census family, got {len(d)}"
    assert d["gene"].notna().all() and d["gene"].is_unique, "regulator identities must be unique"
    n, k = d.n_both.to_numpy(float), d.n_inv.to_numpy(float)
    assert (k <= n).all() and (n > 0).all()
    assert int(n.sum()) == 83_489 and int(k.sum()) == 4_379
    assert int((k > 0).sum()) == 324
    return n, k


def fit(n, k, x, with_slope: bool, with_disp: bool):
    """Beta-binomial (or binomial) with logit mu = b0 + b1*(x - mean x); rho = intra-class corr."""
    def nll(par):
        b1 = par[1] if with_slope else 0.0
        mu = 1 / (1 + np.exp(-(par[0] + b1 * (x - x.mean()))))
        if not with_disp:
            return -stats.binom.logpmf(k, n, mu).sum()
        rho = 1 / (1 + np.exp(-par[-1]))
        a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
        if np.any(a <= 0) or np.any(b <= 0) or not np.all(np.isfinite(a + b)):
            return 1e12
        return -stats.betabinom.logpmf(k, n, a, b).sum()
    p0 = [np.log(0.05 / 0.95)] + ([0.0] if with_slope else []) + ([-3.0] if with_disp else [])
    r = optimize.minimize(nll, p0, method="Nelder-Mead",
                          options=dict(maxiter=40000, fatol=1e-11, xatol=1e-11))
    assert r.success and np.isfinite(r.fun) and np.all(np.isfinite(r.x)), \
        f"burden model optimisation failed: {r.message}"
    return -r.fun, r.x


def gini(v):
    v = np.sort(v)
    i = np.arange(1, len(v) + 1)
    return float((2 * (i * v).sum()) / (len(v) * v.sum()) - (len(v) + 1) / len(v))


def main():
    n, k = load()
    x, N, p = np.log10(n), len(n), k.sum() / n.sum()

    # Tarone Z: classical test of binomial overdispersion
    S = ((k - n * p) ** 2 / (p * (1 - p))).sum()
    Z = float((S - n.sum()) / np.sqrt(2 * (n * (n - 1)).sum()))

    ll_bin, _ = fit(n, k, x, False, False)
    ll_bb, pb = fit(n, k, x, False, True)
    ll_bind, _ = fit(n, k, x, True, False)
    ll_bbd, pd_ = fit(n, k, x, True, True)
    rho_u, rho_a, b1 = 1/(1+np.exp(-pb[-1])), 1/(1+np.exp(-pd_[-1])), float(pd_[1])

    robust = []
    for drop in (0, 5, 10, 20, 50):
        o = np.argsort(k)[::-1][drop:]
        xx = np.log10(n[o])
        a_, _ = fit(n[o], k[o], xx, True, False)
        b_, pp = fit(n[o], k[o], xx, True, True)
        robust.append(dict(dropped=drop, n_regulators=int(len(o)),
                           rho=float(1/(1+np.exp(-pp[-1]))), lrt_chi2=float(2*(b_-a_))))

    # How many regulators sit outside a binomial 95% interval around the fitted degree trend.
    # The naive benchmark 0.05*617 = 31 is wrong: binom.ppf returns a CONSERVATIVE discrete
    # interval (median realised coverage 98.5% here), so the exact independence expectation is
    # 11.5. Even that assumes regulators are independent, which they are not, because their
    # edges share targets. The defensible benchmark is the target-margin-preserving null in
    # scripts/burden_target_null.py.
    # Classify against the BINOMIAL-null fit, not the beta-binomial one: the interval
    # drawn is a binomial interval, so the mean it is centred on must come from the
    # model being tested. Mixing the two gave 147 here against 146 in the figure.
    bin_par = fit(n, k, x, True, False)[1]
    mu_i = 1 / (1 + np.exp(-(bin_par[0] + bin_par[1] * (x - x.mean()))))
    ni = n.astype(int)
    lo_c, hi_c = stats.binom.ppf(0.025, ni, mu_i), stats.binom.ppf(0.975, ni, mu_i)
    outside = int(((k < lo_c) | (k > hi_c)).sum())
    exact_expected = float((stats.binom.cdf(lo_c - 1, ni, mu_i) + stats.binom.sf(hi_c, ni, mu_i)).sum())
    realised_coverage = float(np.median(stats.binom.cdf(hi_c, ni, mu_i) - stats.binom.cdf(lo_c - 1, ni, mu_i)))

    srt = np.sort(k)[::-1]
    cum = srt.cumsum() / k.sum()
    share = lambda f: float(cum[max(1, int(round(f * N))) - 1])

    res = dict(
        _what="Reversal-burden heterogeneity across regulators; replaces the power-law framing.",
        _source=str(SRC.relative_to(ROOT)),
        _script="scripts/burden_heterogeneity.py",
        _source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
        _script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        _caveat=("Technical covariates were tested in burden_covariate_closure.py and do NOT "
                 "explain the dispersion (rho 0.0661 -> 0.0635 adding precision, efficacy and a "
                 "cubic degree term on the same 555 regulators). That excludes the measurable "
                 "technical explanations, not every one, and does not establish biology."),
        _not_new=("This is a single global effect size for the same heterogeneity the "
                  "per-regulator degree-matched Fisher tests already detect (106/617 at 5% FDR)."),
        _inferential_scope=("Tarone and likelihood-ratio quantities use an independent-regulator "
                            "working model and are diagnostics because targets are shared. The "
                            "target-margin permutation calibrates the outside-interval count, not "
                            "those working-model statistics."),
        n_regulators=N, n_edges=int(n.sum()), n_reversals=int(k.sum()), pooled_rate=float(p),
        tarone_Z=Z, tarone_p=float(stats.norm.sf(Z)),
        betabinom_unadjusted=dict(rho=float(rho_u), lrt_chi2=float(2*(ll_bb-ll_bin)),
                                  lrt_p=float(stats.chi2.sf(2*(ll_bb-ll_bin), 1))),
        degree_adjusted=dict(
            slope_per_log10_outdegree=b1,
            slope_lrt_chi2=float(2*(ll_bbd-ll_bb)), slope_p=float(stats.chi2.sf(2*(ll_bbd-ll_bb), 1)),
            residual_rho=float(rho_a), residual_lrt_chi2=float(2*(ll_bbd-ll_bind)),
            residual_p_naive_chi2_1=float(stats.chi2.sf(2*(ll_bbd-ll_bind), 1)),
            # rho = 0 sits on the parameter-space boundary, so the reference is a
            # 0.5*chi2_0 + 0.5*chi2_1 mixture (Self & Liang 1987), which halves the p-value.
            residual_p_boundary_corrected=float(0.5*stats.chi2.sf(2*(ll_bbd-ll_bind), 1)),
            regulators_outside_binomial_95_interval=outside,
            expected_outside_exact_under_independence=exact_expected,
            median_realised_interval_coverage=realised_coverage,
            expected_outside_naive_5pct_WRONG=round(0.05 * N),
            target_preserving_null="see data/burden_heterogeneity/target_null.json"),
        variance_explained=dict(
            spearman_outdegree_vs_count=float(stats.spearmanr(n, k).statistic),
            spearman_outdegree_vs_rate=float(stats.spearmanr(n, k/n).statistic),
            fitted_rate_at_min_outdegree=float(1/(1+np.exp(-(pd_[0]+b1*(x.min()-x.mean()))))),
            fitted_rate_at_max_outdegree=float(1/(1+np.exp(-(pd_[0]+b1*(x.max()-x.mean())))))),
        concentration=dict(gini=gini(k), share_top_5pct=share(.05),
                           share_top_10pct=share(.10), share_top_20pct=share(.20)),
        robustness_drop_top_carriers=robust)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"  Tarone Z={Z:,.1f}  rho_unadj={rho_u:.4f}  rho_adj={rho_a:.4f}  "
          f"slope={b1:+.3f}  Gini={res['concentration']['gini']:.3f}  outside band={outside}/{N}")


if __name__ == "__main__":
    main()
