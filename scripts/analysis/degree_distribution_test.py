#!/usr/bin/env python3
"""Is the regulator out-degree distribution heavy-tailed, and is 'power law' defensible?

WHICH NETWORK. Not the drawn ones. Repository figure R5 is the 532-edge high-magnitude core, whose
degrees are set by the |log2FC| > 0.5 threshold. Figure 1d is the reversing subgraph, where a
regulator's degree is HOW MANY OF ITS EDGES REVERSED, i.e. the outcome variable rather than
its connectivity. The object with actual regulator connectivity is the FULL CENSUS edge table
after the paper-wide Ensembl-ID self-edge mask: 83,489 both-significant trans edges over 617
regulators,
of which both drawn networks are subgraphs. The reversing subgraph is fitted too, but only as
a comparison.

WHAT THIS CAN AND CANNOT CLAIM. The 617 census regulators are not all perturbed genes. An
upstream gate retained 618 perturbations with significant on-target knockdown and at least 50
significant downstream targets in EACH state separately; CARD6 contributed no edge significant
in both endpoint states, leaving 617 regulators represented in the shared-edge census. Thus the
distribution is conditioned on signal-based selection, and a conclusion here is about the shape
over the observed census range, never about the T-cell regulatory network being scale-free.

METHOD (Clauset, Shalizi & Newman 2009), because a straight line on a log-log plot is not
evidence and most "scale-free" claims in biology do not survive this procedure:

  1. discrete power-law MLE for alpha at each candidate xmin, p(k) proportional to
     k^-alpha / zeta(alpha, xmin)
  2. xmin chosen to minimise the KS distance between the empirical and fitted CDFs
  3. GOODNESS OF FIT by semi-parametric bootstrap: synthesise datasets that are power law
     above xmin and resampled-empirical below it, refit each end to end, and take
     p = fraction of synthetic KS >= observed KS. p < 0.1 REJECTS the power law. A large p
     means "not ruled out", never "confirmed".
  4. RELATIVE FIT on the same selected tail against a discretised, tail-conditioned
     lognormal; a shifted geometric (the discrete exponential); and a power law with an
     exponential cutoff. Vuong tests are used only for the two non-nested comparisons.
     The cutoff model contains the pure power law at a boundary, so it is compared
     descriptively by AICc rather than passed through an invalid ordinary Vuong test.

Reads   data/census/census_master_edges.csv.gz  (validated master; 83,489 trans edges)
Writes  data/census/degree_distribution_test.json
        figures/FigS_degree_distribution.{pdf,png}
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp, zeta
from scipy.stats import lognorm, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S

PAPER = Path(__file__).resolve().parents[2]
RD = PAPER / "data"
OUT = PAPER / "figures_final"

SEED = 20260813
N_BOOT = 5000          # synthetic datasets for the goodness-of-fit p
ALPHA_BOUNDS = (1.01, 15.0)
ANALYSIS_VERSION = 4       # v4 reports package-compatible raw R separately from Vuong z


# ----------------------------------------------------------- discrete power law

def _nll_powerlaw(alpha: float, x: np.ndarray, xmin: int) -> float:
    if not (ALPHA_BOUNDS[0] < alpha < ALPHA_BOUNDS[1]):
        return np.inf
    return len(x) * np.log(zeta(alpha, xmin)) + alpha * np.sum(np.log(x))


def fit_alpha(x: np.ndarray, xmin: int) -> float:
    r = minimize_scalar(_nll_powerlaw, bounds=ALPHA_BOUNDS, args=(x, xmin), method="bounded")
    return float(r.x)


def pl_cdf(k: np.ndarray, alpha: float, xmin: int) -> np.ndarray:
    """P(X <= k) for the discrete power law, via the Hurwitz zeta tail."""
    return 1.0 - zeta(alpha, k + 1) / zeta(alpha, xmin)


def ks_distance(x: np.ndarray, alpha: float, xmin: int) -> float:
    """Exact discrete KS distance, evaluated before and after every tied empirical jump."""
    vals, counts = np.unique(x, return_counts=True)
    empirical_after = np.cumsum(counts) / len(x)
    empirical_before = (np.cumsum(counts) - counts) / len(x)
    fitted_after = pl_cdf(vals, alpha, xmin)
    fitted_before = pl_cdf(vals - 1, alpha, xmin)
    return float(max(np.max(np.abs(empirical_after - fitted_after)),
                     np.max(np.abs(empirical_before - fitted_before))))


def choose_xmin(data: np.ndarray) -> tuple[int, float, float, int]:
    cand = np.unique(data)
    cand = cand[cand >= 1]
    cand = cand[:-1]
    best = (None, None, np.inf, 0)
    for xm in cand:
        tail = data[data >= xm]
        if len(tail) < 50:
            continue
        a = fit_alpha(tail, int(xm))
        d = ks_distance(tail, a, int(xm))
        if d < best[2]:
            best = (int(xm), a, d, len(tail))
    return best


def sample_powerlaw(n: int, alpha: float, xmin: int, rng) -> np.ndarray:
    """Exact unbounded discrete-power-law rejection sampler.

    The proposal rounds a continuous Pareto variate with lower bound xmin - 1/2. Convexity
    supplies a finite rejection bound, avoiding the silent finite-tail truncation produced by
    sampling from an arbitrarily capped probability table.
    """
    if n == 0:
        return np.empty(0, dtype=np.int64)
    lower = xmin - 0.5
    normalizer = zeta(alpha, xmin)
    bound = 1.0 / (normalizer * (alpha - 1.0) * lower ** (alpha - 1.0))
    chunks, have = [], 0
    while have < n:
        m = max(64, int((n - have) * 1.05) + 8)
        u = rng.random(m)
        y = lower * np.exp(-np.log1p(-u) / (alpha - 1.0))
        k = np.floor(y + 0.5).astype(np.int64)
        lo = k - 0.5
        q = (lower ** (alpha - 1.0) * lo ** (1.0 - alpha)
             * (-np.expm1((1.0 - alpha) * np.log((k + 0.5) / lo))))
        p = k.astype(float) ** (-alpha) / normalizer
        keep = rng.random(m) < np.minimum(1.0, p / (bound * q))
        accepted = k[keep]
        chunks.append(accepted)
        have += len(accepted)
    return np.concatenate(chunks)[:n]


def gof_pvalue(data: np.ndarray, alpha: float, xmin: int, d_obs: float, rng,
               n_boot: int = N_BOOT) -> tuple[float, int, int]:
    """Semi-parametric bootstrap. Each synthetic set is refitted END TO END, xmin included,
    which is what makes the p-value account for xmin having been chosen from the data."""
    below = data[data < xmin]
    n, n_tail = len(data), int((data >= xmin).sum())
    p_tail = n_tail / n
    worse, valid = 0, 0
    for _ in range(n_boot):
        n_t = rng.binomial(n, p_tail)
        tail = sample_powerlaw(n_t, alpha, xmin, rng)
        rest = rng.choice(below, size=n - n_t, replace=True) if len(below) else np.array([], int)
        syn = np.concatenate([tail, rest]).astype(int)
        xm_s, a_s, d_s, _ = choose_xmin(syn)
        if xm_s is None:
            continue
        valid += 1
        worse += int(d_s >= d_obs)
    # Finite-Monte-Carlo correction: a simulated p-value is never reported as exactly zero.
    return (worse + 1) / (valid + 1), worse, valid


# ------------------------------------------------------- alternative tail models

def loglik_powerlaw(x, alpha, xmin):
    return -alpha * np.log(x) - np.log(zeta(alpha, xmin))


def loglik_lognormal(x, mu, sig, xmin):
    """Integer-bin lognormal mass, conditioned on K >= xmin."""
    upper = lognorm.cdf(x + 0.5, s=sig, scale=np.exp(mu))
    lower = lognorm.cdf(np.maximum(x - 0.5, 1e-12), s=sig, scale=np.exp(mu))
    norm_tail = lognorm.sf(xmin - 0.5, s=sig, scale=np.exp(mu))
    return np.log(np.maximum((upper - lower) / max(norm_tail, 1e-300), 1e-300))


def loglik_exponential(x, q, xmin):
    """Discrete exponential, i.e. geometric mass on K-xmin = 0,1,..."""
    return np.log1p(-q) + (x - xmin) * np.log(q)


def _cutoff_log_normalizer(alpha: float, lam: float, xmin: int, max_observed: int) -> float:
    """Log sum of k^-alpha exp(-lambda*k), k=xmin..infinity, to numerical precision."""
    # exp(-lambda*k) makes a finite numerical grid safe. Thirty-five exponential e-folds
    # beyond xmin leaves < 1e-15 of an already upper-bounded geometric remainder. The
    # observed maximum is always retained, and the cap is merely a guard against pathological
    # optimiser proposals near lambda=0 (the pure PL handles that boundary separately).
    upper = max(max_observed * 2, xmin + 100, int(xmin + 35.0 / max(lam, 1e-9)))
    upper = min(upper, 5_000_000)
    support = np.arange(xmin, upper + 1, dtype=float)
    return float(logsumexp(-alpha * np.log(support) - lam * support))


def loglik_cutoff_powerlaw(x, alpha, lam, xmin):
    log_z = _cutoff_log_normalizer(alpha, lam, xmin, int(np.max(x)))
    return -alpha * np.log(x) - lam * x - log_z


def fit_cutoff_powerlaw(x: np.ndarray, xmin: int) -> tuple[float, float, np.ndarray]:
    """Fit exact discrete p(k) proportional to k^-alpha exp(-lambda*k)."""
    x = np.asarray(x, float)

    def nll(theta):
        alpha, lam = np.exp(theta)
        return -float(loglik_cutoff_powerlaw(x, alpha, lam, xmin).sum())

    starts = [(1.5, 1e-3), (0.05, 0.1), (2.0, 1e-2)]
    fits = [minimize(nll, np.log(start), method="Nelder-Mead",
                     options={"maxiter": 5000, "xatol": 1e-10, "fatol": 1e-9})
            for start in starts]
    opt = min(fits, key=lambda z: z.fun)
    if not opt.success:
        raise RuntimeError(f"cutoff power-law fit failed: {opt.message}")
    alpha, lam = np.exp(opt.x)
    ll = loglik_cutoff_powerlaw(x, float(alpha), float(lam), xmin)
    return float(alpha), float(lam), ll


def aicc(log_likelihood: float, n: int, n_parameters: int) -> float:
    return float(-2 * log_likelihood + 2 * n_parameters
                 + 2 * n_parameters * (n_parameters + 1) / (n - n_parameters - 1))


def fit_tail_models(x: np.ndarray, alpha: float, xmin: int) -> dict:
    """Fit all alternatives on the same integer support and tail conditioning as the PL."""
    x = np.asarray(x, float)
    lx = np.log(x)

    def nll_ln(theta):
        mu, log_sig = theta
        return -float(loglik_lognormal(x, mu, np.exp(log_sig), xmin).sum())

    opt = minimize(nll_ln, [float(lx.mean()), np.log(max(float(lx.std(ddof=1)), 0.1))],
                   method="Nelder-Mead", options={"maxiter": 5000, "xatol": 1e-10,
                                                   "fatol": 1e-10})
    if not opt.success:
        raise RuntimeError(f"discrete lognormal fit failed: {opt.message}")
    mu, sig = float(opt.x[0]), float(np.exp(opt.x[1]))

    mean_excess = float((x - xmin).mean())
    q = mean_excess / (1.0 + mean_excess)

    ll_pl = loglik_powerlaw(x, alpha, xmin)
    ll_ln = loglik_lognormal(x, mu, sig, xmin)
    ll_ex = loglik_exponential(x, q, xmin)
    cutoff_alpha, cutoff_lambda, ll_cutoff = fit_cutoff_powerlaw(x, xmin)
    log_likelihoods = {
        "power_law": float(ll_pl.sum()),
        "lognormal": float(ll_ln.sum()),
        "exponential": float(ll_ex.sum()),
        "cutoff_power_law": float(ll_cutoff.sum()),
    }
    aicc_values = {
        "power_law": aicc(log_likelihoods["power_law"], len(x), 1),
        "lognormal": aicc(log_likelihoods["lognormal"], len(x), 2),
        "exponential": aicc(log_likelihoods["exponential"], len(x), 1),
        "cutoff_power_law": aicc(log_likelihoods["cutoff_power_law"], len(x), 2),
    }
    best_aicc = min(aicc_values.values())
    out = {
        "lognormal_params": {"mu": mu, "sigma": sig},
        "exponential_params": {"q": q, "lambda": -float(np.log(q))},
        "cutoff_power_law_params": {"alpha": cutoff_alpha, "lambda": cutoff_lambda,
                                    "alpha_at_lower_bound": bool(cutoff_alpha < 1e-5)},
        "log_likelihood": log_likelihoods,
        "aicc": aicc_values,
        "delta_aicc": {k: float(v - best_aicc) for k, v in aicc_values.items()},
    }
    for name, ll in [("lognormal", ll_ln), ("exponential", ll_ex)]:
        raw_R, z, p = vuong(ll_pl, ll)
        out[f"power_law_vs_{name}"] = {
            # powerlaw.Fit.distribution_compare() returns this raw summed log-likelihood
            # ratio by default. The normalised statistic is reported separately because it
            # is the quantity used to calculate the non-nested two-sided p-value.
            "R": float(raw_R), "vuong_z": float(z), "p": float(p),
            "favours": ("power law" if raw_R > 0 else name)
            if p < 0.05 else "indistinguishable"}
    out["power_law_vs_cutoff_power_law"] = {
        "delta_log_likelihood": float(ll_pl.sum() - ll_cutoff.sum()),
        "test": ("not tested: cutoff power law contains the pure power law at lambda=0; "
                 "relative support is reported by AICc"),
    }
    return out


def add_global_vuong_fdr(results: dict) -> None:
    """BH-adjust the six non-nested PL-vs-lognormal/geometric comparisons together."""
    tests = []
    keys = ["regulator_outdegree_full_census", "target_indegree_full_census",
            "regulator_outdegree_reversing_subgraph"]
    for result_key in keys:
        for comparison in ("power_law_vs_lognormal", "power_law_vs_exponential"):
            item = results[result_key]["model_comparison"][comparison]
            tests.append((result_key, comparison, float(item["p"])))
    order = np.argsort([x[2] for x in tests])
    adjusted = np.empty(len(tests), float)
    running = 1.0
    for rank_index in range(len(tests) - 1, -1, -1):
        idx = int(order[rank_index])
        rank = rank_index + 1
        running = min(running, tests[idx][2] * len(tests) / rank)
        adjusted[idx] = running
    for (result_key, comparison, _), q in zip(tests, adjusted):
        results[result_key]["model_comparison"][comparison]["q_bh_across_six"] = float(q)


def upgrade_comparison_schema(results: dict) -> None:
    """Upgrade cached v3 comparisons that mislabeled the normalised z statistic as raw R."""
    keys = ["regulator_outdegree_full_census", "target_indegree_full_census",
            "regulator_outdegree_reversing_subgraph"]
    for result_key in keys:
        result = results[result_key]
        comparison = result["model_comparison"]
        n = int(result["n_tail"])
        ll_pl = float(comparison["log_likelihood"]["power_law"])
        for alternative in ("lognormal", "exponential"):
            item = comparison[f"power_law_vs_{alternative}"]
            if "vuong_z" in item:
                continue
            old_sample_sd_z = float(item["R"])
            raw_R = ll_pl - float(comparison["log_likelihood"][alternative])
            z = old_sample_sd_z * np.sqrt(n / (n - 1)) if n > 1 else 0.0
            p = float(2 * norm.sf(abs(z)))
            item.update({
                "R": float(raw_R),
                "vuong_z": float(z),
                "p": p,
                "favours": (("power law" if raw_R > 0 else alternative)
                             if p < 0.05 else "indistinguishable"),
            })
    results["analysis_version"] = ANALYSIS_VERSION


def model_survival(support: np.ndarray, name: str, params: dict, xmin: int) -> np.ndarray:
    """Conditional P(K >= k | K >= xmin) for each fitted discrete-tail model."""
    support = np.asarray(support, int)
    if name == "power_law":
        return zeta(params["alpha"], support) / zeta(params["alpha"], xmin)
    if name == "lognormal":
        return (lognorm.sf(support - 0.5, s=params["sigma"], scale=np.exp(params["mu"])) /
                lognorm.sf(xmin - 0.5, s=params["sigma"], scale=np.exp(params["mu"])))
    if name == "exponential":
        return params["q"] ** (support - xmin)
    if name == "cutoff_power_law":
        alpha, lam = params["alpha"], params["lambda"]
        upper = max(int(support.max()) * 2, xmin + 100,
                    int(xmin + 35.0 / max(lam, 1e-9)))
        upper = min(upper, 5_000_000)
        grid = np.arange(xmin, upper + 1, dtype=float)
        mass = np.exp(-alpha * np.log(grid) - lam * grid
                      - logsumexp(-alpha * np.log(grid) - lam * grid))
        survival = np.cumsum(mass[::-1])[::-1]
        return survival[support - xmin]
    raise ValueError(name)


def vuong(ll1: np.ndarray, ll2: np.ndarray) -> tuple[float, float, float]:
    """Raw log-likelihood ratio, normalised Vuong z, and two-sided p.

    Raw R > 0 favours model 1. This matches the default return value of the Python ``powerlaw``
    package. The p-value is calculated from z = R / sqrt(n * Var[d]); a large p means the two
    models are not distinguishable on this data, which is a real and common answer.
    """
    d = ll1 - ll2
    n = len(d)
    raw_R = float(d.sum())
    # Match powerlaw 1.5/Vuong's population variance convention for the pointwise
    # log-likelihood differences. This changes only the normalising statistic, not raw R.
    variance = float(np.mean((d - d.mean()) ** 2))
    if variance == 0:
        return raw_R, 0.0, 1.0
    z = raw_R / np.sqrt(n * variance)
    return raw_R, float(z), float(2 * norm.sf(abs(z)))


def compare_tail(x: np.ndarray, alpha: float, xmin: int) -> dict:
    return fit_tail_models(x, alpha, xmin)


def analyse(data: np.ndarray, label: str, rng, do_boot: bool = True) -> dict:
    xmin, alpha, d, n_tail = choose_xmin(data)
    at_bound = abs(alpha - ALPHA_BOUNDS[1]) < 1e-3 or abs(alpha - ALPHA_BOUNDS[0]) < 1e-3
    res = {"alpha_at_optimiser_bound": bool(at_bound), "label": label, "n": int(len(data)), "min": int(data.min()),
           "max": int(data.max()), "median": float(np.median(data)),
           "xmin": xmin, "alpha": float(alpha), "n_tail": n_tail,
           "ks_distance": float(d),
           "tail_fraction": round(n_tail / len(data), 4)}
    if do_boot:
        p, exceedances, valid_bootstraps = gof_pvalue(data, alpha, xmin, d, rng)
        bootstrap_mode = "serial; shared deterministic RNG stream"
        res["gof_p"] = round(p, 4)
        res["gof_exceedances"] = int(exceedances)
        res["gof_valid_bootstraps"] = int(valid_bootstraps)
        res["gof_execution"] = bootstrap_mode
        res["power_law_verdict"] = ("REJECTED (p < 0.1)" if p < 0.1
                                    else "not ruled out (this is NOT confirmation)")
    res["model_comparison"] = compare_tail(data[data >= xmin], alpha, xmin)

    # COMBINED VERDICT. Absolute goodness of fit and relative model comparison answer
    # different questions and can disagree, so both are reported.
    if do_boot:
        if res["gof_p"] < 0.1:
            res["verdict"] = "power law REJECTED by goodness of fit"
        else:
            res["verdict"] = "power-law tail compatible by GOF; this is not confirmation"
    return res


def ccdf(x: np.ndarray):
    """Tie-aware empirical survival P(X >= k) at each observed integer value."""
    xs, counts = np.unique(np.asarray(x, int), return_counts=True)
    return xs, np.cumsum(counts[::-1])[::-1] / counts.sum()


def main(figure_only: bool = False) -> None:
    """figure_only reuses data/census/degree_distribution_test.json so that layout
    edits do not re-run the 5,000-draw goodness-of-fit bootstrap."""
    rng = np.random.default_rng(SEED)
    # Everything comes from ONE file now: census_master_edges.csv.gz, built and validated by
    # build_census_master.py against all five published totals (84,044 edges, 617 regulators,
    # 4,379 reversals over 324 regulators, 79,665 same-sign, 2,298/2,081 lost/gained). An
    # earlier version of this figure read a stale kinetics_8hr_edges.csv.gz subset that omitted
    # MCAT and reported 322 reversing regulators. The corrected 8 h subset reports 323 because
    # KAT7 has no finite 8 h value; neither subset is admissible for the Rest--Stim48 topology.
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "direction", "is_self"])
    e = e[~e.is_self].copy()
    assert len(e) == 83_489 and int((e.kind == "rev").sum()) == 4_379
    out_all = e.groupby("regulator").size().to_numpy()
    in_all = e.groupby("target").size().to_numpy()
    rev = e[e.kind == "rev"]
    out_rev = rev.groupby("regulator").size().to_numpy()
    N_REV_REG = int(rev.regulator.nunique())            # 324
    no_outdeg = []                                      # nothing is missing any more

    cache = RD / "census/degree_distribution_test.json"
    if figure_only and cache.exists():
        candidate = json.loads(cache.read_text())
        if (candidate.get("analysis_version") == ANALYSIS_VERSION
                and candidate.get("n_bootstrap") == N_BOOT):
            res = candidate
            print("reusing cached fits (figure-only mode)")
        else:
            res = None
            print("cache predates the corrected discrete-tail analysis; recomputing")
    else:
        res = None
    if res is None:
        parts = {
            "regulator_outdegree_full_census": RD / "census/degree_distribution_test.regulator.json",
            "target_indegree_full_census": RD / "census/degree_distribution_test.target.json",
            "regulator_outdegree_reversing_subgraph": RD / "census/degree_distribution_test.reversal.json",
        }
        if figure_only and all(path.exists() for path in parts.values()):
            res = {"analysis_version": ANALYSIS_VERSION, "seed": SEED,
                   "n_bootstrap": N_BOOT,
                   "_source": ("data/census/census_master_edges.csv.gz, is_self=false "
                               "by endpoint Ensembl identity (83,489 trans census edges)"),
                   "_scope": ("618 perturbations passed the upstream on-target and >=50 "
                              "significant-target gate in each endpoint state; CARD6 contributed "
                              "no shared significant edge, leaving 617 represented census "
                              "regulators. This signal-selected census is not a sample of all "
                              "perturbed genes, so it cannot establish a whole-network "
                              "scale-free claim")}
            res.update({key: json.loads(path.read_text()) for key, path in parts.items()})
            print("assembled corrected single-family bootstrap results")
        else:
            res = {"analysis_version": ANALYSIS_VERSION, "seed": SEED,
                   "n_bootstrap": N_BOOT,
                   "_source": ("data/census/census_master_edges.csv.gz, is_self=false "
                               "by endpoint Ensembl identity (83,489 trans census edges)"),
                   "_scope": ("618 perturbations passed the upstream on-target and >=50 "
                              "significant-target gate in each endpoint state; CARD6 contributed "
                              "no shared significant edge, leaving 617 represented census "
                              "regulators. This signal-selected census is not a sample of all "
                              "perturbed genes, so it cannot establish a whole-network "
                              "scale-free claim"),
                   "regulator_outdegree_full_census": analyse(
                       out_all, "regulator out-degree", rng),
                   "target_indegree_full_census": analyse(
                       in_all, "target in-degree", rng),
                   "regulator_outdegree_reversing_subgraph": analyse(
                       out_rev, "reversals per regulator", rng)}
    upgrade_comparison_schema(res)
    add_global_vuong_fdr(res)
    # ------------------------------------------------------------------ figure
    S.setup()
    # Two plots above a full-width result strip remain legible when embedded at 6.5 inches.
    fig = plt.figure(figsize=(8.2, 7.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.00, 0.82],
                          left=0.085, right=0.985, bottom=0.125, top=0.94,
                          wspace=0.24, hspace=0.34)

    # ---- panel a: which regulators reverse at all, over the same degree axis ----
    # A regulator's out-degree is fixed by the census, so it cannot be split by direction the
    # way targets can. What CAN be shown is which regulators carry any reversal at all.
    ax = fig.add_subplot(gs[0, 0])
    r = res["regulator_outdegree_full_census"]
    rev_regs = set(rev.regulator.unique())     # canonical 324
    od = e.groupby("regulator").size()
    od_rev = od[od.index.isin(rev_regs)].to_numpy()
    od_non = od[~od.index.isin(rev_regs)].to_numpy()
    assert len(od_rev) == N_REV_REG, f"{len(od_rev)} != canonical {N_REV_REG}"
    for dat, col, lab in [(od_non, "#c7ccd1", f"no reversal ({len(od_non)} regulators)"),
                          (od_rev, "#33475b",
                           f"carries $\\geq$1 reversal ({len(od_rev)} regulators)")]:
        xs, cc = ccdf(dat)
        ax.loglog(xs, cc, "o", ms=2.6, color=col, alpha=0.75, label=lab)
    # The reversal-bearing split is descriptive: higher-degree regulators mechanically have
    # more opportunities to carry at least one reversal, so no rank-test p-value is reported.
    ax.set_xlabel("regulator out-degree $k$"); ax.set_ylabel("P(K $\\geq$ k)")
    ax.set_title("Regulator out-degree", fontsize=9.2)
    ax.legend(fontsize=7.0, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, -0.02))
    ax.text(0.02, 0.98,
            f"median $k$: {np.median(od_rev):.0f} (\u22651 reversal) vs "
            f"{np.median(od_non):.0f} (none)",
            transform=ax.transAxes, fontsize=7.4, va="top", color="#3d4348")
    S.panel_tag(ax, "a")

    # ---- panel b: targets split by the DIRECTION of the reversals they receive ----
    ax = fig.add_subplot(gs[0, 1])
    # labels straight from the master file, so the split is the SAME 2,298 / 2,081 shown in
    # Figure 1a rather than a re-derivation that could drift from it
    lost = rev[rev.direction == "lost_on_activation"].groupby("target").size()
    gain = rev[rev.direction == "gained_on_activation"].groupby("target").size()
    xs, cc = ccdf(in_all)
    ax.loglog(xs, cc, "o", ms=2.4, color="#c7ccd1", alpha=0.7,
              label=f"all census edges ({len(in_all):,} targets)")
    # the legend carries the edge counts, which are the SAME 2,298 / 2,081 as Figure 1a, so
    # the match is checkable at a glance without lengthening the title
    for dat, col, lab in [(lost.to_numpy(), S.PAL["LOST"],
                           f"effect LOST on stimulation ({int(lost.sum()):,} edges "
                           f"over {len(lost):,} targets)"),
                          (gain.to_numpy(), S.PAL["GAINED"],
                           f"effect GAINED on stimulation ({int(gain.sum()):,} edges "
                           f"over {len(gain):,} targets)")]:
        xs2, cc2 = ccdf(dat)
        ax.loglog(xs2, cc2, "o", ms=2.8, color=col, alpha=0.8, label=lab)
    both = len(set(lost.index) & set(gain.index))
    ax.set_xlabel("target in-degree $k$"); ax.set_ylabel("P(K $\\geq$ k)")
    ax.set_title("Target in-degree by reversal direction", fontsize=9.2)
    ax.legend(fontsize=6.9, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, -0.02))
    ax.text(0.02, 0.99,
            f"both directions: {both:,}/"
            f"{len(set(lost.index) | set(gain.index)):,} targets "
            f"({both / len(set(lost.index) | set(gain.index)):.1%})",
            transform=ax.transAxes, fontsize=7.4, va="top", color="#3d4348")
    S.panel_tag(ax, "b")

    res["direction_split"] = {
        "n_targets_lost": int(len(lost)), "n_targets_gained": int(len(gain)),
        "n_targets_both": int(both),
        "median_indegree_lost": float(np.median(lost)),
        "median_indegree_gained": float(np.median(gain)),
        "regulator_median_outdeg_reversing": float(np.median(od_rev)),
        "regulator_median_outdeg_nonreversing": float(np.median(od_non)),
        "_caveat": ("the regulator split is guaranteed by sampling: a regulator with more "
                    "edges is more likely to carry at least one reversal, so panel a is "
                    "descriptive only"),
    }

    # ---- panel c: empirical tails and all four fitted discrete candidate models ----
    # This is a graphical fit diagnostic, not a correlation plot. Ordered empirical and
    # theoretical quantiles can correlate near 1 even when a formal GOF test rejects a model.
    carrier = fig.add_subplot(gs[1, :]); carrier.axis("off")
    S.panel_tag(carrier, "c")
    inner = gs[1, :].subgridspec(1, 3, wspace=0.25)
    entries = [("regulator_outdegree_full_census", "regulator out-degree"),
               ("target_indegree_full_census", "target in-degree"),
               ("regulator_outdegree_reversing_subgraph", "reversals per regulator")]
    raw_data = [out_all, in_all, out_rev]
    styles = {
        "power_law": ("#7b3294", "-", "power law"),
        "lognormal": ("#2c7fb8", "--", "lognormal"),
        "exponential": ("#d99000", ":", "geometric"),
        "cutoff_power_law": ("#31966f", "-.", "power law + cutoff"),
    }
    shared_tail_ymin = max(0.4 / max(res[key]["n_tail"] for key, _ in entries), 1e-3)
    for j, (key, nice) in enumerate(entries):
        r = res[key]
        dat = np.asarray(raw_data[j])
        tail_values = dat[dat >= r["xmin"]]
        support = np.arange(r["xmin"], int(tail_values.max()) + 1)
        _, counts = np.unique(tail_values, return_counts=True)
        empirical = np.cumsum(counts[::-1])[::-1] / len(tail_values)
        empirical_x = np.unique(tail_values)
        ax = fig.add_subplot(inner[0, j])
        ax.step(empirical_x, empirical, where="post", color="#222222", lw=1.35,
                label="empirical")
        comp = r["model_comparison"]
        params = {
            "power_law": {"alpha": r["alpha"]},
            "lognormal": comp["lognormal_params"],
            "exponential": comp["exponential_params"],
            "cutoff_power_law": comp["cutoff_power_law_params"],
        }
        for model, (colour, line, label) in styles.items():
            fitted = model_survival(support, model, params[model], r["xmin"])
            ax.plot(support, fitted, color=colour, ls=line, lw=1.2, label=label)
        ax.set_xscale("log"); ax.set_yscale("log")
        tick_values = np.unique(np.rint(np.geomspace(support.min(), support.max(), 3)).astype(int))
        ax.xaxis.set_major_locator(FixedLocator(tick_values))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{int(value):,}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(shared_tail_ymin, 1.08)
        ax.set_xlabel("tail degree $k$" if j < 2 else "reversals $k$")
        if j == 0:
            ax.set_ylabel(r"$P(K\geq k\mid K\geq k_{min})$")
        else:
            ax.tick_params(labelleft=False)
        gof = r["gof_p"]
        verdict = "compatible" if gof >= 0.1 else "rejected"
        ax.set_title(nice, fontsize=8.0, pad=4)
        ax.text(0.03, 0.04,
                f"$k_{{min}}$={r['xmin']}; $n_{{tail}}$={r['n_tail']}\n"
                f"PL GOF $p$={gof:.3f} ({verdict})",
                transform=ax.transAxes, fontsize=6.7, va="bottom", color="#30363b")
    handles, labels = ax.get_legend_handles_labels()
    carrier.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.31),
                   ncol=5, fontsize=6.8, frameon=False, handlelength=2.4,
                   columnspacing=1.2)

    # Include the exact descriptive split behind panels a/b in the machine-readable result.
    cache.write_text(json.dumps(res, indent=2) + "\n")
    part_paths = {
        "regulator_outdegree_full_census": RD / "census/degree_distribution_test.regulator.json",
        "target_indegree_full_census": RD / "census/degree_distribution_test.target.json",
        "regulator_outdegree_reversing_subgraph": RD / "census/degree_distribution_test.reversal.json",
    }
    for result_key, path in part_paths.items():
        path.write_text(json.dumps(res[result_key], indent=2) + "\n")

    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"FigS_degree_distribution.{ext}", dpi=300, facecolor="white")

    for key in res:
        if not isinstance(res[key], dict) or "alpha" not in res[key]:
            continue
        r = res[key]
        print(f"{r['label']:28s} n={r['n']:5d}  kmin={r['xmin']:4d}  alpha={r['alpha']:.3f}  "
              f"KS={r['ks_distance']:.4f}  GOF p={r.get('gof_p', float('nan')):.3f}\n"
              f"      -> {r.get('verdict', '')}")
        for k, v in r["model_comparison"].items():
            if not k.startswith("power_law_vs_"):
                continue
            if "R" in v:
                print(f"      {k:32s} raw R={v['R']:+.3f} z={v['vuong_z']:+.3f} "
                      f"p={v['p']:.3g}  {v['favours']}")
            else:
                print(f"      {k:32s} deltaLL={v['delta_log_likelihood']:+.3f}  "
                      f"{v['test']}")


def run_single_analysis(which: str) -> None:
    """Run one fit/bootstrap family; useful for bounded-memory sequential orchestration."""
    e = pd.read_csv(RD / "census/census_master_edges.csv.gz",
                    usecols=["regulator", "target", "kind", "is_self"])
    e = e[~e.is_self]
    specs = {
        "regulator": (e.groupby("regulator").size().to_numpy(), "regulator out-degree", 0),
        "target": (e.groupby("target").size().to_numpy(), "target in-degree", 1),
        "reversal": (e[e.kind == "rev"].groupby("regulator").size().to_numpy(),
                     "reversals per regulator", 2),
    }
    if which not in specs:
        raise ValueError(f"unknown analysis {which!r}; choose one of {sorted(specs)}")
    data, label, seed_index = specs[which]
    # Advance through earlier families so a standalone family uses the same segment of the
    # deterministic RNG stream as a full run.
    rng = np.random.default_rng(SEED)
    ordered = [specs["regulator"], specs["target"], specs["reversal"]]
    result = None
    for j in range(seed_index + 1):
        dat_j, label_j, _ = ordered[j]
        result = analyse(dat_j, label_j, rng)
    assert result is not None and label_j == label
    out = RD / "census" / f"degree_distribution_test.{which}.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{label}: kmin={result['xmin']}, alpha={result['alpha']}, "
          f"KS={result['ks_distance']}, GOF p={result['gof_p']} -> {out}")


if __name__ == "__main__":
    if "--single" in sys.argv:
        run_single_analysis(sys.argv[sys.argv.index("--single") + 1])
    else:
        main(figure_only="--figure-only" in sys.argv)
