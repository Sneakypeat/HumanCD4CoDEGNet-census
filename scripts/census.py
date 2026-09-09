import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats, optimize

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

EXPECT = {"edges": 84044, "regulators": 617, "reversing": 4379,
          "reversing_regulators": 324, "same_sign": 79665}
EXPECT_TRANS = {"edges": 83489, "regulators": 617, "reversing": 4379,
                "reversing_regulators": 324, "same_sign": 79110,
                "targets": 9554, "reversing_targets": 2577,
                "self_edges_excluded": 555}
EXPECT_ELIGIBILITY = {"state_specific_gate_passed": 618,
                      "represented_in_inclusive_census": 617,
                      "represented_in_trans_census": 617,
                      "gate_passed_without_a_both_significant_edge": 1}
ZERO_EDGE_GATE_PASSERS = ["CARD6"]


def build_census_master() -> None:
    feasibility = json.loads((DATA / "census/signflip_feasibility.json").read_text())
    assert feasibility["n_candidate_perturbations"] == \
        EXPECT_ELIGIBILITY["state_specific_gate_passed"]
    f = pd.read_csv(DATA / "paper1_flips/frozen_edges.tsv", sep="\t")
    f["regulator"] = f.regulator.astype(str)
    f["target"] = f.target.astype(str)
    assert len(f) == len(f.drop_duplicates(["regulator", "target"])), "duplicate edge in census"

    frozen_ids = pd.read_csv(
        DATA / "paper1_flips/frozen_edges_ensembl.tsv", sep="\t",
        usecols=["regulator", "regulator_ensg", "target", "target_ensg"])
    assert frozen_ids[["regulator", "regulator_ensg", "target", "target_ensg"]].notna().all().all()
    for col in ("regulator", "regulator_ensg", "target", "target_ensg"):
        frozen_ids[col] = frozen_ids[col].astype(str)
    assert len(frozen_ids) == EXPECT["edges"]
    assert len(frozen_ids) == len(frozen_ids.drop_duplicates(["regulator", "target"])), \
        "duplicate edge in Ensembl-annotated census"
    frozen_ids = frozen_ids.rename(columns={
        "regulator_ensg": "frozen_regulator_ensg", "target_ensg": "frozen_target_ensg"})
    f = f.merge(frozen_ids, on=["regulator", "target"], how="left", validate="one_to_one")

    raw = pd.read_csv(DATA / "census/atlas_endpoint_identity.tsv", sep="\t")
    required = {"endpoint_role", "symbol", "ensembl_gene_id"}
    assert required <= set(raw.columns)
    assert not raw.duplicated(["endpoint_role", "symbol"]).any()
    raw_reg = (raw[raw.endpoint_role.eq("regulator")][["symbol", "ensembl_gene_id"]]
               .rename(columns={"symbol": "regulator",
                                "ensembl_gene_id": "regulator_ensg"}))
    raw_tgt = (raw[raw.endpoint_role.eq("target")][["symbol", "ensembl_gene_id"]]
               .rename(columns={"symbol": "target", "ensembl_gene_id": "target_ensg"}))
    f = f.merge(raw_reg, on="regulator", how="left", validate="many_to_one").merge(
        raw_tgt, on="target", how="left", validate="many_to_one")
    assert f[["regulator_ensg", "target_ensg"]].notna().all().all(), \
        "an inclusive census edge lacks a raw-atlas Ensembl identity"
    assert f.regulator_ensg.str.startswith("ENSG").all()
    assert f.target_ensg.str.startswith("ENSG").all()
    f["frozen_is_self"] = f.frozen_regulator_ensg.eq(f.frozen_target_ensg)
    f["is_self"] = f.regulator_ensg.eq(f.target_ensg)
    assert f.is_self.equals(f.frozen_is_self), \
        "raw-atlas and frozen ID maps disagree on self-versus-trans membership"
    reg_id_mismatch = int((f.regulator_ensg != f.frozen_regulator_ensg).sum())
    tgt_id_mismatch = int((f.target_ensg != f.frozen_target_ensg).sum())
    reg_symbol_mismatch = int(f.loc[
        f.regulator_ensg != f.frozen_regulator_ensg, "regulator"].nunique())
    tgt_symbol_mismatch = int(f.loc[
        f.target_ensg != f.frozen_target_ensg, "target"].nunique())

    k = pd.read_csv(DATA / "kinetics/kinetics_8hr_edges.csv.gz",
                    usecols=["regulator", "target", "lfc_Stim8hr", "z_Stim8hr",
                             "padj_Stim8hr", "baseMean_Rest", "baseMean_Stim48hr", "kind"])
    k["regulator"] = k.regulator.astype(str)
    k["target"] = k.target.astype(str)
    k = k.rename(columns={"kind": "kind_8h_table"})

    m = f.merge(k, on=["regulator", "target"], how="left", validate="one_to_one")
    m["has_8hr"] = m.lfc_Stim8hr.notna()

    m["kind"] = np.where(m.opposite_sign, "rev", "stable")
    m["direction"] = np.where(~m.opposite_sign, "",
                              np.where(m.lfc_Rest > 0, "lost_on_activation",
                                       "gained_on_activation"))

    both = m[m.kind_8h_table.notna()]
    disagree = int((both.kind != both.kind_8h_table).sum())

    m = m[["regulator", "regulator_ensg", "target", "target_ensg", "is_self",
           "lfc_Rest", "se_Rest", "adjp_Rest",
           "lfc_Stim8hr", "z_Stim8hr", "padj_Stim8hr",
           "lfc_Stim48", "se_Stim48", "adjp_Stim48",
           "baseMean_Rest", "baseMean_Stim48hr",
           "kind", "direction", "has_8hr"]]

    got = {"edges": len(m), "regulators": m.regulator.nunique(),
           "reversing": int((m.kind == "rev").sum()),
           "reversing_regulators": int(m[m.kind == "rev"].regulator.nunique()),
           "same_sign": int((m.kind == "stable").sum())}
    bad = {k2: (v, EXPECT[k2]) for k2, v in got.items() if v != EXPECT[k2]}
    assert not bad, f"master table does not match published totals: {bad}"
    assert disagree == 0, f"{disagree} edges labelled differently by the 8 h table"

    trans = m[~m.is_self]
    trans_got = {
        "edges": len(trans),
        "regulators": int(trans.regulator.nunique()),
        "reversing": int((trans.kind == "rev").sum()),
        "reversing_regulators": int(trans[trans.kind == "rev"].regulator.nunique()),
        "same_sign": int((trans.kind == "stable").sum()),
        "targets": int(trans.target.nunique()),
        "reversing_targets": int(trans[trans.kind == "rev"].target.nunique()),
        "self_edges_excluded": int(m.is_self.sum()),
    }
    trans_bad = {k2: (v, EXPECT_TRANS[k2]) for k2, v in trans_got.items()
                 if v != EXPECT_TRANS[k2]}
    assert not trans_bad, f"trans estimand does not match authoritative totals: {trans_bad}"
    assert int((m.is_self & m.regulator.ne(m.target)).sum()) == 1
    assert m.loc[m.is_self & m.regulator.ne(m.target),
                 ["regulator", "target"]].values.tolist() == [["QNG1", "C9orf64"]]
    assert not (m.loc[m.is_self, "kind"] == "rev").any(), \
        "a self edge was labelled reversing"

    out = DATA / "census/census_master_edges.csv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    # Figures 2, 3 and 4 pin this file's sha256. gzip records an mtime in its
    # header, so a plain rewrite changes the hash even when the table is
    # byte-identical, and simply re-running this script would break all three.
    # Compare the uncompressed contents and rewrite only on a real change.
    fresh = m.to_csv(index=False).encode()
    current = gzip.decompress(out.read_bytes()) if out.is_file() else None
    if current != fresh:
        with gzip.GzipFile(out, "wb", mtime=0) as handle:
            handle.write(fresh)
        print(f"[census] rewrote {out.name}")
    else:
        print(f"[census] {out.name} unchanged")

    lost = int((m.direction == "lost_on_activation").sum())
    gained = int((m.direction == "gained_on_activation").sum())
    manifest = {
        "file": str(out.relative_to(ROOT)),
        "what": ("every inclusive census edge, all three states; the preserved source "
                 "authority. Manuscript and main figures use trans_totals."),
        "sources": {
            "membership_rest_and_48h": "data/paper1_flips/frozen_edges.tsv",
            "endpoint_ensembl_identity": "data/census/atlas_endpoint_identity.tsv",
            "self_mask_crosscheck": "data/paper1_flips/frozen_edges_ensembl.tsv",
            "8h_statistics": "data/kinetics/kinetics_8hr_edges.csv.gz"},
        "totals": got,
        "trans_totals": trans_got,
        "regulator_scope": {
            **EXPECT_ELIGIBILITY,
            "gate_passed_without_a_both_significant_edge_ids": ZERO_EDGE_GATE_PASSERS,
            "definition": ("the upstream gate required on-target significance and at least "
                           "50 significant downstream targets in Rest and Stim48hr; one of "
                           "the 618 gate-passing perturbations (CARD6) contributed no target edge "
                           "significant in both states, leaving 617 represented regulators"),
        },
        "direction_split": {"lost_on_activation": lost, "gained_on_activation": gained},
        "identity_rule": {
            "self_definition": "regulator_ensg == target_ensg",
            "identity_authority": "raw atlas obs/target_contrast and var/gene_ids",
            "symbol_alias_self_edges": ["QNG1->C9orf64"],
            "symbol_equality_self_edges": int(m.regulator.eq(m.target).sum()),
            "ensembl_equality_self_edges": int(m.is_self.sum()),
            "frozen_self_mask_agrees": True,
            "frozen_historical_id_mismatch_rows": {
                "regulator": reg_id_mismatch, "target": tgt_id_mismatch},
            "frozen_historical_id_mismatch_symbols": {
                "regulator": reg_symbol_mismatch, "target": tgt_symbol_mismatch},
        },
        "targets": int(m.target.nunique()),
        "reversing_targets": int(m[m.kind == "rev"].target.nunique()),
        "inclusive_edges_without_8h": int((~m.has_8hr).sum()),
        "inclusive_reversing_edges_without_8h": int(
            ((~m.has_8hr) & (m.kind == "rev")).sum()),
        "trans_edges_without_8h": int((~trans.has_8hr).sum()),
        "trans_reversing_edges_without_8h": int(
            ((~trans.has_8hr) & (trans.kind == "rev")).sum()),
        "regulators_absent_from_8h_table": sorted(
            set(m.regulator) - set(m[m.has_8hr].regulator)),
        "label_disagreements_with_8h_table": disagree,
        "supersedes": ("kinetics_8hr_edges.csv.gz as a census source; that file is the 8 h "
                       "subset and undercounts reversing regulators 324 -> 323 because KAT7's "
                       "two reversing edges have no finite 8 h value"),
    }
    (DATA / "census/census_master_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[census] built {out.name}")


def load_burden():
    d = pd.read_csv(DATA / "census/per_regulator_enrichment_all.csv")
    assert len(d) == 617, f"expected the 617-regulator census family, got {len(d)}"
    n, k = d.n_both.to_numpy(float), d.n_inv.to_numpy(float)
    assert (k <= n).all() and (n > 0).all()
    return n, k


def fit(n, k, x, with_slope: bool, with_disp: bool):
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
    return -r.fun, r.x


def gini(v):
    v = np.sort(v)
    i = np.arange(1, len(v) + 1)
    return float((2 * (i * v).sum()) / (len(v) * v.sum()) - (len(v) + 1) / len(v))


def burden_heterogeneity():
    n, k = load_burden()
    x, N, p = np.log10(n), len(n), k.sum() / n.sum()

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
        _source="data/census/per_regulator_enrichment_all.csv",
        _script="scripts/census.py",
        _caveat=("Technical covariates were tested in burden_covariate_closure.py and do NOT "
                 "explain the dispersion (rho 0.0661 -> 0.0635 adding precision, efficacy and a "
                 "cubic degree term on the same 555 regulators). That excludes the measurable "
                 "technical explanations, not every one, and does not establish biology."),
        _not_new=("This is a single global effect size for the same heterogeneity the "
                  "per-regulator degree-matched Fisher tests already detect (106/617 at 5% FDR)."),
        n_regulators=N, n_edges=int(n.sum()), n_reversals=int(k.sum()), pooled_rate=float(p),
        tarone_Z=Z, tarone_p=float(stats.norm.sf(Z)),
        betabinom_unadjusted=dict(rho=float(rho_u), lrt_chi2=float(2*(ll_bb-ll_bin)),
                                  lrt_p=float(stats.chi2.sf(2*(ll_bb-ll_bin), 1))),
        degree_adjusted=dict(
            slope_per_log10_outdegree=b1,
            slope_lrt_chi2=float(2*(ll_bbd-ll_bb)), slope_p=float(stats.chi2.sf(2*(ll_bbd-ll_bb), 1)),
            residual_rho=float(rho_a), residual_lrt_chi2=float(2*(ll_bbd-ll_bind)),
            residual_p_naive_chi2_1=float(stats.chi2.sf(2*(ll_bbd-ll_bind), 1)),
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

    out = DATA / "burden_heterogeneity/burden_heterogeneity.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2) + "\n")
    print(f"[census] built {out.name}")


def main() -> int:
    build_census_master()
    burden_heterogeneity()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
