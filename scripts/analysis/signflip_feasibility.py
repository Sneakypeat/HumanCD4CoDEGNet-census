"""GO/NO-GO: do trans-regulatory edges INVERT SIGN between Rest and Stim48hr,
beyond what measurement noise can manufacture?

Conservative definition of an inverting edge (perturbation p, target gene g):
  1. significant in BOTH states           adj_p_R < 0.10  AND  adj_p_S < 0.10
  2. OPPOSITE sign of log_fc
  3. the two effects genuinely DIFFER     z = (lfc_R - lfc_S)/sqrt(seR^2 + seS^2), BH-FDR < 0.05
Rule 3 is the one that kills threshold-wobble artefacts: it uses the atlas's own standard errors.
"""
import numpy as np, h5py, fsspec, time, json, os
from scipy import stats
rng = np.random.default_rng(20260715)
URL = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad"
OUT = str(Path(__file__).resolve().parents[2]); os.makedirs(OUT, exist_ok=True)

def open_h5(n=5):
    for a in range(n):
        try: return h5py.File(fsspec.open(URL, block_size=4*1024*1024).open(), "r")
        except Exception as e:
            if a == n-1: raise
            print(f"  retry {a+1} ({type(e).__name__})", flush=True); time.sleep(6*(a+1))
def dec(a): return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
def obs_col(o, k):
    x = o[k]
    if isinstance(x, h5py.Group): return dec(x["categories"][:])[x["codes"][:]]
    return x[:]

h = open_h5(); o = h["obs"]
cond = obs_col(o, "culture_condition")
targ = obs_col(o, "target_contrast_gene_name").astype(str)
onsig = obs_col(o, "ontarget_significant").astype(bool)
ndown = np.nan_to_num(obs_col(o, "n_downstream").astype(float))
genes = dec(h["var"]["gene_name"][:]) if not isinstance(h["var"]["gene_name"], h5py.Group) else None
if genes is None: genes = obs_col(h["var"], "gene_name")
print(f"obs rows {len(cond)} | genes {len(genes)}", flush=True)

# --- candidate perturbations: KD-validated in BOTH states, with real signal in BOTH ---
MIN_OUT = 50
idxR = {t: i for i, (t, c, s, n) in enumerate(zip(targ, cond, onsig, ndown)) if c == "Rest"     and s and n >= MIN_OUT}
idxS = {t: i for i, (t, c, s, n) in enumerate(zip(targ, cond, onsig, ndown)) if c == "Stim48hr" and s and n >= MIN_OUT}
shared = sorted(set(idxR) & set(idxS))
print(f"candidate perturbations (onsig + n_downstream>={MIN_OUT} in BOTH Rest & Stim48hr): {len(shared)}", flush=True)

rows_R = np.array([idxR[t] for t in shared]); rows_S = np.array([idxS[t] for t in shared])
LFC, SE, ADJ = h["layers"]["log_fc"], h["layers"]["lfcSE"], h["layers"]["adj_p_value"]

tot_sig_both = 0; tot_same = 0; tot_opp = 0
flip_p = []; flip_meta = []          # meta: (pert, gene, lfcR, lfcS, seR, seS)
CH = 150
for s0 in range(0, len(shared), CH):
    sl = slice(s0, min(s0+CH, len(shared)))
    rR = rows_R[sl]; rS = rows_S[sl]; names = shared[sl]
    oR = np.argsort(rR); oS = np.argsort(rS)                       # h5py needs increasing indices
    lfcR = LFC[np.sort(rR), :][np.argsort(oR)]; lfcS = LFC[np.sort(rS), :][np.argsort(oS)]
    seR  = SE [np.sort(rR), :][np.argsort(oR)]; seS  = SE [np.sort(rS), :][np.argsort(oS)]
    aR   = ADJ[np.sort(rR), :][np.argsort(oR)]; aS   = ADJ[np.sort(rS), :][np.argsort(oS)]

    sig_both = (aR < 0.10) & (aS < 0.10) & np.isfinite(lfcR) & np.isfinite(lfcS) \
               & (seR > 0) & (seS > 0) & np.isfinite(seR) & np.isfinite(seS)
    opp = sig_both & (np.sign(lfcR) != np.sign(lfcS)) & (lfcR != 0) & (lfcS != 0)
    tot_sig_both += int(sig_both.sum()); tot_opp += int(opp.sum())
    tot_same += int((sig_both & ~opp).sum())

    pi, gi = np.where(opp)
    if len(pi):
        z = (lfcR[pi, gi] - lfcS[pi, gi]) / np.sqrt(seR[pi, gi]**2 + seS[pi, gi]**2)
        p = 2 * stats.norm.sf(np.abs(z))
        flip_p.append(p)
        for k in range(len(pi)):
            flip_meta.append((names[pi[k]], genes[gi[k]], float(lfcR[pi[k], gi[k]]), float(lfcS[pi[k], gi[k]]),
                              float(seR[pi[k], gi[k]]), float(seS[pi[k], gi[k]])))
    print(f"  chunk {s0//CH+1}/{-(-len(shared)//CH)}: sig_both={tot_sig_both:,} opp={tot_opp:,}", flush=True)

flip_p = np.concatenate(flip_p) if flip_p else np.array([])
print("\n" + "="*80)
print(f"edges significant in BOTH states : {tot_sig_both:,}")
print(f"  same sign                      : {tot_same:,}  ({100*tot_same/max(tot_sig_both,1):.2f}%)")
print(f"  OPPOSITE sign (raw candidates) : {tot_opp:,}  ({100*tot_opp/max(tot_sig_both,1):.2f}%)")

res = {"seed": 20260715, "min_outdeg": MIN_OUT, "n_candidate_perturbations": len(shared),
       "edges_sig_both": tot_sig_both, "same_sign": tot_same, "opp_sign_raw": tot_opp}

if len(flip_p):
    # BH-FDR on the SE-aware difference test, among opposite-sign candidates
    m = len(flip_p); order = np.argsort(flip_p); ranked = flip_p[order]
    q = ranked * m / (np.arange(1, m+1)); q = np.minimum.accumulate(q[::-1])[::-1]
    keep = np.zeros(m, bool); keep[order] = q < 0.05
    n_rob = int(keep.sum())
    print(f"  ROBUST inverting edges (BH-FDR<0.05 on the SE-aware difference test): {n_rob:,}")
    print(f"    = {100*n_rob/max(tot_sig_both,1):.3f}% of all edges significant in both states")
    regs = {}
    for k in np.where(keep)[0]:
        regs[flip_meta[k][0]] = regs.get(flip_meta[k][0], 0) + 1
    top = sorted(regs.items(), key=lambda x: -x[1])[:25]
    print(f"  regulators with >=1 robust inverting edge: {len(regs):,}")
    print("\n  TOP INVERTING REGULATORS (n robust sign-inverted targets):")
    for g, c in top: print(f"    {g:12s} {c:6d}")
    # magnitude sanity: are survivors real effects or boundary wobble?
    mags = np.array([min(abs(flip_meta[k][2]), abs(flip_meta[k][3])) for k in np.where(keep)[0]])
    print(f"\n  survivor |log_fc| (the SMALLER of the two states): median {np.median(mags):.3f}, "
          f"25th {np.percentile(mags,25):.3f}, 75th {np.percentile(mags,75):.3f}")
    res.update({"robust_inverting_edges": n_rob, "n_inverting_regulators": len(regs),
                "top_regulators": [{"gene": g, "n_inverted_targets": c} for g, c in top],
                "survivor_min_abs_lfc_median": float(np.median(mags))})
    import csv
    with open(f"{OUT}/data/census/sign_inverting_edges.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["regulator","target","lfc_Rest","lfc_Stim48","se_Rest","se_Stim48"])
        for k in np.where(keep)[0]: w.writerow(flip_meta[k])
    print(f"\n  wrote {OUT}/data/census/sign_inverting_edges.csv")
else:
    print("  NO opposite-sign candidates at all.")
    res["robust_inverting_edges"] = 0
json.dump(res, open(f"{OUT}/data/census/signflip_feasibility.json","w"), indent=2)
print("="*80)
