#!/usr/bin/env python
"""UMAP of the Arce Perturb-CITE-seq SURFACE-PROTEIN (ADT) channel -- the substrate panels of Fig5.

This is the ARCE two-donor pilot data (GSE278572), NOT the 1.13 TB Marson atlas (that pass ran on
HPC and is not local). Deliberately the ADT/protein channel, matching Fig5's other panels, which
read state from the ADT capture library rather than from the RNA being decomposed.

CIRCULARITY GUARD (load-bearing): the CellRanger "Antibody Capture" block contains the four HTO
hashtags (Resting_Teff, Stimulated_Teff, Resting_Treg, Stimulated_Treg) alongside the real ADTs,
and the state label HTO_maxID is by definition the argmax over exactly those four. Leaving them in
makes the embedding and the AUROC circular (label == feature; AUROC == 1.000 exactly). They are
dropped here, as are the Isotype_* nonspecific-binding controls, leaving 130 real surface proteins.

Writes data/arce/arce_adt_umap_coords.csv.gz and data/arce/arce_adt_umap_stats.json into the paper root; the figure
script reads both and hardcodes nothing.
"""
from __future__ import annotations
import json, os
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
import umap

S = str(Path(__file__).resolve().parents[2] / "scripts")
P = str(Path(__file__).resolve().parents[2])
SEED = 0

cells = pd.read_csv(f"{S}/aim3_panel_cells.csv").sort_values("dense_col").reset_index(drop=True)
feat = pd.read_csv(f"{S}/aim3_adt_features.csv").sort_values("dense_row")
NC, NF = len(cells), len(feat)
print(f"cells={NC:,}  antibody-capture features={NF}")

trip = pd.read_csv(f"{S}/aim3_adt_triplets.txt.gz", sep=r"\s+", header=None,
                   names=["f", "c", "v"], dtype={"f": np.int32, "c": np.int32, "v": np.float32})
X = np.zeros((NC, NF), dtype=np.float32)
X[trip.c.values, trip.f.values] = trip.v.values
del trip

HTO = {"Resting_Teff", "Stimulated_Teff", "Resting_Treg", "Stimulated_Treg"}
names = feat.name.values
keep = np.array([(n not in HTO) and (not n.startswith("Isotype_")) for n in names])
dropped_hto = sorted(set(names[~keep]) & HTO)
dropped_iso = sorted(n for n in names[~keep] if n.startswith("Isotype_"))
X = X[:, keep]
print(f"dropped {len(dropped_hto)} hashtags + {len(dropped_iso)} isotype controls "
      f"-> {int(keep.sum())} real surface proteins")

L = np.log1p(X)
X = L - L.mean(axis=1, keepdims=True)          # CLR, standard for ADT
del L

act = np.where(cells.HTO_maxID.str.startswith("Resting"), "Rest", "Stim")
sub = np.where(cells.HTO_maxID.str.contains("Treg"), "Treg", "Teff")
don = cells.donor.values

pcs = PCA(n_components=30, random_state=SEED).fit_transform(X)
emb = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="euclidean",
                random_state=SEED).fit_transform(pcs)

# Activation axis, quantified. PC1 of the CLR ADT space; AUROC is orientation-free.
y = (act == "Stim").astype(int)
auc = lambda m: float(max(roc_auc_score(y[m], pcs[m, 0]), 1 - roc_auc_score(y[m], pcs[m, 0])))
allm = np.ones(NC, bool)
stats = {
    "n_cells": int(NC),
    "n_features_used": int(keep.sum()),
    "dropped_hashtags": dropped_hto,
    "dropped_isotypes": dropped_iso,
    "auroc_activation_all": auc(allm),
    "auroc_activation_by_donor": {str(d): auc(don == d) for d in np.unique(don)},
    "n_by_subset": {k: int(v) for k, v in pd.Series(sub).value_counts().items()},
    "n_by_state": {k: int(v) for k, v in pd.Series(act).value_counts().items()},
    "umap_params": {"n_neighbors": 30, "min_dist": 0.3, "n_pcs": 30, "seed": SEED},
}
print("  AUROC all={auroc_activation_all:.3f}".format(**stats),
      {k: round(v, 3) for k, v in stats["auroc_activation_by_donor"].items()})

pd.DataFrame({"act": act, "donor": don, "subset": sub,
              "u1": emb[:, 0], "u2": emb[:, 1]}).to_csv(
    f"{P}/data/arce/arce_adt_umap_coords.csv.gz", index=False)
json.dump(stats, open(f"{P}/data/arce/arce_adt_umap_stats.json", "w"), indent=2)
print(f"wrote {P}/data/arce/arce_adt_umap_coords.csv.gz + data/arce/arce_adt_umap_stats.json")
