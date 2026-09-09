# Cluster jobs

Three analyses were run on the UAEU HPC because they do not fit on a
workstation. **Their outputs are committed**, so every figure in this
repository rebuilds without access to a cluster. These files document how those
cached outputs were produced; they are not part of the figure build.

| Job | What it does | Cost | Feeds |
|---|---|---|---|
| `atlas_state_decomposition.sbatch` | Per-cell Shapley decomposition over the activation axis, 8 shards of ~130 GB | 24 h, 64 GB | Supplementary Figure S2c,d - the bound showing measured cell-state redistribution is too small to account for the reversals |
| `enrichment_permutation_banks.sbatch` | Network-adjusted rank enrichment: two independent 10,000-draw permutation banks over 2,593 KEGG and GO terms | hours | Figure 3e and Supplementary Results R7 |
| response-shape production (`scripts/response_shape/fig4_production_v5_1.py`) | Donor-balanced two-part decomposition of 1,899 reversals into target detection and positive-cell intensity | hours | Figure 4 |

The source atlas (39 GB) is not redistributed - see Data availability in the
manuscript. Paths inside these files point at the cluster account they ran on
and will not resolve elsewhere; adapt them if you re-run.
