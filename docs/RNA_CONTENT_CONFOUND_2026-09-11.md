# Global RNA content as a confound for sign reversal

Written 2026-09-11, after submission to PLOS Computational Biology, and
incorporated into the manuscript on 2026-09-18 as Supplementary Methods M21
and Supplementary Figure S5. The post-hoc audit was prompted by CIPHER
(Kuznets-Speck et al., bioRxiv 2025.06.27.661814); its status and claim
boundary remain explicit.

## The objection

CIPHER decomposes a perturbation response into a component parallel to the mean
expression vector and one orthogonal to it, and argues that library-size
normalisation removes the parallel component. Applied to this census:

Activated CD4+ T cells blast and carry more RNA per cell than resting ones, so
the two census endpoints have different global means. If a perturbation shifted
total RNA in **opposite** directions at Rest and at 48 h, the compositional
residual for its targets could change sign after normalisation without the
underlying regulation reversing at all.

This is distinct from the cell-state composition confound already tested and
refuted (Supplementary Figure S2c, `atlas_edge_composition`). That one concerns
which cells are in which state; this one concerns how much RNA each cell holds.

The premise is real in this atlas: median total UMI in non-targeting cells is
7,674 at Rest and 11,105 at 48 h, a 1.45-fold increase, measured over all
14,853,023 analysed cells.

## The standing guard

Every log2 fold change in the census is a within-state contrast against that
state's own non-targeting controls. The reference is state-matched, so the
global shift largely cancels. The residual risk is narrower: a perturbation
whose *parallel* response differs in sign between the two states.

## The discriminating test

A global RNA-content artifact acts on a **perturbation**. It moves all of one
regulator's targets together, so reversal direction would be organised at the
regulator level. It cannot make **different** regulators agree on the direction
of reversal at a shared target.

The null shuffles reversal direction within each regulator, holding that
regulator's lost/gained mix exactly fixed. That is the structure the artifact
would impose. Target-level agreement above this null cannot be a per-regulator
normalisation effect.

## Result

| quantity | value |
|---|---|
| targets with at least 2 reversals | 913 |
| of those, single-orientation | 633 (69.3%) |
| regulators with at least 2 reversals | 227 |
| of those, single-orientation | 58 (25.6%) |
| within-regulator permutation null | mean 37.0%, 20,000 draws |
| observed against that null | *p* = 5e-05 |

Direction is coherent within targets and not within regulators. The artifact
predicts the opposite asymmetry. It is not the source of the organisation.

## The direct test, on the full atlas

The test above is structural. The direct form asks the confound's own question:
per perturbed gene, the log2 ratio of median total UMI in perturbed versus
non-targeting cells, computed within each state. It is underpowered on the
local QC subsample, where only 8 census-matched regulators clear 25 cells in
both states, so it ran on the cluster against all eight atlas shards
(`scripts/hpc/rna_content_confound.sbatch`, SLURM 486540 and 486554, UAEU HPC).

Every column the test needs sits in `obs`, so the expression matrix is never
read and the 1.13 TB atlas streams in under four minutes.
14,853,023 cells survive the standard selection, 11,818 of 12,585 perturbations
clear the cell floor, and 616 map to a census regulator, against 8 locally.

### The confound's sharp prediction fails

A manufactured reversal requires the perturbation's global RNA effect to change
**sign** between states. Regulators where it does are not enriched among
reversal carriers. They are depleted.

| | opposite-sign shift | same-sign shift |
|---|---|---|
| regulators | 209 | 407 |
| median reversal rate | 0.0000 | 0.0145 |

Mann-Whitney *p* = 0.0044, Cliff's delta = **-0.132**. Reversal carriage against
opposite-sign shift gives an odds ratio of **0.599** (*p* = 0.0029), and 0.601
after matching on census edge-count quartile. The association runs opposite to
the artifact throughout.

### Reversal direction is unrelated to the RNA shift

The artifact does not merely predict more reversals. It predicts *which way*
they run: a regulator whose perturbed cells gain global RNA on activation
should push its edges one way, and one whose cells lose it the other. Across
the 187 regulators carrying at least three reversals, the change in global RNA
output and the fraction of reversals running "gained on activation" are
unrelated: Spearman rho = **+0.046**, *p* = 0.53.

The edge-level version of the same comparison returns *p* = 0.002 across 4,376
reversal edges, but the effect is negligible (Cliff's delta +0.053, medians
-0.115 against -0.130 log2 units) and the edges are clustered inside 323
regulators, so that test counts the same perturbation up to hundreds of times.
The regulator-level result is the one that holds.

### One association that does survive

The **magnitude** of the between-state shift does track reversal burden:
Spearman rho = +0.182 (*p* = 5.8e-06, n = 616). It is not a power artifact.
Partial Spearman controlling for the regulator's census edge count gives
+0.161, for cell count +0.176, and for both +0.146 (*p* = 2.8e-04), and the
shift gap barely tracks edge count (rho 0.09) and not cell count (rho -0.05).

This is a magnitude association, not a sign association, and it does not
support the artifact. A regulator whose knockdown changes global RNA output
differently at rest and at 48 h is, by that fact, state-dependently active,
which is also what produces state-dependent target effects. The two readouts
share a cause. What the artifact requires, and what is absent, is that the
**direction** of the RNA change set the **direction** of the reversal.

State this association rather than omitting it. It is the honest version of
the result, and the direction tests are what answer the objection.

## Claim boundary

This rules out a per-regulator global RNA-content artifact as the source of
reversal **organisation** and as the source of reversal **direction**. It does
not measure the parallel response component itself, and it is not a general
endorsement or refutation of CIPHER. The magnitude association above stands
and is reported.

## Released files

- `data/census/rna_content_confound.json` contains the structural test
- `data/census/rna_content_umi_shift_full_atlas.json` and
  `data/census/rna_content_per_regulator.tsv` contain the direct full-atlas test
- `scripts/analysis/rna_content_full_atlas.py` and
  `scripts/analysis/rna_content_direction_coherence.py` implement the analyses
- `scripts/hpc/rna_content_normalisation.sbatch` records the production Slurm
  resources and commands
- `scripts/figureS5.py` renders `figures/FigS5.pdf` from the compact outputs
- `docs/Extended_Protocol_M1-M21.pdf` contains the released M21 method
