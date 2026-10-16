# Pre-specification: cell-cycle phase as a candidate axis behind a subset of reversals

Written 2026-09-04, **before any of the analyses below were run**. Companion to
`CELLCYCLE_AXIS_PLAN.md`. Nothing here is a claim; this file fixes what will be computed and how
it will be read, so that the reading cannot be chosen after seeing the numbers.

## Background in one paragraph

A sign reversal can arise because the regulator's relationship to the target changed, or because
the perturbation moved cells along a latent axis so that a population-level coefficient changed
sign with no per-cell change. The paper already bounds the second explanation on the **activation**
axis: measured redistribution was too small for any of 4,102 covered edges (Supplementary Figure
S2c,d). **Cell-cycle phase is a second latent axis, never measured, and it exists only after
stimulation because resting cells are quiescent.** Supplementary Results R6 reports that
cell-cycle genes carry gained effects almost exclusively and that this survives matching on target
induction. That is the pattern the phase-axis explanation predicts, and it is also the pattern a
genuine regulatory change would produce. This file specifies how to tell them apart.

## Stage 1. The pre-registered gene set

The R6 proliferation set was assembled after inspecting the data and is therefore not evidence.
It is replaced here by a set built from a **mechanical rule over an external versioned file**,
`data/genesets/GO_Biological_Process_2023.gmt`, fixed now:

> Take every GO Biological Process 2023 term whose name contains any of
> `Mitotic`, `Cell Cycle`, `DNA Replication`, `Chromosome Segregation`, `Nuclear Division`,
> case-insensitive. The gene set is the union of their members.

No term is selected individually, so no term can be selected because it performed well. The rule
is applied verbatim and the resulting gene count is recorded before any test is read. Genes are
then intersected with the 279 single-orientation targets to give the tested set.

**Primary Stage 1 test.** The induction-matched orientation test already implemented in
`scripts/orientation_function_matched.py`: orientation labels permuted only within quartiles of
each target's own Rest-to-48 h expression change, 10,000 draws, two-sided, count of lost-oriented
members as the statistic. One test, one family of one, no correction needed and none claimed.

**Stage 1 decision.** If the pre-registered set does not reproduce the direction and rough
magnitude of the post-hoc set, the finding does not survive contact with a gene set we did not
choose, and Stage 2 does not run. This is a stopping rule, not a guideline.

## Pre-checks, run alongside Stage 1

Three checks on data already held. They are **descriptive and directional**: each states in
advance what the phase-composition explanation predicts, so that agreement is informative and
disagreement is informative. They carry no multiplicity correction and no p-value is treated as
confirmatory. They exist to stop an expensive Stage 2 that is already contradicted.

**P1. Response-shape signature.** Every supported reversal is already decomposed into a detection
component (how many cells carry detectable target RNA) and a positive-intensity component (how
much, among those cells). A shift in phase composition changes *which cells are cycling*, so it
changes how many cells carry detectable cell-cycle transcripts.
*Prediction if compositional:* reversals on cell-cycle targets are **detection-dominant** more
often than reversals on other targets. Statistic: Fisher exact on detection-dominant versus
intensity-dominant, cell-cycle targets versus the rest.

**P2. Timing fingerprint at 8 h.** Cells enter cycle between 8 h and 48 h, later than the effector
programme fires.
*Prediction if compositional:* reversals on cell-cycle targets are **less often resolved at 8 h**
than reversals on other targets, because the axis has not yet come online. Statistic: Fisher exact
on resolved versus unresolved, and, among resolved edges, on Stim48-sign-matching versus
Rest-sign-retaining.

**P3. Regulator character.** If the effect is a plain proliferation-rate change rather than an
axis-composition change, the regulators carrying it should themselves be growth machinery.
*Statistic:* Fisher exact on whether regulators carrying cell-cycle-target reversals are annotated
to a declared housekeeping/growth family (KEGG ribosome, oxidative phosphorylation, spliceosome,
proteasome, DNA replication) more often than regulators carrying only other reversals.
**Caveat fixed in advance:** this is annotation-based, and annotation density has already produced
one artifact in this project (GO density rho = -0.448, Supplementary Results R3). A positive P3 is
treated as a warning, never as support.

## What Stage 2 would be, if Stage 1 passes

Per-cell S and G2M scores from the Stage 1 gene set; then (2a) whether these knockdowns shift
phase composition at 48 h, (2b) whether measured phase-axis redistribution is large enough to
cross each edge's own sign-change threshold, using the same estimator as the activation-axis bound,
and (2c) whether the orientation contrast survives within phase strata. Families declared
separately, calibration run signal-free first, control-versus-control null mandatory. Stage 2 is
not specified further here because Stage 1 may stop it.

## Reporting rule

Whatever the outcome, it is reported. A Stage 1 failure is written into Supplementary Results R6
as a resolved negative and the Discussion clause is deleted. A Stage 1 pass is reported as a pass
of a pre-registered set and nothing more until Stage 2 runs.

---

## Addendum, 2026-09-04: Stage 1b, written before running

Stage 1 passed and both compositional pre-checks came back flat. Stage 2 requires per-cell phase
scores and may not be affordable. Rather than default to it, two tests are declared here that the
data already in hand can answer well. Both are chosen because they **discriminate** between the
two readings rather than merely describing the pattern again.

### Test A. Continuous response-shape magnitude, not the dominance call

P1 asked only which component was larger. That throws away magnitude, and a binary call has poor
resolution. Test A compares the **oriented continuous** detection and positive-intensity
components between reversals on cell-cycle targets and reversals on all other targets.

*Prediction if compositional:* the detection component is larger in magnitude for cell-cycle
targets, because a shift in which cells are cycling changes how many cells carry detectable
cell-cycle transcript.
*Statistic:* Mann-Whitney two-sided plus Cliff's delta on each component separately, reported with
n for each group. Two tests, declared as one family of two, Holm-corrected.

### Test B. Within-regulator coherence, Rest versus 48 h

This is the discriminating test. If a single latent phase axis drives these effects, then a
perturbation that moves cells along it moves **all** cell-cycle genes together, because they all
read the same axis. Independent regulatory events do not produce that. Cell-cycle genes are
co-regulated at baseline as well, so raw coherence is not informative; the phase axis, however,
exists only after stimulation, since resting cells are quiescent and sit at one point on it.

*Statistic:* for each regulator, over all its comparable edges onto cell-cycle targets, coherence
is the majority-sign fraction of its knockdown coefficients, computed separately at Rest and at
48 h. The quantity of interest is the change, $\Delta C = C_{48} - C_{Rest}$, averaged over
qualifying regulators.
*Control:* for the same regulator, the same statistic on an equally sized set drawn from its own
non-cell-cycle targets, so regulator identity, edge count and out-degree are all held fixed.
*Null:* 10,000 draws of those size-matched control sets; two-sided plus-one Monte Carlo p.
*Inclusion:* regulators with at least five comparable cell-cycle targets and at least five
non-cell-cycle targets in both states.

*Prediction if compositional:* $\Delta C$ is clearly positive for cell-cycle target sets and
larger than for the size-matched control sets, because the axis switches on between the two states.
*Prediction if regulatory:* no excess, since independent regulatory changes carry no reason to
align.

### Reporting rule for Stage 1b

Both tests are reported whatever they show. Neither is a mechanism claim. A flat Test B is
reported as absence of the coherence a single shared axis would produce, which weakens the
compositional reading without establishing the regulatory one. Power and group sizes are reported
alongside every result so that a null can be read as informative or uninformative rather than
being assumed to be either.

---

## Addendum 2, 2026-09-04: Stage 2a, a depth-controlled phase variable

The Stage 2 pre-flight found that the existing per-cell phase call is largely a sequencing-depth
readout: non-G1 cells carry 1.9x the UMI of G1 cells at Rest and 2.3x at Stim48hr, Spearman
between depth and the maximum of the two scores is 0.53 at Rest, and the Rest non-G1 fraction runs
from 0.144 to 0.968 across UMI deciles. The condition difference nonetheless survives inside every
decile. The axis is therefore real and the variable is not usable as it stands. Stage 2 does not
run on the raw call.

### The correction

Fit, per donor and condition, an ordinary least-squares regression of each score on
$\log_{10}(\text{n\_umi})$ **using non-targeting cells only**, then take residuals for every cell.
Fitting on NTC cells matters: it means the correction is estimated where no perturbation effect
exists, so a real perturbation-driven change in phase cannot be absorbed into the depth term.
Phase is re-called from the residual scores by the same rule as the original (G1 when both are at
or below zero, otherwise whichever is larger).

### Acceptance criteria, fixed before the correction is run

Stage 2 proceeds only if all four hold. Any failure means phase is not measurable in this data and
that becomes the reported result.

1. **Depth is removed.** $|\rho(\text{n\_umi}, \max(S', G2M'))| < 0.10$ in both conditions.
2. **The depth gradient is flattened.** The non-G1 fraction across UMI deciles spans less than
   0.25, against 0.82 uncorrected.
3. **The biology survives.** The Stim48hr minus Rest difference in non-G1 fraction remains
   positive and at least 0.20, and remains positive in at least 8 of 10 deciles.
4. **The null is clean.** Non-targeting cells resampled as pseudo-perturbations at matched cell
   counts show a median composition shift indistinguishable from zero, and real perturbations are
   read only as a contrast against that null centre, never against zero.

### If accepted

The same fitted correction is applied genome-scale across the atlas shards, because the 160,000
cell sample carries a median of six cells per perturbed gene per condition and cannot support
per-perturbation estimates. The magnitude convention of the existing atlas composition analysis
carries over unchanged: report the composition magnitude against its random-state zero-information
floor and against the endpoint it would need to reach, never a retention rate.
