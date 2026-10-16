# Pre-specification: co-target structure among reversal-bearing regulators

Written before any statistic was computed. Registered 2026-09-03.

## Why this test and not the one that was abandoned

Supplementary Results R1 withdraws the complex-level enrichment claim: the degree-matched
complex permutation was not uniformly calibrated (false-positive rates 2.0% to 12.5%). That
test asked a **burden** question, does a named complex reverse more often than its degree
predicts. This pre-specification does not revive it.

The question here is **structural**, and it is the one Figure 3d and 3e need. Those panels show
exhaustive intersections, 16 reversal targets shared by CCNC, TADA2B and SUPT7L and 17 shared by
AHR and ARNT, and the legend states they are descriptive. A referee can reply that any two
high-burden regulators share targets by arithmetic. The estimand below asks whether reversal
sharing is concentrated into specific regulator pairs beyond what both margins already force.

## Universe, fixed

The 324 regulators carrying at least one reversal in the frozen census, their non-self
both-significant trans edges, and the 4,379 reversals on them. A regulator/target pair enters the
denominator only if it is comparable, meaning significant in both states and therefore eligible to
reverse. Columns comparable to a single regulator are inert: they cannot contribute to any overlap
and no permutation below can move a reversal into or out of them, so they are dropped and the drop
is recorded.

## Statistic

For regulators i and j, o_ij is the number of targets carrying a reversal from both.

Global statistic S = sum over unordered pairs of o_ij squared. S rises when the same amount of
sharing is concentrated into fewer regulator pairs, which is what a module is. Counts of pairs at
o_ij >= 5 and >= 10 are reported alongside as descriptive readouts, and no threshold is used to
define the test.

## Two nulls, reported separately, never combined

Both reassign reversal labels inside the fixed comparable edge set, so opportunity is preserved by
construction in each.

- **N1 regulator-preserving.** Each regulator keeps its own reversal count, reassigned uniformly
  among its own comparable edges. Target burden is free. Figure 3b already established that
  reversals concentrate on fewer targets than this null predicts, so N1 is expected to be exceeded
  and is reported only to show that the present statistic sits downstream of that result.
- **N2 fixed-margin, decisive.** Curveball trades restricted to the comparable mask. Each
  regulator keeps its reversal count and each target keeps its reversal count. Target
  concentration is held exactly, so the sum over targets of C(k_t, 2), the total quantity of
  sharing, is constant across all N2 draws. Only the allocation of that sharing to regulator pairs
  varies. An excess under N2 is module structure that is not target concentration restated.

10,000 draws for each null, with 5 row-pair trades per row between recorded N2 draws.
Seed 20260903. One-sided upper with the plus-one convention.

## Named modules, declared circular

The Figure 3d and 3e sets were chosen because their members share targets, so their own overlaps
cannot be evidence for the global claim. They are reported as located values against both nulls,
and the global S is recomputed with those five regulators removed as the leave-out control.

## Calibration, run before any p-value is read

The whole pipeline is rerun once with the reversal labels permuted at the top. If the global p is
not consistent with uniform and the named-module p-values do not lose their signal, the design is
discarded rather than reported. This is the check the abandoned complex test failed.

## Decision rule, fixed in advance

Positive only if S exceeds its N2 null at the plus-one one-sided p < 0.05 and the calibration run
is clean. If positive, an unsupervised decomposition follows, and the family choice, decomposition
rather than hard clustering, is taken from Saelens et al. 2018, which benchmarked module detection
methods and found decomposition best for overlapping modules. Overlap is the operative property
here because a regulator may share targets with several others. If negative, no module panel goes
into Figure 3 and the panels stay descriptive as they are.

Output: results_data/census/module_cotarget_structure.json
