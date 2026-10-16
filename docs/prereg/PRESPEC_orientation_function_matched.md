# Pre-specification: does reversal orientation track target function, once activation induction is held?

Written before any matched statistic was computed. Registered 2026-09-03.

## The question and why the earlier tests do not answer it

The network-calibrated test in `fig3b_compound_target_pathway_randomization_freeze` already asked
whether the 279 single-orientation targets are enriched for any KEGG term against a null that
rebuilds the selection rule. It supported no term (leading term regulation of autophagy,
maxT p = 0.481). That closes enrichment of the set.

This asks a different thing: given that a target reverses, does the DIRECTION it reverses in track
its function? The 161 lost-on-activation and 118 gained-on-activation targets serve as each other's
background, so "is a reversal target" cancels by construction rather than by modelling.

## The confound this design exists to hold

Orientation is very strongly associated with the target's own expression trajectory: median log2
change from Rest to 48 h is +1.253 for lost targets and -0.044 for gained targets, Mann-Whitney
p = 1.5e-27 (`results_data/census/orientation_stratified_function.json`). Cytokine-receptor genes
are among the most induced, so an unmatched contrast cannot separate "lost effects sit on cytokine
genes" from "lost effects sit on strongly induced genes".

Every target is therefore assigned to one of **four quartile bins** of its own log2 expression
change, and orientation labels are permuted **within bins**. A term can only win by beating genes
with a comparable induction profile. Four bins rather than more, because orientation is so strongly
predicted by trajectory that finer bins risk becoming single-orientation and carrying no
information; the realised per-bin composition is reported whatever it shows.

## Statistic and null

Per term, the statistic is the number of its member targets that are lost-on-activation. The null
permutes the orientation label within each trajectory bin, holding fixed each bin's lost/gained
composition, each term's membership, and the 279-target universe. 10,000 draws, seed 20260903.
Two-sided p as twice the smaller tail, capped at 1, with the plus-one convention.

## Families, declared in advance

- **A, the pinned 19-term dictionary** already used in the paper. Small, so multiplicity is light.
  Declared, not pre-registered: those 19 terms were grouped after the census was examined, and that
  is stated wherever the family is reported.
- **B, GO BP 2023**, every term with at least 5 of the 279. The honest wide comparison.

Benjamini-Hochberg within each family separately. Terms need at least 5 of the 279 and at least one
target of each orientation somewhere in the set.

The 16 proliferation and 12 cytokine genes named by the earlier exploratory pass are post-hoc. They
are reported as a located observation and belong to neither family.

## Decision rule, fixed before the run

Positive if any term reaches q < 0.10 in family A or family B after matching. Otherwise negative,
and the reportable finding becomes that reversal orientation encodes activation induction, which
explains why the pathway tests keep failing.

## Calibration

The pipeline is rerun with orientation permuted once at the top. If terms still survive, the design
is discarded rather than reported.

Output: `results_data/census/orientation_function_matched.json`
