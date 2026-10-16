# Pre-specification: continuous-score GSEA on excess reversal burden

Written 2026-09-04, before running. Author's own framing, quoted: "a defensible GSEA-like
analysis would rank all 9,554 census targets, including targets with no reversal; use a
network-adjusted target score, such as excess reversal burden or convergence relative to its edge
opportunities; obtain significance from the same regulator-margin-preserving network
permutations." This is a genuinely different question from every prior functional test in this
paper, all of which compared a SELECTED subset (110 convergent targets, 279 single-orientation
targets) against a background. This one ranks the WHOLE 9,554-target universe and asks whether a
gene set trends toward high burden, the way standard GSEA asks whether a set trends toward one end
of a differential-expression ranking.

## Universe

All 9,554 targets carrying at least one both-state-significant (comparable) trans edge in
`census_master_edges.csv.gz`, `is_self == False`. 6,977 of these carry zero reversals and are
KEPT IN, not filtered out, because a continuous score is defined for them too.

## Score

For target $t$, let $\text{obs}(t)$ be its observed reversal count. Under the regulator-margin
null already used for `convergence_test.py` and `module_cotarget_structure.py` (each regulator
$r$ contributes exactly its observed reversal count $k_r$, drawn uniformly without replacement
from its own $n_r$ comparable edges), the expected reversal count is analytic:
$$\text{exp}(t) = \sum_{r \to t} \frac{k_r}{n_r}$$
summed over regulators with a comparable edge into $t$. The score is excess burden,
$\text{score}(t) = \text{obs}(t) - \text{exp}(t)$, defined and generally non-zero even for
zero-reversal targets (a target regulated mostly by high-rate regulators has negative excess even
at zero observed reversals).

## Test

Classic weighted (p=1) Subramanian running-sum enrichment statistic, computed once on the observed
ranking (targets sorted by score, descending) for every gene set passing the size filter, in two
declared families tested and BH-corrected separately: the FULL KEGG Legacy v2026.1 library and the
FULL GO Biological Process 2023 library (not the 19-term pinned subset used in earlier sections,
because this test's premise is to give every process a fair, unselected look).

**Filter, fixed before any ES is read:** term must have between 10 and 300 members inside the
9,554-target universe.

## Null

1,000 draws. In each draw, every regulator's comparable edges are re-sampled with the exact
`convergence_test.py` mechanism (`rng.choice(pool, kobs[r], replace=False)`), giving a permuted
observed-count vector; the analytic expected vector is UNCHANGED (it depends only on the fixed
regulator margins), so permuted score $=$ permuted observed $-$ fixed expected. The full weighted
ES is recomputed on the re-ranked permuted score for every term. Per-term normalized enrichment
score (NES) is ES divided by the mean of the null ES magnitudes for that term (sign-preserving,
the standard GSEA normalization), and the two-sided plus-one p-value is read against the
per-term null NES distribution. BH q within each family.

Seed 20260904.

## Decision rule

Positive for a family if any term reaches BH $q < 0.10$. Both families reported regardless of
outcome. This is exploratory: it asks whether burden trends by process, not whether any individual
edge or target is a validated finding.

## Outputs

`results_data/census/gsea_reversal_burden.json`, one figure summarising both families.
