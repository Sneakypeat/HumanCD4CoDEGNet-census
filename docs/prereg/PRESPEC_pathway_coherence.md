# Pre-specification: pathway-level orientation coherence of reversals

Written 2026-08-24, BEFORE any coherence statistic was computed. Frozen: the filters and the
decision rule below were fixed on term metadata (category, size, edge count) only. No coherence
value was inspected before this file was written.

## Question

Reversals converge on shared targets (convergence_test_results.json: 110 targets receive >=5
reversals vs 16.3 expected, p = 4e-4). The identity of those targets is NOT the effector
programme (10.0% annotated vs 9.74% null, p = 0.463). That closes the *enrichment* question.

The open question is different: **when several reversals land on members of the same pathway, do
they run the same way?** Orientation, not membership.

Orientation is the `direction` column of the trans census: `lost_on_activation` (lfc_Rest > 0,
lfc_Stim48 < 0; n = 2,298) or `gained_on_activation` (lfc_Rest < 0, lfc_Stim48 > 0; n = 2,081).

## Why this survives the annotation-density confound

Annotation density decides how MANY reversals land in a term. It does not decide WHICH WAY they
run. The statistic is computed within a term, conditional on the reversals already there, so a
well-studied gene appearing in more pathways cannot bias the direction. This is the confound that
killed the GO enrichment attempt (rho = -0.448 timing artefact) and it does not apply here.

## Universe

`results_data/census/census_master_edges.csv.gz`, trans only (`is_self == False`):
83,489 edges, 79,110 stable, 4,379 reversing, 324 reversal-bearing regulators, 9,554 targets.

## Annotation

**Primary: MSigDB `c2.cp.kegg_legacy` v2026.1.Hs** (186 sets), fetched 2026-08-24 to
`data/genesets/`. Categories come from KEGG BRITE `br08901` (`data/genesets/kegg_br08901.txt`)
joined through each set's `exactSource` KEGG map id (`KEGG_ALPHA_LINOLENIC_ACID_METABOLISM` ->
`hsa00592` -> BRITE C-level `00592`). The join is exact and numeric: 186/186 sets resolve, zero
unresolved, and no name matching is involved. The script aborts if any set fails to resolve.

Stated plainly in Methods: MSigDB froze this collection at KEGG's 2011 content when KEGG moved to
a subscription licence. It is therefore a current *release* wrapping older *biology* than the
Enrichr `KEGG_2021_Human.gmt` also held in `data/genesets/`. It is used as primary because it is
versioned and citable. The Enrichr 2021 set is carried as a sensitivity precisely to show what the
coverage difference costs (240 terms and 2,043 reversal edges against 118 and 1,450).

For the Enrichr sensitivity only, whose GMT carries no map ids, categories are assigned by BRITE
display name: 305/320 exact, 8 by unique prefix, 7 by a hand-written map recorded in the script.
None of the 7 is a disease term, so the disease exclusion does not depend on any of those calls.

## Term filters (applied before any coherence value is read)

1. **Exclude BRITE category `Human Diseases`.** KEGG's disease entries ("Pathways in cancer",
   "Pathways of neurodegeneration", Alzheimer / ALS / Parkinson / prion disease) are unions of
   signalling modules plus ribosomal, mitochondrial and proteasomal genes. They are the largest
   terms in this dataset by reversal-target count and would top any size-driven ranking, and they
   are not mechanisms in CD4+ T cells.
2. **Exclude BRITE category `Drug Development`.** Not a cellular pathway.
3. **Term size 15-200 eligible targets**, eligible = present as a target in the trans census.
   Lower bound for power, upper bound to drop supersets.
4. **>= 10 reversal edges** in the term.
5. **>= 3 distinct reversal-bearing regulators** in the term.

## Statistic

For term P, let the reversal edges with target in P be E_P. Coherence is the majority-direction
fraction, `C_P = max(n_lost, n_gained) / |E_P|`, bounded below by 0.5.

## Null (this is the whole design)

The dominant confound is REGULATOR IDENTITY, not annotation. A regulator whose reversals are
mostly `lost` will make every pathway it touches look coherent, with no pathway biology involved.

So: **within-regulator permutation of the direction label.** For each regulator independently,
its own direction labels are permuted across its own reversal edges. This holds exactly fixed:
(i) every regulator's count of reversals, (ii) every regulator's lost/gained composition,
(iii) every target's reversal count, (iv) which regulator-target pairs exist. The only thing
destroyed is the association between direction and WHICH target received it. Excess coherence
over this null is therefore pathway structure and nothing else.

10,000 permutations, seed 20260824. Per-term p is one-sided upper (coherence is directional by
construction), (#null >= observed + 1) / (nperm + 1). Multiplicity: Benjamini-Hochberg across
surviving terms only.

## Global test

One headline number in addition to the per-term tests: the edge-weighted mean coherence across
surviving terms, against the same permutation null. This is a single test and does not depend on
the per-term multiplicity.

## Decision rule, fixed in advance

- **Positive** if the global edge-weighted mean coherence exceeds its null 95th percentile
  (p < 0.05) OR at least one term survives BH q < 0.10.
- **Negative** otherwise, and it is reported as a negative in two sentences. A negative is a
  publishable outcome here: it says reversals concentrate on shared targets without imposing a
  shared direction on the pathways those targets sit in.
- Sensitivity, run regardless of outcome and reported regardless of outcome: GO BP 2023 with the
  identical procedure, and CORUM at complex level. These are reported as agreement / disagreement,
  not used to rescue a negative primary.

## Follow-up (conditional, descriptive, and labelled as such)

If the primary is positive, the targets carried by the surviving terms are submitted to ORA
against **KEGG modules** (`rest.kegg.jp`), with background = all reversal targets, NOT the genome.
Testing those genes against KEGG pathways would be circular because the pathways defined them;
modules are a finer partition and the background restriction removes the "is a reversal target"
signal. Even so this step is descriptive, is reported as descriptive, and carries no claim of
independent confirmation.

## Output

- `results_data/census/pathway_coherence.json`
- one panel appended to Figure 3
- 2-3 sentences in Results
