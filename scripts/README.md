# Scripts

## What to run

Nine entry points. Each is self-documenting; read its header for the analysis
it performs, its inputs, and any cost warning.

```bash
./scripts/setup_after_clone.sh     # once, restores file modes Git cannot store

python scripts/census.py           # shared census, needed by Figures 1-4
python scripts/figure1.py          # Figure 1
python scripts/figure2.py          # Figure 2
python scripts/figure3.py          # Figure 3
python scripts/figure4.py          # Figure 4
python scripts/figureS1.py         # Supplementary Figure S1
python scripts/figureS2.py         # Supplementary Figure S2
python scripts/figureS3.py         # Supplementary Figure S3
```

Supplementary Figure S4 has no renderer of its own: it is the published output
of the enrichment job in `hpc/`, whose result is committed under
`results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/`.

`hpc/` documents the three analyses that ran on a cluster. Their outputs are
committed, so nothing here needs one.

## Everything else in this directory

The remaining files are the **renderer history**, and they are here because
they are what produced the published figures, not because they are meant to be
read.

Each figure was revised many times, and every revision was kept as a new file
that calls its predecessor and adjusts the result. Figure 3 reaches eight such
modules across roughly 8,300 lines; Figure 2 reaches ten. The entry points
above drive the final version of each chain.

Flattening those chains into one file per figure would mean re-deriving about
25,000 lines of layout logic, and any error would change a figure that is
already published. The chains are therefore kept intact and the entry points
give a clean interface over them.

Two consequences worth knowing if you modify anything here:

- Several renderers verify the SHA-256 of their dependencies, so editing or
  moving a chain file breaks the figure that pins it.
- Several refuse to overwrite an existing output directory. The entry points
  clear their own output first.
