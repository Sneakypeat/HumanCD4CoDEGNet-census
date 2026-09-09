# HumanCD4CoDEGNet — census

Sign-reversal census of trans-regulatory effects during human CD4+ T-cell activation.

A regulator can raise a target's expression in one cellular state and lower it in another. This repository holds the census that measures how often that happens genome-wide: of 83,489 regulator–target interactions significant in both the resting and 48 h stimulated states, 4,379 (5.2%) reverse sign.

## Data Availability
The released census is provided as a signed, directed graph in `data/census/census_master_edges.csv.gz`. It carries the per-state log2 fold change, standard error and adjusted p-value of every edge.

*   **Source atlas:** genome-scale CRISPRi Perturb-seq of primary human CD4+ T cells,
    Zhu et al., *Cell* (2026), [doi:10.1016/j.cell.2026.08.002](https://doi.org/10.1016/j.cell.2026.08.002).
    The processed differential-expression statistics this work reads are distributed
    through the Virtual Cells Platform; the 39 GB atlas is not redistributed here.
*   **Third-party annotations:** not redistributed. The frozen tables derived from them
    are committed under `data/fig3c_fixed_annotation_stack_v4/source_reference/`, and
    that directory's `go_slim_source_manifest.json` records the release and checksum of
    each download. The releases used were:
    *   GO ontology `go-basic.obo`, release 2026-07-26, from [current.geneontology.org/ontology/go-basic.obo](https://current.geneontology.org/ontology/go-basic.obo)
    *   GO eukaryotic cellular-process ribbon slim, same release
    *   UniProt human GOA annotations `goa_human.gaf`, GAF 2.2, generated 2026-07-28, from [current.geneontology.org/annotations/goa_human.gaf.gz](https://current.geneontology.org/annotations/goa_human.gaf.gz)
    *   UniProt human `goa_human.gpi`, GPI 2.0, generated 2026-07-28
    *   Reactome pathways GMT, v97, from [reactome.org/download-data](https://reactome.org/download-data)

## Prior work

The co-differential-expression approach this census builds on was developed in yeast:

> Nasar MI, Rehman SSU, Ott S, Alam MT. Uncovering coordinated pathway interactions
> through gene co-differential expression in yeast. *Nucleic Acids Research* 54(1),
> gkaf1410 (2026). [doi:10.1093/nar/gkaf1410](https://doi.org/10.1093/nar/gkaf1410)

## Reproducing the Results

**Environment Setup:**
Install the required dependencies using the provided `requirements.txt`:
```bash
pip install -r requirements.txt
```

**Run Figures:**
Execute the polished figure scripts to reproduce the manuscript's figures:
```bash
python scripts/01_process_census.py
python scripts/figure1.py
python scripts/figure2.py
python scripts/figure3.py
python scripts/figure4.py
python scripts/figureS1.py
python scripts/figureS2.py
python scripts/figureS3.py
```

## Compute

The permutation, bootstrap and resampling analyses behind the census ran on the
UAEU HPC cluster under Slurm. Their outputs are committed, so the figure scripts
here run on a laptop.

## Networks

The networks are deposited in NDEx: the full trans-regulatory census
([9c91cf3e-ac18-11f1-a428-005056ae3c32](https://www.ndexbio.org/viewer/networks/9c91cf3e-ac18-11f1-a428-005056ae3c32))
and the sign-reversal subnetwork
([9982abce-9586-11f1-bdc8-005056ae3c32](https://www.ndexbio.org/viewer/networks/9982abce-9586-11f1-bdc8-005056ae3c32)).

## Licence

Released under the MIT Licence. See [LICENSE](LICENSE) for the full text.
