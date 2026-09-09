# HumanCD4CoDEGNet — census

Sign-reversal census of trans-regulatory effects during human CD4+ T-cell activation.

A regulator can raise a target's expression in one cellular state and lower it in another. This repository holds the census that measures how often that happens genome-wide: of 83,489 regulator–target interactions significant in both the resting and 48 h stimulated states, 4,379 (5.2%) reverse sign.

## Data Availability
The released census is provided as a signed, directed graph in `data/census/census_master_edges.csv.gz`. It carries the per-state log2 fold change, standard error and adjusted p-value of every edge.

*   **Raw Source Atlas:** (Zhu et al., genome-scale CRISPRi Perturb-seq of primary human CD4+ T cells) - Available at GEO/Zenodo [Insert Link Here]
*   **Third-party annotations:** The GO ontology, UniProt human GAF/GPI, and Reactome pathway GMT are available via Zenodo [Insert Link Here] or can be downloaded from their respective primary sources.

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

## Citation

If you use this code or the released census, please cite:

> Rehman SS, Muhammad K, Alam MT. Direction of trans-regulatory effects reverses
> genome-wide during human CD4+ T-cell activation. *Manuscript in preparation*.

The networks are also deposited in NDEx: the full trans-regulatory census
([9c91cf3e-ac18-11f1-a428-005056ae3c32](https://www.ndexbio.org/viewer/networks/9c91cf3e-ac18-11f1-a428-005056ae3c32))
and the sign-reversal subnetwork
([9982abce-9586-11f1-bdc8-005056ae3c32](https://www.ndexbio.org/viewer/networks/9982abce-9586-11f1-bdc8-005056ae3c32)).

## Licence

Released under the MIT Licence. See [LICENSE](LICENSE) for the full text.
