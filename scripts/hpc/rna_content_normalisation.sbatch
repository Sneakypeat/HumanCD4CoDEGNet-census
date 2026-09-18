#!/bin/bash
#SBATCH --job-name=rna_content_normalisation
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/rna_content_normalisation.%j.out
#SBATCH --error=logs/rna_content_normalisation.%j.err

set -euo pipefail

: "${PERCELL_UMI_TABLE:?Set PERCELL_UMI_TABLE to the full-atlas per-cell total-UMI CSV.GZ}"

mkdir -p logs
python scripts/analysis/rna_content_full_atlas.py --percell "$PERCELL_UMI_TABLE"
python scripts/analysis/rna_content_direction_coherence.py --permutations 20000 --seed 20260911
