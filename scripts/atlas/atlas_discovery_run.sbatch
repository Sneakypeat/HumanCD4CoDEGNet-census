#!/bin/bash
# =============================================================================
# CORRECTED discovery arm re-run -- UAEU HPC (SLURM).
#
# Regenerates the withdrawn genome-scale state-shift map with the corrected null
# (per-(donor,condition,cell-count) floor + BH-FDR; atlas_discovery_corrected.py).
# This is CHEAP: it reuses the cached composition panel + score model, so there is
# NO re-fit and NO 1.13 TB Pass-B re-run. The occupancy pass streams each shard
# ONCE for per-guide state counts only (~27 MB/shard, tens of minutes each).
#
# PREREQUISITES (already produced by the composition run atlas_run.sbatch):
#   $WORK/atlas_state_panel.txt      (STAGE=gate)
#   $WORK/atlas_score_model.npz      (STAGE=fit; LODO-separation gate passed)
#   $DATA/*.assigned_guide.h5ad      (the 8 per-cell shards)
#
# STAGES (submit in order; each prints the next command):
#   occ        : array 0-7, lightweight state-count pass  -> $WORK/$NAME.occ.npz
#   discovery  : single job, corrected floor+FDR+overlap  -> $WORK/atlas_discovery_corrected.json
#
#   sbatch --array=0-7 --export=STAGE=occ       atlas_discovery_run.sbatch
#   sbatch             --export=STAGE=discovery atlas_discovery_run.sbatch
#
# Then pull $WORK/atlas_discovery_corrected.json + atlas_state_shift_map_corrected.csv.gz
# local and regenerate Fig 6d-f only if n_shifters>0 AND the overlap clears
# (cell_count_matched_permutation_p < 0.05). See DISCOVERY_ARM_RULES.md.
# =============================================================================
#SBATCH --job-name=atlas_discovery
#SBATCH --partition=mediumq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/<hpc-user>/projects/codeg_atlas/logs/%x_%A_%a.out
#SBATCH --error=/home/<hpc-user>/projects/codeg_atlas/logs/%x_%A_%a.err

set -euo pipefail

ROOT=/home/<hpc-user>/projects/codeg_atlas
DATA=$ROOT/data
WORK=$ROOT/work
LOGS=$ROOT/logs
SCRIPTS=$HOME/CoDEG_Tcell/paper/scripts/atlas
CENSUS=$HOME/CoDEG_Tcell/paper/kinetics_8hr_edges.csv.gz
PANEL=$WORK/atlas_state_panel.txt
MODEL=$WORK/atlas_score_model.npz
mkdir -p "$WORK" "$LOGS"

module load anaconda3
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SHARDS=(D1_Rest D1_Stim48hr D2_Rest D2_Stim48hr D3_Rest D3_Stim48hr D4_Rest D4_Stim48hr)

case "${STAGE:?set STAGE=occ|discovery}" in

occ)
  NAME=${SHARDS[${SLURM_ARRAY_TASK_ID:?occ needs --array=0-7}]}
  [[ -s $PANEL ]] || { echo "panel missing: run composition STAGE=gate first"; exit 1; }
  [[ -s $MODEL ]] || { echo "score model missing: run composition STAGE=fit first"; exit 1; }
  [[ -f $DATA/$NAME.assigned_guide.h5ad ]] || { echo "MISSING shard: $NAME"; exit 1; }
  exec 9>"$WORK/.lock.$NAME.occ"; flock -n 9 || { echo "$NAME already running"; exit 0; }
  python3 "$SCRIPTS/atlas_occupancy.py" \
    --h5ad        "$DATA/$NAME.assigned_guide.h5ad" \
    --panel       "$PANEL" --census "$CENSUS" \
    --score-model "$MODEL" \
    --out         "$WORK/$NAME.occ.npz"
  echo "occ done: $NAME  (when all 8 exist -> sbatch --export=STAGE=discovery $0)"
  ;;

discovery)
  n=$(ls "$WORK"/*.occ.npz 2>/dev/null | wc -l | tr -d ' ')
  [[ "$n" == "8" ]] || { echo "need all 8 occ.npz, have $n -- run STAGE=occ (array 0-7) first"; exit 1; }
  python3 "$SCRIPTS/atlas_discovery_corrected.py" \
    --occ     "$WORK"/*.occ.npz \
    --census  "$CENSUS" \
    --out     "$WORK/atlas_discovery_corrected.json" \
    --map-out "$WORK/atlas_state_shift_map_corrected.csv.gz"
  echo "CORRECTED DISCOVERY -> $WORK/atlas_discovery_corrected.json"
  echo "pull local; restore Fig 6d-f ONLY if overlap cell-count-matched p<0.05"
  ;;

*) echo "unknown STAGE=$STAGE"; exit 1;;
esac
