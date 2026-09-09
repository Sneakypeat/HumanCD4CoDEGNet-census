#!/bin/bash
# =============================================================================
# Atlas state-decomposition run -- UAEU HPC (SLURM).
#
# STAGES (submit in order; each prints the next command):
#   0  gate     : atlas_state_panel.py            -- HALTS if the panel is unusable
#   1  passA    : array 0-7, NTC cache per shard  (~2% of bytes)
#   2  fit      : cross-fit the LODO score        (single job, needs all 8 caches)
#   3  passB    : array 0-7, THE full pass        (~1.13 TB, one stream)
#   4  gather   : decompose + inference           (single job)
#
#   sbatch --export=STAGE=gate    atlas_run.sbatch
#   sbatch --export=STAGE=passA   atlas_run.sbatch
#   sbatch --export=STAGE=fit     atlas_run.sbatch
#   sbatch --export=STAGE=passB   atlas_run.sbatch
#   sbatch --export=STAGE=gather  atlas_run.sbatch
#
# Sizing is measured, not assumed: peak ~15 GB resident vs 503 GiB/node;
# <250 GB intermediates vs 136 TB free; ~30-60 min/shard vs mediumq's 10 d.
# mediumq (99 idle nodes) is used for the heavy step -- there is no queue cost
# and it removes any walltime cliff if the chunk layout surprises us.
# =============================================================================
#SBATCH --job-name=atlas_decomp
#SBATCH --partition=mediumq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
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
# Byte sizes VERIFIED against the HDF5 superblock (stored_eof), not the manifest.
declare -A SIZES=(
  [D1_Rest]=153299423232      [D1_Stim48hr]=167712189440
  [D2_Rest]=140574804978      [D2_Stim48hr]=147749289984
  [D3_Rest]=148263075840      [D3_Stim48hr]=147963084800
  [D4_Rest]=127330631680      [D4_Stim48hr]=170928046080
)

# -----------------------------------------------------------------------------
verify_shard () {
  # Size + HDF5 readability. NOTE: the download script's check was `... || echo
  # "WARNING: not readable HDF5"` -- it PRINTS and exits 0, so a right-size but
  # corrupt file reports success. Here a bad file EXITS 1 (failure mode G).
  local name=$1 file=$DATA/$name.assigned_guide.h5ad
  [[ -f $file ]] || { echo "MISSING: $file"; exit 1; }
  local actual expected
  actual=$(stat -c %s "$file")
  expected=${SIZES[$name]:-0}
  if [[ $expected -gt 0 && $actual -ne $expected ]]; then
    echo "SIZE MISMATCH $name: got $actual expected $expected"; exit 1
  fi
  python3 - "$file" <<'PY' || { echo "NOT READABLE HDF5: $file"; exit 1; }
import sys, h5py
with h5py.File(sys.argv[1], "r") as h:
    enc = dict(h["X"].attrs).get("encoding-type", b"")
    enc = enc.decode() if isinstance(enc, bytes) else enc
    assert enc == "csr_matrix", f"X encoding {enc!r}"
    assert h["X/indptr"].shape[0] > 1
    n_var = h["var"].attrs.get("_index", "_index")
    n_var = n_var.decode() if isinstance(n_var, bytes) else n_var
    assert h[f"var/{n_var}"].shape[0] == 18130, "n_var != 18130"
print("HDF5 OK")
PY
  echo "VERIFIED $name ($actual bytes)"
}

# -----------------------------------------------------------------------------
case "${STAGE:?set STAGE=gate|passA|fit|passB|gather}" in

gate)
  verify_shard D1_Rest
  verify_shard D1_Stim48hr
  python3 "$SCRIPTS/atlas_state_panel.py" \
    --rest    "$DATA/D1_Rest.assigned_guide.h5ad" \
    --stim    "$DATA/D1_Stim48hr.assigned_guide.h5ad" \
    --census  "$CENSUS" \
    --h5mu    "$HOME/projects/codeg_donorsplit/data/GWCD4i.DE_stats.by_guide.h5mu" \
              "$HOME/projects/codeg_donorsplit/data/GWCD4i.DE_stats.by_donors.h5mu" \
    --out       "$WORK/atlas_state_panel.json" \
    --panel-out "$PANEL"
  echo "GATE PASSED -> sbatch --array=0-7 --export=STAGE=passA $0"
  ;;

passA)
  #SBATCH --array=0-7
  NAME=${SHARDS[${SLURM_ARRAY_TASK_ID:?passA needs --array=0-7}]}
  [[ -s $PANEL ]] || { echo "panel missing: run STAGE=gate first"; exit 1; }
  verify_shard "$NAME"
  exec 9>"$WORK/.lock.$NAME.A"; flock -n 9 || { echo "$NAME already running"; exit 0; }
  python3 "$SCRIPTS/atlas_extract.py" --pass A \
    --h5ad   "$DATA/$NAME.assigned_guide.h5ad" \
    --panel  "$PANEL" --census "$CENSUS" \
    --out    "$WORK/$NAME.ntcA.npz"
  echo "passA done: $NAME"
  ;;

fit)
  ls "$WORK"/*.ntcA.npz | wc -l | grep -qx 8 || { echo "need all 8 Pass-A caches"; exit 1; }
  python3 "$SCRIPTS/atlas_decompose.py" fit \
    --ntc    "$WORK"/*.ntcA.npz \
    --out    "$MODEL" \
    --report "$WORK/atlas_score_model.json"
  echo "fit done -> sbatch --array=0-7 --export=STAGE=passB $0"
  ;;

passB)
  #SBATCH --array=0-7
  NAME=${SHARDS[${SLURM_ARRAY_TASK_ID:?passB needs --array=0-7}]}
  # FATAL 3b fix: gate the 1.13 TB pass on the LODO separation RESULT, not just on
  # the model file existing. The model .npz is now written only on gate pass, and
  # the JSON report only exists on a successful fit -- check both, and the delta.
  python3 -c "import json,sys; d=json.load(open('$WORK/atlas_score_model.json')); sys.exit(0 if d.get('worst_donor_cliffs_delta',-1)>=0.60 else 1)" \
    || { echo "LODO separation gate NOT passed (delta<0.60) or fit did not complete: rerun STAGE=fit"; exit 1; }
  [[ -s $MODEL ]] || { echo "score model .npz missing despite gate pass -- inconsistent state, halt"; exit 1; }
  verify_shard "$NAME"
  exec 9>"$WORK/.lock.$NAME.B"; flock -n 9 || { echo "$NAME already running"; exit 0; }
  python3 "$SCRIPTS/atlas_extract.py" --pass B \
    --h5ad        "$DATA/$NAME.assigned_guide.h5ad" \
    --panel       "$PANEL" --census "$CENSUS" \
    --score-model "$MODEL" \
    --out         "$WORK/$NAME.sufstat.npz"
  echo "passB done: $NAME"
  ;;

gather)
  ls "$WORK"/*.sufstat.npz | wc -l | grep -qx 8 || { echo "need all 8 Pass-B outputs"; exit 1; }

  # NEGATIVE CONTROL FIRST. A random score preserves every input size and
  # destroys all state meaning. If it still finds state shifts, the pipeline is
  # measuring input size and NOTHING downstream is interpretable.
  python3 "$SCRIPTS/atlas_decompose.py" decompose \
    --shards "$WORK"/*.sufstat.npz --census "$CENSUS" \
    --out    "$WORK/atlas_zeroinfo.npz" \
    --report "$WORK/atlas_zeroinfo.json" \
    --zero-information-state

  python3 "$SCRIPTS/atlas_decompose.py" decompose \
    --shards "$WORK"/*.sufstat.npz --census "$CENSUS" \
    --out    "$WORK/atlas_decomposition.npz" \
    --report "$WORK/atlas_decompose.json" \
    --verify-exactness 50

  python3 "$SCRIPTS/atlas_inference.py" \
    --units   "$WORK/atlas_decomposition_units.csv.gz" \
    --arrays  "$WORK/atlas_decomposition.npz" \
    --census  "$CENSUS" \
    --mu0     "$WORK"/atlas_decomposition.mu0_*.npz \
    --zero-information-run "$WORK/atlas_zeroinfo_units.csv.gz" \
    --out     "$WORK/atlas_results.json" \
    --map-out "$WORK/atlas_state_shift_map.csv.gz"
  echo "RESULTS -> $WORK/atlas_results.json"
  ;;

*) echo "unknown STAGE=$STAGE"; exit 1;;
esac
