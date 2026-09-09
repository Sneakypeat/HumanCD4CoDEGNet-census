#!/usr/bin/env bash
#SBATCH --job-name=broad_gsea_v2
#SBATCH --account=hpcuser
#SBATCH --partition=mediumq
#SBATCH --qos=limit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/home/<hpc-user>/projects/codeg_atlas/logs/broad-gsea-v2-%j.out
#SBATCH --error=/home/<hpc-user>/projects/codeg_atlas/logs/broad-gsea-v2-%j.err

set -euo pipefail
umask 022

: "${RUN_ROOT:?RUN_ROOT must name the checksum-sealed staged CoDEG_Tcell directory}"
: "${BROAD_GSEA_V2_EXPECTED_WRAPPER_SHA256:?expected wrapper SHA-256 is required}"
PROJECT_ROOT="$(realpath "$RUN_ROOT")"
PAPER_ROOT="$PROJECT_ROOT/paper"
PYTHON_BIN="${PYTHON_BIN:-/cm/shared/apps/anaconda3/bin/python3}"
MAIN="$PAPER_ROOT/scripts/broad_alltarget_kegg_go_rank_enrichment_v2.py"
AUDIT="$PAPER_ROOT/scripts/audit_broad_alltarget_kegg_go_rank_enrichment_v2.py"
TESTS="$PAPER_ROOT/tests/test_broad_alltarget_kegg_go_rank_enrichment_v2.py"
WRAPPER="$PAPER_ROOT/scripts/broad_alltarget_kegg_go_rank_enrichment_v2.sbatch"
DESIGN_MD="$PAPER_ROOT/docs/BROAD_ALLTARGET_KEGG_GO_RANK_ENRICHMENT_PREEXECUTION_V2_2026-09-04.md"
DESIGN_JSON="$PAPER_ROOT/results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/preexecution_v2.json"
NETWORK_SCRIPT="$PAPER_ROOT/scripts/fig3b_compound_target_pathway_randomization_v1.py"
SCORE_SCRIPT="$PAPER_ROOT/scripts/fig3e_alltarget_rank_enrichment_v1.py"
OUTPUT="$PAPER_ROOT/results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/production_v2"
CERTIFICATE="$PAPER_ROOT/results_data/census/broad_alltarget_kegg_go_rank_enrichment_v2/independent_audit_v2.json"

EXPECTED_MAIN_SHA256="e709551a1e46701b8ea927856de45afb1c4a9fe356f16b700aca0f175adbd851"
EXPECTED_AUDIT_SHA256="62539cd8541541ac543415f19ee6190a1a5bbce2bbe4ad429dc52bfd2e86d274"
EXPECTED_TESTS_SHA256="e6ba0fd6d9ed8625041eb841dd7a414c2f251cf3ce182588ce1ae14072c0565b"
EXPECTED_DESIGN_MD_SHA256="3456e8d0adbe165a585f0c2df05ad8b68edef923bf24f1cc2f88dc531c8723af"
EXPECTED_DESIGN_JSON_SHA256="070bb493128d366865cca40662812de46aa6aebe1f2d01b081ba69ba19558972"
EXPECTED_NETWORK_SCRIPT_SHA256="4a2f5710f4763f6bc6b7b3b8dcd14309bed56a65c8c0dff9652cdc661688c297"
EXPECTED_SCORE_SCRIPT_SHA256="bebc048e8a33c3dd8ef6ccc2df7af4a3f6f6896a7a60f2b620594581254fc0fe"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export PYTHONHASHSEED=0

test -n "${SLURM_JOB_ID:-}"
test "${SLURM_CPUS_PER_TASK:-}" = "1"
test "${SLURM_MEM_PER_NODE:-0}" -ge 16384
test -x "$PYTHON_BIN"
"$PYTHON_BIN" -c 'import sys,numpy,pandas,scipy; assert sys.version_info >= (3,10)'
test ! -e "$OUTPUT"
test ! -e "$CERTIFICATE"

assert_sha256() {
  local path="$1"
  local expected="$2"
  local observed
  observed="$("$PYTHON_BIN" -c 'import hashlib,sys; h=hashlib.sha256(); f=open(sys.argv[1], "rb"); [h.update(block) for block in iter(lambda: f.read(1048576), b"")]; print(h.hexdigest())' "$path")"
  if [[ "$observed" != "$expected" ]]; then
    echo "staged hash mismatch: $path: $observed != $expected" >&2
    exit 78
  fi
}

assert_sha256 "$MAIN" "$EXPECTED_MAIN_SHA256"
assert_sha256 "$AUDIT" "$EXPECTED_AUDIT_SHA256"
assert_sha256 "$TESTS" "$EXPECTED_TESTS_SHA256"
assert_sha256 "$DESIGN_MD" "$EXPECTED_DESIGN_MD_SHA256"
assert_sha256 "$DESIGN_JSON" "$EXPECTED_DESIGN_JSON_SHA256"
assert_sha256 "$NETWORK_SCRIPT" "$EXPECTED_NETWORK_SCRIPT_SHA256"
assert_sha256 "$SCORE_SCRIPT" "$EXPECTED_SCORE_SCRIPT_SHA256"
assert_sha256 "$PAPER_ROOT/results_data/census/census_master_edges.csv.gz" "0d753121d0e6e69313acd20b29e3431f01ad083a9f3515598f8e2d83775298ad"
assert_sha256 "$PAPER_ROOT/results_data/census/census_master_manifest.json" "1b4b6b7a74c255faf3ba3184b9f7fadf9171dc57167dfde4d9dd539b51dd0af5"
assert_sha256 "$PROJECT_ROOT/data/genesets/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt" "4a87c8d1260e637ec174f5960885519d46fa29da650e70d3a3641a4d41c4fa5e"
assert_sha256 "$PROJECT_ROOT/data/genesets/c2.cp.kegg_legacy.v2026.1.Hs.json" "44d32fd800957f30f8089740fc99c056867d6ffa1862a5e726c37eb9d4b9819f"
assert_sha256 "$PROJECT_ROOT/data/genesets/GO_Biological_Process_2023.gmt" "4be7e1b5ed62ab3a8a4837743fa8c9520fc572c5d28dcad8478ed136189908b9"
assert_sha256 "$PAPER_ROOT/results_data/census/fig3b_compound_target_pathway_randomization_freeze/production_v1/04_null_term_counts.tsv.gz" "28bb1388cf1de70bce8e4e58ee16b4be593560affd938349a9d9dde5ecc23bc7"

export BROAD_GSEA_V2_WRAPPER_PATH="$(realpath "$WRAPPER")"
BROAD_GSEA_V2_WRAPPER_SHA256="$("$PYTHON_BIN" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$WRAPPER")"
if [[ "$BROAD_GSEA_V2_WRAPPER_SHA256" != "$BROAD_GSEA_V2_EXPECTED_WRAPPER_SHA256" ]]; then
  echo "staged wrapper hash mismatch" >&2
  exit 78
fi
export BROAD_GSEA_V2_WRAPPER_SHA256

"$PYTHON_BIN" "$TESTS"
"$PYTHON_BIN" "$MAIN" --mode production --draws-per-bank 10000 --output "$OUTPUT"
"$PYTHON_BIN" "$AUDIT" --output "$OUTPUT" --certificate "$CERTIFICATE"
