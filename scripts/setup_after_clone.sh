#!/bin/bash
# Restore filesystem permissions that Git cannot record.
#
# Figure 4's builder (scripts/response_shape/fig4_production_v5_1.py) treats its
# frozen production bundle as immutable and asserts the exact modes below before
# it will run:
#
#   Figure4ContractError: production result stage must be a non-symlink 0550 directory
#
# Git stores only the executable bit, so a fresh clone gets 0755/0644 and the
# contract fails. Run this once after cloning.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE="$ROOT/results_data/response_shape/production_v5_1"

if [ -d "$STAGE" ]; then
  chmod 0440 "$STAGE"/* 2>/dev/null || true
  chmod 0550 "$STAGE"
  echo "restored 0550/0440 on $STAGE"
else
  echo "note: $STAGE not present; skipping" >&2
fi
