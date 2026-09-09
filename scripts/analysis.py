#!/usr/bin/env python3
"""Recompute the frozen analysis results the figures read.

Each section below is one previously separate analysis script, kept in the
order the results depend on each other. Running this regenerates the JSON
files under data/ that the figure scripts consume, so a reviewer can verify
the numbers rather than take the committed results on trust.

Permutation-based sections carry their own seed and permutation count; those
are recorded in the JSON each writes. Results derived from the 39 GB source
atlas or the cluster runs are NOT recomputed here - they cannot be, from the
data this repository redistributes - and remain as committed inputs.

Run before the figure scripts, or not at all: the committed results are
identical to what this produces.
"""
from __future__ import annotations

import runpy
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECTIONS = [
    ("signflip_feasibility", "regulator eligibility gate"),
    ("prep_paper_data", "figure_data.json, the shared per-figure inputs"),
    ("burden_heterogeneity", "reversal-burden dispersion across regulators"),
    ("degree_distribution_test", "out-degree distribution and heavy-tail fits"),
    ("convergence_test", "target convergence against the degree-preserving null"),
    ("regulator_target_structure", "regulator/target structure statistics"),
    ("target_theory_falsification", "target-conditioned falsification tests"),
    ("target_theory_decisive", "decisive target-orientation test"),
    ("arce_validation_export", "Arce et al. recurrence validation"),
    ("arce_adt_umap", "Arce protein-channel UMAP statistics"),
]


def main() -> int:
    failed = []
    for name, blurb in SECTIONS:
        script = ROOT / "scripts/analysis" / f"{name}.py"
        if not script.is_file():
            print(f"[analysis] MISSING {name}.py")
            failed.append(name)
            continue
        print(f"[analysis] {name}: {blurb}")
        started = time.time()
        argv, cwd = sys.argv[:], Path.cwd()
        try:
            sys.argv = [str(script)]
            runpy.run_path(str(script), run_name="__main__")
        except SystemExit as exc:
            if exc.code:
                failed.append(name)
        except Exception as exc:  # report and continue; one failure is not fatal
            print(f"[analysis]   FAILED: {type(exc).__name__}: {exc}")
            failed.append(name)
        finally:
            sys.argv = argv
        print(f"[analysis]   done in {time.time() - started:.1f}s")
    if failed:
        print(f"[analysis] {len(failed)} section(s) failed: {', '.join(failed)}")
        return 1
    print(f"[analysis] all {len(SECTIONS)} sections completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
