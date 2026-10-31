"""Diagnose the 4DFR (DHFR / methotrexate) redocking failure.

Runs the same redocking with different ligand protonation states and search
depths, over several seeds, and prints a comparison table.

    python scripts/diag_4dfr.py
    python scripts/diag_4dfr.py --seeds 42 7 123 --only A C

Variants
    A  neutral ligand, exhaustiveness 16   (current baseline)
    B  carboxylates charged, exhaustiveness 16   (what a pH 7.4 model gives)
    C  B plus pteridine N1-H+, exhaustiveness 16 (DHFR-bound methotrexate)
    D  neutral ligand, exhaustiveness 32   (tests search depth alone)
"""
from __future__ import annotations

import argparse
import json
import logging
import statistics
import time
from pathlib import Path

from biosim.workflows.redock import SUCCESS_RMSD, run_redock

PDB_ID = "4DFR"
LIGAND = "MTX"

# key -> (description, protonation state, exhaustiveness)
VARIANTS = {
    "A": ("neutral, exh 16 (baseline)", "neutral", 16),
    "B": ("carboxylates charged, exh 16", "anion", 16),
    "C": ("B + N1-H+ (bound state), exh 16", "bound", 16),
    "D": ("neutral, exh 32", "neutral", 32),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 7, 123, 2024, 99])
    parser.add_argument("--only", nargs="+", choices=sorted(VARIANTS), default=sorted(VARIANTS))
    parser.add_argument("--out", default="runs/diag_4dfr")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING)
    out_root = Path(args.out)
    summary = []

    for key in args.only:
        description, protonation, exhaustiveness = VARIANTS[key]
        print(f"\n=== Variant {key}: {description} ===", flush=True)
        tops, bests, passes = [], [], 0
        start = time.time()
        for seed in args.seeds:
            result = run_redock(
                PDB_ID,
                out_root / key / f"seed_{seed}",
                ligand_resname=LIGAND,
                exhaustiveness=exhaustiveness,
                n_poses=9,
                seed=seed,
                protonation=protonation,
            )
            tops.append(result["top_pose_rmsd"])
            bests.append(result["best_rmsd_any_pose"])
            passes += int(result["success_top_pose"])
            print(
                f"  seed {seed:>5}: top RMSD {result['top_pose_rmsd']:6.3f}  "
                f"best any pose {result['best_rmsd_any_pose']:6.3f}  "
                f"{'PASS' if result['success_top_pose'] else 'fail'}",
                flush=True,
            )
        summary.append(
            {
                "variant": key,
                "description": description,
                "seeds_passed": f"{passes}/{len(args.seeds)}",
                "median_top_rmsd": round(statistics.median(tops), 3),
                "best_any_pose": round(min(bests), 3),
                "seconds": round(time.time() - start, 1),
            }
        )

    print(f"\nSummary (pass = top pose within {SUCCESS_RMSD} A of the crystal pose)")
    print(f"{'Var':<4}{'Description':<36}{'Passed':<8}{'Median top':<12}{'Best any':<10}{'Time (s)'}")
    for row in summary:
        print(
            f"{row['variant']:<4}{row['description']:<36}{row['seeds_passed']:<8}"
            f"{row['median_top_rmsd']:<12}{row['best_any_pose']:<10}{row['seconds']}"
        )

    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "diag_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
