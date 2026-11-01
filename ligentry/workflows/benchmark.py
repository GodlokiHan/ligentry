"""Redocking benchmark suite: re-dock crystal ligands across diverse targets.

Success criterion per run: top-ranked pose within 2.0 A heavy-atom RMSD of the
crystal pose. Each complex is docked with one or more seeds; a complex passes
when the top pose passes in a majority of seeds. To extend the suite, add a row
to BENCHMARK_SET.
"""
from __future__ import annotations

import csv
import json
import logging
import statistics
import time
from pathlib import Path

from rdkit import RDLogger

from ligentry.workflows.redock import SUCCESS_RMSD, run_redock

RDLogger.DisableLog("rdApp.warning")  # silence "More than one matching pattern" noise

log = logging.getLogger("ligentry.benchmark")

# (PDB ID, ligand residue name, description)
BENCHMARK_SET: list[tuple[str, str, str]] = [
    ("1HSG", "MK1", "HIV-1 protease / indinavir"),
    ("3PTB", "BEN", "Trypsin / benzamidine"),
    ("1STP", "BTN", "Streptavidin / biotin"),
    ("1M17", "AQ4", "EGFR kinase / erlotinib"),
    ("3ERT", "OHT", "Estrogen receptor alpha / 4-hydroxytamoxifen"),
    ("4DFR", "MTX", "DHFR / methotrexate"),
    ("1KV1", "BMU", "p38 MAP kinase / pyrazole urea (BMU)"),
    ("1EVE", "E20", "Acetylcholinesterase / donepezil"),
]

CSV_FIELDS = [
    "pdb", "ligand", "target", "status", "seeds", "seed_passes", "top_pose_rmsd",
    "best_rmsd_any_pose", "top_vina_score", "passed", "seconds",
]


def run_benchmark(
    out_dir: str | Path,
    only: list[str] | None = None,
    exhaustiveness: int = 16,
    n_poses: int = 9,
    seeds: list[int] | tuple[int, ...] = (42,),
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seeds = list(seeds) or [42]

    selected = [e for e in BENCHMARK_SET if not only or e[0] in {p.upper() for p in only}]
    if not selected:
        raise ValueError(f"No benchmark entries match {only}. Available: "
                         f"{[e[0] for e in BENCHMARK_SET]}")

    rows: list[dict] = []
    for i, (pdb_id, resname, target) in enumerate(selected, start=1):
        log.info("[%d/%d] %s  %s  (%s)", i, len(selected), pdb_id, resname, target)
        row = {
            "pdb": pdb_id, "ligand": resname, "target": target, "status": "ok",
            "seeds": len(seeds), "seed_passes": 0, "top_pose_rmsd": None,
            "best_rmsd_any_pose": None, "top_vina_score": None, "passed": False,
            "seconds": None,
        }
        start = time.time()
        try:
            results = []
            for seed in seeds:
                results.append(
                    run_redock(
                        pdb_id,
                        out_dir / pdb_id / f"seed_{seed}",
                        ligand_resname=resname,
                        exhaustiveness=exhaustiveness,
                        n_poses=n_poses,
                        seed=seed,
                    )
                )
            tops = [r["top_pose_rmsd"] for r in results]
            top_scores = [
                r["poses"][0]["vina_score"] for r in results
                if r["poses"][0]["vina_score"] is not None
            ]
            row["top_pose_rmsd"] = round(statistics.median(tops), 3)
            row["best_rmsd_any_pose"] = min(r["best_rmsd_any_pose"] for r in results)
            row["top_vina_score"] = round(statistics.median(top_scores), 3) if top_scores else None
            row["seed_passes"] = sum(1 for r in results if r["success_top_pose"])
            row["passed"] = row["seed_passes"] * 2 > len(seeds)
        except Exception as exc:  # record and keep going
            row["status"] = f"error: {exc}"
            log.warning("  %s failed: %s", pdb_id, exc)
        row["seconds"] = round(time.time() - start, 1)
        rows.append(row)

    completed = [r for r in rows if r["status"] == "ok"]
    top_rmsds = [r["top_pose_rmsd"] for r in completed]
    n_runs = len(completed) * len(seeds)
    summary = {
        "threshold_angstrom": SUCCESS_RMSD,
        "exhaustiveness": exhaustiveness,
        "seeds": seeds,
        "n_total": len(rows),
        "n_completed": len(completed),
        "n_errors": len(rows) - len(completed),
        "n_passed_top_pose": sum(1 for r in rows if r["passed"]),
        "n_passed_any_pose": sum(
            1 for r in completed if r["best_rmsd_any_pose"] <= SUCCESS_RMSD
        ),
        "seed_pass_rate": (
            round(sum(r["seed_passes"] for r in completed) / n_runs, 3) if n_runs else None
        ),
        "median_top_pose_rmsd": round(statistics.median(top_rmsds), 3) if top_rmsds else None,
        "results": rows,
    }

    with open(out_dir / "benchmark_results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    (out_dir / "benchmark_summary.json").write_text(json.dumps(summary, indent=2))

    return summary
