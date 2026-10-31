"""Command-line interface: `biosim dock`, `redock`, `benchmark`, `plugins`."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import biosim.plugins  # noqa: F401
from biosim.core.plugin import list_plugins
from biosim.data.chembl import fetch_actives
from biosim.workflows.benchmark import run_benchmark
from biosim.workflows.docking import run_docking
from biosim.workflows.redock import run_redock


def _read_smiles_file(path: Path) -> list[tuple[str, str]]:
    ligands = []
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        ligands.append((parts[1] if len(parts) > 1 else f"lig{i}", parts[0]))
    return ligands


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="biosim")
    sub = parser.add_subparsers(dest="command", required=True)

    dock = sub.add_parser("dock", help="Dock a ligand set into a PDB structure")
    dock.add_argument("--pdb", required=True, help="PDB ID with a bound ligand, e.g. 1HSG")
    src = dock.add_mutually_exclusive_group(required=True)
    src.add_argument("--chembl", help="ChEMBL target ID to pull actives from, e.g. CHEMBL243")
    src.add_argument("--smiles-file", type=Path, help="File with lines: SMILES [name]")
    dock.add_argument("--n", type=int, default=20, help="Max ligands from ChEMBL")
    dock.add_argument("--ligand-resname", help="Native ligand residue name (default: largest)")
    dock.add_argument("--exhaustiveness", type=int, default=8)
    dock.add_argument("--out", type=Path, default=Path("runs/dock"))

    redock = sub.add_parser("redock", help="Re-dock one crystal ligand, report RMSD")
    redock.add_argument("--pdb", required=True, help="PDB ID with a bound ligand, e.g. 1HSG")
    redock.add_argument("--ligand-resname", help="Native ligand residue name (default: largest)")
    redock.add_argument(
        "--template-smiles",
        help="SMILES giving the ligand's bond orders (default: fetched from RCSB)",
    )
    redock.add_argument("--exhaustiveness", type=int, default=16)
    redock.add_argument("--n-poses", type=int, default=9)
    redock.add_argument("--out", type=Path, default=Path("runs/redock"))

    bench = sub.add_parser("benchmark", help="Run the redocking benchmark suite")
    bench.add_argument("--only", nargs="*", help="Run only these PDB IDs, e.g. --only 3PTB 1STP")
    bench.add_argument("--exhaustiveness", type=int, default=16)
    bench.add_argument("--out", type=Path, default=Path("runs/benchmark"))

    sub.add_parser("plugins", help="List registered plugins and their availability")

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.command == "plugins":
        for name, desc, missing in list_plugins():
            status = "OK" if not missing else f"missing: {', '.join(missing)}"
            print(f"{name:15s} {status:25s} {desc}")
        return 0

    if args.command == "benchmark":
        summary = run_benchmark(
            args.out, only=args.only, exhaustiveness=args.exhaustiveness
        )
        print(f"\n{'PDB':<5} {'lig':<4} {'top RMSD':>9} {'best RMSD':>10} {'result':<7} target")
        for r in summary["results"]:
            if r["status"] != "ok":
                result, top, best = "ERROR", "-", "-"
            else:
                result = "PASS" if r["passed"] else "FAIL"
                top, best = r["top_pose_rmsd"], r["best_rmsd_any_pose"]
            print(f"{r['pdb']:<5} {r['ligand']:<4} {top:>9} {best:>10} {result:<7} {r['target']}")
            if r["status"] != "ok":
                print(f"      {r['status']}")
        n = summary["n_total"]
        print(
            f"\nTop-pose success (<= {summary['threshold_angstrom']} A): "
            f"{summary['n_passed_top_pose']}/{n}"
        )
        print(f"Success with any of the poses: {summary['n_passed_any_pose']}/{n}")
        if summary["n_errors"]:
            print(f"Errors (counted as failures): {summary['n_errors']}")
        print(f"Median top-pose RMSD: {summary['median_top_pose_rmsd']} A")
        print(f"Details: {args.out / 'benchmark_results.csv'}")
        return 0

    if args.command == "redock":
        result = run_redock(
            args.pdb,
            args.out,
            ligand_resname=args.ligand_resname,
            template_smiles=args.template_smiles,
            exhaustiveness=args.exhaustiveness,
            n_poses=args.n_poses,
        )
        print(f"\nRedocking {result['ligand']} into {result['pdb']}")
        print(f"{'pose':>4}  {'vina':>8}  {'RMSD (A)':>9}")
        for p in result["poses"]:
            score = p["vina_score"] if p["vina_score"] is not None else "n/a"
            print(f"{p['pose']:>4}  {score:>8}  {p['rmsd']:>9}")
        verdict = "PASS" if result["success_top_pose"] else "FAIL"
        print(
            f"\nTop-pose RMSD {result['top_pose_rmsd']} A "
            f"(threshold {result['threshold']} A): {verdict}"
        )
        print(f"Best RMSD among all poses: {result['best_rmsd_any_pose']} A")
        print(f"Details: {args.out / 'redock_results.json'}")
        return 0

    if args.chembl:
        ligands = fetch_actives(args.chembl, limit=args.n)
    else:
        ligands = _read_smiles_file(args.smiles_file)
    if not ligands:
        print("No ligands to dock.", file=sys.stderr)
        return 1

    rows = run_docking(
        args.pdb, ligands, args.out,
        ligand_resname=args.ligand_resname,
        exhaustiveness=args.exhaustiveness,
    )

    print(f"\nTop results (full table: {args.out / 'results.csv'})")
    print(f"{'rank':>4}  {'score':>8}  name")
    for rank, row in enumerate(rows[:10], start=1):
        score = row["vina_score_kcal_mol"]
        print(f"{rank:>4}  {score if score is not None else 'n/a':>8}  {row['name']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
