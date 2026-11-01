"""Break the 4DFR redocking error into methotrexate fragments.

Reads the runs written by scripts/diag_4dfr.py (no docking is repeated) and
reports, for the top-ranked pose and the pose closest to the crystal, the
heavy-atom RMSD of each fragment against the crystal ligand, in place (no
alignment), with symmetric atoms (carboxylate O, phenyl ring) handled.

    python scripts/analyze_4dfr_poses.py
    python scripts/analyze_4dfr_poses.py --variants B --root runs/diag_4dfr

Reading the table
    pteridine small, glutamate large -> head is right, only the floppy tail is off
    pteridine large                  -> the whole ring is misplaced or flipped
"""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

from ligentry.workflows.redock import _crystal_ligand_atoms, _pdb_block, _pose_models

RDLogger.DisableLog("rdApp.warning")

# Bond-order-agnostic (~) so connectivity from proximity bonding is enough.
FRAGMENTS = {
    "pteridine": "[#6]1~[#6]~[#7]~[#6]2~[#7]~[#6](~[#7])~[#7]~[#6](~[#7])~[#6]2~[#7]~1",
    "benzoyl": "[#7]~[#6]1~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[#6](~[#8])~[#7]",
    "glutamate": "[#8]~[#6](~[#8])~[#6](~[#7])~[#6]~[#6]~[#6](~[#8])~[#8]",
}
_QUERIES = {name: Chem.MolFromSmarts(smarts) for name, smarts in FRAGMENTS.items()}


def build_mol(atoms):
    """Connectivity-only molecule from (element, x, y, z) tuples."""
    mol = Chem.MolFromPDBBlock(
        _pdb_block(atoms), proximityBonding=True, removeHs=False, sanitize=False
    )
    if mol is None:
        raise ValueError("RDKit could not build a molecule from the coordinates")
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    return mol


def fragment_rmsd(pose_mol, crystal_mol, name: str) -> float | None:
    query = _QUERIES[name]
    ref = crystal_mol.GetSubstructMatch(query)
    matches = pose_mol.GetSubstructMatches(query, uniquify=False, maxMatches=1000)
    if not ref or not matches:
        return None
    ref_xyz = crystal_mol.GetConformer().GetPositions()[list(ref)]
    pose_xyz = pose_mol.GetConformer().GetPositions()
    return float(
        min(np.sqrt(((pose_xyz[list(m)] - ref_xyz) ** 2).sum(axis=1).mean()) for m in matches)
    )


def analyse_seed(seed_dir: Path, ligand: str) -> dict:
    result = json.loads((seed_dir / "redock_results.json").read_text())
    pdb_files = sorted(seed_dir.glob("*.pdb"))
    pdb_file = next(p for p in pdb_files if not p.stem.endswith("_protein"))
    crystal = build_mol(_crystal_ligand_atoms(pdb_file, ligand))
    models, _scores = _pose_models(str(seed_dir / "poses" / f"{ligand}_redock.pdbqt"))

    whole = [p["rmsd"] for p in result["poses"]]
    closest = int(np.argmin(whole))
    out = {}
    for label, idx in (("top", 0), ("closest", closest)):
        pose = build_mol(models[idx])
        row = {"pose": idx + 1, "whole": whole[idx]}
        for name in FRAGMENTS:
            row[name] = fragment_rmsd(pose, crystal, name)
        out[label] = row
    return out


def fmt(value) -> str:
    return "  n/a" if value is None else f"{value:5.2f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default="runs/diag_4dfr")
    parser.add_argument("--variants", nargs="+", default=["A", "B", "C", "D"])
    parser.add_argument("--ligand", default="MTX")
    args = parser.parse_args()

    for variant in args.variants:
        seed_dirs = sorted((Path(args.root) / variant).glob("seed_*"))
        if not seed_dirs:
            print(f"\nVariant {variant}: no runs found under {args.root}/{variant}")
            continue
        print(f"\n=== Variant {variant} ===")
        print(f"{'seed':<8}{'pose':<9}{'#':<4}{'whole':>6}{'pteridine':>11}{'benzoyl':>9}{'glutamate':>11}")
        collected = {"top": [], "closest": []}
        for seed_dir in seed_dirs:
            rows = analyse_seed(seed_dir, args.ligand)
            for label, row in rows.items():
                collected[label].append(row)
                print(
                    f"{seed_dir.name[5:]:<8}{label:<9}{row['pose']:<4}{row['whole']:>6.2f}"
                    f"{fmt(row['pteridine']):>11}{fmt(row['benzoyl']):>9}{fmt(row['glutamate']):>11}"
                )
        for label, rows in collected.items():
            med = {
                key: statistics.median(r[key] for r in rows if r[key] is not None)
                for key in ("whole", *FRAGMENTS)
            }
            print(
                f"{'median':<8}{label:<9}{'':<4}{med['whole']:>6.2f}"
                f"{med['pteridine']:>11.2f}{med['benzoyl']:>9.2f}{med['glutamate']:>11.2f}"
            )


if __name__ == "__main__":
    main()
