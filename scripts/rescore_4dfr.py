"""Is the 4DFR failure a scoring problem or a sampling problem?

For each saved run this puts the crystal-pose ligand into Vina, scores it and
locally minimises it, then does the same for the docked top-ranked pose and the
docked pose closest to the crystal. Everything is scored with the same Vina
call, so the numbers are directly comparable.

    python scripts/rescore_4dfr.py
    python scripts/rescore_4dfr.py --variants A B --seeds 42 7

Reading the result (energies in kcal/mol, more negative = better)
    native minimised is WORSE than the docked top pose
        -> Vina prefers the wrong pose: a scoring / receptor-model limit, and more
           sampling cannot fix it.
    native minimised is BETTER than the docked top pose
        -> the right basin is more favourable but the search misses it: look at
           the box, ligand preparation and torsion handling.
    "moved" is how far minimisation shifts a pose from where it started (RMSD, A).
    For the crystal pose a large value means Vina does not treat the experimental
    pose as a stable minimum.
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

from ligentry.data.pdb import extract_pocket
from ligentry.workflows.redock import _crystal_ligand_atoms, _mol_from_atoms, _pose_models

RDLogger.DisableLog("rdApp.warning")


def split_models(pdbqt_path: Path) -> list[str]:
    """Each MODEL of a Vina output file as a standalone ligand PDBQT string."""
    models: list[str] = []
    current: list[str] | None = None
    for line in Path(pdbqt_path).read_text().splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            if current:
                models.append("\n".join(current) + "\n")
            current = None
        elif current is not None:
            current.append(line)
    return models


def native_pdbqt(crystal_atoms, template) -> str:
    """Crystal-pose ligand with the docking chemistry, hydrogens placed by geometry."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = _mol_from_atoms(crystal_atoms, template)
    mol = Chem.AddHs(mol, addCoords=True)
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"Meeko failed on the crystal ligand: {err}")
    return pdbqt


def score_and_minimise(vina, pdbqt: str, center, size, template) -> dict:
    """Vina score of the pose as given, then after local minimisation."""
    vina.set_ligand_from_string(pdbqt)
    vina.compute_vina_maps(center=list(center), box_size=list(size))
    scored = vina.score()
    minimised = vina.optimize()

    with tempfile.TemporaryDirectory() as tmp:
        start_file = Path(tmp) / "start.pdbqt"
        start_file.write_text(pdbqt)
        end_file = Path(tmp) / "min.pdbqt"
        vina.write_pose(str(end_file), overwrite=True)
        start_atoms = _pose_models(str(start_file))[0][0]
        end_atoms = _pose_models(str(end_file))[0][0]
    # displacement caused by minimisation (symmetry-aware, in place)
    moved = float(
        rdMolAlign.CalcRMS(
            _mol_from_atoms(end_atoms, template), _mol_from_atoms(start_atoms, template)
        )
    )
    return {
        "score": float(scored[0]),
        "minimised": float(minimised[0]),
        "inter": float(minimised[1]),
        "moved": moved,
    }


def analyse_seed(seed_dir: Path, ligand: str, seed: int) -> list[tuple[str, dict]]:
    from vina import Vina

    result = json.loads((seed_dir / "redock_results.json").read_text())
    template = Chem.MolFromSmiles(result["ligand_smiles"])
    pdb_file = next(p for p in sorted(seed_dir.glob("*.pdb")) if not p.stem.endswith("_protein"))
    pocket = extract_pocket(pdb_file, seed_dir, ligand)
    receptor = next(seed_dir.glob("*_protein.pdbqt"))

    crystal_atoms = _crystal_ligand_atoms(pdb_file, ligand)
    models = split_models(seed_dir / "poses" / f"{ligand}_redock.pdbqt")
    rmsds = [p["rmsd"] for p in result["poses"]]
    closest = int(np.argmin(rmsds))

    vina = Vina(sf_name="vina", seed=seed, verbosity=0)
    vina.set_receptor(str(receptor))

    jobs = [
        ("crystal pose", native_pdbqt(crystal_atoms, template), 0.0),
        ("docked top (#1)", models[0], rmsds[0]),
    ]
    if closest != 0:
        jobs.append((f"docked closest (#{closest + 1})", models[closest], rmsds[closest]))

    rows = []
    for label, pdbqt, rmsd in jobs:
        info = score_and_minimise(vina, pdbqt, pocket.center, pocket.size, template)
        info["rmsd"] = rmsd
        rows.append((label, info))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--root", default="runs/diag_4dfr")
    parser.add_argument("--variants", nargs="+", default=["A", "B"])
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 7, 123, 2024, 99])
    parser.add_argument("--ligand", default="MTX")
    args = parser.parse_args()

    for variant in args.variants:
        print(f"\n=== Variant {variant} ===")
        gaps = []
        for seed in args.seeds:
            seed_dir = Path(args.root) / variant / f"seed_{seed}"
            if not (seed_dir / "redock_results.json").exists():
                print(f"seed {seed}: no run found in {seed_dir}")
                continue
            rows = analyse_seed(seed_dir, args.ligand, seed)
            print(f"\nseed {seed}")
            print(f"  {'pose':<22}{'RMSD':>6}{'scored':>9}{'minimised':>11}{'inter':>8}{'moved':>7}")
            for label, r in rows:
                print(
                    f"  {label:<22}{r['rmsd']:>6.2f}{r['score']:>9.2f}"
                    f"{r['minimised']:>11.2f}{r['inter']:>8.2f}{r['moved']:>7.2f}"
                )
            native, top = rows[0][1], rows[1][1]
            gap = native["minimised"] - top["minimised"]  # > 0: native is worse
            gaps.append(gap)
            verdict = "native WORSE (scoring prefers wrong pose)" if gap > 0 else "native BETTER"
            print(f"  native minus docked top, minimised: {gap:+.2f} kcal/mol  -> {verdict}")
        if gaps:
            worse = sum(g > 0 for g in gaps)
            print(f"\nVariant {variant}: native pose scores worse than the docked top pose "
                  f"in {worse}/{len(gaps)} seeds (median gap {np.median(gaps):+.2f} kcal/mol)")


if __name__ == "__main__":
    main()
