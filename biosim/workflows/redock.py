"""Redocking benchmark: re-dock a crystal ligand and measure RMSD to the crystal pose.

The ligand is rebuilt from scratch (random 3D conformer), so the crystal
coordinates never leak into the docking. RMSD <= 2.0 A for the top-ranked pose
is the usual success criterion.
"""
from __future__ import annotations

import json
import logging
import re
from collections import Counter
from pathlib import Path

import requests
from rdkit import RDLogger

import biosim.plugins  # noqa: F401  (registers plugins)
from biosim.core.plugin import get_plugin
from biosim.data.pdb import download_pdb, extract_pocket

# Best-effort silencing of RDKit's harmless "More than one matching pattern" warnings.
RDLogger.DisableLog("rdApp.warning")
logging.getLogger("rdkit").setLevel(logging.ERROR)

log = logging.getLogger("biosim.redock")

SUCCESS_RMSD = 2.0  # Angstrom

# AutoDock atom types -> chemical elements
_AD_TO_ELEMENT = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "O": "O", "OA": "O",
    "OS": "O", "S": "S", "SA": "S", "F": "F", "Cl": "Cl", "CL": "Cl",
    "Br": "Br", "BR": "Br", "I": "I", "P": "P", "B": "B", "Si": "Si",
    "HD": "H", "H": "H", "HS": "H",
}


def _fetch_ccd_smiles(resname: str) -> str:
    """Get the ideal SMILES for a PDB ligand code from the RCSB component dictionary."""
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    descriptors = data.get("rcsb_chem_comp_descriptor") or {}
    lowered = {k.lower(): v for k, v in descriptors.items()}
    for key in ("smiles_stereo", "smiles"):
        if lowered.get(key):
            return lowered[key]

    for item in data.get("pdbx_chem_comp_descriptor") or []:
        if "SMILES" in str(item.get("type", "")).upper() and item.get("descriptor"):
            return item["descriptor"]

    raise ValueError(
        f"No SMILES found for ligand '{resname}'. Pass one explicitly with --template-smiles."
    )


def _template_from_smiles(smiles: str):
    """Parse a template SMILES and strip every hydrogen, including stereo-defining ones.

    RCSB SMILES such as [H]/N=C(/N)c1ccccc1 keep an explicit H to encode C=N
    geometry. RDKit preserves such an H as a real atom by default, which breaks
    heavy-atom matching against the crystal coordinates.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse template SMILES: {smiles}")
    params = Chem.RemoveHsParameters()
    params.removeDefiningBondStereo = True
    return Chem.RemoveHs(mol, params)


def _crystal_ligand_atoms(pdb_path: Path, resname: str) -> list[tuple[str, float, float, float]]:
    """Heavy atoms of the first copy of the named ligand (first alternate location only)."""
    groups: dict[tuple[str, str], list[tuple[str, float, float, float]]] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != resname.upper():
            continue
        if line[16] not in (" ", "A"):
            continue
        element = line[76:78].strip().capitalize() if len(line) >= 78 else ""
        if not element:
            element = re.sub(r"[^A-Za-z]", "", line[12:16])[:1].upper()
        if element == "H":
            continue
        key = (line[21], line[22:26].strip())
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        groups.setdefault(key, []).append((element, *xyz))
    if not groups:
        raise ValueError(f"Ligand '{resname}' not found in {pdb_path}")
    return max(groups.values(), key=len)


def _pdb_block(atoms: list[tuple[str, float, float, float]]) -> str:
    lines = []
    for i, (el, x, y, z) in enumerate(atoms, start=1):
        name = f"{el}{i}"[:4]
        lines.append(
            f"HETATM{i:5d} {name:<4s} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}"
        )
    lines.append("END")
    return "\n".join(lines)


def _mol_from_atoms_geometric(atoms, template, min_bond: float = 1.0, max_bond: float = 1.95):
    """Fallback: map template atoms onto coordinates by element and bond length.

    Does not rely on RDKit's inferred bonds. Returns a copy of the template with
    the mapped coordinates as its conformer.
    """
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    n = template.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"{len(atoms)} coordinate atoms vs {n} template atoms")

    elements = [a[0] for a in atoms]
    coords = [(a[1], a[2], a[3]) for a in atoms]

    def dist(i: int, j: int) -> float:
        return sum((coords[i][k] - coords[j][k]) ** 2 for k in range(3)) ** 0.5

    # Breadth-first order so every atom (after the first of each fragment) has a
    # neighbour that is already placed.
    order: list[int] = []
    seen: set[int] = set()
    for start in range(n):
        if start in seen:
            continue
        queue = [start]
        seen.add(start)
        while queue:
            a = queue.pop(0)
            order.append(a)
            for nb in template.GetAtomWithIdx(a).GetNeighbors():
                j = nb.GetIdx()
                if j not in seen:
                    seen.add(j)
                    queue.append(j)

    mapping: dict[int, int] = {}
    used: set[int] = set()

    def place(pos: int) -> bool:
        if pos == len(order):
            return True
        t = order[pos]
        t_atom = template.GetAtomWithIdx(t)
        placed_nbrs = [nb.GetIdx() for nb in t_atom.GetNeighbors() if nb.GetIdx() in mapping]
        for c in range(len(atoms)):
            if c in used or elements[c] != t_atom.GetSymbol():
                continue
            if all(min_bond <= dist(c, mapping[nb]) <= max_bond for nb in placed_nbrs):
                mapping[t] = c
                used.add(c)
                if place(pos + 1):
                    return True
                del mapping[t]
                used.discard(c)
        return False

    if not place(0):
        raise ValueError("no geometric mapping of template atoms onto the coordinates")

    mol = Chem.Mol(template)
    mol.RemoveAllConformers()
    conf = Chem.Conformer(n)
    for t, c in mapping.items():
        conf.SetAtomPosition(t, Point3D(*coords[c]))
    mol.AddConformer(conf, assignId=True)
    return mol


def _mol_from_atoms(atoms, template):
    """Build an RDKit molecule from bare coordinates, with bond orders from a template."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    try:
        mol = Chem.MolFromPDBBlock(_pdb_block(atoms), proximityBonding=True, removeHs=False)
        if mol is None:
            raise ValueError("RDKit could not build a molecule from the coordinates")
        return AllChem.AssignBondOrdersFromTemplate(template, mol)
    except ValueError as first_error:
        try:
            return _mol_from_atoms_geometric(atoms, template)
        except ValueError as second_error:
            def formula(elements) -> str:
                return " ".join(f"{el}{n}" for el, n in sorted(Counter(elements).items()))

            raise ValueError(
                f"could not match ligand atoms to the template "
                f"(RDKit: {first_error}; geometric fallback: {second_error}). "
                f"Coordinate atoms: {formula(a[0] for a in atoms)}; "
                f"template atoms: {formula(a.GetSymbol() for a in template.GetAtoms())}; "
                f"template SMILES: {Chem.MolToSmiles(template)}"
            ) from second_error


def _pose_models(pdbqt_path: str) -> tuple[list[list[tuple[str, float, float, float]]], list[float]]:
    """Parse a Vina output PDBQT into per-pose heavy-atom lists and Vina scores."""
    models: list[list[tuple[str, float, float, float]]] = []
    scores: list[float] = []
    current: list[tuple[str, float, float, float]] | None = None

    for line in Path(pdbqt_path).read_text().splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("REMARK VINA RESULT"):
            scores.append(float(line.split()[3]))
        elif line.startswith(("ATOM", "HETATM")):
            if current is None:
                current = []
            ad_type = line.split()[-1]
            element = _AD_TO_ELEMENT.get(ad_type, ad_type[:1].upper())
            if element == "H":
                continue
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            current.append((element, *xyz))
        elif line.startswith("ENDMDL"):
            if current:
                models.append(current)
            current = None
    if current:
        models.append(current)
    return models, scores


def run_redock(
    pdb_id: str,
    workdir: str | Path,
    ligand_resname: str | None = None,
    template_smiles: str | None = None,
    exhaustiveness: int = 16,
    n_poses: int = 9,
    seed: int = 42,
) -> dict:
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from rdkit.Chem.MolStandardize import rdMolStandardize

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    pdb_file = download_pdb(pdb_id, workdir)
    pocket = extract_pocket(pdb_file, workdir, ligand_resname)
    resname = pocket.ligand_resname
    log.info("Native ligand: %s  box center=%s size=%s", resname, pocket.center, pocket.size)

    crystal_atoms = _crystal_ligand_atoms(pdb_file, resname)
    smiles_for_template = template_smiles or _fetch_ccd_smiles(resname)
    template = _template_from_smiles(smiles_for_template)
    if template.GetNumHeavyAtoms() != len(crystal_atoms):
        raise ValueError(
            f"Template has {template.GetNumHeavyAtoms()} heavy atoms but the crystal "
            f"ligand has {len(crystal_atoms)}. Pass the right SMILES with --template-smiles."
        )

    crystal = _mol_from_atoms(crystal_atoms, template)
    Chem.AssignStereochemistryFrom3D(crystal)
    neutral = rdMolStandardize.Uncharger().uncharge(Chem.Mol(crystal))
    dock_smiles = Chem.MolToSmiles(neutral)
    log.info("Ligand SMILES (stereo from crystal): %s", dock_smiles)

    receptor = get_plugin("receptor_prep").run({"receptor_pdb": str(pocket.receptor_pdb)}, workdir)
    prepared = get_plugin("ligand_prep").run({"smiles": dock_smiles, "seed": seed}, workdir)

    poses_dir = workdir / "poses"
    poses_dir.mkdir(exist_ok=True)
    poses_path = poses_dir / f"{resname}_redock.pdbqt"

    log.info(
        "Docking (exhaustiveness=%d, n_poses=%d, seed=%d) ...", exhaustiveness, n_poses, seed
    )
    get_plugin("vina").run(
        {
            "receptor_pdbqt": receptor["receptor_pdbqt"],
            "ligand_pdbqt": prepared["pdbqt"],
            "center": pocket.center,
            "size": pocket.size,
            "exhaustiveness": exhaustiveness,
            "n_poses": n_poses,
            "seed": seed,
            "out_pdbqt": str(poses_path),
        },
        workdir,
    )

    models, scores = _pose_models(str(poses_path))
    if not models:
        raise RuntimeError(f"No poses found in {poses_path}")

    poses = []
    for idx, atoms in enumerate(models):
        pose_mol = _mol_from_atoms(atoms, template)
        rmsd = rdMolAlign.CalcRMS(pose_mol, crystal)  # in place, symmetry-aware, no alignment
        poses.append(
            {
                "pose": idx + 1,
                "vina_score": scores[idx] if idx < len(scores) else None,
                "rmsd": round(rmsd, 3),
            }
        )

    result = {
        "pdb": pdb_id.upper(),
        "ligand": resname,
        "ligand_smiles": dock_smiles,
        "seed": seed,
        "top_pose_rmsd": poses[0]["rmsd"],
        "best_rmsd_any_pose": min(p["rmsd"] for p in poses),
        "success_top_pose": poses[0]["rmsd"] <= SUCCESS_RMSD,
        "success_any_pose": min(p["rmsd"] for p in poses) <= SUCCESS_RMSD,
        "threshold": SUCCESS_RMSD,
        "poses": poses,
    }
    (workdir / "redock_results.json").write_text(json.dumps(result, indent=2))
    return result
