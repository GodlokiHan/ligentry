"""Download PDB structures and derive a docking pocket from the bound ligand."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import requests

# Water, ions, and common crystallization additives are never the ligand of interest.
IGNORED_HET = {
    "HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CU",
    "CD", "NI", "CO", "IOD", "BR", "SO4", "PO4", "NO3", "GOL", "EDO", "PEG",
    "ACT", "DMS", "FMT", "TRS", "MES", "BME", "MSE",
}


@dataclass
class Pocket:
    receptor_pdb: Path
    ligand_resname: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]


def download_pdb(pdb_id: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.upper()
    dest = out_dir / f"{pdb_id}.pdb"
    if dest.exists():
        return dest
    resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    resp.raise_for_status()
    dest.write_text(resp.text)
    return dest


def extract_pocket(
    pdb_path: Path,
    out_dir: Path,
    ligand_resname: str | None = None,
    padding: float = 5.0,
    min_box: float = 18.0,
) -> Pocket:
    """Write a protein-only PDB and compute a docking box around the native ligand.

    If ligand_resname is None, the largest non-ignored HETATM group is used.
    """
    protein_lines: list[str] = []
    ligand_groups: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}

    for line in pdb_path.read_text().splitlines():
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        if line[16] not in (" ", "A"):  # keep first alternate location only
            continue
        if record == "ATOM  ":
            protein_lines.append(line)
            continue
        resname = line[17:20].strip()
        if resname in IGNORED_HET:
            continue
        if ligand_resname and resname != ligand_resname.upper():
            continue
        element = line[76:78].strip() if len(line) >= 78 else ""
        if element == "H":
            continue
        key = (resname, line[21], line[22:26].strip())
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        ligand_groups.setdefault(key, []).append(xyz)

    if not protein_lines:
        raise ValueError(f"No protein atoms found in {pdb_path}")
    if not ligand_groups:
        raise ValueError(
            f"No ligand found in {pdb_path}"
            + (f" with residue name '{ligand_resname}'" if ligand_resname else "")
        )

    (resname, _chain, _resi), coords = max(ligand_groups.items(), key=lambda kv: len(kv[1]))

    mins = [min(c[i] for c in coords) for i in range(3)]
    maxs = [max(c[i] for c in coords) for i in range(3)]
    center = tuple(round((mins[i] + maxs[i]) / 2.0, 3) for i in range(3))
    size = tuple(round(max(maxs[i] - mins[i] + 2 * padding, min_box), 3) for i in range(3))

    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = out_dir / f"{pdb_path.stem}_protein.pdb"
    receptor_pdb.write_text("\n".join(protein_lines) + "\nEND\n")

    return Pocket(receptor_pdb=receptor_pdb, ligand_resname=resname, center=center, size=size)
