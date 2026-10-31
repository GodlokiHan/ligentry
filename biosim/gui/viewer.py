"""3D viewer: protein CA trace plus docked ligand pose, embedded in Qt via PyVista."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor

# AutoDock atom types -> chemical elements
_AD_TO_ELEMENT = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "O": "O", "OA": "O",
    "OS": "O", "S": "S", "SA": "S", "F": "F", "Cl": "Cl", "CL": "Cl",
    "Br": "Br", "BR": "Br", "I": "I", "P": "P", "B": "B", "Si": "Si",
    "HD": "H", "H": "H", "HS": "H",
}
_ELEMENT_COLORS = {
    "C": (140, 200, 90), "N": (60, 90, 230), "O": (230, 60, 60),
    "S": (230, 200, 50), "F": (100, 220, 100), "Cl": (50, 200, 50),
    "Br": (160, 40, 40), "P": (255, 140, 0),
}
_DEFAULT_COLOR = (200, 120, 200)
_BOND_MAX = 1.9  # Angstrom, heavy-atom bond cutoff


def read_ca_traces(pdb_path: str | Path) -> dict[str, list[tuple[float, float, float]]]:
    """C-alpha coordinates per chain, in file order."""
    chains: dict[str, list[tuple[float, float, float]]] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        chains.setdefault(line[21], []).append(xyz)
    return chains


def read_first_pose(pdbqt_path: str | Path) -> list[tuple[str, tuple[float, float, float]]]:
    """Heavy atoms (element, xyz) of the first pose in a Vina output PDBQT."""
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line in Path(pdbqt_path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ad_type = line.split()[-1]
        element = _AD_TO_ELEMENT.get(ad_type, ad_type[:1].upper())
        if element == "H":
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        atoms.append((element, xyz))
    return atoms


class MolViewer(QtInteractor):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.set_background("white")

    def show_complex(self, receptor_pdb: str | Path, pose_pdbqt: str | Path | None = None) -> None:
        self.clear()

        for pts in read_ca_traces(receptor_pdb).values():
            if len(pts) < 2:
                continue
            trace = pv.lines_from_points(np.array(pts))
            self.add_mesh(trace.tube(radius=0.4), color="lightsteelblue", smooth_shading=True)

        centroid = None
        if pose_pdbqt:
            atoms = read_first_pose(pose_pdbqt)
            if atoms:
                centroid = self._add_ligand(atoms)

        self.reset_camera()
        if centroid is not None:
            self.set_focus(centroid.tolist())
        self.render()

    def _add_ligand(self, atoms) -> np.ndarray:
        pts = np.array([a[1] for a in atoms], dtype=float)
        colors = np.array(
            [_ELEMENT_COLORS.get(el, _DEFAULT_COLOR) for el, _ in atoms], dtype=np.uint8
        )

        cells: list[int] = []
        for i in range(len(pts)):
            d = np.linalg.norm(pts[i + 1:] - pts[i], axis=1)
            for k in np.nonzero((d > 0.4) & (d < _BOND_MAX))[0]:
                cells.extend([2, i, i + 1 + int(k)])
        if cells:
            bonds = pv.PolyData(pts)
            bonds.lines = np.array(cells)
            self.add_mesh(bonds.tube(radius=0.12), color="gray")

        cloud = pv.PolyData(pts)
        cloud["rgb"] = colors
        self.add_mesh(
            cloud, scalars="rgb", rgb=True, render_points_as_spheres=True,
            point_size=14, show_scalar_bar=False,
        )
        return pts.mean(axis=0)
