"""Molecular docking with AutoDock Vina (Python bindings)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from ligentry.core.plugin import Plugin, register


@register
class VinaDocking(Plugin):
    name = "vina"
    description = "AutoDock Vina docking of one ligand into a prepared receptor"

    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        from vina import Vina

        v = Vina(sf_name="vina", seed=int(inputs.get("seed", 42)), verbosity=0)
        v.set_receptor(inputs["receptor_pdbqt"])
        v.set_ligand_from_string(inputs["ligand_pdbqt"])
        v.compute_vina_maps(center=list(inputs["center"]), box_size=list(inputs["size"]))
        v.dock(
            exhaustiveness=int(inputs.get("exhaustiveness", 8)),
            n_poses=int(inputs.get("n_poses", 3)),
        )
        score = float(v.energies(n_poses=1)[0][0])  # kcal/mol, more negative = better
        out = Path(inputs["out_pdbqt"])
        out.parent.mkdir(parents=True, exist_ok=True)
        v.write_poses(str(out), n_poses=int(inputs.get("n_poses", 3)), overwrite=True)
        return {"score": score, "poses_pdbqt": str(out)}
