"""Ligand preparation (SMILES -> 3D -> PDBQT) using RDKit and Meeko."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from biosim.core.plugin import Plugin, register

EMBED_SEEDS = (42, 7, 123, 2024, 99)


@register
class LigandPrep(Plugin):
    name = "ligand_prep"
    description = "SMILES -> protonated 3D conformer -> PDBQT string (RDKit + Meeko)"

    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.MolStandardize import rdMolStandardize

        smiles = inputs["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit could not parse SMILES")

        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)  # strip salts
        mol = Chem.AddHs(mol)

        embedded = False
        for seed in EMBED_SEEDS:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True  # far more robust for large/flexible molecules
            if AllChem.EmbedMolecule(mol, params) == 0:
                embedded = True
                break
        if not embedded:
            raise ValueError(f"3D embedding failed after {len(EMBED_SEEDS)} attempts")

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        except Exception:
            pass  # unoptimized geometry is acceptable for docking

        setups = MoleculePreparation().prepare(mol)
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            raise ValueError(f"Meeko failed: {err}")
        return {"pdbqt": pdbqt}
