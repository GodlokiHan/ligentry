"""Ligand preparation (SMILES -> 3D -> PDBQT) using RDKit and Meeko."""
from __future__ import annotations

import random
from pathlib import Path
from typing import Any

from ligentry.core.plugin import Plugin, register

EMBED_SEEDS = (42, 7, 123, 2024, 99)

DEFAULT_MAX_ISOMERS = 8


def enumerate_isomers(
    smiles: str, max_isomers: int = DEFAULT_MAX_ISOMERS, seed: int = 42
) -> tuple[list[str], dict[str, Any]]:
    """Expand unspecified stereocentres / double bonds into explicit stereoisomers.

    Returns (isomer_smiles, info). Only elements that are *unspecified* in the
    input are varied; anything already specified is kept as given.

    - Salts are stripped first (largest fragment), so counterions never count.
    - If the molecule has no unspecified stereo, the input SMILES is returned
      unchanged as the single isomer.
    - If 2**n_unspecified exceeds `max_isomers`, a random sample of `max_isomers`
      is taken (reproducible for a given `seed`) and info["sampled"] is True.
    - max_isomers <= 0 means no cap (use with care: the count doubles per centre).

    info keys: n_unspecified, n_possible, n_returned, sampled.
    Raises ValueError if RDKit cannot parse the SMILES.
    """
    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers,
        StereoEnumerationOptions,
    )
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)

    n_unspecified = sum(
        1
        for si in Chem.FindPotentialStereo(mol)
        if si.specified == Chem.StereoSpecified.Unspecified
    )
    n_possible = 2**n_unspecified
    info: dict[str, Any] = {
        "n_unspecified": n_unspecified,
        "n_possible": n_possible,
        "n_returned": 1,
        "sampled": False,
    }
    if n_unspecified == 0:
        return [smiles], info

    cap = max(int(max_isomers), 0)
    info["sampled"] = cap > 0 and n_possible > cap

    def _enumerate(rand) -> list[str]:
        opts = StereoEnumerationOptions(
            tryEmbedding=False,
            onlyUnassigned=True,
            unique=True,
            maxIsomers=cap,
            rand=rand,
        )
        out: list[str] = []
        for iso in EnumerateStereoisomers(mol, options=opts):
            smi = Chem.MolToSmiles(iso)
            if smi not in out:
                out.append(smi)
        return out

    try:
        isomers = _enumerate(random.Random(seed) if info["sampled"] else None)
    except (TypeError, AttributeError):
        # older RDKit builds that do not accept a Random instance: deterministic sample
        isomers = _enumerate(None)

    info["n_returned"] = len(isomers)
    return isomers, info


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

        seeds = EMBED_SEEDS
        if "seed" in inputs:
            first = int(inputs["seed"])
            seeds = (first,) + tuple(s for s in EMBED_SEEDS if s != first)

        embedded = False
        for seed in seeds:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True  # far more robust for large/flexible molecules
            if AllChem.EmbedMolecule(mol, params) == 0:
                embedded = True
                break
        if not embedded:
            raise ValueError(f"3D embedding failed after {len(seeds)} attempts")

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        except Exception:
            pass  # unoptimized geometry is acceptable for docking

        setups = MoleculePreparation().prepare(mol)
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            raise ValueError(f"Meeko failed: {err}")
        return {"pdbqt": pdbqt}
