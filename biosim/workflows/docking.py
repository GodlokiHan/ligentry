"""End-to-end docking workflow: PDB + ligand list -> ranked results CSV."""
from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

import biosim.plugins  # noqa: F401  (registers plugins)
from biosim.core.plugin import get_plugin
from biosim.data.pdb import download_pdb, extract_pocket

log = logging.getLogger("biosim.docking")


def run_docking(
    pdb_id: str,
    ligands: list[tuple[str, str]],
    workdir: str | Path,
    ligand_resname: str | None = None,
    exhaustiveness: int = 8,
    n_poses: int = 3,
    padding: float = 5.0,
) -> list[dict]:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    log.info("Downloading %s", pdb_id)
    pdb_file = download_pdb(pdb_id, workdir)
    pocket = extract_pocket(pdb_file, workdir, ligand_resname, padding)
    log.info(
        "Pocket from native ligand %s: center=%s size=%s",
        pocket.ligand_resname, pocket.center, pocket.size,
    )

    receptor = get_plugin("receptor_prep").run({"receptor_pdb": str(pocket.receptor_pdb)}, workdir)
    ligand_prep = get_plugin("ligand_prep")
    vina = get_plugin("vina")

    poses_dir = workdir / "poses"
    poses_dir.mkdir(exist_ok=True)

    rows: list[dict] = []
    for i, (name, smiles) in enumerate(ligands, start=1):
        log.info("[%d/%d] %s", i, len(ligands), name)
        row = {"name": name, "smiles": smiles, "vina_score_kcal_mol": None, "status": "ok"}
        try:
            prepared = ligand_prep.run({"smiles": smiles}, workdir)
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            result = vina.run(
                {
                    "receptor_pdbqt": receptor["receptor_pdbqt"],
                    "ligand_pdbqt": prepared["pdbqt"],
                    "center": pocket.center,
                    "size": pocket.size,
                    "exhaustiveness": exhaustiveness,
                    "n_poses": n_poses,
                    "out_pdbqt": str(poses_dir / f"{safe}.pdbqt"),
                },
                workdir,
            )
            row["vina_score_kcal_mol"] = round(result["score"], 3)
        except Exception as exc:  # keep going; record the failure
            row["status"] = f"failed: {exc}"
            log.warning("  %s", row["status"])
        rows.append(row)

    rows.sort(key=lambda r: (r["vina_score_kcal_mol"] is None, r["vina_score_kcal_mol"] or 0.0))

    with open(workdir / "results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["rank", "name", "smiles", "vina_score_kcal_mol", "status"]
        )
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow({"rank": rank, **row})

    return rows
