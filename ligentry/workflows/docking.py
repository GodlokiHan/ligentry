"""End-to-end docking workflow: PDB + ligand list -> ranked results CSV."""
from __future__ import annotations

import csv
import logging
import os
import re
from collections.abc import Callable
from pathlib import Path

import ligentry.plugins  # noqa: F401  (registers plugins)
from ligentry.core.plugin import get_plugin
from ligentry.data.pdb import download_pdb, extract_pocket
from ligentry.plugins.ligand_prep import DEFAULT_MAX_ISOMERS, enumerate_isomers

log = logging.getLogger("ligentry.docking")

CSV_FIELDS = [
    "rank", "name", "smiles", "best_isomer_smiles", "n_isomers",
    "vina_score_kcal_mol", "status",
]


def run_docking(
    pdb_id: str,
    ligands: list[tuple[str, str]],
    workdir: str | Path,
    ligand_resname: str | None = None,
    exhaustiveness: int = 8,
    n_poses: int = 3,
    padding: float = 5.0,
    enumerate_stereo: bool = True,
    max_isomers: int = DEFAULT_MAX_ISOMERS,
    should_stop: Callable[[], bool] | None = None,
) -> list[dict]:
    """Dock every ligand and return rows ranked by Vina score.

    With enumerate_stereo=True, unspecified stereocentres (and double bonds) in a
    ligand's SMILES are expanded into explicit isomers (at most `max_isomers`,
    sampled reproducibly when there are more). Each isomer is docked separately
    and the best-scoring one represents the ligand: its pose is saved as
    poses/<name>.pdbqt and its SMILES goes in `best_isomer_smiles`.
    With enumerate_stereo=False the behaviour is the original one-SMILES-one-dock.

    `should_stop` is polled before each ligand and between isomers (never in the
    middle of a Vina call). When it returns True the run ends early: ligands
    already docked are kept, a ligand interrupted between isomers is kept with
    the isomers docked so far, and results.csv is still written.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    def stop_requested() -> bool:
        return bool(should_stop and should_stop())

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
    stopped = False
    for i, (name, smiles) in enumerate(ligands, start=1):
        if stop_requested():
            stopped = True
            break
        log.info("[%d/%d] %s", i, len(ligands), name)  # the GUI progress bar parses this line
        row = {
            "name": name, "smiles": smiles, "best_isomer_smiles": smiles, "n_isomers": 1,
            "vina_score_kcal_mol": None, "status": "ok",
        }
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
        final_pose = poses_dir / f"{safe}.pdbqt"
        temp_poses: list[Path] = []
        ligand_stopped = False
        try:
            if enumerate_stereo:
                isomers, info = enumerate_isomers(smiles, max_isomers=max_isomers)
                if info["n_unspecified"]:
                    msg = "  %d unspecified stereo element(s) -> %d isomer(s) to dock"
                    log.info(msg, info["n_unspecified"], len(isomers))
                    if info["sampled"]:
                        log.info(
                            "  capped: sampled %d of %d possible isomers",
                            len(isomers), info["n_possible"],
                        )
            else:
                isomers = [smiles]
            row["n_isomers"] = len(isomers)

            results: list[tuple[float, str, Path]] = []
            last_error: Exception | None = None
            for k, iso in enumerate(isomers, start=1):
                if k > 1 and stop_requested():
                    ligand_stopped = True
                    break
                if len(isomers) > 1:
                    log.info("  isomer %d/%d", k, len(isomers))
                    pose_path = poses_dir / f"{safe}__iso{k}.pdbqt"
                    temp_poses.append(pose_path)
                else:
                    pose_path = final_pose
                try:
                    prepared = ligand_prep.run({"smiles": iso}, workdir)
                    result = vina.run(
                        {
                            "receptor_pdbqt": receptor["receptor_pdbqt"],
                            "ligand_pdbqt": prepared["pdbqt"],
                            "center": pocket.center,
                            "size": pocket.size,
                            "exhaustiveness": exhaustiveness,
                            "n_poses": n_poses,
                            "out_pdbqt": str(pose_path),
                        },
                        workdir,
                    )
                    results.append((float(result["score"]), iso, pose_path))
                except Exception as exc:  # one bad isomer must not sink the ligand
                    last_error = exc
                    log.warning("  isomer %d failed: %s", k, exc)

            if not results:
                raise last_error if last_error else RuntimeError("no isomers to dock")

            best_score, best_smiles, best_path = min(results, key=lambda r: r[0])
            if best_path != final_pose:
                os.replace(best_path, final_pose)
            row["vina_score_kcal_mol"] = round(best_score, 3)
            row["best_isomer_smiles"] = best_smiles
            if ligand_stopped:
                row["status"] = f"ok (stopped: {len(results)}/{len(isomers)} isomers docked)"
            elif len(results) < len(isomers):
                row["status"] = f"ok ({len(results)}/{len(isomers)} isomers docked)"
        except Exception as exc:  # keep going; record the failure
            row["status"] = f"failed: {exc}"
            log.warning("  %s", row["status"])
        finally:
            for tmp in temp_poses:
                tmp.unlink(missing_ok=True)
        rows.append(row)
        if ligand_stopped:
            stopped = True
            break

    if stopped:
        log.info("Stopped by user: %d of %d ligand(s) docked", len(rows), len(ligands))

    rows.sort(key=lambda r: (r["vina_score_kcal_mol"] is None, r["vina_score_kcal_mol"] or 0.0))

    with open(workdir / "results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow({"rank": rank, **row})

    return rows
