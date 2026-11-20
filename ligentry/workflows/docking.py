"""End-to-end docking workflow: PDB + ligand list -> ranked results CSV.

Supports GPU batch docking via Uni-Dock with multi-core CPU ligand preparation
for high-throughput libraries (40k-50k+ ligands), as well as CPU AutoDock Vina.
"""
from __future__ import annotations

import csv
import logging
import os
import re
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import ligentry.plugins  # noqa: F401  (registers plugins)
from ligentry.core.plugin import get_plugin
from ligentry.data.pdb import download_pdb, extract_pocket
from ligentry.plugins.ligand_prep import DEFAULT_MAX_ISOMERS, enumerate_isomers

log = logging.getLogger("ligentry.docking")

CSV_FIELDS = [
    "rank", "name", "smiles", "best_isomer_smiles", "n_isomers",
    "ai_pose_score", "pred_pkd", "vina_score_kcal_mol", "status",
]


def _resolve_docking_plugin(engine: str = "auto"):
    pref = engine.lower().strip()
    if pref == "vina":
        return get_plugin("vina"), "vina"
    if pref == "unidock":
        return get_plugin("unidock"), "unidock"
    try:
        return get_plugin("unidock"), "unidock"
    except Exception:
        return get_plugin("vina"), "vina"


def _prep_one_task(task: tuple[int, int, str, str]) -> tuple[int, int, str, str | None, str | None]:
    """Worker function for parallel 3D ligand preparation (RDKit + Meeko)."""
    lig_idx, iso_idx, iso_smiles, out_pdbqt_path = task
    try:
        import ligentry.plugins  # noqa: F401
        from ligentry.core.plugin import get_plugin as _gp

        res = _gp("ligand_prep").run({"smiles": iso_smiles}, Path(out_pdbqt_path).parent)
        Path(out_pdbqt_path).write_text(res["pdbqt"])
        return lig_idx, iso_idx, iso_smiles, out_pdbqt_path, None
    except Exception as exc:
        return lig_idx, iso_idx, iso_smiles, None, str(exc)


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
    rescore: bool = True,
    rescore_device: str = "auto",
    engine: str = "auto",
    batch_size: int = 256,
) -> list[dict]:
    """Dock every ligand and return rows ranked by AI pose score (or Vina score)."""
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
    docker, active_engine = _resolve_docking_plugin(engine)
    log.info("Docking engine: %s", active_engine)

    rescorer = None
    if rescore:
        try:
            rescorer = get_plugin("pose_rescorer")
        except Exception as exc:
            log.warning("AI pose rescorer unavailable, falling back to Vina ranking: %s", exc)

    effective_n_poses = max(int(n_poses), 9) if rescorer is not None else int(n_poses)
    poses_dir = workdir / "poses"
    poses_dir.mkdir(exist_ok=True)

    rows: list[dict] = []
    stopped = False

    if active_engine == "unidock":
        # High-throughput parallel CPU conformer prep + GPU batch docking
        prep_dir = (workdir / "_prep_pdbqt").resolve()
        batch_out_dir = (workdir / "_gpu_batch_out").resolve()
        prep_dir.mkdir(exist_ok=True)
        batch_out_dir.mkdir(exist_ok=True)

        n_workers = max(1, min((os.cpu_count() or 4) - 2, 10))
        total_ligs = len(ligands)

        for chunk_start in range(0, total_ligs, batch_size):
            if stop_requested():
                stopped = True
                break
            chunk = ligands[chunk_start: chunk_start + batch_size]

            # Build isomer tasks for this chunk
            chunk_meta: list[dict] = []
            tasks: list[tuple[int, int, str, str]] = []
            for local_idx, (name, smiles) in enumerate(chunk):
                global_idx = chunk_start + local_idx + 1
                if enumerate_stereo:
                    isomers, _info = enumerate_isomers(smiles, max_isomers=max_isomers)
                else:
                    isomers = [smiles]
                safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
                chunk_meta.append({
                    "global_idx": global_idx,
                    "name": name,
                    "safe": safe,
                    "smiles": smiles,
                    "isomers": isomers,
                    "errors": [],
                    "prepared": [],
                })
                for k, iso in enumerate(isomers, start=1):
                    in_pdbqt = (prep_dir / f"b{chunk_start}_{local_idx}_{k}_{safe}.pdbqt").resolve()
                    tasks.append((local_idx, k, iso, str(in_pdbqt)))

            # Parallel CPU 3D conformer generation
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = [pool.submit(_prep_one_task, t) for t in tasks]
                for fut in as_completed(futures):
                    l_idx, k, iso_smi, pdbqt_file, err = fut.result()
                    if err or not pdbqt_file:
                        chunk_meta[l_idx]["errors"].append(err or "prep failed")
                    else:
                        chunk_meta[l_idx]["prepared"].append((k, iso_smi, Path(pdbqt_file).resolve()))

            if stop_requested():
                stopped = True
                break

            all_in_paths = [
                p.resolve() for m in chunk_meta for (_k, _iso, p) in m["prepared"]
            ]
            batch_res: dict[str, dict] = {}
            if all_in_paths:
                log.info(
                    "GPU batch docking %d conformers (%d-%d of %d ligands) ...",
                    len(all_in_paths), chunk_start + 1, chunk_start + len(chunk), total_ligs,
                )
                batch_out = docker.run(
                    {
                        "receptor_pdbqt": receptor["receptor_pdbqt"],
                        "ligand_pdbqts": all_in_paths,
                        "out_dir": str(batch_out_dir),
                        "center": pocket.center,
                        "size": pocket.size,
                        "exhaustiveness": exhaustiveness,
                        "n_poses": effective_n_poses,
                    },
                    workdir,
                )
                batch_res = batch_out["results"]

            # Rescore and collect best isomer per ligand
            for m in chunk_meta:
                log.info("[%d/%d] %s", m["global_idx"], total_ligs, m["name"])
                row = {
                    "name": m["name"],
                    "smiles": m["smiles"],
                    "best_isomer_smiles": m["smiles"],
                    "n_isomers": len(m["isomers"]),
                    "ai_pose_score": None,
                    "pred_pkd": None,
                    "vina_score_kcal_mol": None,
                    "status": "ok",
                }
                final_pose = poses_dir / f"{m['safe']}.pdbqt"
                scored_isos: list[tuple[float | None, float | None, float, str, Path]] = []

                for _k, iso_smi, in_path in m["prepared"]:
                    hit = batch_res.get(str(in_path.resolve()))
                    if not hit:
                        continue
                    out_p = Path(hit["poses_pdbqt"])
                    vina_s = float(hit["score"])
                    ai_s: float | None = None
                    pkd_s: float | None = None

                    if rescorer is not None:
                        try:
                            rs = rescorer.run(
                                {
                                    "receptor_pdbqt": receptor["receptor_pdbqt"],
                                    "poses_pdbqt": str(out_p),
                                    "out_pdbqt": str(out_p),
                                    "device": rescore_device,
                                },
                                workdir,
                            )
                            ai_s = float(rs["top_ai_score"])
                            pkd_s = float(rs["top_pred_pkd"])
                            if rs.get("top_vina_score") is not None:
                                vina_s = float(rs["top_vina_score"])
                        except Exception as rs_exc:
                            log.warning("  AI rescoring failed for %s: %s", m["name"], rs_exc)

                    scored_isos.append((ai_s, pkd_s, vina_s, iso_smi, out_p))

                if not scored_isos:
                    err_msg = m["errors"][0] if m["errors"] else "docking produced no poses"
                    row["status"] = f"failed: {err_msg}"
                else:
                    best_ai, best_pkd, best_vina, best_smi, best_p = min(
                        scored_isos,
                        key=lambda r: (r[0] is None, -(r[0] or 0.0), r[2]),
                    )
                    os.replace(best_p, final_pose)
                    row["ai_pose_score"] = round(best_ai, 4) if best_ai is not None else None
                    row["pred_pkd"] = round(best_pkd, 2) if best_pkd is not None else None
                    row["vina_score_kcal_mol"] = round(best_vina, 3)
                    row["best_isomer_smiles"] = best_smi
                    if len(scored_isos) < len(m["isomers"]):
                        row["status"] = f"ok ({len(scored_isos)}/{len(m['isomers'])} isomers docked)"

                # Clean up intermediate files for this ligand
                for _k, _iso, in_path in m["prepared"]:
                    in_path.unlink(missing_ok=True)
                    (batch_out_dir / f"{in_path.stem}_out.pdbqt").unlink(missing_ok=True)

                rows.append(row)
    else:
        # Sequential CPU AutoDock Vina fallback
        ligand_prep = get_plugin("ligand_prep")
        for i, (name, smiles) in enumerate(ligands, start=1):
            if stop_requested():
                stopped = True
                break
            log.info("[%d/%d] %s", i, len(ligands), name)
            row = {
                "name": name, "smiles": smiles, "best_isomer_smiles": smiles, "n_isomers": 1,
                "ai_pose_score": None, "pred_pkd": None, "vina_score_kcal_mol": None, "status": "ok",
            }
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            final_pose = poses_dir / f"{safe}.pdbqt"
            temp_poses: list[Path] = []
            ligand_stopped = False
            try:
                if enumerate_stereo:
                    isomers, info = enumerate_isomers(smiles, max_isomers=max_isomers)
                    if info["n_unspecified"]:
                        log.info(
                            "  %d unspecified stereo element(s) -> %d isomer(s) to dock",
                            info["n_unspecified"], len(isomers),
                        )
                else:
                    isomers = [smiles]
                row["n_isomers"] = len(isomers)

                results: list[tuple[float | None, float | None, float, str, Path]] = []
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
                        result = docker.run(
                            {
                                "receptor_pdbqt": receptor["receptor_pdbqt"],
                                "ligand_pdbqt": prepared["pdbqt"],
                                "center": pocket.center,
                                "size": pocket.size,
                                "exhaustiveness": exhaustiveness,
                                "n_poses": effective_n_poses,
                                "out_pdbqt": str(pose_path),
                            },
                            workdir,
                        )
                        vina_score = float(result["score"])
                        ai_score: float | None = None
                        pred_pkd: float | None = None

                        if rescorer is not None:
                            try:
                                rs = rescorer.run(
                                    {
                                        "receptor_pdbqt": receptor["receptor_pdbqt"],
                                        "poses_pdbqt": str(pose_path),
                                        "out_pdbqt": str(pose_path),
                                        "device": rescore_device,
                                    },
                                    workdir,
                                )
                                ai_score = float(rs["top_ai_score"])
                                pred_pkd = float(rs["top_pred_pkd"])
                                if rs.get("top_vina_score") is not None:
                                    vina_score = float(rs["top_vina_score"])
                            except Exception as rs_exc:
                                log.warning("  AI rescoring failed for %s: %s", name, rs_exc)

                        results.append((ai_score, pred_pkd, vina_score, iso, pose_path))
                    except Exception as exc:
                        last_error = exc
                        log.warning("  isomer %d failed: %s", k, exc)

                if not results:
                    raise last_error if last_error else RuntimeError("no isomers to dock")

                best_ai, best_pkd, best_vina, best_smiles, best_path = min(
                    results,
                    key=lambda r: (r[0] is None, -(r[0] or 0.0), r[2]),
                )
                if best_path != final_pose:
                    os.replace(best_path, final_pose)
                row["ai_pose_score"] = round(best_ai, 4) if best_ai is not None else None
                row["pred_pkd"] = round(best_pkd, 2) if best_pkd is not None else None
                row["vina_score_kcal_mol"] = round(best_vina, 3)
                row["best_isomer_smiles"] = best_smiles
                if ligand_stopped:
                    row["status"] = f"ok (stopped: {len(results)}/{len(isomers)} isomers docked)"
                elif len(results) < len(isomers):
                    row["status"] = f"ok ({len(results)}/{len(isomers)} isomers docked)"
            except Exception as exc:
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

    rows.sort(
        key=lambda r: (
            r["vina_score_kcal_mol"] is None,
            r["ai_pose_score"] is None,
            -(r["ai_pose_score"] or 0.0),
            r["vina_score_kcal_mol"] or 0.0,
        )
    )

    with open(workdir / "results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow({"rank": rank, **row})

    return rows
