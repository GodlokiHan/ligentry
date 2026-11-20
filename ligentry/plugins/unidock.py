"""GPU-accelerated molecular docking plugin wrapping Uni-Dock (CUDA Vina engine)."""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from ligentry.core.plugin import Plugin, register

log = logging.getLogger("ligentry.unidock")


def _unidock_cmd_and_env() -> tuple[list[str], dict[str, str]]:
    """Resolve the Uni-Dock executable and prepend its conda env `lib/` to LD_LIBRARY_PATH."""
    env = os.environ.copy()
    unidock_env = Path.home() / "micromamba" / "envs" / "unidock"
    unidock_lib = unidock_env / "lib"
    if unidock_lib.is_dir():
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = f"{unidock_lib}:{existing}" if existing else str(unidock_lib)

    unidock_bin = unidock_env / "bin" / "unidock"
    if unidock_bin.is_file():
        return [str(unidock_bin)], env

    if shutil.which("micromamba"):
        return ["micromamba", "run", "-n", "unidock", "unidock"], env

    return ["unidock"], env


def _parse_unidock_scores(pdbqt_path: Path) -> list[float]:
    """Extract Vina/Uni-Dock binding energies (kcal/mol) from an output PDBQT."""
    scores: list[float] = []
    if not pdbqt_path.is_file():
        return scores
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    scores.append(float(parts[3]))
                except ValueError:
                    pass
    return scores


def _materialize_ligand(item: str | Path, scratch_dir: Path, idx: int) -> Path:
    """Ensure `item` (either a file path or raw PDBQT text) is a real `.pdbqt` file on disk."""
    if isinstance(item, Path):
        return item.resolve()
    text = str(item)
    if "\n" not in text and len(text) < 512 and Path(text).is_file():
        return Path(text).resolve()
    target = scratch_dir / f"lig_{idx}.pdbqt"
    target.write_text(text)
    return target.resolve()


def _run_gpu_batch(
    receptor: Path,
    lig_paths: list[Path],
    out_dir: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int,
    n_poses: int,
    seed: int,
    scoring: str,
) -> dict[str, dict[str, Any]]:
    """Invoke Uni-Dock `--gpu_batch` in chunks safe for OS ARG_MAX limits."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cx, cy, cz = center
    sx, sy, sz = size
    base_cmd, run_env = _unidock_cmd_and_env()

    results: dict[str, dict[str, Any]] = {}
    sub_chunk = 500
    for start in range(0, len(lig_paths), sub_chunk):
        chunk = lig_paths[start: start + sub_chunk]
        cmd = [
            *base_cmd,
            "--receptor", str(receptor.resolve()),
            "--gpu_batch", *[str(p.resolve()) for p in chunk],
            "--dir", str(out_dir.resolve()),
            "--center_x", f"{cx:.3f}",
            "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", f"{sx:.3f}",
            "--size_y", f"{sy:.3f}",
            "--size_z", f"{sz:.3f}",
            "--scoring", scoring,
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses),
            "--seed", str(seed),
            "--verbosity", "1",
        ]
        proc = subprocess.run(cmd, env=run_env, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"Uni-Dock GPU batch failed (code {proc.returncode}): "
                f"{(proc.stderr or proc.stdout).strip()}"
            )

        for lig in chunk:
            out_pdbqt = out_dir / f"{lig.stem}_out.pdbqt"
            scores = _parse_unidock_scores(out_pdbqt)
            if scores:
                results[str(lig)] = {
                    "score": scores[0],
                    "scores": scores,
                    "poses_pdbqt": str(out_pdbqt),
                }
    return results


@register
class UniDocking(Plugin):
    name = "unidock"
    description = "GPU-accelerated AutoDock Vina docking via Uni-Dock (CUDA)"
    required_executables = ("unidock",)

    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        receptor = Path(inputs["receptor_pdbqt"])
        center = tuple(float(v) for v in inputs["center"])
        size = tuple(float(v) for v in inputs["size"])
        exhaustiveness = int(inputs.get("exhaustiveness", 8))
        n_poses = int(inputs.get("n_poses", 9))
        seed = int(inputs.get("seed", 42))
        scoring = str(inputs.get("scoring", "vina"))

        batch_ligands: list[str | Path] | None = inputs.get("ligand_pdbqts")
        if batch_ligands is not None:
            out_dir = Path(inputs["out_dir"])
            if not batch_ligands:
                return {"results": {}}
            with tempfile.TemporaryDirectory(prefix="unidock_in_", dir=workdir) as in_tmp:
                in_dir = Path(in_tmp)
                lig_paths = [
                    _materialize_ligand(item, in_dir, i) for i, item in enumerate(batch_ligands)
                ]
                results = _run_gpu_batch(
                    receptor=receptor,
                    lig_paths=lig_paths,
                    out_dir=out_dir,
                    center=center,  # type: ignore[arg-type]
                    size=size,      # type: ignore[arg-type]
                    exhaustiveness=exhaustiveness,
                    n_poses=n_poses,
                    seed=seed,
                    scoring=scoring,
                )
            return {"results": results}

        # Single-ligand mode (accepts either a PDBQT string from `ligand_prep` or a file path)
        out_pdbqt = Path(inputs.get("out_pdbqt", workdir / "poses.pdbqt")).resolve()
        out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix="unidock_single_", dir=out_pdbqt.parent) as tmp:
            tmp_dir = Path(tmp)
            lig_file = _materialize_ligand(inputs["ligand_pdbqt"], tmp_dir, 0)
            results = _run_gpu_batch(
                receptor=receptor,
                lig_paths=[lig_file],
                out_dir=tmp_dir,
                center=center,  # type: ignore[arg-type]
                size=size,      # type: ignore[arg-type]
                exhaustiveness=exhaustiveness,
                n_poses=n_poses,
                seed=seed,
                scoring=scoring,
            )
            hit = results.get(str(lig_file))
            if not hit:
                raise RuntimeError("Uni-Dock produced no poses for the ligand")
            os.replace(hit["poses_pdbqt"], out_pdbqt)

        scores = _parse_unidock_scores(out_pdbqt)
        return {
            "score": scores[0],
            "scores": scores,
            "poses_pdbqt": str(out_pdbqt),
        }
