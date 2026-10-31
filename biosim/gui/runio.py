"""Load, save, copy and export docking run folders (no Qt, so it is easy to test).

A run folder is whatever `run_docking` writes: results.csv, poses/<name>.pdbqt and
<PDB>_protein.pdb. Opening one needs nothing else, so a finished run can be
reloaded without docking again.
"""
from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path

RESULT_FIELDS = [
    "rank", "name", "smiles", "best_isomer_smiles", "n_isomers",
    "vina_score_kcal_mol", "status", "notes",
]


class RunFolderError(Exception):
    """The folder is not a usable run folder, or a copy/save could not be done."""


def find_receptor(folder: Path) -> tuple[str, Path] | None:
    """(pdb_id, path) of the first <PDB>_protein.pdb in the folder, if any."""
    for path in sorted(Path(folder).glob("*_protein.pdb")):
        return path.name[: -len("_protein.pdb")], path
    return None


def _to_float(text: str | None) -> float | None:
    try:
        return float((text or "").strip())
    except ValueError:
        return None


def _to_int(text: str | None, default: int = 1) -> int:
    try:
        return int((text or "").strip())
    except ValueError:
        return default


def load_run(folder: str | Path) -> tuple[str, list[dict]]:
    """Read results.csv from a run folder. Returns (pdb_id or "", rows).

    Works with CSVs from before stereoisomer enumeration (missing columns get
    defaults) and with CSVs that already carry a `notes` column.
    """
    folder = Path(folder)
    csv_path = folder / "results.csv"
    if not folder.is_dir():
        raise RunFolderError(f"{folder} is not a folder")
    if not csv_path.is_file():
        raise RunFolderError(f"No results.csv in {folder}. Pick a folder from a finished run.")

    rows: list[dict] = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or "name" not in reader.fieldnames:
            raise RunFolderError(f"{csv_path} does not look like a biosim results file")
        for raw in reader:
            smiles = (raw.get("smiles") or "").strip()
            rows.append(
                {
                    "name": (raw.get("name") or "").strip(),
                    "smiles": smiles,
                    "best_isomer_smiles": (raw.get("best_isomer_smiles") or "").strip() or smiles,
                    "n_isomers": _to_int(raw.get("n_isomers")),
                    "vina_score_kcal_mol": _to_float(raw.get("vina_score_kcal_mol")),
                    "status": (raw.get("status") or "").strip(),
                    "notes": raw.get("notes") or "",
                }
            )
    receptor = find_receptor(folder)
    return (receptor[0] if receptor else ""), rows


def _write_rows(path: Path, rows: list[dict]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            score = row.get("vina_score_kcal_mol")
            writer.writerow(
                {
                    "rank": rank,
                    "name": row.get("name", ""),
                    "smiles": row.get("smiles", ""),
                    "best_isomer_smiles": row.get("best_isomer_smiles", ""),
                    "n_isomers": row.get("n_isomers", 1),
                    "vina_score_kcal_mol": "" if score is None else score,
                    "status": row.get("status", ""),
                    "notes": row.get("notes", ""),
                }
            )
    os.replace(tmp, path)  # never leave a half-written results.csv behind


def save_results(folder: str | Path, rows: list[dict]) -> Path:
    """Write results.csv (with notes) into the run folder, atomically."""
    folder = Path(folder)
    if not folder.is_dir():
        raise RunFolderError(f"{folder} does not exist. Use Save As to choose a new folder.")
    path = folder / "results.csv"
    _write_rows(path, rows)
    return path


def export_csv(path: str | Path, rows: list[dict]) -> Path:
    """Write the table to any CSV path (does not touch the run folder)."""
    path = Path(path)
    _write_rows(path, rows)
    return path


def copy_run(src: str | Path, dest: str | Path) -> Path:
    """Copy a whole run folder (receptor, poses, results) to a new, non-existing folder."""
    src, dest = Path(src), Path(dest)
    if not src.is_dir():
        raise RunFolderError(f"{src} does not exist, nothing to copy")
    if dest.exists():
        raise RunFolderError(f"{dest} already exists. Choose a new folder name.")
    if src.resolve() == dest.resolve() or src.resolve() in dest.resolve().parents:
        raise RunFolderError("The new folder cannot be inside the current run folder")
    shutil.copytree(src, dest)
    return dest
