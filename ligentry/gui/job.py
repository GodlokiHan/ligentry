"""A docking job that runs in its own process."""
from __future__ import annotations

import csv
import json
import logging
import os
import sys
import traceback
from pathlib import Path

STOP_FILE = "STOP"


def parse_smiles_text(text: str) -> list[tuple[str, str]]:
    """Lines of 'SMILES [name]' -> [(name, smiles)]."""
    ligands = []
    for i, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        ligands.append((parts[1] if len(parts) > 1 else f"lig{i}", parts[0]))
    return ligands


def load_ligands_file(path: str | Path, limit: int = 100_000) -> list[tuple[str, str]]:
    """Read ligands from a `.smi`, `.txt`, `.tsv`, or `.csv` file up to `limit`."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Ligand library file not found: {p}")

    suffix = p.suffix.lower()
    if suffix in (".csv", ".tsv"):
        delim = "\t" if suffix == ".tsv" else ","
        ligands: list[tuple[str, str]] = []
        with open(p, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            if reader.fieldnames:
                lower_map = {k.strip().lower(): k for k in reader.fieldnames}
                smi_col = next(
                    (lower_map[k] for k in ("smiles", "canonical_smiles", "smi") if k in lower_map),
                    None,
                )
                name_col = next(
                    (
                        lower_map[k]
                        for k in ("name", "id", "title", "molecule_chembl_id", "compound_id")
                        if k in lower_map
                    ),
                    None,
                )
                if smi_col:
                    for idx, row in enumerate(reader, start=1):
                        smi = (row.get(smi_col) or "").strip()
                        if not smi:
                            continue
                        nm = (row.get(name_col) or "").strip() if name_col else ""
                        ligands.append((nm or f"lig{idx}", smi))
                        if len(ligands) >= limit:
                            break
                    return ligands

    # Fallback: plain whitespace-separated `SMILES [name]` lines
    ligands = []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for i, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            ligands.append((parts[1] if len(parts) > 1 else f"lig{i}", parts[0]))
            if len(ligands) >= limit:
                break
    return ligands


def main(argv: list[str] | None = None, *, run=None, fetch=None) -> int:
    """Run one job. `run` and `fetch` can be injected for testing."""
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("ERROR: usage: python -m ligentry.gui.job <job.json>", flush=True)
        return 2
    try:
        spec = json.loads(Path(argv[0]).read_text())
        workdir = Path(spec["workdir"])
        pdb_id = spec["pdb_id"]
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: bad job file: {exc!r}", flush=True)
        return 2

    try:
        os.nice(10)
    except (AttributeError, OSError):
        pass

    logger = logging.getLogger("ligentry")
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    old_level = logger.level
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        limit = int(spec.get("n_ligands", 20))
        chembl = (spec.get("chembl_target") or "").strip()
        smiles_file = (spec.get("smiles_file") or "").strip()

        if chembl:
            print(f"Fetching actives for {chembl} (up to {limit}) ...", flush=True)
            if fetch is None:
                from ligentry.data.chembl import fetch_actives as fetch
            ligands = fetch(chembl, limit=limit)
        elif smiles_file:
            print(f"Loading ligands from {smiles_file} ...", flush=True)
            ligands = load_ligands_file(smiles_file, limit=limit)
        else:
            ligands = parse_smiles_text(spec.get("smiles_text") or "")[:limit]

        if not ligands:
            raise ValueError("No ligands to dock.")

        if run is None:
            from ligentry.workflows.docking import run_docking as run
        stop_path = workdir / STOP_FILE
        run(
            pdb_id,
            ligands,
            workdir,
            exhaustiveness=int(spec.get("exhaustiveness", 8)),
            enumerate_stereo=bool(spec.get("enumerate_stereo", True)),
            max_isomers=int(spec.get("max_isomers", 8)),
            should_stop=stop_path.exists,
            rescore=bool(spec.get("rescore", True)),
            rescore_device=str(spec.get("rescore_device", "auto")),
            engine=str(spec.get("engine", "auto")),
        )
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", flush=True)
        traceback.print_exc()
        return 1
    finally:
        logger.removeHandler(handler)
        logger.setLevel(old_level)


if __name__ == "__main__":
    sys.exit(main())
