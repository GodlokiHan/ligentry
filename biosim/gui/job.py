"""A docking job that runs in its own process.

The GUI starts this with QProcess instead of running docking in a thread. Vina's
C++ docking call can hold Python's global interpreter lock and saturate every
CPU core, which froze the Qt window ("not responding") whenever the desktop
needed it. In a separate process the GUI is never blocked, a crash in a native
library cannot take the window down, and the job runs at a lower priority.

    python -m biosim.gui.job <job.json>

job.json keys: pdb_id, chembl_target, smiles_text, n_ligands, exhaustiveness,
workdir, enumerate_stereo, max_isomers.

Progress goes to stdout as log lines; failures print "ERROR: <message>" and exit
1. To stop early, create the file <workdir>/STOP: the job finishes the current
Vina call, keeps what is done, writes results.csv and exits 0.
"""
from __future__ import annotations

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


def main(argv: list[str] | None = None, *, run=None, fetch=None) -> int:
    """Run one job. `run` and `fetch` can be injected for testing."""
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("ERROR: usage: python -m biosim.gui.job <job.json>", flush=True)
        return 2
    try:
        spec = json.loads(Path(argv[0]).read_text())
        workdir = Path(spec["workdir"])
        pdb_id = spec["pdb_id"]
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: bad job file: {exc!r}", flush=True)
        return 2

    try:
        os.nice(10)  # keep the desktop responsive while all cores are busy
    except (AttributeError, OSError):
        pass

    logger = logging.getLogger("biosim")
    handler = logging.StreamHandler(sys.stdout)  # flushes after every record
    handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    old_level = logger.level
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        chembl = (spec.get("chembl_target") or "").strip()
        if chembl:
            print(f"Fetching actives for {chembl} ...", flush=True)
            if fetch is None:
                from biosim.data.chembl import fetch_actives as fetch
            ligands = fetch(chembl, limit=int(spec.get("n_ligands", 20)))
        else:
            ligands = parse_smiles_text(spec.get("smiles_text") or "")
        if not ligands:
            raise ValueError("No ligands to dock.")

        if run is None:
            from biosim.workflows.docking import run_docking as run
        stop_path = workdir / STOP_FILE
        run(
            pdb_id, ligands, workdir,
            exhaustiveness=int(spec.get("exhaustiveness", 8)),
            enumerate_stereo=bool(spec.get("enumerate_stereo", True)),
            max_isomers=int(spec.get("max_isomers", 8)),
            should_stop=stop_path.exists,
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
