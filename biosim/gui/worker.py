"""Background thread that runs a docking job without freezing the window."""
from __future__ import annotations

import logging
import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal

from biosim.data.chembl import fetch_actives
from biosim.workflows.docking import run_docking


class _SignalHandler(logging.Handler):
    def __init__(self, signal) -> None:
        super().__init__()
        self._signal = signal

    def emit(self, record: logging.LogRecord) -> None:
        self._signal.emit(self.format(record))


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


class DockWorker(QThread):
    log = Signal(str)
    docked = Signal(list)
    failed = Signal(str)

    def __init__(
        self,
        pdb_id: str,
        chembl_target: str,
        smiles_text: str,
        n_ligands: int,
        exhaustiveness: int,
        workdir: Path,
        enumerate_stereo: bool = True,
        max_isomers: int = 8,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.pdb_id = pdb_id
        self.chembl_target = chembl_target
        self.smiles_text = smiles_text
        self.n_ligands = n_ligands
        self.exhaustiveness = exhaustiveness
        self.workdir = workdir
        self.enumerate_stereo = enumerate_stereo
        self.max_isomers = max_isomers
        self._stop = threading.Event()

    def request_stop(self) -> None:
        """Ask the run to end early. It stops between Vina calls, keeping results so far."""
        self._stop.set()

    def run(self) -> None:
        logger = logging.getLogger("biosim")
        handler = _SignalHandler(self.log)
        handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
        old_level = logger.level
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        try:
            if self.chembl_target:
                self.log.emit(f"Fetching actives for {self.chembl_target} ...")
                ligands = fetch_actives(self.chembl_target, limit=self.n_ligands)
            else:
                ligands = parse_smiles_text(self.smiles_text)
            if not ligands:
                raise ValueError("No ligands to dock.")
            rows = run_docking(
                self.pdb_id, ligands, self.workdir, exhaustiveness=self.exhaustiveness,
                enumerate_stereo=self.enumerate_stereo, max_isomers=self.max_isomers,
                should_stop=self._stop.is_set,
            )
            self.docked.emit(rows)
        except Exception as exc:
            self.failed.emit(str(exc))
        finally:
            logger.removeHandler(handler)
            logger.setLevel(old_level)
