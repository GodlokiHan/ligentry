"""Receptor preparation (protein PDB -> rigid PDBQT) using Open Babel."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from ligentry.core.plugin import Plugin, register


@register
class ReceptorPrep(Plugin):
    name = "receptor_prep"
    description = "Protein PDB -> rigid PDBQT with hydrogens at pH 7.4 (Open Babel)"
    required_executables = ("obabel",)

    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        src = Path(inputs["receptor_pdb"])
        out = Path(workdir) / f"{src.stem}.pdbqt"
        cmd = ["obabel", str(src), "-O", str(out), "-xr", "-p", "7.4"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not out.exists():
            raise RuntimeError(f"obabel failed: {proc.stderr.strip()}")
        return {"receptor_pdbqt": str(out)}
